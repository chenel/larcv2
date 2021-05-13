#ifndef __SUPERAG4HITSEGMENT_CXX__
#define __SUPERAG4HITSEGMENT_CXX__

#include <numeric>
#include <unordered_map>

#include "Voxelize.h"

#include "geometry.h"
#include "GenRandom.h"
#include "raybox.h"
#include "SuperaG4HitSegment.h"

namespace larcv {

  static SuperaG4HitSegmentProcessFactory __global_SuperaG4HitSegmentProcessFactory__;

  // ------------------------------------------------------

  SuperaG4HitSegment::SuperaG4HitSegment(const std::string name) : SuperaBase(name)
  { }

  // ------------------------------------------------------

  void SuperaG4HitSegment::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);
    _sparsetensor3d_producer=cfg.get<std::string>("HitTensorProducer");
    _particle_producer =cfg.get<std::string>("ParticleProducer");
    _active_volume_filter = cfg.get<std::vector<std::string> >("ActiveVolumeFilterList");
  }

  // ------------------------------------------------------

  void SuperaG4HitSegment::initialize()
  {
    SuperaBase::initialize();
  }

  // ------------------------------------------------------

  bool SuperaG4HitSegment::process(IOManager& mgr)
  {
    SuperaBase::process(mgr);

    auto evt = this->GetEvent();

    auto ev_hittensor = (EventSparseTensor3D*)(mgr.get_data("sparse3d", _sparsetensor3d_producer));
    auto ev_particles = (EventParticle*)(mgr.get_data("particle",_particle_producer));

    auto meta = ev_hittensor->meta();
    larcv::AABBox<double> box(meta);

    // build a map of which primary particles go with which interaction
    std::unordered_map<std::size_t, std::size_t> primToInt;
    for (std::size_t intIdx = 0; intIdx < evt->Primaries.size(); intIdx++)
    {
      const TG4PrimaryVertex & intxn = evt->Primaries[intIdx];
      for (const TG4PrimaryParticle & part : intxn.Particles)
        primToInt[part.GetTrackId()] = intIdx;
    }

    std::vector<larcv::Particle> particles;
    LARCV_INFO() << "There are " << evt->Trajectories.size() << " TG4Trajectories in this event" << std::endl;
    particles.reserve(evt->Trajectories.size());
    std::size_t partCounter = 0;
    for(auto const& traj : evt->Trajectories)
    {
      ++partCounter;
      const TG4Trajectory * parentTraj = nullptr;

      // work back up to the top of the hierarchy to find the primary this particle came from
      // (we need that info to determine which interaction it's part of)
      const TG4Trajectory * t = &traj;
      while (t->GetParentId() >= 0 && t->GetParentId() < static_cast<int>(evt->Trajectories.size()))
      {
        const TG4Trajectory * p = &evt->Trajectories[t->GetParentId()];
        if (!parentTraj)
          parentTraj = p;
        t = p;
      }
      std::size_t intId = primToInt[t->GetTrackId()];


      larcv::Particle part = this->MakeParticle(traj, parentTraj, box);
      part.id(partCounter);
      part.interaction_id(intId);
      LARCV_INFO() << "Made particle ID= " << part.id() << " from trajectory:" << std::endl
                   << "   In interaction=" << intId << " Track ID=" << part.track_id() << " PDG=" << part.pdg_code()
                        << " KE=" << (part.energy_init() - sqrt(part.energy_init()*part.energy_init() - part.p()*part.p()))
                        << " MeV" << std::endl
                   << "     full extent: "
                        << "[" << part.position().x() << "," << part.position().y() << "," << part.position().z() << "] => "
                        << "[" << part.end_position().x() << "," << part.end_position().y() << "," << part.end_position().z() << "]" << std::endl
                   << "     inside bounding box:"
                        << "[" << part.first_step().x() << "," << part.first_step().y() << "," << part.first_step().z() << "] => "
                        << "[" << part.last_step().x() << "," << part.last_step().y() << "," << part.last_step().z() << "]"
                   << std::endl;
      particles.emplace_back(std::move(part));
    }

    // particle parent IDs can't be assigned until the whole set of particles are assembled...
    for (std::size_t partIdx = 0; partIdx < particles.size(); partIdx++)
    {
      if (particles[partIdx].parent_id() != larcv::kINVALID_INSTANCEID)
        continue;

      // these are primaries and won't have parents anyway
      if (particles[partIdx].parent_track_id() == larcv::kINVALID_UINT)
        continue;

      // parents have to have smaller indices, so walk backwards...
      // (note that this will cause an integer underflow when it crosses 0, hence the odd-looking stop condition)
      for (std::size_t parentIdx = partIdx; parentIdx <= partIdx; parentIdx--)
      {
        if (particles[parentIdx].track_id() == particles[partIdx].parent_track_id())
        {
          particles[partIdx].parent_id(particles[parentIdx].id());
          break;
        }
      }
    }

    larcv::VoxelSet voxelSet;
    for (auto const & detPair : evt->SegmentDetectors)
    {
      LARCV_DEBUG() << "Sensitive detector: " << detPair.first << std::endl;
      LARCV_DEBUG() << "  There are " << detPair.second.size() << " hits" << std::endl;
      if (!_active_volume_filter.empty()
          && std::find(_active_volume_filter.begin(), _active_volume_filter.end(), detPair.first) == _active_volume_filter.end())
      {
        LARCV_DEBUG() << "This volume is not in the accept list (cfg param 'ActiveVolumeFilterList').  Skipping its hits." << std::endl;
        continue;
      }

      for (std::size_t hitIdx = 0; hitIdx < detPair.second.size(); hitIdx++)
      {
        LARCV_DEBUG() << "   hit number: " << hitIdx << std::endl;
        auto const & hitSeg = detPair.second[hitIdx];
        auto newVoxels = MakeVoxels(hitSeg, ev_hittensor->meta(), particles);
        LARCV_DEBUG() << "    made " << newVoxels.size() << " voxels " << std::endl;

        FluctuateEnergy(newVoxels);

        for (std::size_t voxIdx = 0; voxIdx < newVoxels.size(); voxIdx++)
          voxelSet.emplace(std::move(newVoxels[voxIdx]), true);
      }
    }

    LARCV_DEBUG() << "Made " << particles.size() << " true particles" << std::endl;
    LARCV_DEBUG() << "Made " << voxelSet.size() << " voxels, with " << voxelSet.sum() << " total" << std::endl;
    ev_particles->emplace(std::move(particles));
    ev_hittensor->emplace(std::move(voxelSet), ev_hittensor->meta());
    LARCV_DEBUG() << "Done with event." << std::endl;

    return true;
  }

  // ------------------------------------------------------

  void SuperaG4HitSegment::finalize()
  {}

  // ------------------------------------------------------

  larcv::Particle SuperaG4HitSegment::MakeParticle(const TG4Trajectory &traj,
                                                   const TG4Trajectory *parentTraj,
                                                   const larcv::AABBox<double> &bbox)
  {
    larcv::Particle res;
    res.track_id(traj.GetTrackId());
    res.pdg_code(traj.GetPDGCode());
    auto const& mom = traj.GetInitialMomentum();
    res.momentum(mom.Px(),mom.Py(),mom.Pz());
    TLorentzVector start = traj.Points.front().GetPosition();
    start.SetVect(start.Vect() * 0.1);  // convert to cm
    res.position(start.X(), start.Y(), start.Z(), start.T());
    TLorentzVector end = traj.Points.back().GetPosition();
    end.SetVect(end.Vect() * 0.1); // convert to cm
    res.end_position(end.X(), end.Y(), end.Z(), end.T());

    int creationCode = parentTraj ? traj.Points.front().GetSubprocess() : 0;
    int parentPdg = parentTraj ? parentTraj->GetPDGCode() : -1;
    res.creation_process(GetTrajCreationProc(creationCode, parentPdg));
    res.parent_track_id(traj.GetParentId());
    if (parentTraj)
      res.parent_pdg_code(parentTraj->GetPDGCode());

    LARCV_DEBUG() << "Creating particle for TG4Trajectory:" << std::endl;
    LARCV_DEBUG() << "   track id=" << traj.GetTrackId() << std::endl;
    LARCV_DEBUG() << "   parent track id=" << traj.GetParentId() << std::endl;
    LARCV_DEBUG() << "   PDG=" << traj.GetPDGCode() << std::endl;
    LARCV_DEBUG() << "   parent PDG=" << res.parent_pdg_code() << std::endl;
    LARCV_DEBUG() << "   Process/subprocess of initial step:" << traj.Points[0].GetProcess()
                  << "/" << traj.Points[0].GetSubprocess() << std::endl;
    LARCV_DEBUG() << "     --> assigned name:" << res.creation_process() << std::endl;
    LARCV_DEBUG() << "   Trajectory steps:" << std::endl;
    for (std::size_t idx = 1; idx < traj.Points.size(); idx ++)
    {
      TLorentzVector trajP1 = traj.Points[idx - 1].GetPosition();
      TLorentzVector trajP2 = traj.Points[idx].GetPosition();
      trajP1.SetVect(trajP1.Vect() * 0.1);  // convert units to cm
      trajP2.SetVect(trajP2.Vect() * 0.1);  // convert units to cm

      LARCV_DEBUG() << "      (" << trajP1.X() << "," << trajP1.Y() << "," << trajP1.Z() << ")"
                    << " -> (" << trajP2.X() << "," << trajP2.Y() << "," << trajP2.Z() << ")" << std::endl;
    }

      // this will be the first point the trajectory crosses into/out from the bounding box
    const auto FindVertex = [&](int step) -> std::unique_ptr<larcv::Vertex>
    {
      std::unique_ptr<larcv::Vertex> vtx;

      LARCV_DEBUG() << "Finding " << (step > 0 ? "entry" : "exit") << " point of trajectory " << traj.GetTrackId()
                    << " " << (step > 0 ? "into" : "out of") << " bounding box:" << std::endl;

      // note that when step < 0, "idx += step" will eventually cause an integer underflow because idx's type (size_t) is unsigned.
      // this is ok because traj.Points.size() will always be < the maximum integer after that underflow
      // and the loop will be cut off at the right place anyway.
      for (std::size_t idx = (step > 0) ? step : (traj.Points.size()-1 + step); idx < traj.Points.size(); idx += step)
      {
        TLorentzVector trajP1 = traj.Points[idx - step].GetPosition();
        TLorentzVector trajP2 = traj.Points[idx].GetPosition();
        trajP1.SetVect(trajP1.Vect() * 0.1);  // convert units to cm
        trajP2.SetVect(trajP2.Vect() * 0.1);  // convert units to cm

        LARCV_DEBUG() << "  Considering step (" << trajP1.X() << "," << trajP1.Y() << "," << trajP1.Z() << ")"
                      << " -> (" << trajP2.X() << "," << trajP2.Y() << "," << trajP2.Z() << ")" << std::endl;

        larcv::Vec3d p1, p2;
        if (!Intersections(bbox, trajP1.Vect(), trajP2.Vect(),p1, p2))
        {
          LARCV_DEBUG() << "   --> segment does not intersect bounding box." << std::endl;
          continue;
        }

        LARCV_DEBUG() << "   --> intersects bounding box at (" << p1.x << "," << p1.y << "," << p1.z << ")" << std::endl;
        vtx.reset(new larcv::Vertex(p1.x, p1.y, p1.z,
                                    trajP1.T() + (p2 - p1).length() / (trajP2 - trajP1).Vect().Mag() * (trajP2.T() - trajP1.T())) );
        break;
      }
      return vtx;
    };


    if (bbox.contain(start))
      res.first_step(res.position());
    else
    {
      auto entry = FindVertex(1);  // entry point -- step forward
      if (entry)
        res.first_step(*entry);
    }

    if (bbox.contain(end))
      res.last_step(res.end_position());
    else
    {
      auto exit = FindVertex(-1);  // exit point -- step backwards
      if (exit)
        res.last_step(*exit);
    }

    res.energy_init(mom.E());

    double distTraveled = 0;
    for (std::size_t idx = 1; idx < traj.Points.size(); idx++)
      distTraveled += (traj.Points[idx].GetPosition().Vect() - traj.Points[idx-1].GetPosition().Vect()).Mag() * 0.1;
    res.distance_travel(distTraveled);

    LARCV_DEBUG() << "   Final track startpoint: " << res.first_step().dump() << std::endl;
    LARCV_DEBUG() << "   Final track endpoint: " << res.last_step().dump() << std::endl;

    return res;
  }


  // ------------------------------------------------------

  void SuperaG4HitSegment::FluctuateEnergy(std::vector<Voxel> &)
  {
    // just placeholder for now.
    return;
  }

  // ------------------------------------------------------

  std::string SuperaG4HitSegment::GetTrajCreationProc(int creationCode, int parentPdgCode)
  {
    std::string ret;
    switch (creationCode)
    {
      case 0:
        ret = "primary";
        break;

      case TG4TrajectoryPoint::kSubtypeEMIonization:
      {
        if (abs(parentPdgCode) == 11)
          ret = "eIoni";
        else if (abs(parentPdgCode) == 13)
          ret = "muIoni";
        else if (abs(parentPdgCode) > 100)
          ret = "hIoni";
        break;
      }

      case TG4TrajectoryPoint::kSubtypeEMBremsstrahlung:
        ret = "brem";
        break;

        // G4EmProcessSubType::fAnnihilation in G4EmProcessSubType.hh
      case 5:
        ret = "annih";
        break;

      case TG4TrajectoryPoint::kSubtypeEMPhotoelectric:
        ret = "phot";
        break;

      case TG4TrajectoryPoint::kSubtypeEMComptonScattering:
        ret = "compt";
        break;

        // 'EMGammaConversion' is "pair production via linearly polarized gammas on electrons", a.k.a "triple gamma conversion"
      case TG4TrajectoryPoint::kSubtypeEMPairProdByCharged:
      case TG4TrajectoryPoint::kSubtypeEMGammaConversion:
        ret = "conv";
        break;

        // G4EmProcessSubType::fScintillation in G4EmProcessSubType.hh
      case 22:
        ret = "scint";
        break;

      case TG4TrajectoryPoint::kSubtypeHadronElastic:
        ret = "hadElas";
        break;

      case TG4TrajectoryPoint::kSubtypeHadronInelastic:
        ret = "hadInelas";
        break;

        // G4HadronicProcessType::fHadronAtRest in G4HadronicProcessType.hh
      case 151:
      {
        if (parentPdgCode == 211)
          ret = "piPlusCaptureAtRest";
        else if (parentPdgCode == -211)
          ret = "piMinusCaptureAtRest";
        break;
        // don't use a placeholder; just do "unknown" until we decide we need more
      }

        // G4DecayProcessType::DECAY in G4DecayProcessType.hh
      case 201:
        ret = "Decay";
        break;

        //if (prc == "muPairProd")
        //  grp.type = supera::kDelta;

      default:
        ret = "unknown";
        break;
    }

    return ret;
  }

  // ------------------------------------------------------


}

#endif
