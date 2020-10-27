#ifndef __SUPERAG4HITSEGMENT_CXX__
#define __SUPERAG4HITSEGMENT_CXX__

#include <numeric>

#include "Voxelize.h"

#include "geometry.h"
#include "GenRandom.h"
#include "raybox.h"
#include "SuperaG4HitSegment.h"

namespace larcv {

  static SuperaG4HitSegmentProcessFactory __global_SuperaG4HitSegmentProcessFactory__;

  SuperaG4HitSegment::SuperaG4HitSegment(const std::string name) : SuperaBase(name)
  { }

  void SuperaG4HitSegment::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);
    _sparsetensor3d_producer=cfg.get<std::string>("HitTensorProducer");
    _particle_producer =cfg.get<std::string>("ParticleProducer");
  }

  void SuperaG4HitSegment::initialize()
  {
    SuperaBase::initialize();
  }

  bool SuperaG4HitSegment::process(IOManager& mgr)
  {
    SuperaBase::process(mgr);

    auto evt = this->GetEvent();

    auto ev_hittensor = (EventSparseTensor3D*)(mgr.get_data("sparse3d", _sparsetensor3d_producer));
    auto ev_particles = (EventParticle*)(mgr.get_data("particle",_particle_producer));

    auto meta = ev_hittensor->meta();
    larcv::AABBox<double> box(meta);

    std::vector<larcv::Particle> particles;
    particles.reserve(evt->Trajectories.size());
    for(auto const& traj : evt->Trajectories) {
      larcv::Particle part = this->MakeParticle(traj, box);
      part.id(ev_particles->as_vector().size());
      LARCV_INFO() << "Track ID " << part.track_id() << " PDG " << part.pdg_code() << std::endl
                   << "     full extent: "
                        << "[" << part.position().x() << "," << part.position().y() << "," << part.position().z() << "] => "
                        << "[" << part.end_position().x() << "," << part.end_position().y() << "," << part.end_position().z() << "]" << std::endl
                   << "     inside bounding box:"
                        << "[" << part.first_step().x() << "," << part.first_step().y() << "," << part.first_step().z() << "] => "
                        << "[" << part.last_step().x() << "," << part.last_step().y() << "," << part.last_step().z() << "]"
                   << std::endl;
      particles.emplace_back(std::move(part));
      break;
    }

    larcv::VoxelSet voxelSet;
    for (auto const & detPair : evt->SegmentDetectors)
    {
      LARCV_DEBUG() << "Sensitive detector: " << detPair.first << std::endl;
      LARCV_DEBUG() << "  There are " << detPair.second.size() << " hits" << std::endl;
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

  void SuperaG4HitSegment::finalize()
  {}

  larcv::Particle SuperaG4HitSegment::MakeParticle(const TG4Trajectory& traj, const larcv::AABBox<double> & bbox)
  {
    larcv::Particle res;
    res.track_id(traj.GetTrackId());
    res.pdg_code(traj.GetPDGCode());
    auto const& mom = traj.GetInitialMomentum();
    res.momentum(mom.Px(),mom.Py(),mom.Pz());
    auto const& start = traj.Points.front().GetPosition();
    res.position(start.X()/10.,start.Y()/10.,start.Z()/10.,start.T());
    auto const& end = traj.Points.back().GetPosition();
    res.end_position(end.X()/10.,end.Y()/10.,end.Z()/10.,end.T());

    // this will be the first point the trajectory crosses into the bounding box
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
        trajP1.SetRho(trajP1.Mag() * 0.1);  // convert units to cm
        trajP2.SetRho(trajP2.Mag() * 0.1);  // convert units to cm

        LARCV_DEBUG() << "  Considering step (" << trajP1.X() << "," << trajP1.Y() << "," << trajP1.Z() << ")"
                      << " -> " << trajP2.X() << "," << trajP2.Y() << "," << trajP2.Z() << ")" << std::endl;

        larcv::Vec3d p1, p2;
        if (!Intersections(bbox, trajP1.Vect(), trajP2.Vect(),p1, p2))
        {
          LARCV_DEBUG() << "   --> does not intersect bounding box." << std::endl;
          continue;
        }

        LARCV_DEBUG() << "   --> intersects bounding box at (" << p1.x << "," << p1.y << "," << p1.z << ")" << std::endl;
        vtx.reset(new larcv::Vertex(p1.x, p1.y, p1.z,
                                    trajP1.T() + (p2 - p1).length() / (trajP2 - trajP1).Vect().Mag() * (trajP2.T() - trajP1.T())) );
        break;
      }
      return vtx;
    };

    auto entry = FindVertex(1);  // entry point -- step forward
    if (entry)
      res.first_step(*entry);

    auto exit = FindVertex(-1);  // exit point -- step backwards
    if (exit)
      res.last_step(*exit);

    res.energy_init(mom.E());
    res.parent_track_id(traj.GetParentId());

    double distTraveled = 0;
    for (std::size_t idx = 1; idx < traj.Points.size(); idx++)
      distTraveled += (traj.Points[idx].GetPosition().Vect() - traj.Points[idx-1].GetPosition().Vect()).Mag();
    res.distance_travel(distTraveled);

    return res;
  }

  void SuperaG4HitSegment::FluctuateEnergy(std::vector<Voxel> &)
  {
    // just placeholder for now.
    return;
  }


}

#endif
