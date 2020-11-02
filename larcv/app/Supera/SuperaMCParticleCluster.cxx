#include "SuperaMCParticleCluster.h"

#include <numeric>

#include "larcv/core/DataFormat/EventParticle.h"
#include "larcv/core/DataFormat/EventVoxel3D.h"
#include "larcv/core/DataFormat/Particle.h"
#include "Voxelize.h"

namespace larcv
{
  static SuperaMCParticleClusterProcessFactory __global_SuperaMCParticleClusterProcessFactory__;


  // =============================================================================================
  //  public interface
  // =============================================================================================

  // ------------------------------------------------------
  SuperaMCParticleCluster::SuperaMCParticleCluster(const std::string name)
    : SuperaBase(name)
  {}

  // ------------------------------------------------------
  void SuperaMCParticleCluster::configure(const PSet &cfg)
  {
    SuperaBase::configure(cfg);

    _input_particle_label = cfg.get<std::string>("InputParticleLabel");
    _output_label = cfg.get<std::string>("OutputLabel");
    _ref_meta3d_cluster3d = cfg.get<std::string>("Meta3DFromCluster3D", "mcst");
    _ref_meta3d_tensor3d = cfg.get<std::string>("Meta3DFromTensor3D", "");
    //_masked_true_tensor3d = cfg.get<std::string>("MaskedTrueTensor3D","");
    //_masked_true2reco_cluster3d = cfg.get<std::string>("MaskedTrue2RecoCluster3D","");
    _semantic_priority.resize((size_t) (larcv::kShapeUnknown));
    for (size_t i = 0; i < _semantic_priority.size(); ++i)
      _semantic_priority[i] = i;
    _semantic_priority = cfg.get<std::vector<size_t> >("SemanticPriority", _semantic_priority);
    _delta_size = cfg.get<size_t>("DeltaSize", 10);
    _eioni_size = cfg.get<size_t>("IonizationSize", 5);
    _compton_size = cfg.get<size_t>("ComptonSize", 10);
    _edep_threshold = cfg.get<double>("EnergyDepositThreshold", 0.01);
    _projection_id = cfg.get<int>("ProjectionID", -1);
    _use_sed = cfg.get<bool>("UseSimEnergyDeposit");
    _use_sed_points = cfg.get<bool>("UseSimEnergyDepositPoints");
    _use_true_pos = cfg.get<bool>("UseTruePosition", true);
    _check_particle_validity = cfg.get<bool>("CheckParticleValidity", true);
  }

  // ------------------------------------------------------
   void SuperaMCParticleCluster::initialize()
  {
    SuperaBase::initialize();
  }


  // ------------------------------------------------------
  bool SuperaMCParticleCluster::process(IOManager& mgr)
  {
    LARCV_INFO() << "Start processing..." << std::endl;
    SuperaBase::process(mgr);


    larcv::Voxel3DMeta meta3d;

    // load the voxel metadata
    if(!_ref_meta3d_cluster3d.empty()) {
      auto const& ev_cluster3d = mgr.get_data<larcv::EventClusterVoxel3D>(_ref_meta3d_cluster3d);
      meta3d = ev_cluster3d.meta();
    }
    else if(!_ref_meta3d_tensor3d.empty()) {
      auto const& ev_tensor3d = mgr.get_data<larcv::EventSparseTensor3D>(_ref_meta3d_tensor3d);
      meta3d = ev_tensor3d.meta();
    }

    // Build MCParticle List.
    // Note that we made Particles in SuperaG4HitSegment.  We'll use those here.
    TG4Event const *ev = GetEvent();
    const auto ev_particles = dynamic_cast<EventParticle *>(mgr.get_data("particle", _input_particle_label));
    const std::vector<larcv::Particle> & particles = ev_particles->as_vector();
    _mc_part_list.Update(particles, ev->RunId, ev->EventId);

    auto const& trackid2index = _mc_part_list.TrackIdToIndex();

    // Create ParticleGroup
    LARCV_INFO() << "Creating ParticleGroups" << std::endl;
    std::vector<supera::ParticleGroup> part_grp_v = this->CreateParticleGroups(particles);

    // Fill Voxel Information
    LARCV_INFO() << "Analyzing energy deposits" << std::endl;
    this->AnalyzeSimEnergyDeposit(meta3d, part_grp_v);

    // Merge fragments of showers
    LARCV_INFO() << "Merging: shower ionization" << std::endl;
    this->MergeShowerIonizations(part_grp_v);

    // Merge touching LEScatter showers
    LARCV_INFO() << "Merging: touching LEScatters" << std::endl;
    this->MergeShowerTouchingLEScatter(meta3d,part_grp_v, particles);

    // Apply energy threshold
    LARCV_INFO() << "Applying energy threshold" << std::endl;
    this->ApplyEnergyThreshold(part_grp_v);

    // Merge fragments of showers
    LARCV_INFO() << "Merging: shower conversions" << std::endl;
    this->MergeShowerConversion(part_grp_v);

    // Merge touching shower fragments
    // Direct parentage between kShapeShower => kShapeShower/kShapeDelta/kShapeMichel
    LARCV_INFO() << "Merging: shower family touching" << std::endl;
    this->MergeShowerFamilyTouching(meta3d,part_grp_v);

    // Merge touching shower fragments (in a family)
    LARCV_INFO() << "Merging: shower touching" << std::endl;
    this->MergeShowerTouching(meta3d,part_grp_v, particles);

    // merge too small deltas into tracks
    LARCV_INFO() << "Merging: delta rays" << std::endl;
    this->MergeShowerDeltas(part_grp_v);

    // Merge touching LEScatter showers
    LARCV_INFO() << "Merging: touching LEScatters" << std::endl;
    this->MergeShowerTouchingLEScatter(meta3d,part_grp_v, particles);

    // output containers
    std::vector<int> trackid2output(trackid2index.size(), -1);
    std::vector<int> output2trackid;

    // Assign output IDs and relationships
    this->AssignParticleGroupIDs(trackid2index, output2trackid, part_grp_v, trackid2output);

    // At this point, count total number of voxels (will be used for x-check later)
    size_t total_vs_size = std::accumulate(std::begin(part_grp_v),
                                           std::end(part_grp_v),
                                           std::size_t(0),
                                           [](std::size_t s, const supera::ParticleGroup & grp)
                                           {
                                             if (grp.valid && grp.size_all() >= 1)
                                               s += grp.vs.size();
                                             return s;
                                           });

    // For shower orphans, we need to register the most base shower particle in the output (for group)
    LARCV_INFO() << "Searching the root (group) for kShapeShower particles w/ invalid group id... ("
                 << output2trackid.size() << " particles total)" << std::endl;
    this->FixOrphanShowerGroups(particles, output2trackid, part_grp_v, trackid2output);


    // For LEScatter orphans, we need to register the immediate valid (=to be stored) particle
    LARCV_INFO() << "Searching the root (group) for kShapeLEScatter particles w/ invalid group id... ("
                 << output2trackid.size() << " particles total)" << std::endl;
    this->FixOrphanLEScatterGroups(particles, output2trackid, part_grp_v, trackid2output);


    // for shower particles with invalid parent ID, attempt a search
    LARCV_INFO() << "Searching parents for shower particles w/ invalid parent id... ("
                 << output2trackid.size() << " particles total)" << std::endl;
    this->FixInvalidParentShowerGroups(particles, part_grp_v, trackid2output, output2trackid);


    // Now sort out all parent IDs where it's simply not assigned
    // (it's ok to have invalid parent id if parent track id is not stored)
    LARCV_INFO() << "Check all output particle's parent id... ("
                 << output2trackid.size() << " particles total)" << std::endl;
    this->FixUnassignedParentGroups(part_grp_v, trackid2output, output2trackid);

    // Now loop over otuput particle list and check if any remaining group id needs to be assigned
    // Use its parent to group...
    LARCV_INFO() << "Check all output particles for their group id... ("
                 << output2trackid.size() << " particles total)" << std::endl;
    this->FixUnassignedGroups(part_grp_v, output2trackid);

    // Next handle LEScatter group id if not assigned yet
    LARCV_INFO() << "Check LEScatter particle's group id ... ("
                 << output2trackid.size() << " particles total)" << std::endl;
    this->FixUnassignedLEScatterGroups(part_grp_v, output2trackid);

    // Next loop over to find any particle for which first_step is not defined
    LARCV_INFO() << "Check any particle's first step ... ("
                 << output2trackid.size() << " particles total)" << std::endl;
    this->FixFirstStepInfo(part_grp_v, meta3d, output2trackid);

    LARCV_INFO() << "Start storing " << output2trackid.size() << " particles ..." << std::endl;

    // now loop over to create VoxelSet for compton/photoelectron
    std::vector<larcv::Particle> part_v;
    part_v.resize(output2trackid.size());
    auto event_cluster = (EventClusterVoxel3D *) (mgr.get_data("cluster3d", _output_label));
    auto event_cluster_he = (EventClusterVoxel3D *) (mgr.get_data("cluster3d", _output_label + "_highE"));
    auto event_cluster_le = (EventClusterVoxel3D *) (mgr.get_data("cluster3d", _output_label + "_lowE"));
    auto event_leftover = (EventSparseTensor3D *) (mgr.get_data("sparse3d", _output_label + "_leftover"));
    // set meta for all
    event_cluster->meta(meta3d);
    event_cluster_he->meta(meta3d);
    event_cluster_le->meta(meta3d);
    event_leftover->meta(meta3d);


    // Create cluster index tensor to help back-track semantic source particle
    auto event_cindex = (EventSparseTensor3D *) (mgr.get_data("sparse3d", _output_label + "_index"));
    event_cindex->meta(meta3d);
    larcv::VoxelSet cid_vs;
    cid_vs.reserve(total_vs_size);

    this->BuildOutputClusters(output2trackid, event_cluster, event_cluster_he, event_cluster_le, cid_vs, part_grp_v, part_v);


    // Loop to store output cluster/semantic: low energy depositions
    this->BuildLEScatterClusters(particles,
                                 trackid2index,
                                 part_grp_v,
                                 trackid2output,
                                 event_cluster,
                                 event_cluster_le,
                                 cid_vs);

    // create particle ID vs ... overlapped voxel gets higher id number
    auto const &main_vs = event_cluster_he->as_vector();
    auto const &lowe_vs = event_cluster_le->as_vector();

    // Count output voxel count and x-check
    size_t output_vs_size = 0;
    for (auto const &vs : main_vs)
      output_vs_size += vs.size();
    for (auto const &vs : lowe_vs)
      output_vs_size += vs.size();
    LARCV_INFO() << "Voxel count x-check: output = " << output_vs_size << " ... total = " << total_vs_size << std::endl;

    LARCV_INFO() << "Combining remainder into 'leftover' voxels..." << std::endl;
    ConsolidateLeftoverVoxels(meta3d, part_grp_v, total_vs_size, event_leftover, output_vs_size);

    // Loop over to find any "stil valid" supera::kIonization supera::kConversion supera::kComptonHE
    LARCV_INFO() << "Particle list" << std::endl;
    CheckForUnassignedGroups(part_grp_v, output2trackid, part_v, main_vs, lowe_vs);


    // Check the validity of particle id, group id
    if (_check_particle_validity)
      this->CheckParticleValidity(part_grp_v, part_v);

    // create semantic output in 3d
    auto event_segment = (EventSparseTensor3D *) (mgr.get_data("sparse3d", _output_label + "_semantics"));
    event_segment->meta(meta3d);
    larcv::VoxelSet semantic_vs;
    semantic_vs.reserve(total_vs_size);

    // Comptons in 3d
    for (auto const &vs : event_cluster_le->as_vector())
    {
      for (auto const &vox : vs.as_vector())
        semantic_vs.emplace(vox.id(), (float) (larcv::kShapeLEScatter), false);
    }

    // Loop over "high energy" depositions, set semantic labels
    for (size_t index = 0; index < output2trackid.size(); ++index)
    {
      auto const &vs = event_cluster_he->as_vector()[index];
      auto semantic = static_cast<size_t>(part_v[index].shape());
      for (auto const &vox : vs.as_vector())
      {
        auto const &prev = semantic_vs.find(vox.id());
        if (prev.id() == larcv::kINVALID_VOXELID)
        {
          semantic_vs.emplace(vox.id(), semantic, false);
          cid_vs.emplace(vox.id(), index, false);
        }
        else
        {
          size_t prioritized_semantic = this->SemanticPriority(prev.value(), semantic);
          if (prioritized_semantic != prev.value())
          {
            semantic_vs.emplace(vox.id(), semantic, false);
            cid_vs.emplace(vox.id(), index, false);
          }
        }
      }
    }
    // store
    assert(semantic_vs.size() == cid_vs.size());
    event_segment->emplace(std::move(semantic_vs), meta3d);
    event_cindex->emplace(std::move(cid_vs), meta3d);

    // Store output
    auto event_mcp = (EventParticle *) (mgr.get_data("particle", _output_label));
    event_mcp->emplace(std::move(part_v));

    return true;
  }

  void SuperaMCParticleCluster::CheckForUnassignedGroups(const std::vector<supera::ParticleGroup> &part_grp_v,
                                                         const std::vector<int> &output2trackid,
                                                         const std::vector<larcv::Particle> &part_v,
                                                         const std::vector<larcv::VoxelSet> &main_vs,
                                                         const std::vector<larcv::VoxelSet> &lowe_vs) const
  {
    for (size_t index = 0; index < part_v.size(); ++index)
    {
      int trackid = output2trackid[index];
      auto const &grp = part_grp_v[trackid];
      auto const &part = part_v[index];
      auto const &vs0 = main_vs[index];
      auto const &vs1 = lowe_vs[index];
      LARCV_INFO() << "Particle ID " << part.id() << " Track ID " << part.track_id() << " PDG " << part.pdg_code()
                   << " " << part.creation_process() << " ... " << part.energy_init() << " MeV => "
                   << part.energy_deposit() << " MeV "
                   << grp.trackid_v.size() << " children " << vs0.size() << " (" << vs1.size() << ") voxels"
                   << std::endl;
      LARCV_INFO() << "  Parent TrackID " << part.parent_track_id() << " PartID " << part.parent_id()
                   << " PDG " << part.parent_pdg_code() << " " << part.parent_creation_process()
                   << " Ancestor TrackID " << part.ancestor_track_id()
                   << " PDG " << part.ancestor_pdg_code() << " " << part.ancestor_creation_process() << std::endl;
      LARCV_INFO() << "  Group ID: " << part.group_id() << std::endl;
      std::stringstream ss1, ss2;

      ss1 << "  Children particle IDs: " << std::flush;
      for (auto const &child : part.children_id())
        ss1 << child << " " << std::flush;
      ss1 << std::endl;
      LARCV_INFO() << ss1.str();

      ss2 << "  Children track IDs: " << std::flush;
      for (auto const &child : grp.trackid_v)
        ss2 << child << " " << std::flush;
      ss2 << std::endl;
      LARCV_INFO() << ss2.str();

      LARCV_INFO() << "  Start: " << part.first_step().x() << " " << part.first_step().y() << " "
                   << part.first_step().z() << std::endl;
    }
    LARCV_INFO() << "... " << part_v.size() << " particles" << std::endl;
  }
  // SuperaMCParticleCluster::process()

  // =============================================================================================
  //  private methods
  // =============================================================================================

  // ------------------------------------------------------
  void SuperaMCParticleCluster::AnalyzeSimEnergyDeposit(const larcv::Voxel3DMeta& meta,
                                                        std::vector<supera::ParticleGroup>& part_grp_v)
  {
    const auto & ev = this->GetEvent();

    auto const &trackid2index = _mc_part_list.TrackIdToIndex();

    LARCV_INFO() << "Processing TG4HitSegments: " << std::endl;
    std::size_t sedep_counter = 0;
    std::size_t bad_sedep_counter = 0;
    std::set<size_t> missing_trackid;
    std::vector<larcv::Particle> noparticles;  // used in call to MakeVoxelsBelow.  Don't need to update their energy deposits...
    for (const auto & sensitiveDetPair : ev->SegmentDetectors)
    {
      for (const auto & sedep : sensitiveDetPair.second)
      {
        std::stringstream trks;
        std::for_each(std::begin(sedep.Contrib), std::end(sedep.Contrib),
                      [&trks](const int trk) { trks << " " << trk; });
        LARCV_DEBUG() << "Recording edep from tracks " << trks.str()
                      << ", total Edep=" << sedep.EnergyDeposit
                      << ", start pos=(" << sedep.Start.Vect().X() << "," << sedep.Start.Vect().Y() << "," << sedep.Start.Vect().Z() << ")"
                      << std::endl;
        sedep_counter++;


        int track_id = abs(sedep.GetPrimaryId());
        if (track_id >= ((int) (trackid2index.size())))
        {
          bad_sedep_counter++;
          missing_trackid.insert(track_id);
          continue;
        }

        for (const auto &vox : MakeVoxels(sedep, meta, noparticles))
        {
          if (vox.id() == larcv::kINVALID_VOXELID)
          {
            LARCV_DEBUG() << "Skipping edep in invalid voxel: " << vox.id()
                          << ", Edep =" << vox.value()
                          << std::endl;
            continue;
          }
          //ctr_a.insert(vox_id);

          supera::EDep pt;
          pt.x = meta.pos_x(vox.id());
          pt.y = meta.pos_y(vox.id());
          pt.z = meta.pos_z(vox.id());
          pt.t = sedep.Start.T();
          pt.e = vox.value();


          auto &grp = part_grp_v[track_id];
          if (!grp.valid) continue;

          grp.vs.emplace(vox.id(), vox.value(), true);
          grp.AddEDep(pt);
         } // for (vox)
      } // for (sedep)
    } // for (sensitiveDetPair)
    if (bad_sedep_counter)
    {
      LARCV_WARNING() << bad_sedep_counter << " / " << sedep_counter << " Edep-sim TG4HitSegments "
                      << "(from " << missing_trackid.size() << " particles) did not find corresponding MCParticle!"
                      << std::endl;
    }
  } // SuperaMCParticleCluster::AnalyzeSimEnergyDeposit()

  // ------------------------------------------------------
  void SuperaMCParticleCluster::AssignParticleGroupIDs(const std::vector<int> &trackid2index, std::vector<int> &output2trackid,
                                                       std::vector<supera::ParticleGroup> &part_grp_v,
                                                       std::vector<int> &trackid2output) const
  {
    // first create the track id list
    output2trackid.resize(trackid2index.size());
    output2trackid.clear();

    // assign particle group ID numbers and make sure they have all info set
    for (auto & grp : part_grp_v)
    {
      grp.part.energy_deposit((grp.vs.size() ? grp.vs.sum() : 0.));
      size_t output_counter = output2trackid.size();
      if (!grp.valid)
        continue;
      if (grp.size_all() < 1)
        continue;

      // Also define particle "first step" and "last step"
      auto &part = grp.part;
      auto const &first_pt = grp.first_pt;
      auto const &last_pt = grp.last_pt;
      //std::cout<<first_pt.x<< " " << first_pt.y << " " << first_pt.z << std::endl;
      if (first_pt.t != kINVALID_DOUBLE)
        part.first_step(first_pt.x, first_pt.y, first_pt.z, first_pt.t);
      if (last_pt.t != kINVALID_DOUBLE)
        part.last_step(last_pt.x, last_pt.y, last_pt.z, last_pt.t);


      if (grp.shape() == kShapeLEScatter)
        continue;

      grp.part.id(output_counter);
      trackid2output[grp.part.track_id()] = output_counter;
      for (auto const &child : grp.trackid_v)
        trackid2output[child] = output_counter;
      output2trackid.push_back(grp.part.track_id());
      ++output_counter;
    }

    // now assign relationships
    for (auto const &trackid : output2trackid)
    {
      auto &grp = part_grp_v[trackid];
      if (abs(grp.part.pdg_code()) != 11 && abs(grp.part.pdg_code()) != 22)
        continue;
      unsigned int parent_trackid = grp.part.parent_track_id();
      if (trackid2output[parent_trackid] >= 0)
      {
        /*
        if(trackid2output[parent_trackid] < 0)
    grp.part.parent_id(grp.part.id());
        else {
        */
        grp.part.parent_id(trackid2output[parent_trackid]);
        int parent_output_id = trackid2output[parent_trackid];
        int parent_id = output2trackid[parent_output_id];
        if (part_grp_v[parent_id].valid)
          part_grp_v[parent_id].part.children_id(grp.part.id());
      }
    }

    // make sure the primary particle's parent and group id are set (they are themselves)
    for (auto &grp : part_grp_v)
    {
      auto &part = grp.part;
      if (part.track_id() != part.parent_track_id())
        continue;
      part.group_id(part.id());
      part.parent_id(part.id());
    }

  } // SuperaMCParticleCluster::AssignParticleGroupIDs()

  // ------------------------------------------------------
  void SuperaMCParticleCluster::ApplyEnergyThreshold(std::vector<supera::ParticleGroup>& part_grp_v) const
  {
    // Loop again and eliminate voxels that has energy below threshold
    for (auto &grp : part_grp_v)
    {
      larcv::VoxelSet vs;
      vs.reserve(grp.vs.size());
      for (auto const &vox : grp.vs.as_vector())
      {
        if (vox.value() < _edep_threshold) continue;
        vs.emplace(vox.id(), vox.value(), true);
      }
      grp.vs = vs;
      // If compton, here decide whether it should be supera::kComptonHE (high energy)
      if (grp.type == supera::kCompton && grp.vs.size() > _compton_size)
      {
        //std::cout<<"Track ID "<<grp.part.track_id()<<" high energy compton"<<std::endl;
        grp.type = supera::kComptonHE;
      } else if (grp.type == supera::kOtherShower && grp.vs.size() > _compton_size)
      {
        //std::cout<<"Track ID "<<grp.part.track_id()<<" high energy compton"<<std::endl;
        grp.type = supera::kOtherShowerHE;
      }
    }
  } // SuperaMCParticleCluster::ApplyEnergyThreshold()

  // ------------------------------------------------------
  void SuperaMCParticleCluster::CheckParticleValidity(const std::vector<supera::ParticleGroup> &part_grp_v,
                                                      std::vector<larcv::Particle> &part_v) const
  {
    for (size_t part_index = 0; part_index < part_v.size(); ++part_index)
    {
      auto const &part = part_v[part_index];
      if (part.id() != part_index)
      {
        LARCV_CRITICAL() << "Particle index " << part_index << " holds particle w/ an ID " << part.id() << std::endl
                         << part.dump() << std::endl;
        throw std::exception();
      }
      if (part.parent_id() != kINVALID_INSTANCEID && part.parent_id() >= part_v.size())
      {
        LARCV_CRITICAL() << "Particle index " << part_index
                         << " holds particle w/ an invalid parent ID " << part.parent_id() << std::endl
                         << part.dump() << std::endl
                         << part_grp_v[part.parent_track_id()].part.dump() << std::endl;
        throw std::exception();
      }
      if (part.group_id() == kINVALID_INSTANCEID || part.group_id() >= part_v.size())
      {
        LARCV_CRITICAL() << "Particle index " << part_index
                         << " holds particle w/ an invalid group ID " << part.group_id() << std::endl
                         << part.dump() << std::endl;
        throw std::exception();
      }
      if (part.parent_id() == part.id() || part.parent_id() == kINVALID_INSTANCEID)
        continue;
      bool found = false;
      for (auto const &child : part_v[part.parent_id()].children_id())
      {
        if (child == part_index)
        {
          found = true;
          break;
        }
      }
      if (!found)
      {
        LARCV_WARNING() << "Particle index " << part_index
                        << " not in the list of the parent's children (fixing)" << std::endl;
        auto children = part_v[part.parent_id()].children_id();
        children.push_back(part_index);
        part_v[part.parent_id()].children_id(children);
      }
    }
  }

  // ------------------------------------------------------
  void SuperaMCParticleCluster::BuildLEScatterClusters(const std::vector<larcv::Particle> &particles,
                                                       const std::vector<int> &trackid2index,
                                                       std::vector<supera::ParticleGroup> &part_grp_v,
                                                       const std::vector<int> &trackid2output,
                                                       EventClusterVoxel3D *event_cluster,
                                                       EventClusterVoxel3D *event_cluster_le, VoxelSet &cid_vs) const
  {
    for (auto & grp : part_grp_v)
    {
      if (!grp.valid)
        continue;
      if (grp.size_all() < 1)
        continue;
      auto semantic = grp.shape();
      if (semantic != kShapeLEScatter)
      {
        LARCV_CRITICAL() << "Unexpected, valid, >1 pixel count particle type "
                         << semantic << " pixel count " << grp.vs.size()
                         << " (not kShapeLEScatter) at line " << __LINE__
                         << std::endl;
        std::cout << grp.part.dump() << std::endl;
        throw std::exception();
      }
      unsigned int trackid = grp.part.parent_track_id();
      int output_index = -1;
      if (trackid < trackid2output.size())
        output_index = trackid2output[trackid];
      if (output_index < 0)
      {
        // search the first direct, valid parent
        // START_HERE 2020-04-23 ... use ParentTrackIDs(), make sure output_index is not int max
        while (true)
        {
          if (trackid >= trackid2index.size() || trackid2index[trackid] < 0)
            break;
          trackid = particles[trackid2index[trackid]].parent_track_id();
          if (trackid < trackid2output.size())
          {
            output_index = trackid2output[trackid];
            break;
          }
        }
      }
      //std::cout<<"Inspecting Track ID " << grp.part.track_id() << " found output ID " << output_index << std::endl;
      if (output_index < 0)
        continue;

      /*
      # cross-check that touching LEScatter is not the same as found output_index
      for(size_t out_index=0; out_index < output2trackid.size(); ++out_index) {
	bool touching = this->IsTouching(meta3d,grp.vs,part_grp_v[output2trackid[out_index]].vs);
	if(touching && output_index == (int)(out_index)) {
	  auto const& parent = part_grp_v[output2trackid[out_index]];
	  std::cout<< "  Parent Track ID " << parent.part.track_id()
		   << " PDG " << parent.part.pdg_code()
		   << " " << parent.part.creation_process() << std::endl;
	  throw std::exception();
	}
      }
      */

      // fill 3D cluster
      auto &vs_le = event_cluster_le->writeable_voxel_set(output_index);
      auto &vs = event_cluster->writeable_voxel_set(output_index);
      for (auto const &vox : grp.vs.as_vector())
      {
        vs_le.emplace(vox.id(), vox.value(), true);
        vs.emplace(vox.id(), vox.value(), true);
        cid_vs.emplace(vox.id(), output_index, false);
      }
      grp.vs.clear_data();
    }
  }

  // ------------------------------------------------------
  void
  SuperaMCParticleCluster::BuildOutputClusters(const std::vector<int> &output2trackid,
                                               EventClusterVoxel3D *event_cluster,
                                               EventClusterVoxel3D *event_cluster_he,
                                               EventClusterVoxel3D *event_cluster_le, VoxelSet &cid_vs,
                                               std::vector<supera::ParticleGroup> &part_grp_v,
                                               std::vector<larcv::Particle> &part_v) const
  {
    event_cluster->resize(output2trackid.size());
    event_cluster_he->resize(output2trackid.size());
    event_cluster_le->resize(output2trackid.size());
    for (size_t index = 0; index < output2trackid.size(); ++index)
    {
      int trackid = output2trackid[index];
      auto &grp = part_grp_v[trackid];
      // set semantic type
      ShapeType_t semantic = grp.shape();
      if (semantic == kShapeUnknown)
      {
        LARCV_CRITICAL() << "Unexpected type while assigning semantic class: " << grp.type << std::endl;
        auto const &part = grp.part;
        LARCV_CRITICAL() << "Particle ID " << part.id() << " Type " << grp.type << " Valid " << grp.valid
                         << " Track ID " << part.track_id() << " PDG " << part.pdg_code()
                         << " " << part.creation_process() << " ... " << part.energy_init() << " MeV => "
                         << part.energy_deposit() << " MeV "
                         << grp.trackid_v.size() << " children " << grp.vs.size() << " voxels " << grp.vs.sum()
                         << " MeV" << std::endl;
        LARCV_CRITICAL() << "  Parent " << part.parent_track_id() << " PDG " << part.parent_pdg_code()
                         << " " << part.parent_creation_process() << " Ancestor " << part.ancestor_track_id()
                         << " PDG " << part.ancestor_pdg_code() << " " << part.ancestor_creation_process() << std::endl;


        throw std::exception();
      }
      // todo: how should this be fixed?  we don't have the LArSoft MCShowers
      //if(semantic == larcv::kShapeLEScatter && mcs_trackid_s.find(trackid) == mcs_trackid_s.end()) {
      //  LARCV_CRITICAL() << "Unexpected particle to be stored with a shape kShapeLEScatter!" << std::endl;
      //  this->DumpHierarchy(grp.part.track_id(),part_grp_v);
      //  throw std::exception();
      //}
      // Now, we will re-classify some of non LEScatter showers based on pixel count
      // (BUG FIX: use pca or something better)
      if (grp.vs.size() < _compton_size)
      {
        LARCV_INFO() << "Particle ID " << grp.part.id() << " PDG " << grp.part.pdg_code() << " "
                     << grp.part.creation_process() << std::endl
                     << "  ... type switching " << grp.part.shape() << " => " << kShapeLEScatter
                     << " (voxel count " << grp.vs.size() << " < " << _compton_size << ")" << std::endl;
        semantic = kShapeLEScatter;
      }
      // store the shape (semantic) type in particle
      grp.part.shape(semantic);
      // store the voxel count and energy deposit
      grp.part.num_voxels(grp.vs.size());
      grp.part.energy_deposit(grp.vs.sum());
      // set particle
      std::swap(grp.part, part_v[index]);
      grp.part = part_v[index];
      // fill 3d cluster
      event_cluster->writeable_voxel_set(index) = grp.vs;
      if (semantic != kShapeLEScatter)
        event_cluster_he->writeable_voxel_set(index) = grp.vs;
      else
      {
        event_cluster_le->writeable_voxel_set(index) = grp.vs;
        for (auto const &vox : grp.vs.as_vector())
          cid_vs.emplace(vox.id(), index, false);
      }
      //grp.vs.clear_data();
      grp.valid = false;
    }
  }

  // ------------------------------------------------------
  void SuperaMCParticleCluster::ConsolidateLeftoverVoxels(const Voxel3DMeta &meta3d,
                                                          const std::vector<supera::ParticleGroup> &part_grp_v,
                                                          size_t total_vs_size, EventSparseTensor3D *event_leftover,
                                                          size_t output_vs_size) const
  {
    VoxelSet leftover_vs;
    leftover_vs.reserve(total_vs_size - output_vs_size);

    if (total_vs_size > output_vs_size)
    {
      int ctr = 0;
      for (auto &grp : part_grp_v)
      {
        if (grp.size_all() < 1)
          continue;
        for (auto const &vox : grp.vs.as_vector())
          leftover_vs.emplace(vox.id(), vox.value(), true);
        ctr++;
        auto const &part = grp.part;
        LARCV_INFO() << "Particle ID " << part.id() << " Type " << grp.type << " Valid " << grp.valid << " Track ID "
                     << part.track_id() << " PDG " << part.pdg_code()
                     << " " << part.creation_process() << " ... " << part.energy_init() << " MeV => "
                     << part.energy_deposit() << " MeV "
                     << grp.trackid_v.size() << " children " << grp.vs.size() << " voxels " << grp.vs.sum() << " MeV"
                     << std::endl;
        LARCV_INFO() << "  Parent " << part.parent_track_id() << " PDG " << part.parent_pdg_code() << " "
                     << part.parent_creation_process()
                     << " Ancestor " << part.ancestor_track_id() << " PDG " << part.ancestor_pdg_code() << " "
                     << part.ancestor_creation_process() << std::endl;
        LARCV_INFO() << "  Group ID: " << part.group_id() << std::endl;
        std::stringstream ss1, ss2;

        ss1 << "  Children particle IDs: " << std::flush;
        for (auto const &child : part.children_id())
          ss1 << child << " " << std::flush;
        ss1 << std::endl;
        LARCV_INFO() << ss1.str();

        ss2 << "  Children track IDs: " << std::flush;
        for (auto const &child : grp.trackid_v)
          ss2 << child << " " << std::flush;
        ss2 << std::endl;
        LARCV_INFO() << ss2.str();
        LARCV_INFO() << "Above was supposed to be merged..." << std::endl;

      }
      LARCV_INFO() << "... " << ctr << " particles" << std::endl;
    }
    event_leftover->emplace(std::move(leftover_vs), meta3d);
  }

  // ------------------------------------------------------
  std::vector<supera::ParticleGroup>
  SuperaMCParticleCluster::CreateParticleGroups(const std::vector<larcv::Particle>& particles)
  {
    const larcv::Particle invalid_part;

    auto const &parent_pdg_v = _mc_part_list.ParentPdgCode();
    auto const &trackid2index = _mc_part_list.TrackIdToIndex();
    std::vector<supera::ParticleGroup> result(trackid2index.size());
    for (size_t index = 0; index < particles.size(); ++index)
    {

      auto const &mcpart = particles[index];
      int pdg_code = abs(mcpart.pdg_code());
      int mother_index = -1;
      unsigned int track_id = mcpart.track_id();
      if (mcpart.ancestor_track_id() < trackid2index.size())
        mother_index = trackid2index[mcpart.ancestor_track_id()];

      //if(pdg_code != -11 && pdg_code != 11 && pdg_code != 22) continue;
      if (pdg_code > 1000000) continue;

      supera::ParticleGroup grp;
      grp.part = mcpart;

      if (mother_index >= 0)
        grp.part.parent_pdg_code(parent_pdg_v[index]);
      grp.valid = true;

      if (pdg_code == 22 || pdg_code == 11)
      {
        if (pdg_code == 22)
        {
          // photon ... reset first, last, and end position
          grp.type = supera::kPhoton;
          grp.part.first_step(invalid_part.first_step());
          grp.part.last_step(invalid_part.last_step());
          grp.part.end_position(invalid_part.end_position());
        }
        else if (pdg_code == 11)
        {

          const std::string & prc = mcpart.creation_process();
          if (prc == "muIoni" || prc == "hIoni" || prc == "muPairProd")
            grp.type = supera::kDelta;
          else if (prc == "muMinusCaptureAtRest" || prc == "muPlusCaptureAtRest" || prc == "Decay")
            grp.type = supera::kDecay;
          else if (prc == "compt")
            grp.type = supera::kCompton;
          else if (prc == "phot")
            grp.type = supera::kPhotoElectron;
          else if (prc == "eIoni")
            grp.type = supera::kIonization;
          else if (prc == "conv")
            grp.type = supera::kConversion;
          else if (prc == "primary")
            grp.type = supera::kPrimary;
          else
            grp.type = supera::kOtherShower;
        }
        result[track_id] = grp;
      }
      else
      {
        grp.type = supera::kTrack;
        if (grp.part.pdg_code() == 2112)
          grp.type = supera::kNeutron;
        result[track_id] = grp;
      }
      /*
      std::cout<<"Track ID " << grp.part.track_id() << " PDG " << grp.part.pdg_code() << " " << grp.part.creation_process()
	       <<" ... parent Track ID " << grp.part.parent_track_id() << " PDG " << grp.part.parent_pdg_code() << std::endl;
      */
    }

    // fill parentage information

    return result;
  } // SuperaMCParticleCluster::CreateParticleGroups()

  // ------------------------------------------------------
  void SuperaMCParticleCluster::DumpHierarchy(size_t trackid,
                                              const std::vector<supera::ParticleGroup>& part_grp_v) const
  {
    assert(trackid < part_grp_v.size());

    auto const &grp = part_grp_v[trackid];
    LARCV_DEBUG() << std::endl << "#### Dumping particle record for track id "
              << grp.part.track_id() << " ####" << std::endl;
    LARCV_DEBUG() << "id " << grp.part.id() << " from " << grp.part.parent_id() << std::endl
              << "children: " << std::flush;
    for (auto const &child : grp.part.children_id())
      LARCV_DEBUG() << child << " " << std::flush;
    LARCV_DEBUG() << std::endl;
    LARCV_DEBUG() << grp.part.dump() << std::endl;

    size_t parent_trackid = grp.part.parent_track_id();
    while (parent_trackid < part_grp_v.size())
    {

      auto const &parent = part_grp_v[parent_trackid];
      LARCV_DEBUG() << "Parent's group id: " << parent.part.group_id() << " valid? " << parent.valid << std::endl;
      LARCV_DEBUG() << "Parent's children: " << std::flush;
      for (auto const &child : parent.part.children_id())
        LARCV_DEBUG() << child << " " << std::flush;
      LARCV_DEBUG() << std::endl;
      LARCV_DEBUG() << parent.part.dump() << std::endl;
      if (parent_trackid == parent.part.parent_track_id())
        break;
      if (parent_trackid == larcv::kINVALID_UINT)
        break;
      parent_trackid = parent.part.parent_track_id();
    }
    LARCV_DEBUG() << std::endl << std::endl << "#### Dump done ####" << std::endl;
  } // SuperaMCParticleCluster::DumpHierarchy()

  // ------------------------------------------------------
  void SuperaMCParticleCluster::FixFirstStepInfo(std::vector<supera::ParticleGroup> &part_grp_v,
                                                 const Voxel3DMeta &meta3d,
                                                 const std::vector<int> &output2trackid) const
  {
    for (int output_index : output2trackid)
    {
      auto &grp = part_grp_v[output_index];
      auto const &fs = grp.part.first_step();
      if (fs.x() != 0. || fs.y() != 0. || fs.z() != 0. || fs.t() != 0.)
        continue;
      auto const vtx = grp.part.position().as_point3d();
      double min_dist = fabs(kINVALID_DOUBLE);
      Point3D min_pt;
      for (auto const &vox : grp.vs.as_vector())
      {
        auto const pt = meta3d.position(vox.id());
        double dist = pt.squared_distance(vtx);
        if (dist > min_dist)
          continue;
        min_dist = dist;
        min_pt = pt;
      }
      if (min_dist > (sqrt(3.) + 1.e-3))
        grp.part.first_step(kINVALID_DOUBLE, kINVALID_DOUBLE, kINVALID_DOUBLE,
                            kINVALID_DOUBLE);
      else
        grp.part.first_step(min_pt.x, min_pt.y, min_pt.z, grp.part.position().t());

    }
  }
  // ------------------------------------------------------

  void SuperaMCParticleCluster::FixInvalidParentShowerGroups(const std::vector<larcv::Particle> &particles,
                                                             std::vector<supera::ParticleGroup> &part_grp_v,
                                                             std::vector<int> &trackid2output,
                                                             std::vector<int> &output2trackid) const
  {
    for (size_t out_index = 0; out_index < output2trackid.size(); ++out_index)
    {
      int trackid = output2trackid[out_index];
      auto &grp = part_grp_v[trackid];
      if (!grp.valid)
        continue;
      if (grp.part.parent_id() != kINVALID_INSTANCEID)
        continue;
      if (grp.shape() != kShapeShower)
        continue;
      LARCV_INFO() << "Analyzing particle id " << out_index << " trackid " << trackid << std::endl
                   << grp.part.dump() << std::endl;
      int parent_partid = -1;
      unsigned int parent_trackid;
      auto parent_trackid_v = ParentTrackIDs(trackid, particles);
      for (unsigned int idx : parent_trackid_v)
      {
        parent_trackid = idx;
        if (trackid2output[parent_trackid] < 0 || !part_grp_v[parent_trackid].valid)
          continue;
        auto const &parent = part_grp_v[parent_trackid].part;
        // shower parent can be either shower, michel, or delta
        if (parent.shape() == kShapeMichel ||
            parent.shape() == kShapeDelta ||
            parent.shape() == kShapeShower)
          parent_partid = parent.id();
        break;
      }
      /*
      int own_partid = grp.part.id();
      // initiate a search of parent in the valid output particle
      int parent_trackid = grp.part.parent_track_id();
      int parent_partid  = -1;
      while(1) {
	if(parent_trackid >= ((int)(trackid2index.size())) || trackid2index[parent_trackid] <0)
	  break;
	if(parent_trackid < ((int)(trackid2output.size())) &&
	   trackid2output[parent_trackid] >= 0 &&
	   part_grp_v[parent_trackid].valid ) {
	  //parent_partid = trackid2output[parent_trackid];
	  parent_partid = part_grp_v[parent_trackid].part.id();
	  break;
	}
	parent_trackid = larmcp_v[trackid2index[parent_trackid]].Mother();
      }
      */
      if (parent_partid >= 0)
      {
        // assert the group is same
        auto &parent = part_grp_v[output2trackid[parent_partid]];
        if (grp.part.group_id() == kINVALID_INSTANCEID)
        {
          grp.part.group_id(parent.part.group_id());
          for (auto const &child_id : grp.part.children_id())
          {
            auto &child = part_grp_v[output2trackid[child_id]];
            child.part.group_id(parent.part.group_id());
          }
        }
        else
        {
          assert(grp.part.group_id() == part_grp_v[output2trackid[parent_partid]].part.group_id());
        }
        grp.part.parent_id(parent_partid);
        part_grp_v[parent_trackid].part.children_id(grp.part.id());
        LARCV_INFO() << "PartID " << grp.part.id() << " (output index " << out_index << ") assigning parent "
                     << parent_partid << std::endl;
      }
      else
      {
        grp.part.parent_id(grp.part.id());
        if (grp.part.group_id() == kINVALID_INSTANCEID)
          grp.part.group_id(grp.part.id());
        for (auto const &child_id : grp.part.children_id())
        {
          auto &child = part_grp_v[output2trackid[child_id]];
          child.part.group_id(grp.part.id());
        }
        LARCV_INFO() << "PartID " << grp.part.id() << " (output index " << out_index
                     << ") assigning itself as a parent..." << std::endl;
      }
    }
  }


  // ------------------------------------------------------
  void SuperaMCParticleCluster::FixOrphanShowerGroups(const std::vector<larcv::Particle> &particles,
                                                      std::vector<int> &output2trackid,
                                                      std::vector<supera::ParticleGroup> &part_grp_v,
                                                      std::vector<int> &trackid2output) const
  {
    for (size_t out_index = 0; out_index < output2trackid.size(); ++out_index)
    {

      int trackid = output2trackid[out_index];
      auto &grp = part_grp_v[trackid];
      if (!grp.valid)
        continue;
      if (grp.part.group_id() != kINVALID_INSTANCEID)
        continue;
      if (grp.shape() != kShapeShower)
        continue;
      LARCV_INFO() << " #### SHOWER ROOT SEARCH: Analyzing a particle index " << out_index
                   << " track id " << grp.part.track_id() << std::endl
                   << grp.part.dump() << std::endl;

      auto parent_trackid_v = ParentTrackIDs(trackid, particles);
      int root_id = grp.part.id();
      unsigned int root_trackid = grp.part.track_id();
      bool stop = false;
      std::vector<size_t> intermediate_trackid_v;
      intermediate_trackid_v.push_back(trackid);
      for (auto const &parent_trackid : parent_trackid_v)
      {
        auto const &parent = part_grp_v[parent_trackid];
        switch (parent.shape())
        {
          case kShapeShower:
          case kShapeMichel:
          case kShapeDelta:
            // group candidate: check if it is "valid" = exists in the output
            root_trackid = parent_trackid;
            root_id = trackid2output[root_trackid];
            if (root_id >= 0)
            {
              // found the valid group: stop the loop
              stop = true;
              // If not, root_id will be a new output index
            }
            else
            {
              root_id = output2trackid.size();
              // If this particle is invalid, this also needs the group id.
              // Add to intermediate_id_v list so we can set the group id for all of them
              intermediate_trackid_v.push_back(root_trackid);
            }
            stop = (stop || parent.shape() != kShapeShower);
            break;
          case kShapeTrack:
            stop = true;
            break;
          case kShapeUnknown:
            break;
          case kShapeLEScatter:
          case kShapeGhost:
            /*
            LARCV_CRITICAL() << "Unexpected type found while searching for kShapeShower orphans's root!" << std::endl;
            this->DumpHierarchy(trackid,part_grp_v);
            throw std::exception();
            */
            break;
        }
        if (stop)
          break;
      }
      if (root_id < ((int) (output2trackid.size())) && trackid2output[root_trackid] != (int) (root_id))
      {
        LARCV_CRITICAL() << "Logic error for the search of shower root particle for an orphan..." << std::endl
                         << "This particle id=" << out_index << " and track_id=" << trackid << std::endl
                         << "ROOT particle id=" << root_id << " and track_id=" << root_trackid << std::endl;
        DumpHierarchy(trackid, part_grp_v);
        throw std::exception();
      }

      if (((int) (output2trackid.size())) <= root_id)
      {
        output2trackid.push_back(root_trackid);
        // Register the root parent to the output
        LARCV_INFO() << "Adding a new particle to the output to define a group..." << std::endl
                     << "ROOT particle id=" << root_id << " and track_id=" << root_trackid << std::endl
                     << part_grp_v[root_trackid].part.dump() << std::endl;
      }
      assert((size_t) (root_id) < output2trackid.size());

      auto &root = part_grp_v[root_trackid];
      //root.valid = true;
      assert(root.valid);
      root.part.id(root_id);
      root.part.group_id(root_id);
      trackid2output[root_trackid] = root_id;
      for (auto const &child_id : root.part.children_id())
      {
        auto &child = part_grp_v[output2trackid[child_id]];
        if (!child.valid)
          continue;
        assert(child.part.group_id() == kINVALID_INSTANCEID || child.part.group_id() == root_id);
        child.part.group_id(root_id);
      }
      // Set the group ID for THIS + intermediate particles
      for (auto const &child_trackid : intermediate_trackid_v)
      {
        auto &child = part_grp_v[child_trackid];
        if (!child.valid)
          continue;
        assert(child.part.group_id() == kINVALID_INSTANCEID || child.part.group_id() == root_id);
        child.part.group_id(root_id);
      }

      LARCV_INFO() << "... after update ... " << std::endl
                   << part_grp_v[trackid].part.dump() << std::endl;
    }
  }

  // ------------------------------------------------------
  void SuperaMCParticleCluster::FixUnassignedGroups(std::vector<supera::ParticleGroup> &part_grp_v,
                                                    std::vector<int> &output2trackid) const
  {
    for (size_t output_index = 0; output_index < output2trackid.size(); ++output_index)
    {
      auto &grp = part_grp_v[output2trackid[output_index]];
      if (grp.part.group_id() != kINVALID_INSTANCEID)
        continue;
      auto shape = grp.shape();
      auto parent_shape = kShapeUnknown;
      auto parent_partid = grp.part.parent_id();
      //auto parent_groupid = larcv::kINVALID_INSTANCEID;
      // If delta, its own grouping

      switch (shape)
      {
        case kShapeLEScatter:
          // if LEScatter, we handle later (next loop)
          break;
        case kShapeDelta:
        case kShapeMichel:
        case kShapeTrack:
          // If delta, Michel, or track, it's own group
          grp.part.group_id(output_index);
          for (auto const &child_index : grp.part.children_id())
          {
            part_grp_v[output2trackid[child_index]].part.group_id(output_index);
          }
          break;

        case kShapeShower:
          // If shower && no parent, consider it as a primary = assign group id for all children
          if (parent_partid == kINVALID_INSTANCEID)
          {
            grp.part.group_id(output_index);
            for (auto const &child_index : grp.part.children_id())
              part_grp_v[output2trackid[child_index]].part.group_id(output_index);
            continue;
          }
          parent_shape = part_grp_v[output2trackid[parent_partid]].shape();
          switch (parent_shape)
          {
            case kShapeMichel:
            case kShapeDelta:
              grp.part.group_id(parent_partid);
              for (auto const &child_index : grp.part.children_id())
              {
                part_grp_v[output2trackid[child_index]].part.group_id(parent_partid);
              }
              break;
            case kShapeTrack:
              grp.part.group_id(output_index);
              for (auto const &child_index : grp.part.children_id())
              {
                part_grp_v[output2trackid[child_index]].part.group_id(output_index);
              }
              break;
            case kShapeShower:
              LARCV_CRITICAL() << "Unexpected case: a shower has no group id while being a child of another shower..."
                               << std::endl;
              DumpHierarchy(grp.part.track_id(), part_grp_v);
              throw std::exception();
              /*
              // COMMENTED OUT as this is no longer expected
              parent_groupid = part_grp_v[output2trackid[parent_partid]].part.group_id();
              if(parent_groupid != larcv::kINVALID_INSTANCEID) {
                grp.part.group_id(parent_groupid);
                for(auto const& child_index : grp.part.children_id()) {
                  part_grp_v[output2trackid[child_index]].part.group_id(parent_groupid);
                }
              }
              */
              break;
            case kShapeLEScatter:
              LARCV_ERROR() << "Logic error: shower parent shape cannot be LEScatter!" << std::endl;
              throw std::exception();
            default:
              LARCV_ERROR() << "Unexpected larcv::ShapeType_t encountered at " << __LINE__ << std::endl;
              throw std::exception();
          }
          break;
        case kShapeGhost:
        case kShapeUnknown:
          LARCV_ERROR() << "Unexpected larcv::ShapeType_t encountered at " << __LINE__ << std::endl;
          throw std::exception();
      }
    }
  }

  // ------------------------------------------------------
  void SuperaMCParticleCluster::FixUnassignedLEScatterGroups(std::vector<supera::ParticleGroup> &part_grp_v,
                                                             const std::vector<int> &output2trackid) const
  {
    for (size_t output_index = 0; output_index < output2trackid.size(); ++output_index)
    {
      auto &grp = part_grp_v[output2trackid[output_index]];
      if (grp.part.group_id() != kINVALID_INSTANCEID)
        continue;
      if (grp.shape() == kShapeLEScatter)
      {
        // assign parent's group, otherwise leave as is = kINVALID_INSTANCEID
        auto parent_partid = grp.part.parent_id();
        if (parent_partid == kINVALID_INSTANCEID)
          continue;
        grp.part.group_id(part_grp_v[output2trackid[parent_partid]].part.group_id());
      }
    }
  }


  // ------------------------------------------------------

  void SuperaMCParticleCluster::FixUnassignedParentGroups(std::vector<supera::ParticleGroup> &part_grp_v,
                                                          std::vector<int> &trackid2output,
                                                          std::vector<int> &output2trackid) const
  {
    for (int output_index : output2trackid)
    {
      auto &grp = part_grp_v[output_index];
      auto parent_trackid = grp.part.parent_track_id();
      auto parent_id = grp.part.parent_id();
      auto &parent = part_grp_v[parent_trackid].part;
      // if parent_id is invalid, try if parent_trackid can help out
      if (parent_id == kINVALID_INSTANCEID &&
          trackid2output[parent_trackid] >= 0)
      {
        parent_id = trackid2output[parent_trackid];
        grp.part.parent_id(parent_id);
      }
      if (parent_id == kINVALID_INSTANCEID)
        continue;
      // if parent id is set, make sure this particle is in the children
      auto children = parent.children_id();
      bool add = true;
      for (auto const &child : children)
      {
        if (child != grp.part.id())
          continue;
        add = false;
        break;
      }
      if (add)
      {
        children.push_back(grp.part.id());
        parent.children_id(children);
      }
    }
  }

  // ------------------------------------------------------

 void SuperaMCParticleCluster::FixOrphanLEScatterGroups(const std::vector<larcv::Particle> &particles,
                                                         const std::vector<int> &output2trackid,
                                                         std::vector<supera::ParticleGroup> &part_grp_v,
                                                         std::vector<int> &trackid2output) const
  {
    for (size_t out_index = 0; out_index < output2trackid.size(); ++out_index)
    {
      int trackid = output2trackid[out_index];
      auto &grp = part_grp_v[trackid];
      if (grp.part.group_id() != kINVALID_INSTANCEID)
        continue;
      if (grp.shape() != kShapeLEScatter)
        continue;
      LARCV_INFO() << " #### LEScatter ROOT SEARCH #### " << std::endl
                   << " Analyzing a particle index " << out_index << " id " << grp.part.id() << std::endl
                   << grp.part.dump() << std::endl;

      auto parent_trackid_v = ParentTrackIDs(trackid, particles);
      size_t group_id = kINVALID_INSTANCEID;
      bool stop = false;
      for (auto const &parent_trackid : parent_trackid_v)
      {
        auto const &parent = part_grp_v[parent_trackid];
        switch (parent.shape())
        {
          case kShapeShower:
          case kShapeMichel:
          case kShapeDelta:
          case kShapeTrack:
          case kShapeLEScatter:
            // group candidate: check if it is "valid" = exists in the output
            if (parent.valid)
            {
              group_id = trackid2output[parent_trackid];
              // found the valid group: stop the loop
              stop = true;
            }
            break;
          case kShapeUnknown:
          case kShapeGhost:
            LARCV_CRITICAL() << "Unexpected type found while searching for kShapeLEScatter orphans's root!"
                             << std::endl;
            throw std::exception();
            break;
        }
        if (stop)
          break;
      }
      if (group_id == kINVALID_INSTANCEID)
      {
        LARCV_INFO() << "Ignoring kShapeLEScatter particle as its root particle (for group id) is not to be stored..."
                     << std::endl
                     << grp.part.dump() << std::endl;
        continue;
      }
      LARCV_INFO() << "Assigning a group ID " << group_id << " to kShapeLEScatter orphan" << std::endl
                   << "  Track ID " << grp.part.track_id() << " PDG " << grp.part.pdg_code()
                   << " " << grp.part.creation_process() << std::endl;
      grp.part.group_id(group_id);
      trackid2output[trackid] = group_id;
    }
  }

  // ------------------------------------------------------
  bool SuperaMCParticleCluster::IsTouching(const Voxel3DMeta& meta, const VoxelSet& vs1, const VoxelSet& vs2) const
  {

    bool touching = false;
    size_t ix1, iy1, iz1;
    size_t ix2, iy2, iz2;
    size_t diffx, diffy, diffz;

    for (auto const &vox1 : vs1.as_vector())
    {
      meta.id_to_xyz_index(vox1.id(), ix1, iy1, iz1);
      for (auto const &vox2 : vs2.as_vector())
      {
        meta.id_to_xyz_index(vox2.id(), ix2, iy2, iz2);
        if (ix1 > ix2) diffx = ix1 - ix2; else diffx = ix2 - ix1;
        if (iy1 > iy2) diffy = iy1 - iy2; else diffy = iy2 - iy1;
        if (iz1 > iz2) diffz = iz1 - iz2; else diffz = iz2 - iz1;
        touching = diffx <= 1 && diffy <= 1 && diffz <= 1;
        if (touching)
        {
          //std::cout<<"Touching ("<<ix1<<","<<iy1<<","<<iz1<<") ("<<ix2<<","<<iy2<<","<<iz2<<")"<<std::endl;
          break;
        }
      }
      if (touching) break;
    }

    return touching;
  } // SuperaMCParticleCluster::IsTouching()

  // ------------------------------------------------------
  void SuperaMCParticleCluster::MergeShowerConversion(std::vector<supera::ParticleGroup>& part_grp_v)
  {
    int merge_ctr = 0;
    int invalid_ctr = 0;
    do
    {
      merge_ctr = 0;
      for (auto &grp : part_grp_v)
      {
        if (!grp.valid) continue;
        //if(grp.type != supera::kIonization && grp.type != supera::kConversion) continue;
        if (grp.type != supera::kConversion) continue;
        // merge to a valid "parent"
        bool parent_found = false;
        unsigned int parent_index = grp.part.parent_track_id();
        unsigned int parent_index_before = grp.part.track_id();
        while (true)
        {
          //std::cout<< "Inspecting: " << grp.part.track_id() << " => " << parent_index << std::endl;
          // todo: this is always false since parent_track_id() always returns an unsigned int.  what's the "unknown parent" value?
          if (parent_index < 0)
          {
            LARCV_DEBUG() << "Invalid parent track id " << parent_index
                          << " Could not find a parent for " << grp.part.track_id() << " PDG " << grp.part.pdg_code()
                          << " " << grp.part.creation_process() << " E = " << grp.part.energy_init()
                          << " (" << grp.part.energy_deposit() << ") MeV" << std::endl;
            auto const &parent = part_grp_v[parent_index_before].part;
            LARCV_DEBUG() << "Previous parent: " << parent.track_id() << " PDG " << parent.pdg_code()
                          << " " << parent.creation_process()
                          << std::endl;
            parent_found = false;
            invalid_ctr++;
            break;
            //throw std::exception();
          }
          auto const &parent = part_grp_v[parent_index];
          parent_found = parent.valid;
          if (parent_found) break;
          else
          {
            unsigned int ancestor_index = parent.part.parent_track_id();
            if (ancestor_index == parent_index)
            {
              LARCV_INFO() << "Particle " << parent_index << " is root and invalid particle..." << std::endl
                           << "PDG " << parent.part.pdg_code() << " " << parent.part.creation_process() << std::endl;
              break;
            }
            parent_index_before = parent_index;
            parent_index = ancestor_index;
          }
        }
        // if parent is found, merge
        if (parent_found)
        {
          auto &parent = part_grp_v[parent_index];
          parent.Merge(grp);
          merge_ctr++;
        }
      }
      LARCV_INFO() << "Merge counter: " << merge_ctr << " invalid counter: " << invalid_ctr << std::endl;
    } while (merge_ctr > 0);
  }  // SuperaMCParticleCluster::MergeShowerConversion()

  // ------------------------------------------------------
  void SuperaMCParticleCluster::MergeShowerDeltas(std::vector<supera::ParticleGroup>& part_grp_v) const
  {
    for (auto &grp : part_grp_v)
    {
      //if(grp.type != supera::kDelta) continue;
      if (grp.shape() != larcv::kShapeDelta) continue;
      unsigned int parent_trackid = grp.part.parent_track_id();
      auto &parent = part_grp_v[parent_trackid];
      if (!parent.valid) continue;

      // if voxel count is smaller than delta ray requirement, simply merge
      if (grp.vs.size() < _delta_size)
      {
        // if parent is found, merge
        /*
        if(grp.vs.size()>0) {
          std::cout<<"Merging delta " << grp.part.track_id() << " PDG " << grp.part.pdg_code()
             << " " << grp.part.creation_process() << " vox count " << grp.vs.size() << std::endl
             <<" ... parent found " << parent.part.track_id()
             << " PDG " << parent.part.pdg_code() << " " << parent.part.creation_process()
             << std::endl;
          for(auto const& vs : grp.vs2d_v)
            std::cout<<vs.size() << " " << std::flush;
          std::cout<<std::endl;
        }
        */
        parent.Merge(grp);
      }
      else
      {
        // check unique number of voxels
        size_t unique_voxel_count = 0;
        for (auto const &vox : grp.vs.as_vector())
        {
          if (parent.vs.find(vox.id()).id() == larcv::kINVALID_VOXELID)
            ++unique_voxel_count;
        }
        if (unique_voxel_count < _delta_size)
        {
          // if parent is found, merge
          /*
          if(unique_voxel_count>0) {
            std::cout<<"Merging delta " << grp.part.track_id() << " PDG " << grp.part.pdg_code()
               << " " << grp.part.creation_process() << " vox count " << grp.vs.size() << std::endl
               <<" ... parent found " << parent.part.track_id()
               << " PDG " << parent.part.pdg_code() << " " << parent.part.creation_process()
               << std::endl;
            for(auto const& vs : grp.vs2d_v)
              std::cout<<vs.size() << " " << std::flush;
            std::cout<<std::endl;
          }
          */
          parent.Merge(grp);
        }
        /*
        else{
          std::cout<<"NOT merging delta " << grp.part.track_id() << " PDG " << grp.part.pdg_code()
             << " " << grp.part.creation_process() << " vox count " << grp.vs.size() << std::endl
             <<" ... parent found " << parent.part.track_id()
             << " PDG " << parent.part.pdg_code() << " " << parent.part.creation_process()
             << std::endl;
          for(auto const& vs : grp.vs2d_v)
            std::cout<<vs.size() << " " << std::flush;
          std::cout<<std::endl;

        }
        */
      }
    }
  } // SuperaMCParticleCluster::MergeShowerDeltas()

  // ------------------------------------------------------
  void SuperaMCParticleCluster::MergeShowerFamilyTouching(const larcv::Voxel3DMeta& meta,
                                                          std::vector<supera::ParticleGroup>& part_grp_v)
  {
    // Merge touching shower fragments
    // Direct parentage between kShapeShower => kShapeShower/kShapeDelta/kShapeMichel
    int merge_ctr = 0;
    int invalid_ctr = 0;
    do {
      merge_ctr = 0;
      for(auto& grp : part_grp_v) {
        if(!grp.valid) continue;
        if(grp.shape() != larcv::kShapeShower) continue;
        // search for a possible parent
        int parent_trackid = -1;
        // a direct parent ?
        if(part_grp_v[grp.part.parent_track_id()].valid)
          parent_trackid = grp.part.parent_track_id();
        else {
          for(size_t shower_trackid = 0; shower_trackid<part_grp_v.size(); ++shower_trackid) {
            auto const& candidate_grp = part_grp_v[shower_trackid];
            if(shower_trackid == grp.part.parent_track_id() || !candidate_grp.valid) continue;
            for(auto const& trackid : candidate_grp.trackid_v) {
              if(trackid != grp.part.parent_track_id()) continue;
              parent_trackid = shower_trackid;
              break;
            }
            if(parent_trackid >= 0) break;
          }
        }
        if(parent_trackid < 0 || parent_trackid == (int)(grp.part.track_id())) continue;
        auto& parent = part_grp_v[parent_trackid];
        //auto parent_type = part_grp_v[parent_trackid].type;
        //if(parent_type == supera::kTrack || parent_type == supera::kNeutron) continue;
        if(parent.shape() != larcv::kShapeShower && parent.shape() != larcv::kShapeDelta && parent.shape() != larcv::kShapeMichel) continue;
        if(this->IsTouching(meta,grp.vs,parent.vs)) {
          // if parent is found, merge
          parent.Merge(grp);
          merge_ctr++;
        }
      }
      LARCV_INFO() << "Merge counter: " << merge_ctr << " invalid counter: " << invalid_ctr << std::endl;
    }while(merge_ctr>0);
  }

  // ------------------------------------------------------
  void SuperaMCParticleCluster::MergeShowerIonizations(std::vector<supera::ParticleGroup>& part_grp_v)
  {
    // Loop over particles of a type kIonization (=touching to its parent physically by definition)
    // If a parent is found, merge to the parent
    int merge_ctr = 0;
    int invalid_ctr = 0;
    do
    {
      merge_ctr = 0;
      for (auto &grp : part_grp_v)
      {
        if (!grp.valid) continue;
        if (grp.type != supera::kIonization) continue;
        // merge to a valid "parent"
        bool parent_found = false;
        unsigned int parent_index = grp.part.parent_track_id();
        unsigned int parent_index_before = grp.part.track_id();
        while (true)
        {
          //std::cout<< "Inspecting: " << grp.part.track_id() << " => " << parent_index << std::endl;
          // todo: this test is always false since parent_track_id() returns an unsigned int. what's the "unknown parent" value?
          if (parent_index < 0)
          {
            LARCV_DEBUG() << "Invalid parent track id " << parent_index
                          << " Could not find a parent for " << grp.part.track_id() << " PDG " << grp.part.pdg_code()
                          << " " << grp.part.creation_process() << " E = " << grp.part.energy_init()
                          << " (" << grp.part.energy_deposit() << ") MeV" << std::endl;
            auto const &parent = part_grp_v[parent_index_before].part;
            LARCV_DEBUG() << "Previous parent: " << parent.track_id() << " PDG " << parent.pdg_code()
                          << " " << parent.creation_process()
                          << std::endl;
            parent_found = false;
            invalid_ctr++;
            break;
            //throw std::exception();
          }
          auto const &parent = part_grp_v[parent_index];
          parent_found = parent.valid;
          if (parent_found) break;
          else
          {
            unsigned int ancestor_index = parent.part.parent_track_id();
            if (ancestor_index == parent_index)
            {
              LARCV_INFO() << "Particle " << parent_index << " is root and invalid particle..." << std::endl
                           << "PDG " << parent.part.pdg_code() << " " << parent.part.creation_process() << std::endl;
              break;
            }
            parent_index_before = parent_index;
            parent_index = ancestor_index;
          }
        }
        // if parent is found, merge
        if (parent_found)
        {
          auto &parent = part_grp_v[parent_index];
          parent.Merge(grp);
          merge_ctr++;
        }
      } // for (grp)
      LARCV_INFO() << "Ionization merge counter: " << merge_ctr << " invalid counter: " << invalid_ctr << std::endl;
    } while (merge_ctr > 0);
  } // SuperaMCParticleCluster::MergeShowerIonizations()

  // ------------------------------------------------------
  void SuperaMCParticleCluster::MergeShowerTouching(const larcv::Voxel3DMeta& meta,
                                                    std::vector<supera::ParticleGroup>& part_grp_v,
                                                    const std::vector<larcv::Particle> & particles)
  {
    // Go over all pair-wise combination of two shower instances
    // For each shower, find all consecutive parents of shower/michel/delta type (break if track found)
    // If there is a common parent in two list AND if two showers are physically touching, merge
    int merge_ctr = 0;
    do
    {
      merge_ctr = 0;
      for (size_t i = 0; i < part_grp_v.size(); ++i)
      {
        auto &grp_a = part_grp_v[i];
        if (!grp_a.valid) continue;
        if (grp_a.shape() != larcv::kShapeShower) continue;
        for (size_t j = 0; j < part_grp_v.size(); ++j)
        {
          if (i == j) continue;
          auto &grp_b = part_grp_v[j];
          if (!grp_b.valid) continue;
          if (grp_b.shape() != larcv::kShapeShower) continue;

          // check if these showers share the parentage
          // list a's parents
          size_t trackid = i;
          std::set<size_t> parent_list_a;
          std::set<size_t> parent_list_b;
          /*
          while(1){
            auto const& parent_a = part_grp_v[trackid];
            if(parent_a.part.parent_track_id() >= part_grp_v.size())
              break;
            if(parent_a.part.parent_track_id() == parent_a.part.track_id())
              break;
            trackid = parent_a.part.parent_track_id();
            if(parent_a.shape() == larcv::kShapeMichel ||
               parent_a.shape() == larcv::kShapeShower ||
               parent_a.shape() == larcv::kShapeDelta )
              parent_list_a.insert(trackid);
            else if(parent_a.shape() == larcv::kShapeTrack ||
              parent_a.shape() == larcv::kShapeUnknown)
              break;

            if(trackid < part_grp_v.size() && part_grp_v[trackid].part.parent_track_id() == trackid)
              break;
          }
          */
          auto parents_a = this->ParentShowerTrackIDs(trackid, part_grp_v, particles);
          for (auto const &parent_trackid : parents_a) parent_list_a.insert(parent_trackid);
          parent_list_a.insert(trackid);

          trackid = j;
          /*
          while(1){
            auto const& parent_b = part_grp_v[trackid];
            if(parent_b.part.parent_track_id() >= part_grp_v.size())
              break;
            if(parent_b.part.parent_track_id() == parent_b.part.track_id())
              break;
            trackid = parent_b.part.parent_track_id();
            if(parent_b.shape() == larcv::kShapeMichel ||
               parent_b.shape() == larcv::kShapeShower ||
               parent_b.shape() == larcv::kShapeDelta )
              parent_list_b.insert(trackid);
            else if(parent_b.shape() == larcv::kShapeTrack ||
              parent_b.shape() == larcv::kShapeUnknown)
              break;
            if(trackid < part_grp_v.size() && part_grp_v[trackid].part.parent_track_id() == trackid)
              break;
          }
          */
          auto parents_b = this->ParentShowerTrackIDs(trackid, part_grp_v, particles);
          for (auto const &parent_trackid : parents_b) parent_list_b.insert(parent_trackid);
          parent_list_b.insert(trackid);

          bool merge = false;
          for (auto const &parent_trackid : parent_list_a)
          {
            if (parent_list_b.find(parent_trackid) != parent_list_b.end())
              merge = true;
            if (merge) break;
          }
          for (auto const &parent_trackid : parent_list_b)
          {
            if (parent_list_a.find(parent_trackid) != parent_list_a.end())
              merge = true;
            if (merge) break;
          }

          if (merge && this->IsTouching(meta, grp_a.vs, grp_b.vs))
          {
            if (grp_a.vs.size() < grp_b.vs.size())
              grp_b.Merge(grp_a);
            else
              grp_a.Merge(grp_b);
            merge_ctr++;
          }
        }
      }
      LARCV_INFO() << "Merge counter: " << merge_ctr << std::endl;
    } while (merge_ctr > 0);
  } // SuperaMCParticleCluster::MergeShowerTouching()

  // ------------------------------------------------------
  void SuperaMCParticleCluster::MergeShowerTouchingLEScatter(const larcv::Voxel3DMeta& meta,
                                                             std::vector<supera::ParticleGroup>& part_grp_v,
                                                             const std::vector<larcv::Particle> & particles)
  {
    size_t merge_ctr = 1;
    while (merge_ctr)
    {
      merge_ctr = 0;
      for (auto &grp : part_grp_v)
      {
        if (!grp.valid || grp.vs.size() < 1 || grp.shape() != larcv::kShapeLEScatter) continue;
        // Find all direct shower-type or other LEScatter type parent
        //auto const& parents = this->ParentShowerTrackIDs(grp.part.track_id(), part_grp_v, true);
        auto const &parents = this->ParentTrackIDs(grp.part.track_id(), particles);
        /*
        std::cout<<"Inspecting LEScatter Track ID " << grp.part.track_id()
           << " PDG " << grp.part.pdg_code()
           << " " << grp.part.creation_process() << std::endl;
        std::cout<< "  ... parents:"<<std::flush;
        for(auto const& parent_trackid : parents) std::cout<<" "<<parent_trackid;
        std::cout<<std::endl;
        */
        for (auto const &parent_trackid : parents)
        {
          auto &parent = part_grp_v[parent_trackid];
          if (!parent.valid || parent.vs.size() < 1) continue;
          if (this->IsTouching(meta, grp.vs, parent.vs))
          {
            parent.Merge(grp);
            merge_ctr++;
            break;
          }
        } // for (parent_trackid)
      } // for (grp)
    } // while (merge_ctr)
  } // SuperaMCParticleCluster::MergeShowerTouchingLEScatter()

  // ------------------------------------------------------
  std::vector<unsigned int>
  SuperaMCParticleCluster::ParentShowerTrackIDs(size_t trackid,
                                                const std::vector<supera::ParticleGroup>& part_grp_v,
                                                const std::vector<larcv::Particle> & particles,
                                                bool include_lescatter) const
  {
    auto parents = this->ParentTrackIDs(trackid, particles);
    std::vector<unsigned int> result;
    result.reserve(parents.size());

    for(auto const& parent_id : parents) {

      if(parent_id >= part_grp_v.size()) continue;

      auto const& grp = part_grp_v[parent_id];
      if(!grp.valid) continue;

      if(grp.shape() == larcv::kShapeTrack ||
         grp.shape() == larcv::kShapeUnknown)
        break;

      if(grp.shape() == larcv::kShapeMichel ||
         grp.shape() == larcv::kShapeShower ||
         grp.shape() == larcv::kShapeDelta  ||
         (grp.shape() == larcv::kShapeLEScatter && include_lescatter))
        result.push_back(parent_id);
    }
    return result;
  }

  // ------------------------------------------------------
  std::vector<unsigned int>
  SuperaMCParticleCluster::ParentTrackIDs(size_t trackid,
                                          const std::vector<larcv::Particle> & particles) const
  {
    auto const &trackid2index = _mc_part_list.TrackIdToIndex();
    std::vector<unsigned int> result;

    if (trackid >= trackid2index.size() || trackid2index[trackid] < 0)
      return result;

    unsigned int parent_trackid = particles[trackid2index[trackid]].parent_track_id();
    std::set<int> accessed;
    while (parent_trackid > 0 &&
           (size_t) (parent_trackid) < trackid2index.size() &&
           trackid2index[parent_trackid] >= 0)
    {
      if (accessed.find(parent_trackid) != accessed.end())
      {
        LARCV_CRITICAL() << "LOOP-LOGIC-ERROR for ParentTrackIDs for track id " << trackid << std::endl;
        for (size_t parent_cand_idx = 0; parent_cand_idx < result.size(); ++parent_cand_idx)
        {
          auto const &parent_cand_trackid = result[parent_cand_idx];
          auto const &mcp = particles[trackid2index.at(parent_cand_trackid)];
          LARCV_CRITICAL() << "Parent " << parent_cand_idx
                           << " Track ID " << mcp.track_id() //<< " (" << parent_trackid << ")"
                           << " PDG " << mcp.pdg_code()
                           << " Process " << mcp.creation_process()
                           << " Mother " << mcp.parent_track_id() << std::endl;
        }
        throw std::exception();
      }

      auto const &parent = particles[trackid2index[parent_trackid]];
      result.push_back(parent_trackid);
      accessed.insert(parent_trackid);
      if (parent.parent_track_id() == parent_trackid) break;
      parent_trackid = parent.parent_track_id();
    }
    return result;
  }

  // ------------------------------------------------------
  size_t SuperaMCParticleCluster::SemanticPriority(size_t a, size_t b) const
  {
    if (a == b)
      return a;
    for (auto const &semantic : _semantic_priority)
    {
      if (a == semantic)
        return a;
      if (b == semantic)
        return b;
    }
    return a;
  }

}
