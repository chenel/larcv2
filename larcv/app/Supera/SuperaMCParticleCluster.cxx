#include "SuperaMCParticleCluster.h"

#include <numeric>

#include "larcv/core/DataFormat/EventParticle.h"
#include "larcv/core/DataFormat/EventVoxel3D.h"
#include "larcv/core/DataFormat/Particle.h"
#include "Voxelize.h"

namespace larcv
{
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

    // todo: do we need these?
    auto tpc_v = cfg.get<std::vector<unsigned short> >("TPCList");
    auto cryostat_v = cfg.get<std::vector<unsigned short> >("CryostatList");
    auto plane_v = cfg.get<std::vector<unsigned short> >("PlaneList");

    // todo: how do we establish world bounds without a stored geometry?
//    larcv::Point3D min_pt(1.e9,1.e9,1.e9);
//    larcv::Point3D max_pt(-1.e9,-1.e9,-1.e9);
//    for(auto const& tpc_id : tpc_v) {
//      auto geop = lar::providerFrom<geo::Geometry>();
//      for(size_t c=0; c<geop->Ncryostats(); ++c) {
//        auto const& cryostat = geop->Cryostat(c);
//        if(!cryostat.HasTPC(tpc_id)) continue;
//        auto const& tpcabox = cryostat.TPC(tpc_id).ActiveBoundingBox();
//        if(min_pt.x > tpcabox.MinX()) min_pt.x = tpcabox.MinX();
//        if(min_pt.y > tpcabox.MinY()) min_pt.y = tpcabox.MinY();
//        if(min_pt.z > tpcabox.MinZ()) min_pt.z = tpcabox.MinZ();
//        if(max_pt.x < tpcabox.MaxX()) max_pt.x = tpcabox.MaxX();
//        if(max_pt.y < tpcabox.MaxY()) max_pt.y = tpcabox.MaxY();
//        if(max_pt.z < tpcabox.MaxZ()) max_pt.z = tpcabox.MaxZ();
//        break;
//      }
//    }
//    _world_bounds.update(min_pt,max_pt);


    // todo: do we need this stuff that contains info about the geometry?
//    _scan.clear();
//    auto geop = lar::providerFrom<geo::Geometry>();
//    _scan.resize(geop->Ncryostats());
//    for (size_t cryoid = 0; cryoid < _scan.size(); ++cryoid)
//    {
//      auto const &cryostat = geop->Cryostat(cryoid);
//      auto &scan_cryo = _scan[cryoid];
//      scan_cryo.resize(cryostat.NTPC());
//      for (size_t tpcid = 0; tpcid < scan_cryo.size(); ++tpcid)
//      {
//        auto const &tpc = cryostat.TPC(tpcid);
//        auto &scan_tpc = scan_cryo[tpcid];
//        scan_tpc.resize(tpc.Nplanes(), -1);
//      }
//    }
//    //for(size_t cryo_id=0; cryo_id<_scan.size(); ++cryo_id){
//    _valid_nplanes = 0;
//    for (auto const &cryo_id : cryostat_v)
//    {
//      auto const &cryostat = geop->Cryostat(cryo_id);
//      for (auto const &tpc_id : tpc_v)
//      {
//        if (!cryostat.HasTPC(tpc_id)) continue;
//        auto const &tpc = cryostat.TPC(tpc_id);
//        for (auto const &plane_id : plane_v)
//        {
//          if (!tpc.HasPlane(plane_id)) continue;
//          _scan[cryo_id][tpc_id][plane_id] = _valid_nplanes;
//          //std::cout<<cryo_id<<" "<<tpc_id<<" "<<plane_id<<" ... "<<_valid_nplanes<< std::endl;
//          ++_valid_nplanes;
//        }
//      }
//    }
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
    const auto ev_particles = dynamic_cast<EventParticle *>(mgr.get_data("particle", _particle_producer));
    const std::vector<larcv::Particle> & particles = ev_particles->as_vector();
    _mc_part_list.Update(particles, ev->RunId, ev->EventId);

    auto const& trackid2index = _mc_part_list.TrackIdToIndex();

    // Create ParticleGroup
    LARCV_INFO() << "Creating ParticleGroups" << std::endl;
    std::vector<supera::ParticleGroup> part_grp_v = this->CreateParticleGroups(particles);

    // Fill Voxel Information
    LARCV_INFO() << "Analyzing energy deposits" << std::endl;
    this->AnalyzeSimEnergyDeposit(meta3d, part_grp_v, mgr);

    // Merge fragments of showers
    LARCV_INFO() << "Merging: shower ionization" << std::endl;
    this->MergeShowerIonizations(part_grp_v);

    // Merge touching LEScatter showers
    LARCV_INFO() << "Merging: touching LEScatters" << std::endl;
    this->MergeShowerTouchingLEScatter(meta3d,part_grp_v);

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
    this->MergeShowerTouching(meta3d,part_grp_v);

    // Merge touching showers in 2D
    LARCV_INFO() << "Merging: shower touching (2d)" << std::endl;
    this->MergeShowerTouching2D(part_grp_v);

    // merge too small deltas into tracks
    LARCV_INFO() << "Merging: delta rays" << std::endl;
    this->MergeShowerDeltas(part_grp_v);

    // Merge touching LEScatter showers
    LARCV_INFO() << "Merging: touching LEScatters" << std::endl;
    this->MergeShowerTouchingLEScatter(meta3d,part_grp_v);

    // Assign output IDs
    // For particles in MCShower/MCTrack collection, make sure to keep them
    std::set<unsigned int> mcs_trackid_s;
    auto const& mcs_v = LArData<supera::LArMCShower_t>();
    for(auto const& mcs : mcs_v) {mcs_trackid_s.insert(mcs.TrackID());}
    std::vector<int> trackid2output(trackid2index.size(),-1);
    std::vector<int> output2trackid;
    output2trackid.reserve(trackid2index.size());
    for(size_t trackid=0; trackid<part_grp_v.size(); trackid++) {
      auto& grp = part_grp_v[trackid];
      grp.part.energy_deposit((grp.vs.size() ? grp.vs.sum() : 0.));
      size_t output_counter = output2trackid.size();
      if(mcs_trackid_s.find(trackid) != mcs_trackid_s.end()){
        grp.valid=true;
        grp.part.group_id(output_counter);
      }
      else{
        if(!grp.valid) continue;
        if(grp.size_all()<1) continue;
        if(grp.shape() == larcv::kShapeLEScatter) continue;
      }

      grp.part.id(output_counter);
      trackid2output[grp.part.track_id()] = output_counter;
      for(auto const& child : grp.trackid_v)
        trackid2output[child] = output_counter;
      output2trackid.push_back(grp.part.track_id());
      ++output_counter;
    }


    // Assign relationships
    for(auto const& trackid : output2trackid) {
      auto& grp = part_grp_v[trackid];
      if( abs(grp.part.pdg_code()) != 11 && abs(grp.part.pdg_code()) != 22 ) continue;
      int parent_trackid = grp.part.parent_track_id();
      if(trackid2output[parent_trackid] >= 0) {
        /*
        if(trackid2output[parent_trackid] < 0)
    grp.part.parent_id(grp.part.id());
        else {
        */
        grp.part.parent_id(trackid2output[parent_trackid]);
        int parent_output_id = trackid2output[parent_trackid];
        int parent_id = output2trackid[parent_output_id];
        if(part_grp_v[parent_id].valid)
          part_grp_v[parent_id].part.children_id(grp.part.id());
      }
    }


    // At this point, count total number of voxels (will be used for x-check later)
    size_t total_vs_size = 0;
    for(auto& grp : part_grp_v) {
      if(!grp.valid) continue;
      if(grp.size_all()<1) continue;
      total_vs_size += grp.vs.size();
      // Also define particle "first step" and "last step"
      auto& part = grp.part;
      auto const& first_pt = grp.first_pt;
      auto const& last_pt  = grp.last_pt;
      //std::cout<<first_pt.x<< " " << first_pt.y << " " << first_pt.z << std::endl;
      if(first_pt.t != larcv::kINVALID_DOUBLE)
        part.first_step(first_pt.x,first_pt.y,first_pt.z,first_pt.t);
      if(last_pt.t  != larcv::kINVALID_DOUBLE)
        part.last_step(last_pt.x,last_pt.y,last_pt.z,last_pt.t);
    }


    // loop over MCShower to assign parent/ancestor information
    LARCV_INFO() << "Processing MCShower array: " << mcs_v.size() << std::endl;
    for(auto const& mcs : mcs_v) {
      int track_id = mcs.TrackID();
      if(track_id >= ((int)(trackid2output.size()))) {
        LARCV_INFO() << "MCShower " << track_id << " PDG " << mcs.PdgCode()
                     << " not found in output group..." << std::endl;
        continue;
      }
      int output_id = trackid2output[track_id];
      //int group_id  = -1;
      int group_id  = output_id;
      if(output_id >= 0) {
        auto& grp = part_grp_v[track_id];
        assert(grp.part.group_id() == larcv::kINVALID_INSTANCEID);
        grp.part.group_id(group_id);
        /*
        if(grp.part.group_id() == larcv::kINVALID_INSTANCEID) {
          if(group_id < 0) {
            group_id = group_counter + 1;
            ++group_counter;
          }
          grp.part.group_id(group_id);
        }
        */
        // see if first step is not set yet
        if(grp.first_pt.t == larcv::kINVALID_DOUBLE)
          grp.part.first_step(mcs.DetProfile().X(),mcs.DetProfile().Y(),mcs.DetProfile().Z(),mcs.DetProfile().T());
        grp.part.parent_position(mcs.MotherStart().X(),
                                 mcs.MotherStart().Y(),
                                 mcs.MotherStart().Z(),
                                 mcs.MotherStart().T());
        grp.part.parent_creation_process(mcs.MotherProcess());
        grp.part.ancestor_position(mcs.AncestorStart().X(),
                                   mcs.AncestorStart().Y(),
                                   mcs.AncestorStart().Z(),
                                   mcs.AncestorStart().T());
        grp.part.ancestor_track_id(mcs.AncestorTrackID());
        grp.part.ancestor_pdg_code(mcs.AncestorPdgCode());
        grp.part.ancestor_creation_process(mcs.AncestorProcess());
      }

      for(auto const& child : mcs.DaughterTrackID()) {
        //if(child < trackid2output.size() && trackid2output[child] < 0)
        if(child < trackid2output.size() && trackid2output[child] >= 0) {
          trackid2output[child] = output_id;
          auto& grp = part_grp_v[child];
          assert(grp.part.group_id() == larcv::kINVALID_INSTANCEID);
          grp.part.group_id(group_id);
          /*
          if(grp.part.group_id() == larcv::kINVALID_INSTANCEID) {
            if(group_id < 0) {
              group_id = group_counter + 1;
              ++group_counter;
            }
            grp.part.group_id(group_id);
          }
          */
          grp.part.ancestor_position(mcs.AncestorStart().X(),
                                     mcs.AncestorStart().Y(),
                                     mcs.AncestorStart().Z(),
                                     mcs.AncestorStart().T());
          grp.part.ancestor_track_id(mcs.AncestorTrackID());
          grp.part.ancestor_pdg_code(mcs.AncestorPdgCode());
          grp.part.ancestor_creation_process(mcs.AncestorProcess());
        }
      }
    }


    // loop over MCTrack to assign parent/ancestor information
    auto const& mct_v = LArData<supera::LArMCTrack_t>();
    LARCV_INFO() << "Processing MCTrack array: " << mct_v.size() << std::endl;
    for(auto const& mct : mct_v) {
      int track_id  = mct.TrackID();
      int output_id = trackid2output[track_id];
      //int group_id  = -1;
      int group_id  = output_id;
      if(output_id >= 0) {
        auto& grp = part_grp_v[track_id];
        assert(grp.part.group_id() == larcv::kINVALID_INSTANCEID);
        /*
        if(group_id < 0) {
          group_id = group_counter + 1;
          ++group_counter;
        }
        */
        grp.part.group_id(group_id);
        if(grp.first_pt.t == larcv::kINVALID_DOUBLE && mct.size())
          grp.part.first_step(mct.front().X(),mct.front().Y(),mct.front().Z(),mct.front().T());
        if(grp.last_pt.t == larcv::kINVALID_DOUBLE && mct.size())
          grp.part.last_step(mct.back().X(),mct.back().Y(),mct.back().Z(),mct.back().T());
        grp.part.parent_position(mct.MotherStart().X(),
                                 mct.MotherStart().Y(),
                                 mct.MotherStart().Z(),
                                 mct.MotherStart().T());
        grp.part.parent_creation_process(mct.MotherProcess());
        grp.part.ancestor_position(mct.AncestorStart().X(),
                                   mct.AncestorStart().Y(),
                                   mct.AncestorStart().Z(),
                                   mct.AncestorStart().T());
        grp.part.ancestor_track_id(mct.AncestorTrackID());
        grp.part.ancestor_pdg_code(mct.AncestorPdgCode());
        grp.part.ancestor_creation_process(mct.AncestorProcess());
      }
      for(size_t output_index=0; output_index<output2trackid.size(); ++output_index) {
        int output_trackid = output2trackid[output_index];
        auto& grp = part_grp_v[output_trackid];
        if((int)(grp.part.parent_track_id()) != track_id) continue;
        //group ID should not be distinct for track children
        /*
        if(group_id < 0) {
          group_id = group_counter + 1;
          ++group_counter;
        }
        */
        grp.part.group_id(output_index);
        grp.part.ancestor_position(mct.AncestorStart().X(),
                                   mct.AncestorStart().Y(),
                                   mct.AncestorStart().Z(),
                                   mct.AncestorStart().T());
        grp.part.ancestor_track_id(mct.AncestorTrackID());
        grp.part.ancestor_pdg_code(mct.AncestorPdgCode());
        grp.part.ancestor_creation_process(mct.AncestorProcess());
      }
    }


    // Make sure the primary particle's parent and group id are set (they are themselves)
    for(auto& grp : part_grp_v) {
      auto& part = grp.part;
      if(part.track_id() != part.parent_track_id()) continue;
      part.group_id(part.id());
      part.parent_id(part.id());
    }


    // For shower orphans, we need to register the most base shower particle in the output (for group)
    LARCV_INFO() << "Searching the root (group) for kShapeShower particles w/ invalid group id... ("
                 << output2trackid.size() << " particles total)" << std::endl;
    for(size_t out_index=0; out_index<output2trackid.size(); ++out_index) {

      int trackid = output2trackid[out_index];
      auto& grp = part_grp_v[trackid];
      if(!grp.valid) continue;
      if(grp.part.group_id() != larcv::kINVALID_INSTANCEID) continue;
      if(grp.shape() != larcv::kShapeShower) continue;
      LARCV_INFO() <<" #### SHOWER ROOT SEARCH: Analyzing a particle index " << out_index
                   << " track id " << grp.part.track_id() << std::endl
                   << grp.part.dump() << std::endl;

      auto parent_trackid_v = this->ParentTrackIDs(trackid);
      int root_id = grp.part.id();
      int root_trackid = grp.part.track_id();
      bool stop = false;
      std::vector<size_t> intermediate_trackid_v;
      intermediate_trackid_v.push_back(trackid);
      for(auto const& parent_trackid : parent_trackid_v) {
        auto const& parent = part_grp_v[parent_trackid];
        switch(parent.shape()) {
          case larcv::kShapeShower:
          case larcv::kShapeMichel:
          case larcv::kShapeDelta:
            // group candidate: check if it is "valid" = exists in the output
            root_trackid = parent_trackid;
            root_id = trackid2output[root_trackid];
            if(root_id >= 0) {
              // found the valid group: stop the loop
              stop = true;
              // If not, root_id will be a new output index
            }else{
              root_id = output2trackid.size();
              // If this particle is invalid, this also needs the group id.
              // Add to intermediate_id_v list so we can set the group id for all of them
              intermediate_trackid_v.push_back(root_trackid);
            }
            stop = (stop || parent.shape() != larcv::kShapeShower);
            break;
          case larcv::kShapeTrack:
            stop = true;
            break;
          case larcv::kShapeUnknown:
            break;
          case larcv::kShapeLEScatter:
          case larcv::kShapeGhost:
            /*
            LARCV_CRITICAL() << "Unexpected type found while searching for kShapeShower orphans's root!" << std::endl;
            this->DumpHierarchy(trackid,part_grp_v);
            throw std::exception();
            */
            break;
        };
        if(stop) break;
      }
      if( root_id < ((int)(output2trackid.size())) && trackid2output[root_trackid] != (int)(root_id) ) {
        LARCV_CRITICAL() << "Logic error for the search of shower root particle for an orphan..." << std::endl
                         << "This particle id=" << out_index << " and track_id=" << trackid << std::endl
                         << "ROOT particle id=" << root_id  << " and track_id=" << root_trackid << std::endl;
        this->DumpHierarchy(trackid,part_grp_v);
        throw std::exception();
      }

      if(((int)(output2trackid.size())) <= root_id) {
        output2trackid.push_back(root_trackid);
        // Register the root parent to the output
        LARCV_INFO() << "Adding a new particle to the output to define a group..." << std::endl
                     << "ROOT particle id=" << root_id  << " and track_id=" << root_trackid << std::endl
                     << part_grp_v[root_trackid].part.dump() << std::endl;
      }
      assert((size_t)(root_id) < output2trackid.size());

      auto& root = part_grp_v[root_trackid];
      //root.valid = true;
      assert(root.valid);
      root.part.id(root_id);
      root.part.group_id(root_id);
      trackid2output[root_trackid] = root_id;
      for(auto const& child_id : root.part.children_id()) {
        auto& child = part_grp_v[output2trackid[child_id]];
        if(!child.valid) continue;
        assert(child.part.group_id() == larcv::kINVALID_INSTANCEID || child.part.group_id() == root_id);
        child.part.group_id(root_id);
      }
      // Set the group ID for THIS + intermediate particles
      for(auto const& child_trackid : intermediate_trackid_v) {
        auto& child = part_grp_v[child_trackid];
        if(!child.valid) continue;
        assert(child.part.group_id() == larcv::kINVALID_INSTANCEID || child.part.group_id() == root_id);
        child.part.group_id(root_id);
      }

      LARCV_INFO() << "... after update ... " << std::endl
                   << part_grp_v[trackid].part.dump() << std::endl;
    }


    // For LEScatter orphans, we need to register the immediate valid (=to be stored) particle
    LARCV_INFO() << "Searching the root (group) for kShapeLEScatter particles w/ invalid group id... ("
                 << output2trackid.size() << " particles total)" << std::endl;
    for(size_t out_index=0; out_index<output2trackid.size(); ++out_index) {
      int trackid = output2trackid[out_index];
      auto& grp = part_grp_v[trackid];
      if(grp.part.group_id() != larcv::kINVALID_INSTANCEID) continue;
      if(grp.shape() != larcv::kShapeLEScatter) continue;
      LARCV_INFO() <<" #### LEScatter ROOT SEARCH #### " << std::endl
                   <<" Analyzing a particle index " << out_index << " id " << grp.part.id() << std::endl
                   << grp.part.dump() << std::endl;

      auto parent_trackid_v = this->ParentTrackIDs(trackid);
      size_t group_id = larcv::kINVALID_INSTANCEID;
      bool stop = false;
      for(auto const& parent_trackid : parent_trackid_v) {
        auto const& parent = part_grp_v[parent_trackid];
        switch(parent.shape()) {
          case larcv::kShapeShower:
          case larcv::kShapeMichel:
          case larcv::kShapeDelta:
          case larcv::kShapeTrack:
          case larcv::kShapeLEScatter:
            // group candidate: check if it is "valid" = exists in the output
            if(parent.valid) {
              group_id = trackid2output[parent_trackid];
              // found the valid group: stop the loop
              stop = true;
            }
            break;
          case larcv::kShapeUnknown:
          case larcv::kShapeGhost:
            LARCV_CRITICAL() << "Unexpected type found while searching for kShapeLEScatter orphans's root!" << std::endl;
            throw std::exception();
            break;
        };
        if(stop) break;
      }
      if(group_id == larcv::kINVALID_INSTANCEID) {
        LARCV_INFO() << "Ignoring kShapeLEScatter particle as its root particle (for group id) is not to be stored..." << std::endl
                     << grp.part.dump() << std::endl;
        continue;
      }
      LARCV_INFO() << "Assigning a group ID " << group_id << " to kShapeLEScatter orphan" << std::endl
                   << "  Track ID " << grp.part.track_id() << " PDG " << grp.part.pdg_code()
                   << " " << grp.part.creation_process() << std::endl;
      grp.part.group_id(group_id);
      trackid2output[trackid] = group_id;
    }


    // for shower particles with invalid parent ID, attempt a search
    LARCV_INFO() << "Searching parents for shower particles w/ invalid parent id... ("
                 << output2trackid.size() << " particles total)" << std::endl;
    for(size_t out_index=0; out_index<output2trackid.size(); ++out_index) {
      int trackid = output2trackid[out_index];
      auto& grp = part_grp_v[trackid];
      if(!grp.valid) continue;
      if(grp.part.parent_id() != larcv::kINVALID_INSTANCEID) continue;
      if(grp.shape() != larcv::kShapeShower)
        continue;
      LARCV_INFO() << "Analyzing particle id " << out_index << " trackid " << trackid << std::endl
                   << grp.part.dump() << std::endl;
      int parent_partid = -1;
      unsigned int parent_trackid;
      auto parent_trackid_v = this->ParentTrackIDs(trackid);
      for(size_t idx=0; idx<parent_trackid_v.size(); ++idx) {
        parent_trackid = parent_trackid_v[idx];
        if(trackid2output[parent_trackid] < 0 || !part_grp_v[parent_trackid].valid)
          continue;
        auto const& parent = part_grp_v[parent_trackid].part;
        // shower parent can be either shower, michel, or delta
        if(parent.shape() == larcv::kShapeMichel ||
           parent.shape() == larcv::kShapeDelta  ||
           parent.shape() == larcv::kShapeShower)
          parent_partid  = parent.id();
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
      if(parent_partid >=0) {
        // assert the group is same
        auto& parent = part_grp_v[output2trackid[parent_partid]];
        if(grp.part.group_id() == larcv::kINVALID_INSTANCEID) {
          grp.part.group_id(parent.part.group_id());
          for(auto const& child_id : grp.part.children_id()) {
            auto& child = part_grp_v[output2trackid[child_id]];
            child.part.group_id(parent.part.group_id());
          }
        }else{
          assert(grp.part.group_id() == part_grp_v[output2trackid[parent_partid]].group_id());
        }
        grp.part.parent_id(parent_partid);
        part_grp_v[parent_trackid].part.children_id(grp.part.id());
        LARCV_INFO() << "PartID " << grp.part.id() << " (output index " << out_index << ") assigning parent " << parent_partid << std::endl;
      }else{
        grp.part.parent_id(grp.part.id());
        if(grp.part.group_id() == larcv::kINVALID_INSTANCEID)
          grp.part.group_id(grp.part.id());
        for(auto const& child_id : grp.part.children_id()) {
          auto& child = part_grp_v[output2trackid[child_id]];
          child.part.group_id(grp.part.id());
        }
        LARCV_INFO() << "PartID " << grp.part.id() << " (output index " << out_index << ") assigning itself as a parent..." << std::endl;
      }
    }


    // Now sort out all parent IDs where it's simply not assigned
    // (it's ok to have invalid parent id if parent track id is not stored)
    LARCV_INFO() << "Check all output particle's parent id... ("
                 << output2trackid.size() << " particles total)" << std::endl;
    for(size_t output_index=0; output_index<output2trackid.size(); ++output_index) {
      auto& grp = part_grp_v[output2trackid[output_index]];
      auto parent_trackid = grp.part.parent_track_id();
      auto parent_id = grp.part.parent_id();
      auto& parent = part_grp_v[parent_trackid].part;
      // if parent_id is invalid, try if parent_trackid can help out
      if(parent_id == larcv::kINVALID_INSTANCEID &&
         trackid2output[parent_trackid] >= 0) {
        parent_id = trackid2output[parent_trackid];
        grp.part.parent_id(parent_id);
      }
      if(parent_id == larcv::kINVALID_INSTANCEID) continue;
      // if parent id is set, make sure this particle is in the children
      auto children = parent.children_id();
      bool add=true;
      for(auto const& child : children) {
        if(child != grp.part.id()) continue;
        add = false; break;
      }
      if(add) {children.push_back(grp.part.id()); parent.children_id(children);}
    }

    // Now loop over otuput particle list and check if any remaining group id needs to be assigned
    // Use its parent to group...
    LARCV_INFO() <<  "Check all output particles for their group id... ("
                 << output2trackid.size() << " particles total)" << std::endl;
    for(size_t output_index=0; output_index<output2trackid.size(); ++output_index) {
      auto& grp = part_grp_v[output2trackid[output_index]];
      if(grp.part.group_id() != larcv::kINVALID_INSTANCEID) continue;
      auto shape = grp.shape();
      auto parent_shape  = larcv::kShapeUnknown;
      auto parent_partid = grp.part.parent_id();
      //auto parent_groupid = larcv::kINVALID_INSTANCEID;
      // If delta, its own grouping

      switch(shape) {
        case larcv::kShapeLEScatter:
          // if LEScatter, we handle later (next loop)
          break;
        case larcv::kShapeDelta:
        case larcv::kShapeMichel:
        case larcv::kShapeTrack:
          // If delta, Michel, or track, it's own group
          grp.part.group_id(output_index);
          for(auto const& child_index : grp.part.children_id()) {
            part_grp_v[output2trackid[child_index]].part.group_id(output_index);
          }
          break;

        case larcv::kShapeShower:
          // If shower && no parent, consider it as a primary = assign group id for all children
          if(parent_partid == kINVALID_INSTANCEID) {
            grp.part.group_id(output_index);
            for(auto const& child_index : grp.part.children_id())
              part_grp_v[output2trackid[child_index]].part.group_id(output_index);
            continue;
          }
          parent_shape = part_grp_v[output2trackid[parent_partid]].shape();
          switch(parent_shape) {
            case larcv::kShapeMichel:
            case larcv::kShapeDelta:
              grp.part.group_id(parent_partid);
              for(auto const& child_index : grp.part.children_id()) {
                part_grp_v[output2trackid[child_index]].part.group_id(parent_partid);
              }
              break;
            case larcv::kShapeTrack:
              grp.part.group_id(output_index);
              for(auto const& child_index : grp.part.children_id()) {
                part_grp_v[output2trackid[child_index]].part.group_id(output_index);
              }
              break;
            case larcv::kShapeShower:
              LARCV_CRITICAL() << "Unexpected case: a shower has no group id while being a child of another shower..." << std::endl;
              this->DumpHierarchy(grp.part.track_id(),part_grp_v);
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
            case larcv::kShapeLEScatter:
              LARCV_ERROR() << "Logic error: shower parent shape cannot be LEScatter!" <<std::endl;
              throw std::exception();
            default:
              LARCV_ERROR() << "Unexpected larcv::ShapeType_t encountered at " << __LINE__ <<std::endl;
              throw std::exception();
          };
          break;
        case larcv::kShapeGhost:
        case larcv::kShapeUnknown:
          LARCV_ERROR() << "Unexpected larcv::ShapeType_t encountered at " << __LINE__ <<std::endl;
          throw std::exception();
      }
    };

    /*
    // GroupID / parentage search
    for(size_t trackid=0; trackid<part_grp_v.size(); ++trackid) {
      auto grp& part_grp_v[trackid];
      if(trackid2output[trackid]<0)

      if(grp.shape()==larcv::kShapeTrack)

    }
    */


    // Next handle LEScatter group id if not assigned yet
    LARCV_INFO() <<  "Check LEScatter particle's group id ... ("
                 << output2trackid.size() << " particles total)" << std::endl;
    for(size_t output_index=0; output_index<output2trackid.size(); ++output_index) {
      auto& grp = part_grp_v[output2trackid[output_index]];
      if(grp.part.group_id() != larcv::kINVALID_INSTANCEID) continue;
      if(grp.shape() == larcv::kShapeLEScatter) {
        // assign parent's group, otherwise leave as is = kINVALID_INSTANCEID
        auto parent_partid = grp.part.parent_id();
        if(parent_partid == larcv::kINVALID_INSTANCEID) continue;
        grp.part.group_id( part_grp_v[output2trackid[parent_partid]].part.group_id() );
      }
    }

    // Next loop over to find any particle for which first_step is not defined
    LARCV_INFO() <<  "Check any particle's first step ... ("
                 << output2trackid.size() << " particles total)" << std::endl;
    for(size_t output_index=0; output_index<output2trackid.size(); ++output_index) {
      auto& grp = part_grp_v[output2trackid[output_index]];
      auto const& fs = grp.part.first_step();
      if(fs.x()!=0. || fs.y()!=0. || fs.z()!=0. || fs.t()!=0.) continue;
      auto const vtx = grp.part.position().as_point3d();
      double min_dist = std::fabs(larcv::kINVALID_DOUBLE);
      larcv::Point3D min_pt;
      for(auto const& vox : grp.vs.as_vector()) {
        auto const pt = meta3d.position(vox.id());
        double dist = pt.squared_distance(vtx);
        if(dist > min_dist) continue;
        min_dist = dist;
        min_pt = pt;
      }
      if(min_dist > (sqrt(3.) + 1.e-3)) grp.part.first_step(larcv::kINVALID_DOUBLE,larcv::kINVALID_DOUBLE,larcv::kINVALID_DOUBLE,larcv::kINVALID_DOUBLE);
      else grp.part.first_step(min_pt.x, min_pt.y, min_pt.z, grp.part.position().t());

    }

    LARCV_INFO() <<  "Start storing " << output2trackid.size() << " particles ...y" << std::endl;

    // now loop over to create VoxelSet for compton/photoelectron
    std::vector<larcv::Particle> part_v; part_v.resize(output2trackid.size());
    auto event_cluster    = (EventClusterVoxel3D*)(mgr.get_data("cluster3d",_output_label));
    auto event_cluster_he = (EventClusterVoxel3D*)(mgr.get_data("cluster3d",_output_label + "_highE"));
    auto event_cluster_le = (EventClusterVoxel3D*)(mgr.get_data("cluster3d",_output_label + "_lowE"));
    auto event_leftover   = (EventSparseTensor3D*)(mgr.get_data("sparse3d",_output_label + "_leftover"));
    // set meta for all
    event_cluster->meta(meta3d);
    event_cluster_he->meta(meta3d);
    event_cluster_le->meta(meta3d);
    event_leftover->meta(meta3d);

    std::vector<larcv::ClusterPixel2D> vsa2d_v;          vsa2d_v.resize(_valid_nplanes);
    std::vector<larcv::ClusterPixel2D> vsa2d_he_v;       vsa2d_he_v.resize(_valid_nplanes);
    std::vector<larcv::ClusterPixel2D> vsa2d_le_v;       vsa2d_le_v.resize(_valid_nplanes);
    for(size_t plane_idx=0; plane_idx<_valid_nplanes; ++plane_idx) {
      vsa2d_v[plane_idx].resize(output2trackid.size());
      vsa2d_he_v[plane_idx].resize(output2trackid.size());
      vsa2d_le_v[plane_idx].resize(output2trackid.size());
      vsa2d_v[plane_idx].meta(_meta2d_v[plane_idx]);
      vsa2d_he_v[plane_idx].meta(_meta2d_v[plane_idx]);
      vsa2d_le_v[plane_idx].meta(_meta2d_v[plane_idx]);
    }


    // Create cluster index tensor to help back-track semantic source particle
    auto event_cindex = (EventSparseTensor3D*)(mgr.get_data("sparse3d",_output_label + "_index"));
    event_cindex->meta(meta3d);
    larcv::VoxelSet cid_vs; cid_vs.reserve(total_vs_size);

    event_cluster->resize(output2trackid.size());
    event_cluster_he->resize(output2trackid.size());
    event_cluster_le->resize(output2trackid.size());
    for(size_t index=0; index<output2trackid.size(); ++index) {
      int trackid = output2trackid[index];
      auto& grp   = part_grp_v[trackid];
      // set semantic type
      larcv::ShapeType_t semantic = grp.shape();
      if(semantic == larcv::kShapeUnknown) {
        LARCV_CRITICAL() << "Unexpected type while assigning semantic class: " << grp.type << std::endl;
        auto const& part = grp.part;
        LARCV_CRITICAL() << "Particle ID " << part.id() << " Type " << grp.type << " Valid " << grp.valid
                         << " Track ID " << part.track_id() << " PDG " << part.pdg_code()
                         << " " << part.creation_process() << " ... " << part.energy_init() << " MeV => " << part.energy_deposit() << " MeV "
                         << grp.trackid_v.size() << " children " << grp.vs.size() << " voxels " << grp.vs.sum() << " MeV" << std::endl;
        LARCV_CRITICAL() << "  Parent " << part.parent_track_id() << " PDG " << part.parent_pdg_code()
                         << " " << part.parent_creation_process() << " Ancestor " << part.ancestor_track_id()
                         << " PDG " << part.ancestor_pdg_code() << " " << part.ancestor_creation_process() << std::endl;


        throw std::exception();
      }
      if(semantic == larcv::kShapeLEScatter && mcs_trackid_s.find(trackid) == mcs_trackid_s.end()) {
        LARCV_CRITICAL() << "Unexpected particle to be stored with a shape kShapeLEScatter!" << std::endl;
        this->DumpHierarchy(grp.part.track_id(),part_grp_v);
        throw std::exception();
      }
      // Now, we will re-classify some of non LEScatter showers based on pixel count
      // (BUG FIX: use pca or something better)
      if(grp.vs.size() < _compton_size) {
        LARCV_INFO() << "Particle ID " << grp.part.id() << " PDG " << grp.part.pdg_code() << " " << grp.part.creation_process() << std::endl
                     << "  ... type switching " << grp.part.shape() << " => " << larcv::kShapeLEScatter
                     << " (voxel count " << grp.vs.size() << " < " << _compton_size << ")" << std::endl;
        semantic = larcv::kShapeLEScatter;
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
      if(semantic != kShapeLEScatter)
        event_cluster_he->writeable_voxel_set(index) = grp.vs;
      else {
        event_cluster_le->writeable_voxel_set(index) = grp.vs;
        for(auto const& vox : grp.vs.as_vector())
          cid_vs.emplace(vox.id(),index,false);
      }
      //grp.vs.clear_data();
      // fill 2d cluster
      for(size_t plane_idx=0; plane_idx<_valid_nplanes; ++plane_idx) {
        auto& vsa2d = vsa2d_v[plane_idx];
        vsa2d.writeable_voxel_set(index) = grp.vs2d_v[plane_idx];
        if(semantic != kShapeLEScatter) {
          auto& vsa2d_he = vsa2d_he_v[plane_idx];
          vsa2d_he.writeable_voxel_set(index) = grp.vs2d_v[plane_idx];
        }else{
          auto& vsa2d_le = vsa2d_le_v[plane_idx];
          vsa2d_le.writeable_voxel_set(index) = grp.vs2d_v[plane_idx];
        }
        grp.vs2d_v[plane_idx].clear_data();
      }
      grp.valid=false;
    }


    // Loop to store output cluster/semantic: low energy depositions
    //for(auto& grp : part_grp_v) {
    for(size_t grp_idx=0; grp_idx<part_grp_v.size(); ++grp_idx) {
      auto& grp = part_grp_v[grp_idx];
      if(!grp.valid) continue;
      if(grp.size_all()<1) continue;
      auto semantic = grp.shape();
      if(semantic != larcv::kShapeLEScatter) {
        LARCV_CRITICAL() << "Unexpected, valid, >1 pixel count particle type "
                         << semantic << " pixel count " << grp.vs.size()
                         << " (not kShapeLEScatter) at line " << __LINE__
                         << std::endl;
        std::cout<<grp.part.dump()<<std::endl;
        throw std::exception();
      }
      int trackid = grp.part.parent_track_id();
      int output_index = -1;
      if(trackid < ((int)(trackid2output.size()))) output_index = trackid2output[trackid];
      if(output_index<0) {
        // search the first direct, valid parent
        // START_HERE 2020-04-23 ... use ParentTrackIDs(), make sure output_index is not int max
        while(1) {
          if(trackid >= ((int)(trackid2index.size())) || trackid2index[trackid] < 0) break;
          trackid = larmcp_v[trackid2index[trackid]].Mother();
          if(trackid < ((int)(trackid2output.size()))) {
            output_index = trackid2output[trackid];
            break;
          }
        }
      }
      //std::cout<<"Inspecting Track ID " << grp.part.track_id() << " found output ID " << output_index << std::endl;
      if(output_index<0) continue;

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
      auto& vs_le = event_cluster_le->writeable_voxel_set(output_index);
      auto& vs    = event_cluster->writeable_voxel_set(output_index);
      for(auto const& vox : grp.vs.as_vector()) {
        vs_le.emplace(vox.id(),vox.value(),true);
        vs.emplace(vox.id(),vox.value(),true);
        cid_vs.emplace(vox.id(),output_index,false);
      }
      grp.vs.clear_data();
      // fill 2d cluster
      for(size_t plane_idx=0; plane_idx<_valid_nplanes; ++plane_idx) {
        auto& vs2d = vsa2d_v[plane_idx].writeable_voxel_set(output_index);
        auto& vs2d_le = vsa2d_le_v[plane_idx].writeable_voxel_set(output_index);
        for(auto const& vox : grp.vs2d_v[plane_idx].as_vector()) {
          vs2d.emplace(vox.id(),vox.value(),true);
          vs2d_le.emplace(vox.id(),vox.value(),true);
        }
        grp.vs2d_v[plane_idx].clear_data();
      }
    }

    // create particle ID vs ... overlapped voxel gets higher id number
    auto const& main_vs = event_cluster_he->as_vector();
    auto const& lowe_vs = event_cluster_le->as_vector();

    // Count output voxel count and x-check
    size_t output_vs_size = 0;
    for(auto const& vs : main_vs) output_vs_size += vs.size();
    for(auto const& vs : lowe_vs) output_vs_size += vs.size();
    LARCV_INFO() << "Voxel count x-check: output = " << output_vs_size << " ... total = " << total_vs_size << std::endl;

    LARCV_INFO()<<"Combined reminders..."<<std::endl;
    larcv::VoxelSet leftover_vs; leftover_vs.reserve(total_vs_size - output_vs_size);
    std::vector<larcv::VoxelSet> leftover2d_vs(_valid_nplanes);
    for(auto& vs : leftover2d_vs) vs.reserve(total_vs_size - output_vs_size);

    if(total_vs_size > output_vs_size) {
      int ctr= 0;
      for(auto& grp : part_grp_v) {
        if(grp.size_all()<1) continue;
        for(auto const& vox : grp.vs.as_vector()) leftover_vs.emplace(vox.id(),vox.value(),true);
        for(size_t plane_idx=0; plane_idx<_valid_nplanes; ++plane_idx) {
          for(auto const& vox : grp.vs2d_v[plane_idx].as_vector()) {
            leftover2d_vs[plane_idx].emplace(vox.id(),vox.value(),true);
          }
        }
        ctr++;
        auto const& part = grp.part;
        LARCV_INFO() << "Particle ID " << part.id() << " Type " << grp.type << " Valid " << grp.valid << " Track ID " << part.track_id() << " PDG " << part.pdg_code()
                     << " " << part.creation_process() << " ... " << part.energy_init() << " MeV => " << part.energy_deposit() << " MeV "
                     << grp.trackid_v.size() << " children " << grp.vs.size() << " voxels " << grp.vs.sum() << " MeV" << std::endl;
        LARCV_INFO() << "  Parent " << part.parent_track_id() << " PDG " << part.parent_pdg_code() << " " << part.parent_creation_process()
                     << " Ancestor " << part.ancestor_track_id() << " PDG " << part.ancestor_pdg_code() << " " << part.ancestor_creation_process() << std::endl;
        LARCV_INFO() << "  Group ID: " << part.group_id() << std::endl;
        std::stringstream ss1,ss2;

        ss1 << "  Children particle IDs: " << std::flush;
        for(auto const& child : part.children_id()) ss1 << child << " " << std::flush;
        ss1 << std::endl;
        LARCV_INFO() << ss1.str();

        ss2 << "  Children track IDs: " << std::flush;
        for(auto const& child : grp.trackid_v) ss2 << child << " " << std::flush;
        ss2 << std::endl;
        LARCV_INFO() << ss2.str();
        LARCV_INFO() << "Above was supposed to be merged..." << std::endl;

      }
      LARCV_INFO() << "... " << ctr << " particles" << std::endl;
    }
    event_leftover->emplace(std::move(leftover_vs),meta3d);

    // Loop over to find any "stil valid" supera::kIonization supera::kConversion supera::kComptonHE
    LARCV_INFO() << "Particle list" << std::endl;
    for(size_t index = 0; index < part_v.size(); ++index) {
      int trackid = output2trackid[index];
      auto const& grp  = part_grp_v[trackid];
      auto const& part = part_v[index];
      auto const& vs0  = main_vs[index];
      auto const& vs1  = lowe_vs[index];
      LARCV_INFO() << "Particle ID " << part.id() << " Track ID " << part.track_id() << " PDG " << part.pdg_code()
                   << " " << part.creation_process() << " ... " << part.energy_init() << " MeV => " << part.energy_deposit() << " MeV "
                   << grp.trackid_v.size() << " children " << vs0.size() << " (" << vs1.size() << ") voxels" << std::endl;
      LARCV_INFO() << "  Parent TrackID " << part.parent_track_id() << " PartID " << part.parent_id()
                   << " PDG " << part.parent_pdg_code() << " " << part.parent_creation_process()
                   << " Ancestor TrackID " << part.ancestor_track_id()
                   << " PDG " << part.ancestor_pdg_code() << " " << part.ancestor_creation_process() << std::endl;
      LARCV_INFO() << "  Group ID: " << part.group_id() << std::endl;
      std::stringstream ss1,ss2;

      ss1 << "  Children particle IDs: " << std::flush;
      for(auto const& child : part.children_id()) ss1 << child << " " << std::flush;
      ss1 << std::endl;
      LARCV_INFO() << ss1.str();

      ss2 << "  Children track IDs: " << std::flush;
      for(auto const& child : grp.trackid_v) ss2 << child << " " << std::flush;
      ss2 << std::endl;
      LARCV_INFO() << ss2.str();

      LARCV_INFO() << "  Start: " << part.first_step().x() << " " << part.first_step().y() << " " << part.first_step().z() << std::endl;
    }
    LARCV_INFO() << "... " << part_v.size() << " particles" << std::endl;


    // Check the validity of particle id, group id
    if(_check_particle_validity) {
      for(size_t part_index=0; part_index < part_v.size(); ++part_index) {
        auto const& part = part_v[part_index];
        if(part.id() != part_index) {
          LARCV_CRITICAL() << "Particle index " << part_index << " holds particle w/ an ID " << part.id() << std::endl
                           << part.dump() << std::endl;
          throw std::exception();
        }
        if(part.parent_id() != larcv::kINVALID_INSTANCEID && part.parent_id() >= part_v.size()) {
          LARCV_CRITICAL() << "Particle index " << part_index
                           << " holds particle w/ an invalid parent ID " << part.parent_id() << std::endl
                           << part.dump() << std::endl
                           << part_grp_v[part.parent_track_id()].part.dump() << std::endl;
          throw std::exception();
        }
        if(part.group_id() == larcv::kINVALID_INSTANCEID || part.group_id() >= part_v.size()) {
          LARCV_CRITICAL() << "Particle index " << part_index
                           << " holds particle w/ an invalid group ID " << part.group_id() << std::endl
                           << part.dump() << std::endl;
          throw std::exception();
        }
        if(part.parent_id() == part.id() || part.parent_id() == larcv::kINVALID_INSTANCEID) continue;
        bool found=false;
        for(auto const& child : part_v[part.parent_id()].children_id()) {
          if(child == part_index) { found = true; break; }
        }
        if(!found) {
          LARCV_WARNING() << "Particle index " << part_index
                          << " not in the list of the parent's children (fixing)" << std::endl;
          auto children = part_v[part.parent_id()].children_id();
          children.push_back(part_index);
          part_v[part.parent_id()].children_id(children);
        }
      }
    }


    // create semantic output in 3d
    auto event_segment = (EventSparseTensor3D*)(mgr.get_data("sparse3d",_output_label + "_semantics"));
    event_segment->meta(meta3d);
    larcv::VoxelSet semantic_vs; semantic_vs.reserve(total_vs_size);
    // create semantic output in 2d
    std::vector<larcv::VoxelSet> semantic2d_vs_v;
    semantic2d_vs_v.resize(_valid_nplanes);
    for(auto& vs : semantic2d_vs_v) vs.reserve(total_vs_size);

    /*
    for(size_t index=0; index<main_vs.size(); ++index) {
      auto const& vs = main_vs[index];
      for(auto const& vox : vs.as_vector()) cid_vs.emplace(vox.id(),(float)(index),false);
    }
    */

    // Comptons in 3d
    for(auto const& vs : event_cluster_le->as_vector()) {
      for(auto const& vox : vs.as_vector()) semantic_vs.emplace(vox.id(),(float)(larcv::kShapeLEScatter),false);
    }
    // Comptons in 2d
    for(size_t plane_idx=0; plane_idx<_valid_nplanes; ++plane_idx) {
      auto& semantic2d = semantic2d_vs_v[plane_idx];
      auto& vsa2d_le   = vsa2d_le_v[plane_idx];
      for(auto const& vs : vsa2d_le.as_vector()) {
        for(auto const& vox : vs.as_vector()) semantic2d.emplace(vox.id(),(float)(larcv::kShapeLEScatter),false);
      }
    }

    // Loop over "high energy" depositions, set semantic labels
    for(size_t index=0; index<output2trackid.size(); ++index) {
      auto const& vs  = event_cluster_he->as_vector()[index];
      size_t semantic = (size_t)(part_v[index].shape());
      for(auto const& vox : vs.as_vector()) {
        auto const& prev = semantic_vs.find(vox.id());
        if(prev.id() == larcv::kINVALID_VOXELID) {
          semantic_vs.emplace(vox.id(),semantic,false);
          cid_vs.emplace(vox.id(),index,false);
        }
        else {
          size_t prioritized_semantic = this->SemanticPriority(prev.value(),semantic);
          if(prioritized_semantic != prev.value()) {
            semantic_vs.emplace(vox.id(),semantic,false);
            cid_vs.emplace(vox.id(),index,false);
          }
        }
      }
      for(size_t plane_id=0; plane_id < _valid_nplanes; ++plane_id) {
        auto& semantic2d_vs = semantic2d_vs_v[plane_id];
        auto const& vsa2d_he = vsa2d_he_v[plane_id];
        auto const& vs2d = vsa2d_he.as_vector()[index];
        for(auto const& vox : vs2d.as_vector()) {
          auto const& prev = semantic2d_vs.find(vox.id());
          if(prev.id() == larcv::kINVALID_VOXELID)
            semantic2d_vs.emplace(vox.id(),semantic,false);
          else
            semantic2d_vs.emplace(vox.id(),this->SemanticPriority(((size_t)(prev.value())),semantic),false);
        }
      }
    }
    // store
    assert(semantic_vs.size() == cid_vs.size());
    event_segment->emplace(std::move(semantic_vs),meta3d);
    event_cindex->emplace(std::move(cid_vs),meta3d);

    for(size_t cryo_id=0; cryo_id<_scan.size(); ++cryo_id) {
      auto const& tpcs = _scan[cryo_id];
      for(size_t tpc_id=0; tpc_id<tpcs.size(); ++tpc_id) {
        auto const& planes = tpcs[tpc_id];
        for(size_t plane_id=0; plane_id<planes.size(); ++plane_id) {
          auto const& idx = planes.at(plane_id);
          if(idx < 0) continue;
          if(idx >= (int)(_meta2d_v.size())) {
            LARCV_CRITICAL() <<idx << "Unexpected: " << _meta2d_v.size() << " " << semantic2d_vs_v.size() << std::endl;
            throw larbys();
          }
          auto meta2d = _meta2d_v[idx];
          std::string suffix = "_" + std::to_string(cryo_id) + "_" + std::to_string(tpc_id) + "_" + std::to_string(plane_id);
          std::string output_producer = _output_label + "_semantics2d" + suffix;
          auto semantic2d_output = (larcv::EventSparseTensor2D*)(mgr.get_data("sparse2d",output_producer));
          semantic2d_output->emplace(std::move(semantic2d_vs_v[idx]),std::move(meta2d));
          output_producer = _output_label + suffix;
          auto& cluster2d_output  = mgr.get_data<larcv::EventClusterPixel2D>(output_producer);
          cluster2d_output.emplace(std::move(vsa2d_v[idx]));

          output_producer = _output_label + "_he" + suffix;
          auto& cluster2d_he_output  = mgr.get_data<larcv::EventClusterPixel2D>(output_producer);
          cluster2d_he_output.emplace(std::move(vsa2d_he_v[idx]));

          output_producer = _output_label + "_le" + suffix;
          auto& cluster2d_le_output  = mgr.get_data<larcv::EventClusterPixel2D>(output_producer);
          cluster2d_le_output.emplace(std::move(vsa2d_le_v[idx]));

          output_producer = _output_label + "_leftover" + suffix;
          auto& tensor2d_leftover  = mgr.get_data<larcv::EventSparseTensor2D>(output_producer);
          meta2d = _meta2d_v[idx];
          tensor2d_leftover.emplace(std::move(leftover2d_vs[idx]),std::move(meta2d));
          //std::cout<<cryo_id<<" "<<tpc_id<<" "<<plane_id<<" ... " << semantic2d_output.as_vector().size() << " " << semantic2d_output.as_vector().front().meta().id() <<std::endl;
        }
      }
    }

    // Store output
    auto event_mcp = (EventParticle*)(mgr.get_data("particle",_output_label));
    event_mcp->emplace(std::move(part_v));

    return true;
  }

  // =============================================================================================
  //  private methods
  // =============================================================================================

  // ------------------------------------------------------
  void SuperaMCParticleCluster::AnalyzeSimEnergyDeposit(const larcv::Voxel3DMeta& meta,
                                                        std::vector<supera::ParticleGroup>& part_grp_v,
                                                        larcv::IOManager& mgr)
  {
    const auto & ev = this->GetEvent();
    ;

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
  }

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
      int track_id = mcpart.track_id();
      if (mcpart.ancestor_track_id() < ((int) (trackid2index.size())))
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
        int parent_index = grp.part.parent_track_id();
        int parent_index_before = grp.part.track_id();
        while (true)
        {
          //std::cout<< "Inspecting: " << grp.part.track_id() << " => " << parent_index << std::endl;
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
            int ancestor_index = parent.part.parent_track_id();
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
  void SuperaMCParticleCluster::MergeShowerTouchingLEScatter(const larcv::Voxel3DMeta& meta,
                                                             std::vector<supera::ParticleGroup>& part_grp_v)
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
        auto const &parents = this->ParentTrackIDs(grp.part.track_id());
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
  SuperaMCParticleCluster::ParentTrackIDs(size_t trackid,
                                          const std::vector<larcv::Particle> & particles) const
  {
    auto const &trackid2index = _mc_part_list.TrackIdToIndex();
    std::vector<unsigned int> result;

    if (trackid >= trackid2index.size() || trackid2index[trackid] < 0)
      return result;

    int parent_trackid = particles[trackid2index[trackid]].parent_track_id();
    std::set<int> accessed;
    while (parent_trackid > 0 &&
           (size_t) (parent_trackid) < trackid2index.size() &&
           trackid2index[parent_trackid] >= 0)
    {
      if (accessed.find(parent_trackid) != accessed.end())
      {
        LARCV_CRITICAL() << "LOOP-LOGIC-ERROR for ParentTrackIDs for track id " << trackid << std::endl;
        for (size_t parent_idx = 0; parent_idx < result.size(); ++parent_idx)
        {
          auto const &parent_trackid = result[parent_idx];
          auto const &mcp = particles[trackid2index.at(parent_trackid)];
          LARCV_CRITICAL() << "Parent " << parent_idx
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

}
