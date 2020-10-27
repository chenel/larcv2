/**
 * \file SuperMCParticleCluster.h
 *
 * \brief Algorithm to do clustering of MC particles in groups like EM showers
 *
 * @author J. Wolcott <jwolcott@fnal.gov>
 */

#ifndef LARCV2_SUPERAMCPARTICLECLUSTER_H
#define LARCV2_SUPERAMCPARTICLECLUSTER_H

#include <string>

#include "larcv/core/Base/PSet.h"

#include "MCParticleList.h"
#include "SuperaBase.h"
#include "SuperaMCParticleClusterData.h"

class TG4Trajectory;

namespace larcv
{
  class IOManager;

  class Particle;

  class Voxel3DMeta;
  class VoxelSet;

  /**
   * \class SuperaMCParticleCluster
   *  Clusters together true MC particles that belong together (e.g., EM showers)
   */
  class SuperaMCParticleCluster : public SuperaBase
  {
    public:
      /// Default constructor
      SuperaMCParticleCluster(const std::string name = "SuperaMCParticleCluster");

      /// Default destructor
      ~SuperaMCParticleCluster() override = default;

      void configure(const PSet&) override;
      void initialize() override;
      bool process(IOManager& mgr) override;

    private:
      // todo: this is duplicated from SuperaG4HitSegment
//      larcv::Particle MakeParticle(const TG4Trajectory&) const;

//      bool IsTouching(const Voxel3DMeta& meta, const VoxelSet& vs1, const VoxelSet& vs2) const;
//
//      std::vector<supera::ParticleGroup> CreateParticleGroups();
//
//      void AnalyzeSimEnergyDeposit(const larcv::Voxel3DMeta& meta,
//                                   std::vector<supera::ParticleGroup>& part_grp_v,
//                                   larcv::IOManager& mgr);
//
//      void AnalyzeSimChannel(const larcv::Voxel3DMeta& meta,
//                             std::vector<supera::ParticleGroup>& part_grp_v,
//                             larcv::IOManager& mgr);
//      void AnalyzeFirstLastStep(const larcv::Voxel3DMeta& meta,
//                                std::vector<supera::ParticleGroup>& part_grp_v);
//      void MergeShowerTouchingLEScatter(const larcv::Voxel3DMeta& meta,
//                                        std::vector<supera::ParticleGroup>& part_grp_v);
//      void MergeShowerIonizations(std::vector<supera::ParticleGroup>& part_grp_v);
//      void ApplyEnergyThreshold(std::vector<supera::ParticleGroup>& part_grp_v);
//      void MergeShowerConversion(std::vector<supera::ParticleGroup>& part_grp_v);
//      void MergeShowerFamilyTouching(const larcv::Voxel3DMeta& meta,
//                                     std::vector<supera::ParticleGroup>& part_grp_v);
//      void MergeShowerTouching(const larcv::Voxel3DMeta& meta,
//                               std::vector<supera::ParticleGroup>& part_grp_v);
//      void MergeShowerTouching2D(std::vector<supera::ParticleGroup>& part_grp_v);
//      void MergeShowerDeltas(std::vector<supera::ParticleGroup>& part_grp_v);
//      void DumpHierarchy(size_t trackid,
//                         const std::vector<supera::ParticleGroup>& part_grp_v) const;
//      std::vector<unsigned int> ParentTrackIDs(size_t trackid) const;
//      std::vector<unsigned int> ParentShowerTrackIDs(size_t trackid,
//                                                     const std::vector<supera::ParticleGroup>& part_grp_v,
//                                                     bool include_lescatter=false) const;

      supera::MCParticleList _mc_part_list;
      std::string _ref_meta3d_cluster3d;
      std::string _ref_meta3d_tensor3d;
      std::string _ref_meta2d_tensor2d;
      std::string _output_label;
      std::string _masked_true2reco_cluster3d;
      std::string _masked_true_tensor3d;
      bool _use_sed_points;
      size_t _eioni_size;
      size_t _delta_size;
      size_t _compton_size;
      double _compton_energy;
      double _edep_threshold;
      bool _use_true_pos;
      bool _use_sed;
      bool _check_particle_validity;
      int  _projection_id;

      // todo: we aren't currently setting these in configure(), so don't allow them to be used
//      larcv::BBox3D _world_bounds;
//      std::vector<std::vector<std::vector<int> > > _scan;
//      size_t _valid_nplanes;

      std::vector<size_t> _semantic_priority;

  };

}
#endif //LARCV2_SUPERAMCPARTICLECLUSTER_H
