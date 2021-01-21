/**
 * \file SuperMCParticleClusterData.h
 *
 * \brief Storage classes for use in SuperaMCParticleCluster
 *        Ported from DeepLearnPhysics/Supera
 *
 * @author J. Wolcott <jwolcott@fnal.gov>
 */
#ifndef LARCV2_SUPERAMCPARTICLECLUSTERDATA_H
#define LARCV2_SUPERAMCPARTICLECLUSTERDATA_H

#include "larcv/core/Base/LArCVTypes.h"
#include "larcv/core/DataFormat/Particle.h"
#include "larcv/core/DataFormat/Voxel.h"

namespace supera
{
  enum ProcessType
  {
    kTrack,
    kNeutron,
    kPhoton,
    kPrimary,
    kCompton,       ///< detach low E shower
    kComptonHE,     ///< detach high E shower
    kDelta,         ///< attach-mu low E special
    kConversion,    ///< detach high E gamma
    kIonization,    ///< attach-e low E
    kPhotoElectron, ///< detatch low E
    kDecay,         ///< attach high E
    kOtherShower,   ///< anything else (low E)
    kOtherShowerHE, ///< anything else (high E)
    kInvalidProcess
  };

  // --------------------------------------------------------

  struct EDep
  {
    double x = larcv::kINVALID_DOUBLE;
    double y = larcv::kINVALID_DOUBLE;
    double z = larcv::kINVALID_DOUBLE;
    double t = larcv::kINVALID_DOUBLE;
    double e = larcv::kINVALID_DOUBLE;
  };

  // --------------------------------------------------------

  struct ParticleGroup
  {
    // ----------- data first

    /// the "top" particle
    larcv::Particle part;

    /// all the GEANT4 tracks corresponding to this particle and those grouped with it
    std::vector<size_t> trackid_v;

    /// Should this ParticleGroup be used?  (set to false, e.g., for groups that have been subsumed into others)
    bool valid;

    /// What sort of grouping does this correspond to?
    ProcessType type;

    /// voxels associated with this particle group
    larcv::VoxelSet vs;

    EDep first_pt;
    EDep last_pt;

    // ----------- now methods

    // constructor
    ParticleGroup();

    /// Associate a true energy deposition to this particle group
    void AddEDep(const EDep& pt);

    /// Add another particle group into this one
    ///
    /// \param child         The child particle to be added
    /// \param updatePoints  Update the beginning/end of this particle with the child's if the child extends further in either direction
    void Merge(ParticleGroup &child, bool updatePoints = true);

    /// Total count of voxels (3D + 2D) associated with this particle group
    std::size_t size_all() const;

    /// Semantic classification (see larcv::ShapeType_t)
    larcv::ShapeType_t shape() const;

  };

  // --------------------------------------------------------

  /// Sum up the total number of voxels in a collection of ParticleGroups.
  ///
  /// \param pgs           Vector of ParticleGroups to be considered
  /// \param inclInvalid   Include ParticleGroups for which !group.valid? (i.e. those that have been merged into another)
  /// \return              the count
  std::size_t CountVoxels(const std::vector<ParticleGroup>& pgs, bool inclInvalid=false);

  /// Sum up the total deposited energy of voxels in a collection of ParticleGroups.
  ///
  /// \param pgs           Vector of ParticleGroups to be considered
  /// \param inclInvalid   Include ParticleGroups for which !group.valid? (i.e. those that have been merged into another)
  /// \return              the summed energy
  float SumVoxelsEdep(const std::vector<ParticleGroup>& pgs, bool inclInvalid= false);


}

#endif //LARCV2_SUPERAMCPARTICLECLUSTERDATA_H
