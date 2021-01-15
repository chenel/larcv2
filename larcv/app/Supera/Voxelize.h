/**
 * \file Voxelize.h
 *
 * \brief Utilities for segmenting truth energy deposits into 3D voxels
 *
 * @author J. Wolcott <jwolcott@fnal.gov>
 */

#ifndef LARCV2_VOXELIZE_H
#define LARCV2_VOXELIZE_H

#include <vector>

#include "larcv/core/DataFormat/Particle.h"
#include "larcv/core/DataFormat/Voxel.h"
#include "larcv/core/DataFormat/Voxel3DMeta.h"
#include "EDepSim/TG4Event.h"

namespace larcv
{
  // some forward declarations
  template <typename T>
  class AABBox;

  template <typename T>
  class Vec3;
  typedef Vec3<float> Vec3f;
  typedef Vec3<double> Vec3d;



  /// Split an edep-sim TG4HitSegment true energy deposit into voxels
  ///
  /// \param[in]  hitSegment   The TG4HitSegment to operate on
  /// \param[in]  meta         Metadata about the voxel array
  /// \param[in] particles     Collection of true particles whose recorded energy deposits will be updated.  Pass empty vector to disable
  /// \return                  A vector of voxels
  std::vector <larcv::Voxel>
  MakeVoxels(const ::TG4HitSegment &hitSegment,
             const larcv::Voxel3DMeta &meta,
             std::vector <larcv::Particle> &particles);

  /// Split an edep-sim TG4HitSegment true energy deposit into voxels
  /// (same as above, but when no Particles' energy deposits are to be updated)
  /// \param[in]  hitSegment   The TG4HitSegment to operate on
  /// \param[in]  meta         Metadata about the voxel array
  /// \return                  A vector of voxels
  std::vector <larcv::Voxel>
  MakeVoxels(const ::TG4HitSegment &hitSegment,
             const larcv::Voxel3DMeta &meta);

  /// Where (if anywhere) does a line segment intersect a given bounding box?
  /// (If the entire line segment is contained, the entry and exit points
  ///  will be set to the start and stop points provided.)
  ///
  /// \param bbox        Bounding box in question
  /// \param startPoint  Start point of the line segment in question
  /// \param stopPoint   Stop point of the line segment in question
  /// \param entryPoint  Computed entry point of the line segment into the box, if any
  /// \param exitPoint   Computed exit point of the line segment from the box, if any
  /// \return            Number of intersections (will be 0, 1, or 2)
  template <typename T>
  char Intersections(const larcv::AABBox<T> & bbox,
                     const TVector3 & startPoint,
                     const TVector3 & stopPoint,
                     larcv::Vec3<T> & entryPoint,
                     larcv::Vec3<T> & exitPoint);

}

#endif //LARCV2_VOXELIZE_H
