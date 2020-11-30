/**
 * \file SuperaG4HitSegment.h
 *
 * \ingroup Package_Name
 *
 * \brief Extract truth info & voxelize true hit segments from edepsim
 *
 * @author J. Wolcott <jwolcott@fnal.gov>
 */

/** \addtogroup Package_Name

    @{*/
#ifndef __SUPERAG4HITSEGMENT_H__
#define __SUPERAG4HITSEGMENT_H__
//#ifndef __CINT__
//#ifndef __CLING__
#include "SuperaBase.h"
#include "EDepSim/TG4Event.h"
#include "larcv/core/DataFormat/EventParticle.h"
#include "larcv/core/DataFormat/EventVoxel3D.h"

// forward declarations
class TVector3;

namespace larcv {

  template <typename T>
  class AABBox;

  /**
     \class SuperaG4HitSegment
     Responsible for defining a rectangular volume boundary and voxelization
  */
  class SuperaG4HitSegment : public SuperaBase
  {

    public:

      /// Default constructor
      SuperaG4HitSegment(const std::string name = "SuperaG4HitSegment");

      /// Default destructor
      ~SuperaG4HitSegment() {}

      void configure(const PSet&);

      void initialize();

      bool process(IOManager& mgr);

      void finalize();

    private:
      /// Apply Landau fluctuations to the set of energies in these voxels
      void FluctuateEnergy(std::vector<Voxel> &voxels);

      static std::string GetTrajCreationProc(int creationCode, int parentPdgCode) ;

      larcv::Particle MakeParticle(const TG4Trajectory &traj,
                                   const TG4Trajectory *parentTraj,
                                   const larcv::AABBox<double> &bbox);


      std::string _sparsetensor3d_producer;
      std::string _particle_producer;
      std::vector<std::string> _active_volume_filter;  ///<  should we keep hits from only certain active volumes?  (keeps all if none specified)

  };

  /**
     \class larcv::SuperaBBoxInteractionFactory
     \brief A concrete factory class for larcv::SuperaG4HitSegment
  */
  class SuperaG4HitSegmentProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    SuperaG4HitSegmentProcessFactory() { ProcessFactory::get().add_factory("SuperaG4HitSegment", this); }
    /// dtor
    ~SuperaG4HitSegmentProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) { return new SuperaG4HitSegment(instance_name); }
  };

}
//#endif
//#endif
#endif
/** @} */ // end of doxygen group

