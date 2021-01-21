#include "SuperaMCParticleClusterData.h"

#include <numeric>

#include "larcv/core/Base/larcv_logger.h"

namespace supera
{
  // -------------------------------------

  ParticleGroup::ParticleGroup()
    : valid(false), type(kInvalidProcess)
  {
    last_pt.t = -1.e9;
  }

  // -------------------------------------

  void ParticleGroup::AddEDep(const EDep& pt)
  {
    if(pt.x == larcv::kINVALID_DOUBLE) return;
    //start.AddEDep(pt);
    if(pt.t < first_pt.t)
    {
      first_pt = pt;
      LARCV_SDEBUG() << " Start point of track " << part.track_id()
                     << " now (" << this->first_pt.x << "," << this->first_pt.y << "," << this->first_pt.z << ")"
                     << std::endl;
    }
    if(pt.t > last_pt.t)
    {
      last_pt = pt;
      LARCV_SDEBUG() << " End point of track " << part.track_id()
                     << " now (" << this->last_pt.x << "," << this->last_pt.y << "," << this->last_pt.z << ")"
                     << std::endl;
    }
  }

  // -------------------------------------

  void ParticleGroup::Merge(ParticleGroup &child, bool updatePoints)
  {
    for(auto const& vox : child.vs.as_vector())
      this->vs.emplace(vox.id(),vox.value(),true);

    LARCV_SDEBUG() << "Parent track id " << this->part.track_id()
                   << " PDG " << this->part.pdg_code() << " " << this->part.creation_process() << std::endl
                   << "  ... merging " << child.part.track_id()
                   << " PDG " << child.part.pdg_code() << " " << child.part.creation_process() << std::endl;
    /*
    for(auto const& pt : child.start.pts)
this->AddEDep(pt);
    */
    if(updatePoints)
    {
      this->AddEDep(child.last_pt);
      this->AddEDep(child.first_pt);
    }

    this->trackid_v.push_back(child.part.track_id());
    for(auto const& trackid : child.trackid_v)
      this->trackid_v.push_back(trackid);
    child.vs.clear_data();
    child.valid=false;
  }

  // -------------------------------------

  larcv::ShapeType_t ParticleGroup::shape() const
  {
    // identify delta ray
    if(type == kInvalidProcess)
      return larcv::kShapeUnknown;
    if(type == kDelta)
      return larcv::kShapeDelta;
    if(type == kNeutron) //return larcv::kShapeUnknown;
      return larcv::kShapeLEScatter;
    if(part.pdg_code() == 11 || part.pdg_code() == 22 || part.pdg_code() == -11)
    {
      if(type == kComptonHE || type == kPhoton || type == kPrimary || type == kConversion || type==kOtherShowerHE)
        return larcv::kShapeShower;
      if(type == kDecay)
      {
        if(part.parent_pdg_code() == 13 || part.parent_pdg_code() == -13)
          return larcv::kShapeMichel;
        else
          return larcv::kShapeShower;
      }
      return larcv::kShapeLEScatter;
    }
    else
      return larcv::kShapeTrack;
  }

  // -------------------------------------

  std::size_t ParticleGroup::size_all() const
  {
    std::size_t res=vs.size();
    return res;
  }

  // -------------------------------------

  std::size_t CountVoxels(const std::vector<ParticleGroup> &pgs, bool inclInvalid)
  {
    return std::accumulate(std::begin(pgs),
                           std::end(pgs),
                           std::size_t(0),
                           [inclInvalid](std::size_t s, const supera::ParticleGroup & grp)
                           {
                             if (!inclInvalid || grp.valid)
                               s += grp.vs.size();
                             return s;
                           });
  }

  // -------------------------------------

  float SumVoxelsEdep(const std::vector<ParticleGroup> &pgs, bool inclInvalid)
  {
    return std::accumulate(std::begin(pgs),
                           std::end(pgs),
                           0.0f,
                           [inclInvalid](double s, const supera::ParticleGroup & grp)
                           {
                             if (!inclInvalid || grp.valid)
                               s += grp.vs.sum();
                             return s;
                           });
  }
}