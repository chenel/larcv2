#include "Voxelize.h"

#include <numeric>

#include "SuperaG4HitSegment.h"
#include "GenRandom.h"
#include "geometry.h"
#include "raybox.h"

namespace larcv
{
  // --------------------------------------------------
  std::vector<larcv::Voxel> MakeVoxels(const ::TG4HitSegment &hitSegment,
                                       const larcv::Voxel3DMeta &meta)
  {
    std::vector<larcv::Particle> noparticles;
    return MakeVoxels(hitSegment, meta, noparticles);
  }

    // --------------------------------------------------
  std::vector<larcv::Voxel> MakeVoxels(const ::TG4HitSegment &hitSegment,
                                       const larcv::Voxel3DMeta &meta,
                                       std::vector<larcv::Particle> &particles)
  {
    std::vector<larcv::Voxel> voxels;
    larcv::AABBox<double> box(meta);

    const double epsilon = 1.e-3;
    double dist_travel = 0.;
    double energy_deposit = 0.;
    double smallest_side = std::min(meta.size_voxel_x(),std::min(meta.size_voxel_y(),meta.size_voxel_z()));
    LARCV_SDEBUG() <<"World: " << box.bounds[0] << " => " << box.bounds[1] << std::endl;

    TVector3 start = hitSegment.GetStart().Vect();
    TVector3 end = hitSegment.GetStop().Vect();
    start *= 0.1;  // convert unit to cm
    end *= 0.1;

    const auto FindParticle = [&particles](int trackId) -> larcv::Particle*
    {
      auto itPart = std::find_if(particles.begin(), particles.end(),
                                 [=](const larcv::Particle & p) { return p.track_id() == static_cast<unsigned int>(trackId); });
      if (itPart == particles.end())
        return nullptr;

      return &(*itPart);
    };

    int trackId; // = hitSegment.GetPrimaryId();
    if (hitSegment.Contrib.size() == 1)
      trackId = hitSegment.Contrib.front();
    //else if (std::find(hitSegment.Contrib.begin(), hitSegment.Contrib.end(), hitSegment.GetPrimaryId()) != hitSegment.Contrib.end())
    //  trackId = hitSegment.GetPrimaryId();
    else
    {
      LARCV_SWARNING() << "Could not determine which GEANT track ID to assign edep-sim energy to!" << std::endl;
      std::stringstream trks;
      std::for_each(std::begin(hitSegment.Contrib), std::end(hitSegment.Contrib),
                    [&trks](const int trk) { trks << " " << trk; });
      LARCV_SWARNING() << "  Chose the first of:" << trks.str() << std::endl;
      trackId = hitSegment.Contrib.front();
    }
    auto particle = FindParticle(trackId);

    LARCV_SDEBUG() << "Voxelizing TG4HitSegment for GEANT track " << trackId
                 << " from (" << start.x() << "," << start.y() << "," << start.z() << ")"
                 << " to (" << end.x() << "," << end.y() << "," << end.z() << ")"
                 << ", length = " << (end - start).Mag() << " cm"
                 << std::endl;

    larcv::Vec3d pt0, pt1;
    char crossings = Intersections(box, start, end, pt0, pt1);

    if(crossings == 0) {
      LARCV_SDEBUG() << "No crossing point found..." << std::endl;
      return voxels;
    }

    LARCV_SDEBUG() << "   Intersects with bounding box at"
                  << " (" << pt0.x << "," << pt0.y << "," << pt0.z << ")"
                  << " and (" << pt1.x << "," << pt1.y << "," << pt1.z << ")"
                  << std::endl;

    larcv::Vec3d dir = pt1 - pt0;
    double length = dir.length();
    dir.normalize();
    larcv::Ray<double> ray(pt0, dir);


    voxels.reserve(voxels.size() + (size_t)(length / smallest_side));
    double t0, t1, dist_section;
    dist_section = 0.;
    size_t nx, ny, nz;
    t0=t1=0.;
    //size_t ctr=0;
    while(true) {
      //ctr += 1;
      //if(ctr>=10) break;
      // define the inspection box
      Vec3d pt = pt0 + dir * (t1 + epsilon);
      LARCV_SDEBUG() << "    New point: " << pt << std::endl;
      auto vox_id = meta.id((double)(pt.x), (double)(pt.y), (double)(pt.z));
      if(vox_id==larcv::kINVALID_VOXELID) break;
      meta.id_to_xyz_index(vox_id, nx, ny, nz);
      box.bounds[0].x = meta.min_x() + nx * meta.size_voxel_x();
      box.bounds[0].y = meta.min_y() + ny * meta.size_voxel_y();
      box.bounds[0].z = meta.min_z() + nz * meta.size_voxel_z();
      box.bounds[1].x = box.bounds[0].x + meta.size_voxel_x();
      box.bounds[1].y = box.bounds[0].y + meta.size_voxel_y();
      box.bounds[1].z = box.bounds[0].z + meta.size_voxel_z();
      LARCV_SDEBUG() << "    Inspecting a voxel id " << vox_id << " ... " << box.bounds[0] << " => " << box.bounds[1] << std::endl;
      auto cross = box.intersect(ray,t0,t1);

      // no crossing
      if(cross==0) {
        LARCV_SERROR() << "      No crossing (not expected) ... breaking" << std::endl;
        break;
      }
      double dx;
      if(cross==1) {
        LARCV_SDEBUG() << "      One crossing: " << pt0 + dir * t1 << std::endl;
        dx = std::min(t1,length);
      }else {
        LARCV_SDEBUG() << "      Two crossing" << pt0 + dir * t0 << " => " << pt0 + dir * t1 << std::endl;
        if (t0 > length)
          dx = length;
        else if (t1 > length)
          dx = length - t0;
        else
          dx = t1 - t0;
      }
      /*
      res_pt[0] = (nx+0.5) * meta.size_voxel_x();
      res_pt[1] = (ny+0.5) * meta.size_voxel_y();
      res_pt[2] = (nz+0.5) * meta.size_voxel_z();
      res_pt[3] = dx;
      res.push_back(res_pt);
      */
      double energyInVoxel = dx / length * hitSegment.GetEnergyDeposit();
      if (energyInVoxel < 0)
      {
        LARCV_SCRITICAL() << "Voxel with negative energy deposited!" << std::endl
                          << "  ID = " << vox_id << std::endl
                          << "  edep computed from:" << std::endl
                          << "      dx = " << dx
                              << ", length = " << length
                              << ", TG4HitSegment edep = " << hitSegment.GetEnergyDeposit()
                          << std::endl;
        throw larbys("Voxelize.cxx:MakeVoxels(): Negative energy deposit in voxel");
      }
      else if (energyInVoxel > 0)
        voxels.emplace_back(vox_id, energyInVoxel);
      dist_travel += dx;
      dist_section += dx;
      energy_deposit += energyInVoxel;
      //LARCV_SDEBUG() << "      Registering voxel id " << vox_id << " at distance fraction " << t1/length << std::endl;
      LARCV_SDEBUG() << "      Registering voxel id " << vox_id << " t1 =" << t1 << " (total length = " << length << ")" << std::endl;
      if(t1>length) {
        LARCV_SDEBUG() << "      Reached the segment end (t1 = " << t1 << " fractional length " << t1/length << ") ... breaking" << std::endl;
        break;
      }

      LARCV_SDEBUG() << "      Updated t1 = " << t1 << " (fractional length " << t1/length << ")" << std::endl;
    }
    LARCV_SDEBUG() << "Made " << voxels.size() << " voxels with total Edep = " << energy_deposit << std::endl;

    if (particle)
    {
      particle->energy_deposit(particle->energy_deposit() + energy_deposit);
      particle->num_voxels(particle->num_voxels() + voxels.size());
    }
    return voxels;
  }

  // --------------------------------------------------

  template <typename T>
  char Intersections(const AABBox<T> &bbox,
                     const TVector3 &startPoint,
                     const TVector3 &stopPoint,
                     Vec3<T> &entryPoint,
                     Vec3<T> &exitPoint)
  {
    TVector3 displVec = (stopPoint - startPoint);
    TVector3 dir = displVec.Unit();
    larcv::Ray<T> ray(startPoint, dir);

    bool startContained = bbox.contain(startPoint);
    bool stopContained = bbox.contain(stopPoint);

    if (startContained)
      entryPoint = startPoint;
    if (stopContained)
      exitPoint = stopPoint;

    if(!startContained || !stopContained)
    {
      double t0, t1;
      int cross = bbox.intersect(ray, t0, t1);
      // note that AABBox::intersect() will trace the ray to infinity in both directions,
      // which may result in intersections beyond our segment of interest
      if (cross > 0)
      {
        if ((!startContained && t0 < 0) || t0 > displVec.Mag())
          cross--;
        if (t1 < 0 || t1 > displVec.Mag())
          cross--;
      }

      if (cross > 0)
      {
        const T epsilon = 0.0001;
        if (!startContained)
          entryPoint = startPoint + (t0 + epsilon) * dir;

        if (!stopContained)
          exitPoint = startPoint + (t1 - epsilon) * dir;
      }

      LARCV_SDEBUG() << "Number of crossings=" << cross
                    << " for bounding box " << bbox.bounds[0] << "-" << bbox.bounds[1]
                    << " and ray between " << "(" << startPoint.x() << "," << startPoint.y() << "," << startPoint.z() << ")"
                    << " and (" <<  stopPoint.x() << "," << stopPoint.y() << "," << stopPoint.z() << ")" << std::endl;
      if (cross > 0)
      {
        LARCV_SDEBUG() << "Start point contained?: " << startContained << std::endl;
        if (!startContained)
          LARCV_SDEBUG() << "  entry point: " << entryPoint << "; t0=" << t0 << std::endl;
        LARCV_SDEBUG() << "Stop point contained?: " << stopContained << std::endl;
        if (!stopContained)
          LARCV_SDEBUG() << "  exit point: " << exitPoint << "; t1=" << t1 << std::endl;
      }

      if (cross == 1 && startContained == stopContained)
      {
        LARCV_SERROR() << "Unexpected number of crossings (" << cross << ")"
                      << " for bounding box and ray between "
                      << "(" << startPoint.x() << "," << startPoint.y() << "," << startPoint.z() << ")"
                      << " and (" <<  stopPoint.x() << "," << stopPoint.y() << "," << stopPoint.z() << ")" << std::endl;
        LARCV_SERROR() << "Start point contained?: " << startContained << ".  Stop point contained?: " << stopContained << std::endl;
      }

      return cross;
    }

    return 2;
  }

  // instantiate the template for the type(s) we care about
  template char Intersections(const AABBox<double> &bbox,
                              const TVector3 &startPoint,
                              const TVector3 &stopPoint,
                              Vec3d &entryPoint,
                              Vec3d &exitPoint);


}