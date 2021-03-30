#ifndef __PARTICLECOUNTFILTER_CXX__
#define __PARTICLECOUNTFILTER_CXX__

#include "ParticleCountFilter.h"
#include "larcv/core/DataFormat/EventParticle.h"

namespace larcv {

  static ParticleCountFilterProcessFactory __global_ParticleCountFilterProcessFactory__;

  ParticleCountFilter::ParticleCountFilter(const std::string name)
    : ProcessBase(name)
  {}
    
  void ParticleCountFilter::configure(const PSet& cfg)
  {
    _part_producer = cfg.get<std::string>("ParticleProducer");
    _max_part_count = cfg.get<size_t>("MaxCount", std::numeric_limits<size_t>::max());
    _min_part_count = cfg.get<size_t>("MinCount",0);
    auto shapes = cfg.get<std::vector<int>>("ParticleShapes");
    for (const auto & shape : shapes)
      _part_shapes.insert(static_cast<larcv::ShapeType_t>(shape));
  }

  void ParticleCountFilter::initialize()
  {}

  bool ParticleCountFilter::process(IOManager& mgr)
  {
    auto const& part_v = mgr.get_data<larcv::EventParticle>(_part_producer).as_vector();
    if(part_v.size() >= _part_count_v.size()) {
      _part_count_v.resize(part_v.size()+1,0);
    }
    std::size_t size = 0;
    if (!_part_shapes.empty())
    {
      for (const auto & part : part_v)
      {
        if (_part_shapes.find(part.shape()) != _part_shapes.end())
          size++;
      }
    }
    else
      size = part_v.size();
    _part_count_v[size] += 1;
    
    return (_min_part_count <= part_v.size() && part_v.size() <= _max_part_count);
  }

  void ParticleCountFilter::finalize()
  {
    double total_count = 0;
    for(auto const& v : _part_count_v) total_count += v;
    LARCV_NORMAL() << "Reporting Particle counts. Total events processed: " << (int)total_count << std::endl;
    for(size_t i=0; i<_part_count_v.size(); ++i) 
      LARCV_NORMAL() << "    Multi=" << i 
		     << " ... " << _part_count_v[i] << " events (" 
		     << _part_count_v[i] / total_count * 100 << " %)" << std::endl;
  }

}
#endif
