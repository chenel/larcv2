#ifndef SUPERA_MCPARTICLELIST_CXX
#define SUPERA_MCPARTICLELIST_CXX

#include "MCParticleList.h"

#include "larcv/core/DataFormat/Particle.h"

namespace supera{

  void MCParticleList::Update(const std::vector<larcv::Particle> &particles, int run, int event)
  {
    if(run == _run && event == _event) return;

    _run = run;
    _event = event;

    _trackid_v.resize(particles.size());
    _pdgcode_v.resize(particles.size());
    _parent_index_v.resize(particles.size());
    _parent_trackid_v.resize(particles.size());
    _parent_pdg_v.resize(particles.size());

    _trackid2index.clear();
    _trackid2index.resize(std::max(_trackid2index.size(), particles.size()), -1);

    for(size_t index=0; index < particles.size(); ++index) {
      auto const& mcpart = particles[index];
      _trackid_v[index] = mcpart.track_id();
      _pdgcode_v[index] = abs(mcpart.pdg_code());
      _parent_trackid_v[index] = mcpart.parent_id();
      if(mcpart.track_id() >= _trackid2index.size())
        _trackid2index.resize(mcpart.track_id()+1,-1);
      _trackid2index[mcpart.track_id()] = index;
    }

    for(size_t index=0; index < particles.size(); ++index) {
      auto const& mcpart = particles[index];
      int mother_id      = mcpart.parent_id();
      int mother_index   = -1;
      if(mother_id < ((int)(_trackid2index.size()))) {
        mother_index = std::distance(particles.begin(),
                                     std::find_if(particles.begin(), particles.end(),
                                                  [=](const larcv::Particle & part) {return part.id() == mother_id;}));
        if(mother_index >= 0) {
          _parent_pdg_v[index] = particles[mother_index].pdg_code();
          _parent_index_v[index] = mother_index;
        }
      }
    }
  }
}
#endif
