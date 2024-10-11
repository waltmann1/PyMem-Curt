#include "integrator_brownian_vertex.hpp"


void IntegratorBrownianMeshVertex::poststep(void)
{
  for (int id = 0; id < _system.Numvertices; id ++)
  {
    real3 force_rnd;
    force_rnd.x = force_rnd.y = force_rnd.z = 0.0;
    if (this->get_temperature() > 0.0)
    {
      force_rnd.x = B *rng->gauss_rng();
      force_rnd.y = B *rng->gauss_rng();
      force_rnd.z = B *rng->gauss_rng();
    }
    // Update particle position
    _system.vertices[id].r.x += mu * this->get_time_step() * _system.vertices[id].forceC.x + sqrt_dt * force_rnd.x;
    _system.vertices[id].r.y += mu * this->get_time_step() * _system.vertices[id].forceC.y + sqrt_dt * force_rnd.y;
    _system.vertices[id].r.z += mu * this->get_time_step() * _system.vertices[id].forceC.z + sqrt_dt * force_rnd.z;
    
    _system.vertices[id].v.x = _system.vertices[id].forceC.x + sqrt_dt * force_rnd.x;
    _system.vertices[id].v.y = _system.vertices[id].forceC.y + sqrt_dt * force_rnd.y;
    _system.vertices[id].v.z = _system.vertices[id].forceC.z + sqrt_dt * force_rnd.z;
  }
}

