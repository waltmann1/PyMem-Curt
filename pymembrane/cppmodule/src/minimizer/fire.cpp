/*
* @Author: siyu
* @Date:   2020-07-26 11:06:58
* @Last Modified by:   siyu
* @Last Modified time: 2020-08-11 22:55:44
*/
#include "fire.hpp"

/** brief Constant times a vector **/
#define Xvec1(v, a, v1) \
  (v.x = (a)*v1.x),     \
      (v.y = (a)*v1.y), \
      (v.z = (a)*v1.z) /** brief Constant times a vector **/
#define Xvec2(v, a, v1, b, v2)     \
  (v.x = (a)*v1.x + (b)*v2.x),     \
      (v.y = (a)*v1.y + (b)*v2.y), \
      (v.z = (a)*v1.z + (b)*v2.z) /** brief Constant times a vector **/
#define Xvec3(v, a, v1, b, v2, c, v3)         \
  (v.x = (a)*v1.x + (b)*v2.x + (c)*v3.x),     \
      (v.y = (a)*v1.y + (b)*v2.y + (c)*v3.y), \
      (v.z = (a)*v1.z + (b)*v2.z + (c)*v3.z) /** brief Constant times a vector **/
#define Xvec4(v, a, v1, b, v2, c, v3, d, v4)             \
  (v.x = (a)*v1.x + (b)*v2.x + (c)*v3.x + (d)*v4.x),     \
      (v.y = (a)*v1.y + (b)*v2.y + (c)*v3.y + (d)*v4.y), \
      (v.z = (a)*v1.z + (b)*v2.z + (c)*v3.z + (d)*v4.z)

void MinimizerMeshFIRE::reset(void)
{
  m_converged = false;
  m_n_since_negative = m_nmin + 1;
  m_n_since_start = 0;
  m_alpha = m_alpha_start;
  m_energy_total = 0.0;
  m_old_energy = 0.0;
  m_dE = 0.0;
  m_fnorm = 0.0;
  m_dt = 0.005;
  this->update_time_step_parameters();
  for (int id = 0; id < _system.Numvertices; id++)
  {
    _system.vertices[id].v.x = 0.0;
    _system.vertices[id].v.y = 0.0;
    _system.vertices[id].v.z = 0.0;
  }
}

void MinimizerMeshFIRE::poststep1(void)
{
  // NVT, velocity verlet
  for (int id = 0; id < _system.Numvertices; id++)
  {
    Xvec3(_system.vertices[id].r,
          1., _system.vertices[id].r,
          m_dt, _system.vertices[id].v,
          0.5 * m_dt * m_dt, _system.vertices[id].forceC);
    Xvec2(_system.vertices[id].v,
          1., _system.vertices[id].v,
          0.5 * m_dt, _system.vertices[id].forceC);
  }
}
void MinimizerMeshFIRE::poststep2(void)
{
  for (int id = 0; id < _system.Numvertices; id++)
  {
    Xvec2(_system.vertices[id].v,
          1., _system.vertices[id].v,
          0.5 * m_dt, _system.vertices[id].forceC);
  }
}
void MinimizerMeshFIRE::poststep3(void)
{
  //FIRE
  double P = 0.;
  double vnorm = 0.;
  double fnorm = 0.;
  double tnorm = 0.;
  double energy = 0.;

  for (int id = 0; id < _system.Numvertices; id++)
  {
    P += vdot(_system.vertices[id].forceC, _system.vertices[id].v);
    fnorm += vdot(_system.vertices[id].forceC, _system.vertices[id].forceC);
    vnorm += vdot(_system.vertices[id].v, _system.vertices[id].v);
    energy += _system.vertices[id].energy;
  }

  for (int id = 0; id < _system.Numedges; id++)
    energy += _system.edges[id].energy;
  for (int id = 0; id < _system.Numfaces; id++)
    energy += _system.faces[id].energy;

  fnorm = sqrt(fnorm);
  vnorm = sqrt(vnorm);
  m_fnorm = fnorm;
  m_energy_total = energy;
  m_dE = m_energy_total - m_old_energy;
  // energy = energy/(_system.Numvertices);

  // if ((fnorm/sqrt(3.)/(_system.Numvertices) < m_ftol && fabs(energy-m_old_energy) < m_etol) && m_n_since_start >= m_run_minsteps)
  if (fnorm < m_ftol && fabs(energy - m_old_energy) < m_etol && m_n_since_start >= m_run_minsteps)
  {
    py::print(name, " converged, ftol:" , m_ftol , " m_etol:" , m_etol);
    m_converged = true;
    return;
  }

  // double factor_t;
  // if (fabs(fnorm) > EPSILON)
  //     factor_t = m_alpha*vnorm/fnorm;
  // else
  //     factor_t = 1.0;

  for (int id = 0; id < _system.Numvertices; id++)
  {
    Xvec2(_system.vertices[id].v,
          1.0 - m_alpha, _system.vertices[id].v,
          m_alpha * vnorm / fnorm, _system.vertices[id].forceC);
  }

  if (P > 0.0)
  {
    m_n_since_negative++;
    if (m_n_since_negative > m_nmin)
    {
      m_dt = std::min(m_dt, m_dT_max);
      m_alpha *= m_falpha;
    }
  }
  else if (P <= 0.0)
  {
    m_dt = m_dt * m_fdec;
    m_alpha = m_alpha_start;
    m_n_since_negative = 0;

    //py::print( this->name, " zero velocities");
    //py::print("P:" , P , " fnorm:" , fnorm , " vnorm:" , vnorm , " delta energy:" , fabs(energy - m_old_energy) , " energy/N:" , energy / (_system.Numvertices) , "alpha:" , m_alpha);

    for (int id = 0; id < _system.Numvertices; id++)
    {
      _system.vertices[id].v.x = 0.0;
      _system.vertices[id].v.y = 0.0;
      _system.vertices[id].v.z = 0.0;
    }
  }

  m_n_since_start++;
  m_old_energy = energy;
}

void MinimizerMeshFIRE::minimize(void)
{
  for (auto step = 0; step < m_max_iter; step++)
  {
    _evolver.reset_mesh_forces();
    _evolver.compute_mesh_forces();
    poststep1();

    _evolver.reset_mesh_forces();
    _evolver.compute_mesh_forces();
    poststep2();

    //enforce the constraints before calculate the energy
    _evolver.enforce_mesh_constraints();

    _evolver.reset_mesh_energy();
    _evolver.compute_mesh_energy();
    poststep3();

    if (m_converged == true)
      break;
  }
  //py::print(" Minimizer ", this->get_info());
}
