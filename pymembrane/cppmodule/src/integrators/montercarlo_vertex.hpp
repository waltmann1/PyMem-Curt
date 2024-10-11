/*!
 * \file montercarlo_vertex.hpp
 * \author Daniel Matoz , fdamatoz@gmail.com
 * \date 29-Dec-2018
 * \brief Declaration of IntegratorBrownianVertex class
 */

#ifndef __montercarlo_vertex_HPP__
#define __montercarlo_vertex_HPP__

/** @addtogroup integrators Vertex Brownian Integrator
 *  @brief IntegratorBrownianVertex class
 *  @{
 */

#include <iostream>
#include "montecarlointegrator.hpp"
#include "../rng/rng.hpp"
#include "../types/globaltypes.hpp"
#include "../utils/fromstring.hpp"
#include "../system/systemclass.hpp"
/**
 * @class MonteCarloIntegratorVertex
 * @brief Integrator Brownian class implements Brownian dynamics for the vertex position. Particle director will not be integrated
 */
class MonteCarloIntegratorVertex : public MonteCarloIntegrator
{
public:
  /** @brief VertexIntegrator Constructor */
  MonteCarloIntegratorVertex(SystemClass &system, VertexCompute &potentials) : MonteCarloIntegrator(system, potentials)
  {
    type = "monte carlo";
    name = "vertex move";
    this->set_default_properties();
  }
  /** @brief destructor */
  ~MonteCarloIntegratorVertex() {}

  void set_default_properties(void)
  {
    host::vector<bool> _freezed_vertex(NUM_TYPES_ALLOWED, false);
    freezed_vertex = _freezed_vertex;
    this->set_temperature(0.0);
    m_seed = 123456; ///default value
    m_rng = std::make_unique<RNG>(m_seed);
    m_spherical_move = false;
    m_dx = 1e-2;
    m_dy = 1e-2;
    m_dz = 1e-2;
  }

  using MonteCarloIntegrator::set_property;
  void set_property(std::map<std::string, std::string> &value_map) override
  {
    for (const auto &item : value_map)
    {
      if (item.first.compare("T") == 0)
      {
        this->set_temperature(util::from_string_double(item.second));
        this->update_temperature_parameters();
      }
      else if (item.first.compare("dr") == 0)
      {
        m_dx = util::from_string_double(item.second);
        m_dy = util::from_string_double(item.second);
        m_dz = util::from_string_double(item.second);
      }
      else if (item.first.compare("drx") == 0)
      {
        m_dx = util::from_string_double(item.second);
      }
      else if (item.first.compare("dry") == 0)
      {
        m_dy = util::from_string_double(item.second);
      }
      else if (item.first.compare("drz") == 0)
      {
        m_dz = util::from_string_double(item.second);
      }
      else if (item.first.compare("seed") == 0)
      {
        m_seed = uint(util::from_string_int(item.second));
        m_rng = std::make_unique<RNG>(m_seed);
      }
      else if (item.first.compare("every step") == 0)
      {
        this->set_integrate_every(util::from_string_int(item.second));
      }
      else if (item.first.compare("spherical_move") == 0)
      {
        m_spherical_move = util::from_string_bool(item.second);
      }
      else
        this->print_warning_property_name(name);
    }
  }

  void set_property(std::map<std::string, std::map<std::string, std::string>> &value_map) override
  {
    for (const auto &item : value_map)
    {
      if (item.first.compare("freeze") == 0)
      {
        host::vector<bool> _freezed_vertex = util::from_dict_to_vector_types(host::copy(freezed_vertex), item.second);
        freezed_vertex = _freezed_vertex;
      }
      else
        this->print_warning_property_name(item.first);
    }
  }
  int integrate(void);

private:
  double m_dx, m_dy, m_dz;
  unsigned int m_seed; //!< random number seed;
  RNG_ptr m_rng;       //!< Random number generator
  bool m_spherical_move;
  host::vector<bool> freezed_vertex;
};

typedef std::shared_ptr<MonteCarloIntegratorVertex> MonteCarloIntegratorVertex_ptr;

#endif

/** @} */
