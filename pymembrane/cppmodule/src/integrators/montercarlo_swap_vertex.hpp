/*!
 * \file montercarlo_vertex.hpp
 * \author Daniel Matoz , fdamatoz@gmail.com
 * \date 29-Dec-2018
 * \brief Declaration of IntegratorBrownianVertex class
 */

#ifndef __montercarlo_swap_vertex_HPP__
#define __montercarlo_swap_vertex_HPP__

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
 * @class MonteCarloIntegratorSwapVertex
 * @brief Integrator Brownian class implements Brownian dynamics for the vertex position. Particle director will not be integrated
 */
class MonteCarloIntegratorSwapVertex : public MonteCarloIntegrator
{
public:
  MonteCarloIntegratorSwapVertex(SystemClass &system, VertexCompute &potentials) : MonteCarloIntegrator(system, potentials)
  {
    name = "vertex swap";
    type = "monte carlo";
    this->set_default_properties();
  }
  ~MonteCarloIntegratorSwapVertex() {}

  void set_default_properties(void)
  {
    host::vector<bool> _freezed_vertex(NUM_TYPES_ALLOWED, false);
    freezed_vertex = _freezed_vertex;
    this->set_temperature(0.0);
    m_seed = 123456; ///default value
    m_rng = std::make_unique<RNG>(m_seed);
    stochastic_tunnelling = false;
    stochastic_gamma = 0.0;
  }
  using MonteCarloIntegrator::set_property;
  void set_property(std::map<std::string, std::string> &value_map)
  {
    for (const auto &item : value_map)
    {
      if (item.first.compare("T") == 0)
      {
        this->set_temperature(util::from_string_double(item.second));
        this->update_temperature_parameters();
      }
      else if (item.first.compare("gamma") == 0)
      {
        stochastic_tunnelling = true;
        stochastic_gamma = util::from_string_double(item.second);
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
      else if (item.first.compare("stochastic_tunnelling") == 0)
      {
        stochastic_tunnelling = util::from_string_bool(item.second);
      }
      else
        print_warning_property_name(item.first);
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
  bool stochastic_tunnelling;
  double stochastic_gamma;
  host::vector<bool> freezed_vertex;
};

typedef std::shared_ptr<MonteCarloIntegratorSwapVertex> MonteCarloIntegratorSwapVertex_ptr;

#endif

/** @} */
