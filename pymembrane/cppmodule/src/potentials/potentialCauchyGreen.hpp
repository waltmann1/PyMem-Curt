/*!
* @file potentialCauchyGreen.hpp
* @brief ComputeVertexCauchyGreenEnergy class
* @author Daniel Matoz, fdamatoz@gmail.com
* @date 20-Sept-2017
*/

#ifndef __POTENTIALCAUCHYGREEN_HPP__
#define __POTENTIALCAUCHYGREEN_HPP__

/** @addtogroup computeenergy Compute Vertex Cauchy-Green Streaching Energy
 *  @brief ComputeVertexCauchyGreenEnergy definitions
 *  @{
 */

//host
#include "computeforceclass.hpp"

//

#include "../system/systemclass.hpp"
#include "../mesh/computegeometry.hpp"

/**
 * @class ComputeVertexCauchyGreenEnergy
 * @brief ComputeVertexCauchyGreenEnergy Compute Vertex Cauchy-Green Streaching Energy
 */
class ComputeVertexCauchyGreenEnergy : public ComputeForceClass
{
public:
  ComputeVertexCauchyGreenEnergy(SystemClass &system) : ComputeForceClass(system)
  {
    name = "Cauchy-Green"; //!< potential name
    type = "face";        //!< potential type
    this->set_default_properties();
  }

  ~ComputeVertexCauchyGreenEnergy() {}

  void set_default_properties(void) override
  {
    host::vector<real> _E(NUM_TYPES_ALLOWED, 1.0);
    Y = _E;
    host::vector<real> _nu(NUM_TYPES_ALLOWED, (1.0 / 3.0));
    nu = _nu;
    host::vector<real> _h(NUM_TYPES_ALLOWED, 1.0);
    h = _h;
  }

  using ComputeForceClass::set_property;
  void set_property(std::map<std::string, std::map<std::string, std::string>> &region_map) override
  {
    for (const auto &item : region_map)
    {
      if (item.first.compare("E") == 0)
      {
        host::vector<real> _E = util::from_dict_to_vector_types(host::copy(Y), item.second);
        Y = _E;
      }
      else if (item.first.compare("nu") == 0)
      {
        host::vector<real> _nu = util::from_dict_to_vector_types(host::copy(nu), item.second);
        nu = _nu;
      }
      else if (item.first.compare("h") == 0)
      {
        host::vector<real> _h = util::from_dict_to_vector_types(host::copy(h), item.second);
        h = _h;
      }
      else
        this->print_warning_property_name(item.first);
    }
  }
  std::map<std::string, std::string> get_info(void)
  {
    std::map<std::string, std::string> value;
    value["name"] = name;
    value["type"] = type;
    value["E"] = util::to_string_vec(Y);
    value["nu"] = util::to_string_vec(nu);
    value["h"] = util::to_string_vec(h);
    return value;
  }

  void compute_energy(void);
  void compute(void);
  double compute_edge_energy(int);
  double compute_vertex_energy(int);

private:
  host::vector<real> Y;  //!< Young's modulus
  host::vector<real> nu; //!< Poison's ratio
  host::vector<real> h;  //!< Thickness of the film
};

#endif

/** @} */