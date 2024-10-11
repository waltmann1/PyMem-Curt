/*!
* @file potentialbending.hpp
* @brief ComputeVertexBendingEnergy class
* @author Daniel Matoz, fdamatoz@gmail.com
* @date 20-Sept-2017
*/
#ifndef __POTENTIALBENDING_HPP__
#define __POTENTIALBENDING_HPP__

/** @addtogroup computeenergy Compute Vertex Bending Energy
 *  @brief ComputeVertexBendingEnergy definitions
 *  @{
 */
//host
#include "computeforceclass.hpp"

//
#include "../system/systemclass.hpp"
#include "../mesh/computegeometry.hpp"
/**
 * @class ComputeVertexBendingEnergy
 * @brief ComputeVertexBendingEnergy Compute the Seung-Nelson Bending Energy/Force
 */
class ComputeVertexBendingEnergy : public ComputeForceClass
{
public:
  ComputeVertexBendingEnergy(SystemClass &system) : ComputeForceClass(system)
  {
    name = "Bending"; //!< potential name
    type = "edge";    //!< potential type
    this->set_default_properties();
  }

  ~ComputeVertexBendingEnergy() {}

  void set_default_properties(void) override
  {
    host::vector<real> _kappa(NUM_TYPES_ALLOWED, 1.0);
    kappa = _kappa;
  }

  using ComputeForceClass::set_property;
  void set_property(std::map<std::string, std::map<std::string, std::string>> &region_map) override
  {
    for (const auto &item : region_map)
    {
      if (item.first.compare("kappa") == 0)
      {
        host::vector<real> _kappa = util::from_dict_to_vector_types(host::copy(kappa), item.second);
        kappa = _kappa;
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
    value["kappa"] = util::to_string_vec(kappa);
    return value;
  }

  void compute_energy(void);
  void compute(void);
  double compute_edge_energy(int);
  double compute_vertex_energy(int);

private:
  host::vector<real> kappa; //!< bending rigidity
};

#endif

/** @} */