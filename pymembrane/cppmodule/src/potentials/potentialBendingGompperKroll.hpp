/*!
* @file potentialbending.hpp
* @brief ComputeVertexBendingGKEnergy class
* @author Daniel Matoz, fdamatoz@gmail.com
* @date 20-Sept-2017
*/
#ifndef __POTENTIALBENDINGGOMPPERKROLL_HPP__
#define __POTENTIALBENDINGGOMPPERKROLL_HPP__

/** @addtogroup computeenergy Compute Vertex Bending Energy
 *  @brief ComputeVertexBendingGKEnergy definitions
 *  @{
 */
//host
#include "computeforceclass.hpp"

//
#include "../system/systemclass.hpp"
#include "../mesh/computegeometry.hpp"
/**
 * @class ComputeVertexBendingGKEnergy
 * @brief ComputeVertexBendingGKEnergy Compute the Seung-Nelson Bending Energy/Force
 */
class ComputeVertexBendingGKEnergy : public ComputeForceClass
{
public:
  ComputeVertexBendingGKEnergy(SystemClass &system) : ComputeForceClass(system)
  {
    name = "BendingGK"; //!< potential name
    type = "vertex";    //!< potential type
    this->set_default_properties();
  }

  ~ComputeVertexBendingGKEnergy() {}

  void set_default_properties(void) override
  {
    host::vector<real> _kappaH(NUM_TYPES_ALLOWED, 1.0);
    kappaH = _kappaH;
    host::vector<real> _H0(NUM_TYPES_ALLOWED, 0.0);
    H0 = _H0;
    host::vector<real> _kappaG(NUM_TYPES_ALLOWED, 0.0);
    kappaG = _kappaG;
  }

  using ComputeForceClass::set_property;
  void set_property(std::map<std::string, std::map<std::string, std::string>> &region_map) override
  {
    for (const auto &item : region_map)
    {
      if (item.first.compare("kappaH") == 0)
      {
        host::vector<real> _kappaH= util::from_dict_to_vector_types(host::copy(kappaH), item.second);
        kappaH = _kappaH;
      }
      else if (item.first.compare("H0") == 0)
      {
        host::vector<real> _H0 = util::from_dict_to_vector_types(host::copy(H0), item.second);
        H0 = _H0;
      }
      else if (item.first.compare("kappaG") == 0)
      {
        host::vector<real> _kappaG = util::from_dict_to_vector_types(host::copy(kappaG), item.second);
        kappaG = _kappaG;
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
    value["kappaH"] = util::to_string_vec(kappaH);
    value["H0"] = util::to_string_vec(H0);
    value["kappaG"] = util::to_string_vec(kappaG);
    return value;
  }

  void compute_energy(void);
  void compute(void);
  double compute_edge_energy(int);
  double compute_vertex_energy(int);

private:
  host::vector<double> kappaH; 
  host::vector<double> H0;    
  host::vector<double> kappaG;  //!< bending rigidity
};

#endif

/** @} */