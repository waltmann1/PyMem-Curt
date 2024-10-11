/*!
* @file potentiallimit.hpp
* @brief ComputeVertexSubstrateEnergy class
* @author Daniel Matoz, fdamatoz@gmail.com
* @date 20-Sept-2017
*/

#ifndef __potentialsubstrate_HPP__
#define __potentialsubstrate_HPP__

/** @addtogroup computeenergy Compute tethering potential Energy
 *  @brief ComputeVertexSubstrateEnergy definitions
 *  @{
 */
#include <math.h>
//host
#include "computeforceclass.hpp"

//

#include "../system/systemclass.hpp"
#include "../mesh/computegeometry.hpp"

/**
 * @class ComputeVertexSubstrateEnergy
 * @brief ComputeVertexSubstrateEnergy the tethering potential, if the edge length is larger than lmax then the energy is equal to inf
 */
 
/**
 * @brief Example class that acts on the edges 
 * 
 */
class ComputeVertexSubstrateEnergy : public ComputeForceClass
{
public:
  ComputeVertexSubstrateEnergy(SystemClass &system) : ComputeForceClass(system)
  {
    name = "Substrate-harmonic"; //!< potential name
    type = "vertex";  //!< potential type
    this->set_default_properties();
  }
  ~ComputeVertexSubstrateEnergy() {}

  void set_default_properties(void) override
  {
    host::vector<real> _k(NUM_TYPES_ALLOWED, 0.0);
    kz1 = _k;
    kz2 = _k;
  }

  using ComputeForceClass::set_property;
  void set_property(std::map<std::string, std::map<std::string, std::string>> &region_map) override
  {
    for (const auto &item : region_map)
    {
      if (item.first.compare("kz+") == 0)
      {
        host::vector<real> _k = util::from_dict_to_vector_types(host::copy(kz1), item.second);
        kz1 = _k;
      }
      else if (item.first.compare("kz-") == 0)
      {
        host::vector<real> _k = util::from_dict_to_vector_types(host::copy(kz2), item.second);
        kz2 = _k;
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
    value["kz1"] = util::to_string_vec(kz1);
    value["kz2"] = util::to_string_vec(kz2);
    return value;
  }

  void compute_energy(void){}
  void compute(void);
  double compute_edge_energy(int){return 0.0;}
  double compute_vertex_energy(int){return 0.0;}

private:
  host::vector<real> kz1; //!< maximum edge allowed
  host::vector<real> kz2; //!< maximum edge allowed
};

#endif

/** @} */