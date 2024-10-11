/*!
* @file potentiallimit.hpp
* @brief ComputeVertexLimitEnergy class
* @author Daniel Matoz, fdamatoz@gmail.com
* @date 20-Sept-2017
*/

#ifndef __POTENTIALLIMIT_HPP__
#define __POTENTIALLIMIT_HPP__

/** @addtogroup computeenergy Compute tethering potential Energy
 *  @brief ComputeVertexLimitEnergy definitions
 *  @{
 */
#include <math.h>
//host
#include "computeforceclass.hpp"

//

#include "../system/systemclass.hpp"
#include "../mesh/computegeometry.hpp"

/**
 * @class ComputeVertexLimitEnergy
 * @brief ComputeVertexLimitEnergy the tethering potential, if the edge length is larger than lmax then the energy is equal to inf
 */
 
/**
 * @brief Example class that acts on the edges 
 * 
 */
class ComputeVertexLimitEnergy : public ComputeForceClass
{
public:
  ComputeVertexLimitEnergy(SystemClass &system) : ComputeForceClass(system)
  {
    name = "Limit"; //!< potential name
    type = "edge";  //!< potential type
    this->set_default_properties();
  }
  ~ComputeVertexLimitEnergy() {}

  void set_default_properties(void) override
  {
    host::vector<real> _lmax(NUM_TYPES_ALLOWED, 100.0);
    lmax = _lmax;
    host::vector<real> _lmin(NUM_TYPES_ALLOWED, 0.0);
    lmin = _lmin;
  }

  using ComputeForceClass::set_property;
  void set_property(std::map<std::string, std::map<std::string, std::string>> &region_map) override
  {
    for (const auto &item : region_map)
    {
      if (item.first.compare("lmax") == 0)
      {
        host::vector<real> _lmax = util::from_dict_to_vector_types(host::copy(lmax), item.second);
        lmax = _lmax;
      }
      else if (item.first.compare("lmin") == 0)
      {
        host::vector<real> _lmin = util::from_dict_to_vector_types(host::copy(lmin), item.second);
        lmin = _lmin;
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
    value["lmax"] = util::to_string_vec(lmax);
    value["lmin"] = util::to_string_vec(lmin);
    return value;
  }

  void compute_energy(void);
  void compute(void);
  double compute_edge_energy(int);
  double compute_vertex_energy(int);

private:
  host::vector<real> lmax; //!< maximum edge allowed
  host::vector<real> lmin; //!< maximum edge allowed
};

#endif

/** @} */