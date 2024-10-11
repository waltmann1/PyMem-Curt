/*!
* @file potentialharmonic.hpp
* @brief ComputeVertexHarmonicEnergy class
* @author Daniel Matoz, fdamatoz@gmail.com
* @date 20-Sept-2017
*/

#ifndef __POTENTIALHARMONIC_HPP__
#define __POTENTIALHARMONIC_HPP__

/** @addtogroup computeenergy Compute Harmonic Edge Energy
 *  @brief ComputeVertexHarmonicEnergy definitions
 *  @{
 */
//host
#include "computeforceclass.hpp"

//

#include "../system/systemclass.hpp"
#include "../mesh/computegeometry.hpp"

/**
 * @class ComputeVertexHarmonicEnergy
 * @brief ComputeVertexHarmonicEnergy compute Vertex Harmonic (edge) Energy
 */
class ComputeVertexHarmonicEnergy : public ComputeForceClass
{
public:
  ComputeVertexHarmonicEnergy(SystemClass &system) : ComputeForceClass(system)
  {
    name = "Harmonic"; //!< potential name
    type = "edge";     //!< potential type
    this->set_default_properties();
  }

  ~ComputeVertexHarmonicEnergy() {}

  void set_default_properties(void) override
  {
    host::vector<real> _k(NUM_TYPES_ALLOWED, 0.0);
    k = _k;
    host::vector<real> _l0(NUM_TYPES_ALLOWED, 1.0);
    l0 = _l0;
  }

  using ComputeForceClass::set_property;
  void set_property(std::map<std::string, std::map<std::string, std::string>> &region_map) override
  {
    for (const auto &item : region_map)
    {
      if (item.first.compare("k") == 0)
      {
        host::vector<real> _k = util::from_dict_to_vector_types(host::copy(k), item.second);
        k = _k;
      }
      else if (item.first.compare("l0") == 0)
      {
        host::vector<real> _l0 = util::from_dict_to_vector_types(host::copy(l0), item.second);
        l0 = _l0;
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
    value["k"] = util::to_string_vec(k);
    value["l0"] = util::to_string_vec(l0);
    return value;
  }

  void compute_energy(void);
  void compute(void);
  double compute_edge_energy(int);
  double compute_vertex_energy(int);
  void compute_stress(void);
  void compute_atomic_stress(void);

private:
  host::vector<real> k;  //!< spring constant
  host::vector<real> l0; //!< rest length
};

#endif

/** @} */