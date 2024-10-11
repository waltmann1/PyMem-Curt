/*!
* @file potentialharmonic.hpp
* @brief ComputeVertexLineTension class
* @author Daniel Matoz, fdamatoz@gmail.com
* @date 20-Sept-2017
*/

#ifndef __POTENTIALLINETENSION_HPP__
#define __POTENTIALLINETENSION_HPP__

/** @addtogroup computeenergy Compute Harmonic Edge Energy
 *  @brief ComputeVertexLineTension definitions
 *  @{
 */
//host
#include "computeforceclass.hpp"

//
#include "../system/systemclass.hpp"
#include "../mesh/computegeometry.hpp"

/**
 * @class ComputeVertexLineTension
 * @brief Iterate over all neighbors of a given vertex. If a current neighbor is
 * of a different type add gamma to the line tension energy.
 * Since in principle we can associate different gamma to each vertex 
 * we use a simple algebraic men to determine the line tension between 
 * vertices with different local gamma. This will be useful only if we have a 
 * three or more components system.
 */
class ComputeVertexLineTension : public ComputeForceClass
{
public:
  ComputeVertexLineTension(SystemClass &system) : ComputeForceClass(system)
  {
    name = "Line Tension"; //!< potential name
    type = "vertex";       //!< potential type
    this->set_default_properties();
  }

  ~ComputeVertexLineTension() {}

  void set_default_properties(void)
  {
    std::vector<double> _gamma(49, 0.0);
    gamma = _gamma;
    flag_avg = false;
    flag_scale = false;
  }

  using ComputeForceClass::set_property;
  void set_property(std::map<std::string, std::map<std::string, std::string>> &region_map) override
  {
    for (const auto &item : region_map)
    {
      if (item.first.compare("gamma") == 0)
      {
        host::vector<real> _gamma = util::from_dict_to_vector_types(host::copy(gamma), item.second);
        gamma = _gamma;
      }
      else
        this->print_warning_property_name(item.first);
    }
  }

  void set_property(std::map<std::string, std::string> &region_map) override
  {
    for (const auto &item : region_map)
    {
      if (item.first.compare("avg") == 0)
        flag_avg = util::from_string_bool(item.second);
      else if (item.first.compare("scale") == 0)
        flag_scale = util::from_string_bool(item.second);
      else
        this->print_warning_property_name(name);
    }
  }

    int get_type(int t1, int t2)
    {
      return int(t1 * t1 + t2 * t1 + t2 * t2);
    }


    std::map<std::string, std::string> get_info(void)
  {
    std::map<std::string, std::string> value;
    value["name"] = name;
    value["type"] = type;
    value["gamma"] = util::to_string_vec(gamma);
    value["avg"] = util::to_string(flag_avg);
    value["scale"] = util::to_string(flag_scale);
    return value;
  }

  void compute_energy(void);
  double compute_edge_energy(int);
  double compute_vertex_energy(int);

private:
  host::vector<real> gamma; //!< spring constant
  bool flag_avg;            //!< if is true then use a simple algebraic betwen type is use it
  bool flag_scale;
};

#endif

/** @} */