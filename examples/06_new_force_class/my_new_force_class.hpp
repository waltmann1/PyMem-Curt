#ifndef __my_new_force_HPP__
#define __my_new_force_HPP__

#include <math.h>
#include "computeforceclass.hpp"
#include "../system/systemclass.hpp"
#include "../mesh/computegeometry.hpp"


/**
 * @brief Example of force class that acts on the edges 
 */
class MyNewForceClass : public ComputeForceClass
{
public:
  MyNewForceClass(SystemClass &system) : ComputeForceClass(system)
  {
    name = "my new force"; 	//!< potential name
    type = "edge";      	//!< potential type
    this->set_default_properties();
  }
  ~MyNewForceClass() {}
  /* variable initialization */
  void set_default_properties(void) override;
  /* properties definitions to python interface
     see my_new_force_class.cpp*/
  using ComputeForceClass::set_property;
  void set_property(std::map<std::string, std::map<std::string, std::string>> &region_map) override;
  /* output to properties for printing*/
  std::map<std::string, std::string> get_info(void) override;

  /*compute the energy in a given edge */
  real compute_edge_energy(int) override;	
  /*compute the energy in a given vertex */
  real compute_vertex_energy(int) override;	
  /* compute the energy of the whole system */
  void compute_energy(void) override;    	
  /* compute the forces in the whole system */
  void compute(void) override;			

private:
  host::vector<real> my_property; 			//!< user property
};

#endif