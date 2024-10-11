#ifndef __POTENTIALVOLUME_HPP__
#define __POTENTIALVOLUME_HPP__

/** @addtogroup computeenergy Compute tethering potential Energy
 *  @brief ComputeVertexVolumeEnergy definitions
 *  @{
 */
#include <math.h>
//host

#include "computeforceclass.hpp"

//
#include "../system/systemclass.hpp"
#include "../mesh/computegeometry.hpp"

/**
 * @class ComputeVertexVolumeEnergy
 * @brief ComputeVertexVolumeEnergy the tethering potential, if the edge length is larger than lmax then the energy is equal to inf
 */
class ComputeVertexVolumeEnergy : public ComputeForceClass
{
public:
  ComputeVertexVolumeEnergy(SystemClass &system) : ComputeForceClass(system)
  {
    name = "volumeconstraint"; //!< potential name
    type = "face";  //!< potential type
    this->set_defaults_property();
  }
  ~ComputeVertexVolumeEnergy() {}

  void set_defaults_property(void)
  {
      lambda = 0.0;
      vol0 = this->compute_volume();
  }

  void set_property(const std::string& name, std::map<int, double> &region_value_map)
  {
    if (name.compare("lambda") == 0)
    {
        lambda = region_value_map[0];
    }
    else if (name.compare("vol0") == 0)
    {
        vol0 = region_value_map[0];
    }
    else
      this->print_warning_property_name(name);
  }

  std::map<std::string, std::string> get_info(void) 
  {
    std::map<std::string, std::string> value;
    value["name"] = name;
    value["type"] = type;
    value["lambda"] = util::to_string(lambda);
    value["vol0"] = util::to_string(vol0);
    return value;
  }

  void compute_energies(void);
  void compute_vertex_forces(void);
  void compute_vertex_force(int) {}
  double compute_volume(void);


  void compute_energy(void);
  void compute(void);
  double compute_vertex_energy(int);
  //edge_calculation methods
  double compute_edge_energy(int){ return 0.0; }
//  double compute_energy(int){ return 0.0; }

private:
  double lambda, vol0;
  
};

#endif

/** @} */