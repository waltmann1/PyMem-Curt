#ifndef __POTENTIALDIHEDRAL_HPP__
#define __POTENTIALDIHEDRAL_HPP__

//host
#include "computeforceclass.hpp"

//
#include "../system/systemclass.hpp"
#include "../mesh/computegeometry.hpp"
/**
 * @class ComputeVertexDihedralEnergy
 * @brief ComputeVertexDihedralEnergy Compute the Dihedral Harmonics Energy/Force
 */
class ComputeVertexDihedralEnergy : public ComputeForceClass
{
public:
  ComputeVertexDihedralEnergy(SystemClass &system) : ComputeForceClass(system)
  {
    name = "Dihedral"; //!< potential name
    type = "edge";    //!< potential type
    this->set_default_properties();
  }

  ~ComputeVertexDihedralEnergy() {}

  void set_default_properties(void) override
  {
    host::vector<real> _kappa(NUM_TYPES_ALLOWED, 1.0);
    kappa = _kappa;
    host::vector<real> _theta0(NUM_TYPES_ALLOWED, defPI);
    theta0 = _theta0;
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
      else if (item.first.compare("theta0") == 0)
      {
        host::vector<real> _theta0 = util::from_dict_to_vector_types(host::copy(theta0), item.second);
        theta0 = _theta0;
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
    value["theta0"] = util::to_string_vec(theta0);
    return value;
  }

  void compute_energy(void);
  void compute(void);
  double compute_edge_energy(int);
  double compute_vertex_energy(int);
  void compute_stress(void);
  void compute_atomic_stress(void);
private:
  host::vector<real> kappa; //!< bending rigidity
  host::vector<real> theta0; //!< preferred dihedral angle
};

#endif

/** @} */