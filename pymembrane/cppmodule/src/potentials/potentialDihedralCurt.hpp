#ifndef __POTENTIALDIHEDRALCURT_HPP__
#define __POTENTIALDIHEDRALCURT_HPP__

//host
#include "computeforceclass.hpp"

//
#include "../system/systemclass.hpp"
#include "../mesh/computegeometry.hpp"
/**
 * @class ComputeVertexDihedralCurtEnergy
 * @brief ComputeVertexDihedralCurtEnergy Compute the Dihedral Harmonics Energy/Force the way Curt wants to
 */
class ComputeVertexDihedralCurtEnergy : public ComputeForceClass
{
public:
  ComputeVertexDihedralCurtEnergy(SystemClass &system) : ComputeForceClass(system)
  {
    name = "DihedralCurt"; //!< potential name
    type = "vertex";    //!< potential type
    this->set_default_properties();
  }

  ~ComputeVertexDihedralCurtEnergy() {}

  void set_default_properties(void) override
  {
    host::vector<real> _kappa(49, 0.0);
    kappa = _kappa;
    host::vector<real> _theta0(49, 0.0);
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
  void update_vertex_normal(int);
//  void compute_stress(void);
//  void compute_atomic_stress(void);
private:
  host::vector<real> kappa; //!< bending rigidity
  host::vector<real> theta0; //!< preferred dihedral angle
};

#endif

/** @} */