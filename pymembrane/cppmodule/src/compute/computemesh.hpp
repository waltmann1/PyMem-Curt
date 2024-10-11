#ifndef __computemesh_hpp__
#define __computemesh_hpp__

#include <memory>
#include <map>

#include "../types/globaltypes.hpp"
#include "../types/hostvector.hpp"

class SystemClass;  //forward declaration
class EvolverClass;

class ComputeMesh
{
public:
    ComputeMesh(SystemClass& system) : _system(system)
    {
    }
    void compute_vertex_normals(bool vertex_normal_angle_weight = false);
    void compute_face_normals(void);
    host::vector<real> compute_edge_lengths(void);
    host::vector<real> gaussiancurvature(void);
    host::vector<real> meancurvature(void);
    host::vector<host::vector<real>> compute_face_metric(void);
    real compute_mesh_volume(void);
    host::vector<real> compute_face_area(void);
    real compute_mesh_area(void);
    host::vector<real> compute_vertex_area(void);

    host::vector<real3> get_vertex_normals(void);

    /// Energy and force measures
    std::map<std::string, real> compute_mesh_energy(EvolverClass &);
    host::vector<real3> compute_vertex_forces(EvolverClass &);
    realTensor compute_stresses(EvolverClass &, const bool &);
    realTensor get_stresses(EvolverClass &, const bool &);
    std::vector<real> compute_pressure(EvolverClass &);
    realTensor compute_kinetic_energy_tensor(void);
    real compute_kinetic_energy(void);
    real compute_temperature(void);
    realTensor compute_stresses_virial(EvolverClass &, const bool &);
    std::vector<realTensor> compute_stresses_atom(EvolverClass &, const bool &);

private:
    SystemClass& _system;
};


#endif
