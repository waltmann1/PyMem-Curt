#ifndef __potentialConstantAreaTriangle_HPP__
#define __potentialConstantAreaTriangle_HPP__

/** @addtogroup computeenergy Compute Vertex Cauchy-Green Streaching Energy
 *  @brief ComputeVertexConstantAreaTriangleEnergy definitions
 *  @{
 */

//host
#include "computeforceclass.hpp"

//

#include "../system/systemclass.hpp"
#include "../mesh/computegeometry.hpp"

/**
 * @class ComputeVertexConstantAreaTriangleEnergy
 * @brief ComputeVertexConstantAreaTriangleEnergy Compute Vertex Cauchy-Green Streaching Energy
 */

/*! For a given triangle t compute contribution from the 
    surface tension caused by change in the surface area. 
    
    \f$ V(r) = \frac{1}{2} \sigma \left(A - A_0\right)^2 \f$
    
    where \f$ A \f$ is current area of the triangle, \f$ A_0 \f$  is the native (given) triangle area, and 
    \f$ \sigma \f$ if the surface tension.
*/
class ComputeVertexConstantAreaTriangleEnergy : public ComputeForceClass
{
public:
    ComputeVertexConstantAreaTriangleEnergy(SystemClass &system) : ComputeForceClass(system)
    {
        name = "Constant Area"; //!< potential name
        type = "face";                   //!< potential type
        this->set_default_properties();
    }

    ~ComputeVertexConstantAreaTriangleEnergy() {}

    void set_default_properties(void) override
    {
        host::vector<real> _sigma(NUM_TYPES_ALLOWED, 0.0);
        sigma = _sigma;
        host::vector<real> _target_area(NUM_TYPES_ALLOWED, 1.0);
        target_area = _target_area;
    }

    using ComputeForceClass::set_property;
    void set_property(std::map<std::string, std::map<std::string, std::string>> &region_map) override
    {
        for (const auto &item : region_map)
        {
            if (item.first.compare("sigma") == 0)
            {
                host::vector<real> _sigma = util::from_dict_to_vector_types(host::copy(sigma), item.second);
                sigma = _sigma;
            }
            else if (item.first.compare("target_area") == 0)
            {
                host::vector<real> _target_area = util::from_dict_to_vector_types(host::copy(target_area), item.second);
                target_area = _target_area;
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
        value["sigma"] = util::to_string_vec(sigma);
        value["target_area"] = util::to_string_vec(target_area);

        return value;
    }

    void compute_energy(void);
    void compute(void);
    double compute_edge_energy(int);
    double compute_vertex_energy(int);

private:
    host::vector<real> sigma;       //!< strength
    host::vector<real> target_area; //!< target equilibirum area
};

#endif

/** @} */