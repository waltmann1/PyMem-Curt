#include "my_new_force_class.hpp"


void MyNewForceClass::set_default_properties(void)
{
    //!< Default property set to zero
    for (auto &vi : my_property)
        vi = 0.0;
}

void MyNewForceClass::set_property(std::map<std::string, std::map<std::string, std::string>> &region_map)
{
    for (const auto &item : region_map)
    {
        if (item.first.compare("my property") == 0)
            my_property = util::from_dict_to_vector_types(host::copy(my_property), item.second);
        else
            this->print_warning_property_name(item.first);
    }
}

std::map<std::string, std::string> MyNewForceClass::get_info(void)
{
    std::map<std::string, std::string> value;
    value["name"] = name;
    value["type"] = type;
    value["my property"] = util::to_string_vec(my_property);
    return value;
}

real MyNewForceClass::compute_edge_energy(int)
{
    /**
     * USER-CODE Here add the energy definition
     */
    real edge_energy = ... ;
    /**
     * USER-CODE-END
     * 
     */
    return edge_energy;
}

void MyNewForceClass::compute_energy(void)
{
    for (int edge_index = 0; edge_index < Numedges; edge_index++)
        _system.edges[edge_index]+=this->compute_edge_energy(edge_index);
}   


real MyNewForceClass::compute_vertex_energy(int vertex_index)
{
    /**
     * @brief use the halfedge structure to loop over the edges 
     * that emmanates from the vertex_index
     * 
     */
    real energy = 0.0;
    int he = _system.vertices[query_vertex_index]._hedge;
    int first = he;
    do
    {
        int edge_index = _system.halfedges[he].edge;
        int type = _system.edges[edge_index].type;
        int v1 = _system.edges[edge_index].i;
        real3 _r1 = _system.vertices[v1].r;
        int v2 = _system.edges[edge_index].j;
        real3 _r2 = _system.vertices[v2].r;
        energy +=0.5*this->compute_edge_energy(edge_index);
        // MOVE TO THE NEXT edge
        int he_pair = _system.halfedges[he].pair;
        he = _system.halfedges[he_pair].next;
    } while ((he != first));
    return energy;
}


void MyNewForceClass::compute(void)       
{
    for (int edge_index = 0; edge_index < Numedges; edge_index++)
    {
        real3 _rij;
        vsub(_rij, vertices[ edges[edge_index].j].r, vertices[ edges[edge_index].i].r);
        /**
         * USER-CODE Here add the energy definition
         */
        real fval = ... ;
        /**
         * USER-CODE-END
         * 
         */
        vertices[v1].forceC.x+= fval * _rij.x;
        vertices[v1].forceC.y+= fval * _rij.y;
        vertices[v1].forceC.z+= fval * _rij.z;
        vertices[v2].forceC.x+= -1.0 * fval * _rij.x;
        vertices[v2].forceC.y+= -1.0 * fval * _rij.y;
        vertices[v2].forceC.z+= -1.0 * fval * _rij.z;
    }
}