#include "potentialLineTension.hpp"

void ComputeVertexLineTension::compute_energy(void)
{
    for(int vertex_index=0;vertex_index<_system.Numvertices;vertex_index++)
    {
        _system.vertices[vertex_index].energy=this->compute_vertex_energy(vertex_index);
    }
}

double ComputeVertexLineTension::compute_edge_energy(int query_edge_index)
{
    double edge_energy = this->compute_vertex_energy(_system.edges[query_edge_index].i);
    edge_energy+=this->compute_vertex_energy(_system.edges[query_edge_index].j);
    return edge_energy;
}

double ComputeVertexLineTension::compute_vertex_energy(int query_vertex_index)
{
    double energy = 0.0;
    ///< get the triangle that this vertex is part of
    int he = _system.vertices[query_vertex_index]._hedge;
    int first = he;
    int type_0 = _system.vertices[query_vertex_index].type;
    real3 r0 = _system.vertices[query_vertex_index].r;

    double scale_gamma, gamma_type;
    do
    {
        int v1 = _system.halfedges[he].vert_to;
        int type_1 = _system.vertices[v1].type;
//        if(type_0!=type_1)
//        {
//            if(flag_avg==true)
//                gamma_type = 0.5*(gamma[type_0]+gamma[type_1]);
//            else
//                gamma_type = gamma[type_0];
//            if(flag_scale==true)
//            {
//                real3 rij;
//                real3 r1 = _system.vertices[query_vertex_index].r;
//                vsub(rij, r1, r0);
//                scale_gamma = sqrt(vdot(rij, rij));
//            }
//            else
//                scale_gamma = 1.0;
//            energy+= scale_gamma*gamma_type;
//        }
        energy += gamma[get_type(type_0, type_1)];

        // MOVE TO THE NEXT edge
        int he_pair = _system.halfedges[he].pair;
        he = _system.halfedges[he_pair].next;
    }while((he!=first));
    return energy;
}