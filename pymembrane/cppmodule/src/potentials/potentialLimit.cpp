#include "potentialLimit.hpp"

/*! Iterate over all edges that origiante in the given vertex. If the 
    edge is shorter than \a m_min or longer than \a m_max infinity (actually
    a very large number, \f$ 10^{15} \f$ is returned. Otherwise, this method 
    returns zero.
    
    \f$ E_{limit} = \infty \f$ if \f$ r < l_{min} \f$ or \f$ r > l_{max} \f$ and
    \f$ E_{limit} = 0 \f$ otherwise
    \param v Pointer to the vertex
*/

double ComputeEdgeLimitEnergy_dev(const real3 r1,
                                  const real3 r2,
                                  const double lmax,
                                  const double lmin)
{
    real3 rij;
    vsub(rij, r2, r1);

    double dr = sqrt(vdot(rij, rij));

    double energy = 0.0;
    if (dr > lmax)
        energy = BIG_VERTEX_ENERGY_LIMIT;
    else if (dr < lmin)
        energy = BIG_VERTEX_ENERGY_LIMIT;
    return energy;
}

void ComputeVertexLimitEnergy_kernel(int Numedges,
                                     HE_EdgeProp *edges,
                                     HE_VertexProp *vertices,
                                     double *_lmax,
                                     double *_lmin)
{
    for (int edge_index = 0; edge_index < Numedges; edge_index++)
    {
        int type = edges[edge_index].type;

        int v1 = edges[edge_index].i;
        real3 _r1 = vertices[v1].r;

        int v2 = edges[edge_index].j;
        real3 _r2 = vertices[v2].r;

        edges[edge_index].energy += ComputeEdgeLimitEnergy_dev(_r1, _r2, _lmax[type], _lmin[type]);
    }
}

void ComputeVertexLimitEnergy::compute_energy(void)
{

    ComputeVertexLimitEnergy_kernel(_system.Numedges,
                                    &_system.edges[0],
                                    &_system.vertices[0],
                                    &lmax[0],
                                    &lmin[0]);
}

double ComputeVertexLimitEnergy::compute_edge_energy(int query_edge_index)
{
    int type = _system.edges[query_edge_index].type;

    int v1 = _system.edges[query_edge_index].i;
    real3 _r1 = _system.vertices[v1].r;

    int v2 = _system.edges[query_edge_index].j;
    real3 _r2 = _system.vertices[v2].r;

    double edge_energy = ComputeEdgeLimitEnergy_dev(_r1, _r2, lmax[type], lmin[type]);

    return edge_energy;
}

double ComputeVertexLimitEnergy::compute_vertex_energy(int query_vertex_index)
{

    double energy = 0.0;
    ///< get the triangle that this vertex is part of
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
        energy += ComputeEdgeLimitEnergy_dev(_r1, _r2, lmax[type], lmin[type]);
        // MOVE TO THE NEXT edge
        int he_pair = _system.halfedges[he].pair;
        he = _system.halfedges[he_pair].next;
    } while ((he != first));
    return energy;
}

real ComputeVertexLimitForce(real3 rij,
                             real lmax,
                             real lmin)
{
    real dr = sqrt(vdot(rij, rij));
    real fval = 0.0;
    //if(dr > lmax) fval = ( dr - lmax)/dr;
    if (dr > lmax)
        fval = 1.0 / dr;
    else if (dr < lmin)
        fval = -1.0 / dr;
    return fval;
}

void ComputeVertexLimitForce_kernel(int Numedges,
                                    HE_VertexProp *vertices,
                                    HE_EdgeProp *edges,
                                    real *_lmax,
                                    real *_lmin)
{
    for (int edge_index = 0; edge_index < Numedges; edge_index++)
    {
        int type = edges[edge_index].type;
        int v1 = edges[edge_index].i;
        real3 _r1 = vertices[v1].r;
        int v2 = edges[edge_index].j;
        real3 _r2 = vertices[v2].r;
        real3 _rij;
        vsub(_rij, _r2, _r1);

        real fval = ComputeVertexLimitForce(_rij, _lmax[type], _lmin[type]);

        vertices[v1].forceC.x+= fval * _rij.x;
        vertices[v1].forceC.y+= fval * _rij.y;
        vertices[v1].forceC.z+= fval * _rij.z;
        vertices[v2].forceC.x+= -1.0 * fval * _rij.x;
        vertices[v2].forceC.y+= -1.0 * fval * _rij.y;
        vertices[v2].forceC.z+= -1.0 * fval * _rij.z;
    }
}
void ComputeVertexLimitEnergy::compute(void)
{

    ComputeVertexLimitForce_kernel(_system.Numedges,
                                   &_system.vertices[0],
                                   &_system.edges[0],
                                   &lmax[0],
                                   &lmin[0]);
}
