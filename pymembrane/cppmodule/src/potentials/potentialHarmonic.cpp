#include "potentialHarmonic.hpp"

double ComputeEdgeHarmonicEnergy(const real3 r1,
                                 const real3 r2,
                                 const double k,
                                 const double l0)
{
    real3 rij;
    vsub(rij, r2, r1);

    double dr = sqrt(vdot(rij, rij));

    double energy = 0.5 * k * (dr - l0) * (dr - l0);

    return energy;
}

void ComputeVertexHarmonicEnergy_kernel(const int Numedges,
                                        HE_EdgeProp *edges,
                                        const HE_VertexProp *vertices,
                                        const double *__restrict__ _k,
                                        const double *__restrict__ _l0)
{
    for (int edge_index = 0; edge_index < Numedges; edge_index++)
    {
        int type = edges[edge_index].type;

        int v1 = edges[edge_index].i;
        real3 _r1 = vertices[v1].r;

        int v2 = edges[edge_index].j;
        real3 _r2 = vertices[v2].r;

        double energy = ComputeEdgeHarmonicEnergy(_r1, _r2, _k[type], _l0[type]);

        ///ADD ENERGY TO THAT EDGE
        edges[edge_index].energy += energy;
    }
}

double ComputeVertexHarmonicEnergy::compute_edge_energy(int query_edge_index)
{
    int type = _system.edges[query_edge_index].type;
    int v1 = _system.edges[query_edge_index].i;
    real3 _r1 = _system.vertices[v1].r;
    int v2 = _system.edges[query_edge_index].j;
    real3 _r2 = _system.vertices[v2].r;
    double edge_energy = ComputeEdgeHarmonicEnergy(_r1, _r2, k[type], l0[type]);
    return edge_energy;
}

double ComputeVertexHarmonicEnergy::compute_vertex_energy(int query_vertex_index)
{

    double energy = 0.0;
    ///< get the triangle that this vertex is part of
    int he = _system.vertices[query_vertex_index]._hedge;
    int first = he;
    do
    {
        // DO SOMETHING WITH THAT FACE
        int edge_index = _system.halfedges[he].edge;
        int type = _system.edges[edge_index].type;
        int v1 = _system.edges[edge_index].i;
        real3 _r1 = _system.vertices[v1].r;
        int v2 = _system.edges[edge_index].j;
        real3 _r2 = _system.vertices[v2].r;
        energy += ComputeEdgeHarmonicEnergy(_r1, _r2, k[type], l0[type]);
        // MOVE TO THE NEXT FACE
        int he_pair = _system.halfedges[he].pair;
        he = _system.halfedges[he_pair].next;
    } while ((he != first));
    return energy;
}

void ComputeVertexHarmonicEnergy::compute_energy(void)
{

    ComputeVertexHarmonicEnergy_kernel(_system.Numedges,
                                       &_system.edges[0],
                                       &_system.vertices[0],
                                       &k[0],
                                       &l0[0]);
}

double ComputeVertexHarmonicForce(const real3 rij,
                                  const double k,
                                  const double l0)
{
    double dr = sqrt(vdot(rij, rij));
    double fval = k * (dr - l0) / dr;
    return fval;
}

void ComputeVertexHarmonicForce_kernel(const int Numedges,
                                       HE_VertexProp *vertices,
                                       const HE_EdgeProp *__restrict__ edges,
                                       const double *__restrict__ _k,
                                       const double *__restrict__ _l0)
{
    //int i = blockIdx.x*blockDim.x + threadIdx.x;
    //if(  i < Numvertices )
    for (int edge_index = 0; edge_index < Numedges; edge_index++)
    {

        int v1 = edges[edge_index].i;
        real3 _r1 = vertices[v1].r;

        int v2 = edges[edge_index].j;
        real3 _r2 = vertices[v2].r;

        real3 _rij;
        vsub(_rij, _r2, _r1);

        int type = edges[edge_index].type;

        double fval = ComputeVertexHarmonicForce(_rij, _k[type], _l0[type]);

        vertices[v1].forceC.x += fval * _rij.x;
        vertices[v1].forceC.y += fval * _rij.y;
        vertices[v1].forceC.z += fval * _rij.z;

        vertices[v2].forceC.x += -1.0 * fval * _rij.x;
        vertices[v2].forceC.y += -1.0 * fval * _rij.y;
        vertices[v2].forceC.z += -1.0 * fval * _rij.z;
    }
}

void ComputeVertexHarmonicEnergy::compute(void)
{

    ComputeVertexHarmonicForce_kernel(_system.Numedges,
                                      &_system.vertices[0],
                                      &_system.edges[0],
                                      &k[0],
                                      &l0[0]);
}

void ComputeVertexHarmonicStress_kernel(const int Numedges,
                                        HE_VertexProp *vertices,
                                        const HE_EdgeProp *__restrict__ edges,
                                        const double *__restrict__ _k,
                                        const double *__restrict__ _l0,
                                        realTensor *stress_group_edges)
{
    //int i = blockIdx.x*blockDim.x + threadIdx.x;
    //if(  i < Numvertices )
    for (int edge_index = 0; edge_index < Numedges; edge_index++)
    {

        int v1 = edges[edge_index].i;
        real3 r1 = vertices[v1].r;

        int v2 = edges[edge_index].j;
        real3 r2 = vertices[v2].r;

        real3 r12;
        vsub(r12, r2, r1);

        int type = edges[edge_index].type;

        double fval = ComputeVertexHarmonicForce(r12, _k[type], _l0[type]);

        //This might be wrong so have to be checked
        //double check J. Chem. Phys. 131, 154107 (2009) page 4 Eq. 21
        //Asume that v1 is in the local replica then contruct the r1, r2 based on it
        real3 uw_r2;
        uw_r2 = host::vector_sum(r1, r12);

        real3 F2, F1;
        F1.x = fval * r12.x;
        F1.y = fval * r12.y;
        F1.z = fval * r12.z;

        F2.x = -fval * r12.x;
        F2.y = -fval * r12.y;
        F2.z = -fval * r12.z;

        stress_group_edges[edge_index].xx += r1.x * F1.x + uw_r2.x * F2.x;
        stress_group_edges[edge_index].xy += r1.x * F1.y + uw_r2.x * F2.y;
        stress_group_edges[edge_index].xz += r1.x * F1.z + uw_r2.x * F2.z;

        stress_group_edges[edge_index].yx += r1.y * F1.x + uw_r2.y * F2.x;
        stress_group_edges[edge_index].yy += r1.y * F1.y + uw_r2.y * F2.y;
        stress_group_edges[edge_index].yz += r1.y * F1.z + uw_r2.y * F2.z;

        stress_group_edges[edge_index].zx += r1.z * F1.x + uw_r2.z * F2.x;
        stress_group_edges[edge_index].zy += r1.z * F1.y + uw_r2.z * F2.y;
        stress_group_edges[edge_index].zz += r1.z * F1.z + uw_r2.z * F2.z;
    }
}

void ComputeVertexHarmonicEnergy::compute_stress(void)
{
    ComputeVertexHarmonicStress_kernel(_system.Numedges,
                                       &_system.vertices[0],
                                       &_system.edges[0],
                                       &k[0],
                                       &l0[0],
                                       &_system.stress_group_edges[0]);
}

void ComputeVertexHarmonicStressAtom_kernel(const int Numedges,
                                            HE_VertexProp *vertices,
                                            const HE_EdgeProp *__restrict__ edges,
                                            const double *__restrict__ _k,
                                            const double *__restrict__ _l0,
                                            realTensor *stress_virial_atom)
{
    //int i = blockIdx.x*blockDim.x + threadIdx.x;
    //if(  i < Numvertices )
    for (int edge_index = 0; edge_index < Numedges; edge_index++)
    {

        int v1 = edges[edge_index].i;
        real3 r1 = vertices[v1].r;

        int v2 = edges[edge_index].j;
        real3 r2 = vertices[v2].r;

        real3 r12;
        vsub(r12, r2, r1);

        int type = edges[edge_index].type;

        double fval = ComputeVertexHarmonicForce(r12, _k[type], _l0[type]);

        //This might be wrong so have to be checked
        //double check J. Chem. Phys. 131, 154107 (2009) page 4 Eq. 21
        //Asume that v1 is in the local replica then contruct the r1, r2 based on it
        real3 uw_r2;
        uw_r2 = host::vector_sum(r1, r12);

        real3 F2, F1;
        F1.x = fval * r12.x;
        F1.y = fval * r12.y;
        F1.z = fval * r12.z;

        F2.x = -fval * r12.x;
        F2.y = -fval * r12.y;
        F2.z = -fval * r12.z;

        //virial as we know it with PBC
        stress_virial_atom[v1].xx += 0.5 * r12.x * F1.x;
        stress_virial_atom[v1].xy += 0.5 * r12.x * F1.y;
        stress_virial_atom[v1].xz += 0.5 * r12.x * F1.z;
        stress_virial_atom[v1].yx += 0.5 * r12.y * F1.x;
        stress_virial_atom[v1].yy += 0.5 * r12.y * F1.y;
        stress_virial_atom[v1].yz += 0.5 * r12.y * F1.z;
        stress_virial_atom[v1].zx += 0.5 * r12.z * F1.x;
        stress_virial_atom[v1].zy += 0.5 * r12.z * F1.y;
        stress_virial_atom[v1].zz += 0.5 * r12.z * F1.z;

        stress_virial_atom[v2].xx += -0.5 * r12.x * F2.x;
        stress_virial_atom[v2].xy += -0.5 * r12.x * F2.y;
        stress_virial_atom[v2].xz += -0.5 * r12.x * F2.z;
        stress_virial_atom[v2].yx += -0.5 * r12.y * F2.x;
        stress_virial_atom[v2].yy += -0.5 * r12.y * F2.y;
        stress_virial_atom[v2].yz += -0.5 * r12.y * F2.z;
        stress_virial_atom[v2].zx += -0.5 * r12.z * F2.x;
        stress_virial_atom[v2].zy += -0.5 * r12.z * F2.y;
        stress_virial_atom[v2].zz += -0.5 * r12.z * F2.z;
    }
}

void ComputeVertexHarmonicEnergy::compute_atomic_stress(void)
{
    ComputeVertexHarmonicStressAtom_kernel(_system.Numedges,
                                           &_system.vertices[0],
                                           &_system.edges[0],
                                           &k[0],
                                           &l0[0],
                                           &_system.stress_virial_atom[0]);
}
