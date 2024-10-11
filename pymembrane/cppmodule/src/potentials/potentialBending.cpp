#include "potentialBending.hpp"

/*! @fn double ComputeVertexBendingEnergy::compute_energy(int vertex)
    @brief Use normals that belong to the two traingles that share a vertex
    and calculate bengind energy as
    \f$ E_{bend} = \frac 1 2 k_b \left( \hat{n}_i- \hat{n}_j \right)^2 r \f$, where
    \f$ k_b \f$ is the bending rigidity, \f$ r \f$ is edge length, and \f$ \hat{n}_i, \hat{n}_j \f$ are the normal unit vectors of the triangles meeting at the edge.  
    the method used here is the Seung and Nelson bending energy, Seung, H. S. & Nelson, D. R. Defects in flexible membranes with crystalline order. Phys. Rev. A 38, 1005–1018 (1988). 
    @param vertex id
    @return the bending energy at vertex
*/

double ComputeVertexBendingEnergy_dev(const real3 nk,
                                    const real3 nl,
                                    const real3 nk_ref,
                                    const real3 nl_ref,
                                    const double kappa_tilde)
{
    double norm_nk = sqrt(vdot(nk, nk));
    double norm_nl = sqrt(vdot(nl, nl));

    double nknl = vdot(nk, nl) / (norm_nk * norm_nl);

    double norm_nk_ref = sqrt(vdot(nk_ref, nk_ref));
    double norm_nl_ref = sqrt(vdot(nl_ref, nl_ref));

    double nknl_ref = vdot(nk_ref, nl_ref) / (norm_nk_ref * norm_nl_ref);

    double energy = kappa_tilde * (1.0 - nknl * nknl_ref - sqrt(fabs(1.0 - nknl * nknl)) * sqrt((1.0 - nknl_ref * nknl_ref)));
    //double energy = kappa_tilde*(1.0 - nknl);

    return energy;
}

void ComputeVertexBendingEnergy_kernel(int Numedges,
                                       HE_EdgeProp *edges,
                                       const HE_FaceProp *faces,
                                       const HE_VertexProp *vertices,
                                       const double *_kappa)
{
    for (int edge_index = 0; edge_index < Numedges; edge_index++)
    {

        if (edges[edge_index].boundary == false)
        {
            int type = edges[edge_index].type;

            double kappa_tilde = (2.0 / sqrt(3.0)) * (_kappa[type]); //kappa

            int v0 = edges[edge_index].v0;
            int v1 = edges[edge_index].v1;
            int v2 = edges[edge_index].v2;
            int v3 = edges[edge_index].v3;
            /// vertices
            real3 r0 = vertices[v0].r;
            real3 r1 = vertices[v1].r;
            real3 r2 = vertices[v2].r;
            real3 r3 = vertices[v3].r;

            real3 nk, nl;
            nk = host::compute_normal_triangle(r0, r1, r2);
            nl = host::compute_normal_triangle(r0, r2, r3);

            real3 nk_ref, nl_ref;
            nk_ref = faces[edges[edge_index].face_k].normal_reference;
            nl_ref = faces[edges[edge_index].face_l].normal_reference;

            /*<<------------------------------------------>>/
            double norm_nk = sqrt(vdot(nk, nk));
            double norm_nl = sqrt(vdot(nl, nl));            
        
            double nknl = vdot(nk,nl)/(norm_nk*norm_nl);
        
            double norm_nk_ref = sqrt(vdot(nk_ref, nk_ref));
            double norm_nl_ref = sqrt(vdot(nl_ref, nl_ref));
        
            double nknl_ref = vdot(nk_ref, nl_ref)/(norm_nk_ref*norm_nl_ref);

            //if((1.0 - nknl*nknl)<0.0 || (1.0 - nknl_ref*nknl_ref)<0.0) printf(" norm_nk_ref = %f  norm_nl_ref =%f 1.0-nknl_ref*nknl_ref = %f 1.0-nknl*nknl = %f \n",norm_nk_ref, norm_nl_ref, 1.0-(nknl_ref*nknl_ref 1.0-(nknl*nknl));
            //printf("nk_ref = %f %f %f, nl_ref = %f %f %f norm_nk_ref = %f  norm_nl_ref =%f nknl_ref = %f \n",nk_ref.x,nk_ref.y,nk_ref.z,nl_ref.x,nl_ref.y,nl_ref.z,norm_nk_ref, norm_nl_ref, nknl_ref);
            //<<------------------------------------------>>*/
            double energy = ComputeVertexBendingEnergy_dev(nk, nl, nk_ref, nl_ref, kappa_tilde);

            ///Add energy to that edge
            edges[edge_index].energy += energy;
        }
    }
}

void ComputeVertexBendingEnergy::compute_energy(void)
{

    ComputeVertexBendingEnergy_kernel(_system.Numedges,
                                      &_system.edges[0],
                                      &_system.faces[0],
                                      &_system.vertices[0],
                                      &kappa[0]);
}

// this is done in the cpu the user is responsible for calling get_host_mesh()
double ComputeVertexBendingEnergy::compute_edge_energy(int query_edge_index)
{
    //we need to loop the 4 edges that are connected to the edge_index
    auto edge_index_vec = _system.get_edge_neighbours_host(query_edge_index);
    //reset energy
    double edge_energy = 0.0;
    for (auto edge_index : edge_index_vec)
    {
        if (_system.edges[edge_index].boundary == false)
        {
            int type = _system.edges[edge_index].type;

            double kappa_tilde = (2.0 / sqrt(3.0)) * (kappa[type]); //kappa

            int v0 = _system.edges[edge_index].v0;
            int v1 = _system.edges[edge_index].v1;
            int v2 = _system.edges[edge_index].v2;
            int v3 = _system.edges[edge_index].v3;
            /// vertices
            real3 r0 = _system.vertices[v0].r;
            real3 r1 = _system.vertices[v1].r;
            real3 r2 = _system.vertices[v2].r;
            real3 r3 = _system.vertices[v3].r;

            real3 nk, nl;
            nk = host::compute_normal_triangle(r0, r1, r2);
            nl = host::compute_normal_triangle(r0, r2, r3);

            real3 nk_ref, nl_ref;
            nk_ref = _system.faces[_system.edges[edge_index].face_k].normal_reference;
            nl_ref = _system.faces[_system.edges[edge_index].face_l].normal_reference;
            edge_energy += ComputeVertexBendingEnergy_dev(nk, nl, nk_ref, nl_ref, kappa_tilde);
        }
    }
    return edge_energy;
}

// this is done in the cpu the user is responsible for calling get_host_mesh()
double ComputeVertexBendingEnergy::compute_vertex_energy(int query_vertex_index)
{
    double energy = 0.0;
    int he = _system.vertices[query_vertex_index]._hedge;
    int first = he;
    int he_vec[2];
    do
    {
        he_vec[0] = he;
        he_vec[1] = _system.halfedges[he].next;
        for (auto he_index : he_vec)
        {
            int edge_index = _system.halfedges[he_index].edge;
            if (_system.edges[edge_index].boundary == false)
            {
                int type = _system.edges[edge_index].type;

                double kappa_tilde = (2.0 / sqrt(3.0)) * (kappa[type]); //kappa

                int v0 = _system.edges[edge_index].v0;
                int v1 = _system.edges[edge_index].v1;
                int v2 = _system.edges[edge_index].v2;
                int v3 = _system.edges[edge_index].v3;
                /// vertices
                real3 r0 = _system.vertices[v0].r;
                real3 r1 = _system.vertices[v1].r;
                real3 r2 = _system.vertices[v2].r;
                real3 r3 = _system.vertices[v3].r;

                real3 nk, nl;
                nk = host::compute_normal_triangle(r0, r1, r2);
                nl = host::compute_normal_triangle(r0, r2, r3);

                real3 nk_ref, nl_ref;
                nk_ref = _system.faces[_system.edges[edge_index].face_k].normal_reference;
                nl_ref = _system.faces[_system.edges[edge_index].face_l].normal_reference;
                energy += ComputeVertexBendingEnergy_dev(nk, nl, nk_ref, nl_ref, kappa_tilde);
            }
        }
        int he_prev = _system.halfedges[he].prev;
        he = _system.halfedges[he_prev].pair;
    } while (he != first);
    return energy;
}

/*! @fn double ComputeVertexBendingEnergy::compute_vertex_force(int vertex) 
    @brief Use normals that belong to the two traingles that share a vertex
    and calculate bengind force as
    \f$ \vec{F} = -\nabla E_{bend} = \frac 1 2 k_b \left( \hat{n}_i- \hat{n}_j \right)^2 r \f$, where
    \f$ k_b \f$ is the bending rigidity, \f$ r \f$ is edge length, and \f$ \hat{n}_i, \hat{n}_j \f$ are the normal unit vectors of the triangles meeting at the edge. 
    the method used here is the Seung and Nelson bending energy, Seung, H. S. & Nelson, D. R. Defects in flexible membranes with crystalline order. Phys. Rev. A 38, 1005–1018 (1988). 
    @param vertex id
    @return the bending energy at vertex
*/

forceMatrix ComputeVertexBendingForce_dev(const real3 r0,
                                          const real3 r1,
                                          const real3 r2,
                                          const real3 r3,
                                          const double kappa_tilde)
{
    real3 r01, r02, r03;
    real3 nk, nl;
    double s, Ak, Al, AkAl;

    double r01_dot_r02;
    double r01_dot_r03;
    double r02_dot_r03;
    double r01_dot_r01;
    double r02_dot_r02;
    double r03_dot_r03;

    real3 forceM11, forceM12, forceM13;

    nk = host::compute_normal_triangle(r0, r1, r2);
    Ak = 0.5 * sqrt(vdot(nk, nk));
    nl = host::compute_normal_triangle(r0, r2, r3);
    Al = 0.5 * sqrt(vdot(nl, nl));
    s = vdot(nk, nl);
    AkAl = Ak * Al;

    vsub(r01, r1, r0);
    vsub(r02, r2, r0);
    vsub(r03, r3, r0);
    r01_dot_r02 = vdot(r01, r02);
    r01_dot_r03 = vdot(r01, r03);
    r02_dot_r03 = vdot(r02, r03);
    r01_dot_r01 = vdot(r01, r01);
    r02_dot_r02 = vdot(r02, r02);
    r03_dot_r03 = vdot(r03, r03);

    forceM11.x = forceM11.y = forceM11.z = 0.0;
    forceM12.x = forceM12.y = forceM12.z = 0.0;
    forceM13.x = forceM13.y = forceM13.z = 0.0;

    //(r02 · r03) r02 − (r02 · r02) r03
    //(r02 · r03) r01 + (r01 · r02) r03 − 2 (r01 · r03) r02
    //(r01 · r02) r02 − (r02 · r02) r01
    forceM11.x += (r02_dot_r03)*r02.x - (r02_dot_r02)*r03.x;
    forceM11.y += (r02_dot_r03)*r02.y - (r02_dot_r02)*r03.y;
    forceM11.z += (r02_dot_r03)*r02.z - (r02_dot_r02)*r03.z;

    forceM12.x += (r02_dot_r03)*r01.x + (r01_dot_r02)*r03.x - 2.0 * (r01_dot_r03)*r02.x;
    forceM12.y += (r02_dot_r03)*r01.y + (r01_dot_r02)*r03.y - 2.0 * (r01_dot_r03)*r02.y;
    forceM12.z += (r02_dot_r03)*r01.z + (r01_dot_r02)*r03.z - 2.0 * (r01_dot_r03)*r02.z;

    forceM13.x += (r01_dot_r02)*r02.x - (r02_dot_r02)*r01.x;
    forceM13.y += (r01_dot_r02)*r02.y - (r02_dot_r02)*r01.y;
    forceM13.z += (r01_dot_r02)*r02.z - (r02_dot_r02)*r01.z;

    //(-s/(4*Ak*Ak))((r02 · r02) r01 − (r01 · r02) r02)
    //(-s/(4*Ak*Ak))(r01 · r01) r02 − (r01 · r02) r01
    //0
    forceM11.x += (-s / (4.0 * Ak * Ak)) * ((r02_dot_r02)*r01.x - (r01_dot_r02)*r02.x);
    forceM11.y += (-s / (4.0 * Ak * Ak)) * ((r02_dot_r02)*r01.y - (r01_dot_r02)*r02.y);
    forceM11.z += (-s / (4 * Ak * Ak)) * ((r02_dot_r02)*r01.z - (r01_dot_r02)*r02.z);

    forceM12.x += (-s / (4.0 * Ak * Ak)) * ((r01_dot_r01)*r02.x - (r01_dot_r02)*r01.x);
    forceM12.y += (-s / (4.0 * Ak * Ak)) * ((r01_dot_r01)*r02.y - (r01_dot_r02)*r01.y);
    forceM12.z += (-s / (4.0 * Ak * Ak)) * ((r01_dot_r01)*r02.z - (r01_dot_r02)*r01.z);

    //0
    //(-s/(4*Al*Al))((r03 · r03) r02 − (r02 · r03) r03)
    //(-s/(4*Al*Al))((r02 · r02) r03 − (r02 · r03) r02)

    forceM12.x += (-s / (4.0 * Al * Al)) * ((r03_dot_r03)*r02.x - (r02_dot_r03)*r03.x);
    forceM12.y += (-s / (4.0 * Al * Al)) * ((r03_dot_r03)*r02.y - (r02_dot_r03)*r03.y);
    forceM12.z += (-s / (4.0 * Al * Al)) * ((r03_dot_r03)*r02.z - (r02_dot_r03)*r03.z);

    forceM13.x += (-s / (4.0 * Al * Al)) * ((r02_dot_r02)*r03.x - (r02_dot_r03)*r02.x);
    forceM13.y += (-s / (4.0 * Al * Al)) * ((r02_dot_r02)*r03.y - (r02_dot_r03)*r02.y);
    forceM13.z += (-s / (4.0 * Al * Al)) * ((r02_dot_r02)*r03.z - (r02_dot_r03)*r02.z);

    ///Final
    double factor = kappa_tilde / (4.0 * AkAl);
    forceM11.x *= factor;
    forceM12.x *= factor;
    forceM13.x *= factor;
    forceM11.y *= factor;
    forceM12.y *= factor;
    forceM13.y *= factor;
    forceM11.z *= factor;
    forceM12.z *= factor;
    forceM13.z *= factor;

    forceMatrix result;

    result.forceM11 = forceM11;
    result.forceM12 = forceM12;
    result.forceM13 = forceM13;

    return result;
}

forceMatrix scale_BendingForceMatrix_dev(const real3 nk,
                                         const real3 nl,
                                         const real3 nk_ref,
                                         const real3 nl_ref,
                                         forceMatrix fval)
{
    double norm_nk = sqrt(vdot(nk, nk));
    double norm_nl = sqrt(vdot(nl, nl));

    double nknl = vdot(nk, nl) / (norm_nk * norm_nl);
    double fac_nknl = 1.0 - nknl * nknl;

    double norm_nk_ref = sqrt(vdot(nk_ref, nk_ref));
    double norm_nl_ref = sqrt(vdot(nl_ref, nl_ref));

    double nknl_ref = vdot(nk_ref, nl_ref) / (norm_nk_ref * norm_nl_ref);
    double fac_nknl_ref = 1.0 - nknl_ref * nknl_ref;

    double factor = 0.0;
    if (fac_nknl > 0.0 && fac_nknl_ref >= 0.0) //save guard for ridges and weird cases
    {
        factor = nknl_ref * (1.0 - sqrt(fac_nknl_ref / fac_nknl));
    }
    /*else
    {
        printf("err fac_nknl = %f  fac_nknl_ref = %f \n",fac_nknl, fac_nknl_ref);
    }*/

    fval.forceM11.x *= factor;
    fval.forceM12.x *= factor;
    fval.forceM13.x *= factor;
    fval.forceM11.y *= factor;
    fval.forceM12.y *= factor;
    fval.forceM13.y *= factor;
    fval.forceM11.z *= factor;
    fval.forceM12.z *= factor;
    fval.forceM13.z *= factor;

    return fval;
}

void ComputeVertexBendingForce_kernel(int Numedges,
                                      const HE_EdgeProp *edges,
                                      const HE_FaceProp *faces,
                                      HE_VertexProp *vertices,
                                      const double *_kappa)
{
    for (int edge_index = 0; edge_index < Numedges; edge_index++)
    {

        if (edges[edge_index].boundary == false)
        {
            int type = edges[edge_index].type;

            double kappa_tilde = (2.0 / sqrt(3.0)) * (_kappa[type]); //kappa

            int v0 = edges[edge_index].v0;
            int v1 = edges[edge_index].v1;
            int v2 = edges[edge_index].v2;
            int v3 = edges[edge_index].v3;
            /// vertices
            real3 r0 = vertices[v0].r;
            real3 r1 = vertices[v1].r;
            real3 r2 = vertices[v2].r;
            real3 r3 = vertices[v3].r;

            forceMatrix gradEe = ComputeVertexBendingForce_dev(r0, r1, r2, r3, kappa_tilde);

            /* Scale matrix of force by reference configuration*/
            real3 nk, nl;
            nk = host::compute_normal_triangle(r0, r1, r2);
            nl = host::compute_normal_triangle(r0, r2, r3);

            real3 nk_ref, nl_ref;
            nk_ref = faces[edges[edge_index].face_k].normal_reference;
            nl_ref = faces[edges[edge_index].face_l].normal_reference;
            // std::cout<<"nk_ref:"<<nk_ref.x<<","<<nk_ref.y<<","<<nk_ref.z<<"\n"<<" nl_ref:"<<nl_ref.x<<","<<nl_ref.y<<","<<nl_ref.z<<std::endl;
            //forceMatrix fval = gradEe;//scale_BendingForceMatrix(nk ,nl, nk_ref, nl_ref, gradEe);
            forceMatrix fval = scale_BendingForceMatrix_dev(nk, nl, nk_ref, nl_ref, gradEe);

            //v0
            vertices[v0].forceC.x += -fval.forceM11.x - fval.forceM12.x - fval.forceM13.x;
            vertices[v0].forceC.y += -fval.forceM11.y - fval.forceM12.y - fval.forceM13.y;
            vertices[v0].forceC.z += -fval.forceM11.z - fval.forceM12.z - fval.forceM13.z;

            //v1
            vertices[v1].forceC.x += fval.forceM11.x;
            vertices[v1].forceC.y += fval.forceM11.y;
            vertices[v1].forceC.z += fval.forceM11.z;

            //v2
            vertices[v2].forceC.x += fval.forceM12.x;
            vertices[v2].forceC.y += fval.forceM12.y;
            vertices[v2].forceC.z += fval.forceM12.z;

            //v3
            vertices[v3].forceC.x += fval.forceM13.x;
            vertices[v3].forceC.y += fval.forceM13.y;
            vertices[v3].forceC.z += fval.forceM13.z;
        }
    }
}

void ComputeVertexBendingEnergy::compute(void)
{

    ComputeVertexBendingForce_kernel(_system.Numedges,
                                     &_system.edges[0],
                                     &_system.faces[0],
                                     &_system.vertices[0],
                                     &kappa[0]);
}