#include "potentialDihedralCurt.hpp"

/*! @fn double ComputeEdgeDihedralCurtEnergy::compute_energy(int vertex)
    @brief
    calculate bengind energy use dihedral harmonics as
    \f$ E_{bend} = \frac 1 2 k_b (1 - \cos(\theta - \theta_0)) \f$, where
    \f$ k_b \f$ is the bending rigidity, \f$ \theta \f$ and \f$ \theta_0 \f$ are the angle and preferred dihedral angle
    of the two vertices in the edge
    @param vertex id
    @return the bending energy at vertex
*/
#define Xvec1(v, a, v1)   \
    (v.x = (a)*v1.x),     \
        (v.y = (a)*v1.y), \
        (v.z = (a)*v1.z)
#define Xvec2(v, a, v1, b, v2)       \
    (v.x = (a)*v1.x + (b)*v2.x),     \
        (v.y = (a)*v1.y + (b)*v2.y), \
        (v.z = (a)*v1.z + (b)*v2.z)

int get_type(int t1, int t2)
{
    return int(t1 * t1 + t2 * t1 + t2 * t2);
}


real3 average_vertex_normal(int query_vertex_index,
                             HE_EdgeProp *edges,
                             const HE_HalfEdgeProp *halfedges,
                             const HE_FaceProp *faces,
                             const HE_VertexProp *vertices)
{

    real3 nv, nf;
    nv.x = nv.y = nv.z = 0.0;
    int v0=query_vertex_index,v1,v2;
    int he = vertices[query_vertex_index]._hedge;
    int first = he;

    do
    {
        int edge_index = halfedges[he].edge;
        if (halfedges[he].boundary == false)
        {
            v1 = halfedges[he].vert_to;
            int he_next = halfedges[he].next;
            v2 = halfedges[he_next].vert_to;
            nf = host::compute_normal_triangle(vertices[v0].r, vertices[v1].r, vertices[v2].r);
            nv.x += nf.x;
            nv.y += nf.y;
            nv.z += nf.z;
        }
        int he_pair = halfedges[he].pair;
        he = halfedges[he_pair].next;
    } while ((he != first));

    nv = host::unit_vector(nv);
    return nv;
}


double compute_harmonic_dihedral(const real3 nl,
                                 const real3 rl,
                                 const real3 nk,
                                 const real3 rk,
                                 const double kappa,
                                 const double theta0)
{
    //double nksq, nlsq, r02norm;
    //nksq = vdot(nk, nk);
    //nlsq = vdot(nl, nl);
    //r02norm = sqrt(vdot(r02, r02));

    //double r02inv, nk2inv, nl2inv, nknlinv;
    //r02inv = nk2inv = nl2inv = 0.0;
    //if (r02norm > 0.0)
     //   r02inv = 1.0 / r02norm;
    //if (nksq > 0.0)
    //    nk2inv = 1.0 / nksq;
    //if (nlsq > 0.0)
     //   nl2inv = 1.0 / nlsq;
    //nknlinv = sqrt(nk2inv * nl2inv);

    double cos_theta, sin_theta, cos_theta_0, sin_theta_0, theta, lengthr,  lengthpoint;
    int sign = 1;

    real3 pointk, pointl;

    pointk.x = pointk.y = pointk.z = 0.0;
    pointk.x = rk.x + nk.x;
    pointk.y = rk.y + nk.y;
    pointk.z = rk.z + nk.z;

    pointl.x = pointl.y = pointl.z = 0.0;
    pointl.x = rl.x + nl.x;
    pointl.y = rl.y + nl.y;
    pointl.z = rl.z + nl.z;

    lengthr = sqrt( (rk.x - rl.x) *  (rk.x - rl.x) +   (rk.y - rl.y) *  (rk.y - rl.y)  +   (rk.z - rl.z) *  (rk.z - rl.z));

    lengthpoint = sqrt( (pointk.x - pointl.x) *  (pointk.x - pointl.x) +   (pointk.y - pointl.y) *  (pointk.y - pointl.y)
            +   (pointk.z - pointl.z) *  (pointk.z - pointl.z));

    if(lengthr > lengthpoint)
    {
        sign = -1;
    }

    cos_theta = vdot(nk, nl);

    if (cos_theta >= 1.0)
        cos_theta = 0.999999;
    if (cos_theta <= -1.0)
        cos_theta = -0.999999;

    theta = sign * acos(cos_theta);
    sin_theta = sin(theta);

    sin_theta_0 = sin(theta0);
    cos_theta_0 = cos(theta0);

    double energy = 0.5 * kappa * (1.0 - cos_theta * cos_theta_0 - sin_theta * sin_theta_0);
    //double energy = kappa * (1.0 - cos_theta * cos_theta_0 - sin_theta * sin_theta_0);
            //double diff = theta- theta0;
            //double energy = 0.5 * kappa * diff * diff;

    //py::print("running the potential calculator energy of this vertex is ", energy, "because kappa cos_theta", kappa,
    //cos_theta, "sin_theta cos_theta0 sin_theta0  ", sin_theta, cos_theta_0, sin_theta_0, "but nk nl", nk, nl);

    return energy;

}


void compute_harmonic_dihedral_force(const real3 nl,
                                 const real3 nk,
                                 const int vl,
                                 const int vk, HE_VertexProp *vertices,
                                 const double kappa,
                                 const double theta0)
{
    //double nksq, nlsq, r02norm;
    //nksq = vdot(nk, nk);
    //nlsq = vdot(nl, nl);
    //r02norm = sqrt(vdot(r02, r02));

    //double r02inv, nk2inv, nl2inv, nknlinv;
    //r02inv = nk2inv = nl2inv = 0.0;
    //if (r02norm > 0.0)
    //    r02inv = 1.0 / r02norm;
    //if (nksq > 0.0)
    //    nk2inv = 1.0 / nksq;
    //if (nlsq > 0.0)
     //   nl2inv = 1.0 / nlsq;
    //nknlinv = sqrt(nk2inv * nl2inv);

    double cos_theta, sin_theta, cos_theta_0, sin_theta_0, theta;
    //cos_theta = -vdot(nk, nl) * nknlinv;
    //sin_theta = -r02norm * nknlinv * vdot(nk, r32);

    theta = vdot(nk, nl);
    cos_theta = cos(theta);
    sin_theta = sin(theta);
    sin_theta_0 = sin(theta0);
    cos_theta_0 = cos(theta0);

    if (cos_theta > 1.0)
        cos_theta = 1.0;
    if (cos_theta < -1.0)
        cos_theta = -1.0;

    double df = -kappa * (sin_theta * cos_theta_0 - cos_theta * sin_theta_0);

    //how to distribute this force????!!??

    vertices[vl].forceC.x += 0.5 * df * nl.x;
    vertices[vl].forceC.y += 0.5 * df * nl.y;
    vertices[vl].forceC.z += 0.5 * df * nl.z;

    vertices[vk].forceC.x += 0.5 * df * nk.x;
    vertices[vk].forceC.y += 0.5 * df * nk.y;
    vertices[vk].forceC.z += 0.5 * df * nk.z;

    //double dE0k = vdot(r01, r02) * nk2inv * r02inv - r02norm * nk2inv;
    //double dE0l = -vdot(r32, r02) * nl2inv * r02inv;
    //double dE1k = nk2inv * r02norm;
    //double dE2k = -vdot(r01, r02) * nk2inv * r02inv;
    //double dE2l = vdot(r32, r02) * nl2inv * r02inv - r02norm * nl2inv;
    //double dE3l = nl2inv * r02norm;

    //real3 F0, F1, F2, F3;
    //Xvec2(F0, df * dE0k, nk, df * dE0l, nl);
    //Xvec1(F1, df * dE1k, nk);
    //Xvec2(F2, df * dE2k, nk, df * dE2l, nl);
    //Xvec1(F3, df * dE3l, nl);

    //vertices[v0].forceC.x += F0.x;
    //vertices[v0].forceC.y += F0.y;
    //vertices[v0].forceC.z += F0.z;

    //vertices[v1].forceC.x += F1.x;
    //vertices[v1].forceC.y += F1.y;
    //vertices[v1].forceC.z += F1.z;

    //vertices[v2].forceC.x += F2.x;
    //vertices[v2].forceC.y += F2.y;
    //vertices[v2].forceC.z += F2.z;

    //vertices[v3].forceC.x += F3.x;
    //vertices[v3].forceC.y += F3.y;
    //vertices[v3].forceC.z += F3.z;

}


double ComputeVertexDihedralCurtEnergy_quick(int query_vertex_index,
                                           HE_EdgeProp *edges,
                                           const HE_HalfEdgeProp *halfedges,
                                           const HE_FaceProp *faces,
                                           const HE_VertexProp *vertices,
                                             const double *_kappa,
                                             const double *_theta0)
{
    real3 n0, n1, r0, r1;
    double energy = 0;
    int v1;
    int v0= query_vertex_index;
    int he = vertices[query_vertex_index]._hedge;
    n0 = vertices[v0].normal;
    r0 = vertices[v0].r;
    int first = he;
    do
    {
        int edge_index = halfedges[he].edge;
        int type = edges[edge_index].type;
        if (halfedges[he].boundary == false)
        {
            v1 = halfedges[he].vert_to;

            n1 = vertices[v1].normal;
            r1 = vertices[v1].r;
            //type should be calculated from the types of the two vertices
            type = get_type(vertices[v0].type, vertices[v1].type);
            //type should be calculated from the types of the two vertices
            energy = energy + 0.5 * compute_harmonic_dihedral(n0, r0, n1, r1, _kappa[type], _theta0[type]);
        }
        int he_pair = halfedges[he].pair;
        he = halfedges[he_pair].next;
    } while ((he != first));

    //py::print("found the energy of this vertex in _quick, it is ", energy);

    return energy;

}


double ComputeVertexDihedralCurtEnergy_kernel(int query_vertex_index,
                                        HE_EdgeProp *edges,
                                            const HE_HalfEdgeProp *halfedges,
                                        const HE_FaceProp *faces,
                                        const HE_VertexProp *vertices,
                                        const double *_kappa,
                                        const double *_theta0)
{
    //need to get all of the edges corresponding to the hexagons around the vertices

    real3 n0, n1, r0, r1;
    double energy = 0;
    int v1;
    int v0= query_vertex_index;
    int he = vertices[query_vertex_index]._hedge;
    n0 = average_vertex_normal(v0, edges, halfedges, faces, vertices);
    r0 = vertices[v0].r;
    int first = he;
    do
    {
        int edge_index = halfedges[he].edge;
        int type = edges[edge_index].type;
        if (halfedges[he].boundary == false)
        {
            v1 = halfedges[he].vert_to;

            n1 = average_vertex_normal(v1, edges, halfedges, faces, vertices);
            r1 = vertices[v1].r;
            type = get_type(vertices[v0].type, vertices[v1].type);
            //type should be calculated from the types of the two vertices
            energy = energy + 0.5 * compute_harmonic_dihedral(n0, r0, n1, r1, _kappa[type], _theta0[type]);
        }
        int he_pair = halfedges[he].pair;
        he = halfedges[he_pair].next;
    } while ((he != first));

    return energy;
}

void ComputeVertexDihedralCurtEnergy::compute_energy(void)
{

    for (int vertex_index = 0; vertex_index < _system.Numvertices; vertex_index++)
    {
        double energy=ComputeVertexDihedralCurtEnergy_quick(vertex_index, &_system.edges[0], &_system.halfedges[0],
                &_system.faces[0], &_system.vertices[0], &kappa[0], &theta0[0]);
        _system.vertices[vertex_index].energy += energy;
    }


}

void ComputeVertexDihedralCurtEnergy::update_vertex_normal(int query_vertex_index)
{
    real3 n0, n1;
    int v1;
    int v0= query_vertex_index;
    int he = _system.vertices[v0]._hedge;
    n0 = average_vertex_normal(v0, &_system.edges[0], &_system.halfedges[0], &_system.faces[0], &_system.vertices[0]);
    _system.vertices[v0].normal.x = n0.x;
    _system.vertices[v0].normal.y = n0.y;
    _system.vertices[v0].normal.z = n0.z;
    int first = he;
    do
    {
        int edge_index = _system.halfedges[he].edge;
        int type = _system.edges[edge_index].type;
        if (_system.halfedges[he].boundary == false)
        {
            v1 = _system.halfedges[he].vert_to;

            n1 = average_vertex_normal(v1, &_system.edges[0], &_system.halfedges[0], &_system.faces[0], &_system.vertices[0]);
            _system.vertices[v1].normal.x = n1.x;
            _system.vertices[v1].normal.y = n1.y;
            _system.vertices[v1].normal.z = n1.z;
        }
        int he_pair = _system.halfedges[he].pair;
        he = _system.halfedges[he_pair].next;
    } while ((he != first));

}


// this is done in the cpu the user is responsible for calling get_host_mesh()
double ComputeVertexDihedralCurtEnergy::compute_edge_energy(int query_edge_index)
{

    host::vector<int> v_index_vec(4);
    v_index_vec[0] = _system.edges[query_edge_index].v0;
    v_index_vec[1] = _system.edges[query_edge_index].v1;
    v_index_vec[2] = _system.edges[query_edge_index].v2;
    v_index_vec[3] = _system.edges[query_edge_index].v3;
    double edge_energy = 0.0;
    for (auto v_index : v_index_vec)
    {
        edge_energy+=ComputeVertexDihedralCurtEnergy_quick(v_index, &_system.edges[0], &_system.halfedges[0],
                &_system.faces[0], &_system.vertices[0], &kappa[0], &theta0[0]);
    }

    return(edge_energy);
}

// this is done in the cpu the user is responsible for calling get_host_mesh()
double ComputeVertexDihedralCurtEnergy::compute_vertex_energy(int query_vertex_index)
{
    double energy = ComputeVertexDihedralCurtEnergy_quick(query_vertex_index, &_system.edges[0],&_system.halfedges[0],
            &_system.faces[0], &_system.vertices[0], &kappa[0], &theta0[0]);

    int he = _system.vertices[query_vertex_index]._hedge;
    int first = he;
    do
    {
        energy += ComputeVertexDihedralCurtEnergy_quick(_system.halfedges[he].vert_to, &_system.edges[0],&_system.halfedges[0],
                &_system.faces[0], &_system.vertices[0], &kappa[0], &theta0[0]);
        int he_pair = _system.halfedges[he].pair;
        he = _system.halfedges[he_pair].next;
    } while ((he != first));
    return energy;

}




void ComputeVertexCurtBendingForce_kernel(int query_vertex_index,
                                      const HE_EdgeProp *edges,
                                          const HE_HalfEdgeProp *halfedges,
                                      const HE_FaceProp *faces,
                                      HE_VertexProp *vertices,
                                      const double *_kappa,
                                      const double *_theta0)
{

    real3 n0, n1;
    double energy = 0;
    int v1;
    int v0= query_vertex_index;
    int he = vertices[query_vertex_index]._hedge;
    int first = he;
    do
    {
        int edge_index = halfedges[he].edge;
        int type = edges[edge_index].type;
        if (halfedges[he].boundary == false)
        {
            v1 = halfedges[he].vert_to;
            n0 = average_vertex_normal(v0, edges, halfedges, faces, vertices);
            n1 = average_vertex_normal(v1, edges, halfedges, faces, vertices);
            //type should be calculated from the types of the two vertices
            type = get_type(vertices[v0].type, vertices[v1].type);
            compute_harmonic_dihedral_force(n0, n1, v0, v1, vertices, _kappa[type], _theta0[type]);
        }
        int he_pair = halfedges[he].pair;
        he = halfedges[he_pair].next;
    } while ((he != first));

}

void ComputeVertexDihedralCurtEnergy::compute(void)
{

    for (int vertex_index = 0; vertex_index < _system.Numvertices; vertex_index++)
    {
        ComputeVertexCurtBendingForce_kernel(vertex_index, &_system.edges[0], &_system.halfedges[0],
                                             &_system.faces[0], &_system.vertices[0], &kappa[0], &theta0[0]);
    }

}


//I don't think anything below here is meaningful

//void ComputeVertexBendingStress_kernel(int Numedges,
//                                       const HE_EdgeProp *edges,
//                                       const HE_FaceProp *faces,
//                                       HE_VertexProp *vertices,
//                                       const double *_kappa,
//                                       const double *_theta0,
//                                       realTensor *stress_group_edges)
//{
//    for (int edge_index = 0; edge_index < Numedges; edge_index++)
//    {
//
//        if (edges[edge_index].boundary == false)
//        {
//            int type = edges[edge_index].type;
//
//            int v0 = edges[edge_index].v0;
//            int v1 = edges[edge_index].v1;
//            int v2 = edges[edge_index].v2;
//            int v3 = edges[edge_index].v3;
//            /// vertices
//            real3 r0 = vertices[v0].r;
//            real3 r1 = vertices[v1].r;
//            real3 r2 = vertices[v2].r;
//            real3 r3 = vertices[v3].r;
//
//            real3 r01, r20, r02, r32;
//            vsub(r01, r0, r1);
//            vsub(r02, r0, r2);
//            vsub(r20, r2, r0);
//            vsub(r32, r3, r2);
//
//            real3 nk, nl;
//            vcross(nk, r01, r02);
//            vcross(nl, r32, r02);
//
//            double nksq, nlsq, r02norm;
//            nksq = vdot(nk, nk);
//            nlsq = vdot(nl, nl);
//            r02norm = sqrt(vdot(r02, r02));
//
//            double r02inv, nk2inv, nl2inv, nknlinv;
//            r02inv = nk2inv = nl2inv = 0.0;
//            if (r02norm > 0.0)
//                r02inv = 1.0 / r02norm;
//            if (nksq > 0.0)
//                nk2inv = 1.0 / nksq;
//            if (nlsq > 0.0)
//                nl2inv = 1.0 / nlsq;
//            nknlinv = sqrt(nk2inv * nl2inv);
//
//            double cos_theta, sin_theta, cos_theta_0, sin_theta_0;
//            cos_theta = -vdot(nk, nl) * nknlinv;
//            sin_theta = -r02norm * nknlinv * vdot(nk, r32);
//            sin_theta_0 = sin(_theta0[type]);
//            cos_theta_0 = cos(_theta0[type]);
//
//            if (cos_theta > 1.0)
//                cos_theta = 1.0;
//            if (cos_theta < -1.0)
//                cos_theta = -1.0;
//
//            double df = -_kappa[type] * (sin_theta * cos_theta_0 - cos_theta * sin_theta_0);
//
//            double dE0k = vdot(r01, r02) * nk2inv * r02inv - r02norm * nk2inv;
//            double dE0l = -vdot(r32, r02) * nl2inv * r02inv;
//            double dE1k = nk2inv * r02norm;
//            double dE2k = -vdot(r01, r02) * nk2inv * r02inv;
//            double dE2l = vdot(r32, r02) * nl2inv * r02inv - r02norm * nl2inv;
//            double dE3l = nl2inv * r02norm;
//
//            real3 F0, F1, F2, F3;
//            Xvec2(F0, df * dE0k, nk, df * dE0l, nl);
//            Xvec1(F1, df * dE1k, nk);
//            Xvec2(F2, df * dE2k, nk, df * dE2l, nl);
//            Xvec1(F3, df * dE3l, nl);
//
//            //This might be wrong so have to be checked
//            //double check J. Chem. Phys. 131, 154107 (2009) page 4 Eq. 21
//            //Asume that v0 is in the local replica then contruct the r1, r2, r3 based on it
//            //real3 r01, r02,
//            real3 r03;
//            r01.x *= -1.0;
//            r01.y *= -1.0;
//            r01.z *= -1.0; //r01 = host::vector_subtract(r1, r0, _box);
//            r02.x *= -1.0;
//            r02.y *= -1.0;
//            r02.z *= -1.0; //r02 = host::vector_subtract(r2, r0, _box);
//            r03 = host::vector_subtract(r3, r0);
//            real3 uw_r3, uw_r2, uw_r1 /*,uw_r0,*/;
//            //uw_r0 = r0;
//            uw_r1 = host::vector_sum(r0, r01);
//            uw_r2 = host::vector_sum(r0, r02);
//            uw_r3 = host::vector_sum(r0, r03);
//
//            stress_group_edges[edge_index].xx += r0.x * F0.x + uw_r1.x * F1.x + uw_r2.x * F2.x + uw_r3.x * F3.x;
//            stress_group_edges[edge_index].xy += r0.x * F0.y + uw_r1.x * F1.y + uw_r2.x * F2.y + uw_r3.x * F3.y;
//            stress_group_edges[edge_index].xz += r0.x * F0.z + uw_r1.x * F1.z + uw_r2.x * F2.z + uw_r3.x * F3.z;
//
//            stress_group_edges[edge_index].yx += r0.y * F0.x + uw_r1.y * F1.x + uw_r2.y * F2.x + uw_r3.y * F3.x;
//            stress_group_edges[edge_index].yy += r0.y * F0.y + uw_r1.y * F1.y + uw_r2.y * F2.y + uw_r3.y * F3.y;
//            stress_group_edges[edge_index].yz += r0.y * F0.z + uw_r1.y * F1.z + uw_r2.y * F2.z + uw_r3.y * F3.z;
//
//            stress_group_edges[edge_index].zx += r0.z * F0.x + uw_r1.z * F1.x + uw_r2.z * F2.x + uw_r3.z * F3.x;
//            stress_group_edges[edge_index].zy += r0.z * F0.y + uw_r1.z * F1.y + uw_r2.z * F2.y + uw_r3.z * F3.y;
//            stress_group_edges[edge_index].zz += r0.z * F0.z + uw_r1.z * F1.z + uw_r2.z * F2.z + uw_r3.z * F3.z;
//        }
//    }
//}
//
//void ComputeVertexDihedralCurtEnergy::compute_stress(void)
//{
//    ComputeVertexBendingStress_kernel(_system.Numedges,
//                                      &_system.edges[0],
//                                      &_system.faces[0],
//                                      &_system.vertices[0],
//                                      &kappa[0],
//                                      &theta0[0],
//                                      &_system.stress_group_edges[0]);
//}
//
//void ComputeVertexBendingStressAtom_kernel(int Numedges,
//                                           const HE_EdgeProp *edges,
//                                           const HE_FaceProp *faces,
//                                           HE_VertexProp *vertices,
//                                           const double *_kappa,
//                                           const double *_theta0,
//                                           realTensor *stress_virial_atom)
//{
//    for (int edge_index = 0; edge_index < Numedges; edge_index++)
//    {
//
//        if (edges[edge_index].boundary == false)
//        {
//            int type = edges[edge_index].type;
//
//            int v0 = edges[edge_index].v0;
//            int v1 = edges[edge_index].v1;
//            int v2 = edges[edge_index].v2;
//            int v3 = edges[edge_index].v3;
//            /// vertices
//            real3 r0 = vertices[v0].r;
//            real3 r1 = vertices[v1].r;
//            real3 r2 = vertices[v2].r;
//            real3 r3 = vertices[v3].r;
//
//            real3 r01, r20, r02, r32;
//            vsub(r01, r0, r1);
//            vsub(r02, r0, r2);
//            vsub(r20, r2, r0);
//            vsub(r32, r3, r2);
//
//            real3 nk, nl;
//            vcross(nk, r01, r02);
//            vcross(nl, r32, r02);
//
//            double nksq, nlsq, r02norm;
//            nksq = vdot(nk, nk);
//            nlsq = vdot(nl, nl);
//            r02norm = sqrt(vdot(r02, r02));
//
//            double r02inv, nk2inv, nl2inv, nknlinv;
//            r02inv = nk2inv = nl2inv = 0.0;
//            if (r02norm > 0.0)
//                r02inv = 1.0 / r02norm;
//            if (nksq > 0.0)
//                nk2inv = 1.0 / nksq;
//            if (nlsq > 0.0)
//                nl2inv = 1.0 / nlsq;
//            nknlinv = sqrt(nk2inv * nl2inv);
//
//            double cos_theta, sin_theta, cos_theta_0, sin_theta_0;
//            cos_theta = -vdot(nk, nl) * nknlinv;
//            sin_theta = -r02norm * nknlinv * vdot(nk, r32);
//            sin_theta_0 = sin(_theta0[type]);
//            cos_theta_0 = cos(_theta0[type]);
//
//            if (cos_theta > 1.0)
//                cos_theta = 1.0;
//            if (cos_theta < -1.0)
//                cos_theta = -1.0;
//
//            double df = -_kappa[type] * (sin_theta * cos_theta_0 - cos_theta * sin_theta_0);
//
//            double dE0k = vdot(r01, r02) * nk2inv * r02inv - r02norm * nk2inv;
//            double dE0l = -vdot(r32, r02) * nl2inv * r02inv;
//            double dE1k = nk2inv * r02norm;
//            double dE2k = -vdot(r01, r02) * nk2inv * r02inv;
//            double dE2l = vdot(r32, r02) * nl2inv * r02inv - r02norm * nl2inv;
//            double dE3l = nl2inv * r02norm;
//
//            real3 F0, F1, F2, F3;
//            Xvec2(F0, df * dE0k, nk, df * dE0l, nl);
//            Xvec1(F1, df * dE1k, nk);
//            Xvec2(F2, df * dE2k, nk, df * dE2l, nl);
//            Xvec1(F3, df * dE3l, nl);
//
//            /*
//            FOR NOW USE THE GROUP METHOD FOR AN EDGE AND DEVIDE BY 4 IN EACH ATOM
//            */
//            //This might be wrong so have to be checked
//            //double check J. Chem. Phys. 131, 154107 (2009) page 4 Eq. 21
//            //Asume that v0 is in the local replica then contruct the r1, r2, r3 based on it
//            //real3 r01, r02,
//            real3 r03;
//            r01.x *= -1.0;
//            r01.y *= -1.0;
//            r01.z *= -1.0; //r01 = host::vector_subtract(r1, r0, _box);
//            r02.x *= -1.0;
//            r02.y *= -1.0;
//            r02.z *= -1.0; //r02 = host::vector_subtract(r2, r0, _box);
//            r03 = host::vector_subtract(r3, r0);
//            real3 uw_r3, uw_r2, uw_r1 /*,uw_r0,*/;
//            //uw_r0 = r0;
//            uw_r1 = host::vector_sum(r0, r01);
//            uw_r2 = host::vector_sum(r0, r02);
//            uw_r3 = host::vector_sum(r0, r03);
//
//            realTensor stress_group_edge;
//            stress_group_edge.xx = r0.x * F0.x + uw_r1.x * F1.x + uw_r2.x * F2.x + uw_r3.x * F3.x;
//            stress_group_edge.xy = r0.x * F0.y + uw_r1.x * F1.y + uw_r2.x * F2.y + uw_r3.x * F3.y;
//            stress_group_edge.xz = r0.x * F0.z + uw_r1.x * F1.z + uw_r2.x * F2.z + uw_r3.x * F3.z;
//
//            stress_group_edge.yx = r0.y * F0.x + uw_r1.y * F1.x + uw_r2.y * F2.x + uw_r3.y * F3.x;
//            stress_group_edge.yy = r0.y * F0.y + uw_r1.y * F1.y + uw_r2.y * F2.y + uw_r3.y * F3.y;
//            stress_group_edge.yz = r0.y * F0.z + uw_r1.y * F1.z + uw_r2.y * F2.z + uw_r3.y * F3.z;
//
//            stress_group_edge.zx = r0.z * F0.x + uw_r1.z * F1.x + uw_r2.z * F2.x + uw_r3.z * F3.x;
//            stress_group_edge.zy = r0.z * F0.y + uw_r1.z * F1.y + uw_r2.z * F2.y + uw_r3.z * F3.y;
//            stress_group_edge.zz = r0.z * F0.z + uw_r1.z * F1.z + uw_r2.z * F2.z + uw_r3.z * F3.z;
//
//            int vvec[4] = {v0, v1, v2, v3};
//            for (auto v : vvec)
//            {
//                stress_virial_atom[v].xx += 0.25 * stress_group_edge.xx;
//                stress_virial_atom[v].xy += 0.25 * stress_group_edge.xy;
//                stress_virial_atom[v].xz += 0.25 * stress_group_edge.xz;
//
//                stress_virial_atom[v].yx += 0.25 * stress_group_edge.yx;
//                stress_virial_atom[v].yy += 0.25 * stress_group_edge.yy;
//                stress_virial_atom[v].yz += 0.25 * stress_group_edge.yz;
//
//                stress_virial_atom[v].zx += 0.25 * stress_group_edge.zx;
//                stress_virial_atom[v].zy += 0.25 * stress_group_edge.zy;
//                stress_virial_atom[v].zz += 0.25 * stress_group_edge.zz;
//            }
//        }
//    }
//}
//
//void ComputeVertexDihedralCurtEnergy::compute_atomic_stress(void)
//{
//    ComputeVertexBendingStressAtom_kernel(_system.Numedges,
//                                          &_system.edges[0],
//                                          &_system.faces[0],
//                                          &_system.vertices[0],
//                                          &kappa[0],
//                                          &theta0[0],
//                                          &_system.stress_virial_atom[0]);
//}