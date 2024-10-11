#include "potentialCauchyGreen.hpp"

double ComputeVertexCauchyGreenEnergy_dev(const double *__restrict__ g_now,
                                          const double *__restrict__ g_reference,
                                          const double *__restrict__ g_reference_inv,
                                          const double Y,
                                          const double nu,
                                          const double h)
{
    double triangle_reference_area = host::compute_area_triangle_from_metric(g_reference);
    //     _                                         _     _                     _     _      _
    //    |   g_reference_inv[0]  g_reference_inv[1]  |   |   g_now[0]  g_now[1]  |   |  1  0  |
    //F = |                                           | x |                       | - |        |
    //    |_  g_reference_inv[1]  g_reference_inv[2] _|   |_  g_now[1]  g_now[2] _|   |_ 0  1 _|
    double F11 = g_reference_inv[0] * g_now[0] + g_reference_inv[1] * g_now[1] - 1.0;
    double F12 = g_reference_inv[0] * g_now[1] + g_reference_inv[1] * g_now[2];
    double F21 = g_reference_inv[1] * g_now[0] + g_reference_inv[2] * g_now[1];
    double F22 = g_reference_inv[1] * g_now[1] + g_reference_inv[2] * g_now[2] - 1.0;

    ///capture the face type
    double coeff_1 = Y * h * triangle_reference_area / (8.0 * (1.0 + nu));
    double coeff_2 = coeff_1 * (nu / (1.0 - nu));
    double energy = coeff_1 * (F11 * F11 + F12 * F21 + F12 * F21 + F22 * F22) + coeff_2 * (F11 + F22) * (F11 + F22);

    return energy;
}

void ComputeVertexCauchyGreenEnergy_kernel(const int Numfaces,
                                           HE_FaceProp *faces,
                                           const HE_VertexProp *vertices,
                                           const double *__restrict__ _Y,
                                           const double *__restrict__ _nu,
                                           const double *__restrict__ _h)
{
    for (int face_index = 0; face_index < Numfaces; face_index++)
    {
        int type = faces[face_index].type;
        int v1 = faces[face_index].v1;
        int v2 = faces[face_index].v2;
        int v3 = faces[face_index].v3;

        double Ydev = _Y[type];
        double nudev = _nu[type];
        double hdev = _h[type];

        // compute
        double g_now[3];
        host::compute_form_factor_triangle(g_now, vertices[v1].r, vertices[v2].r, vertices[v3].r);
        double energy = ComputeVertexCauchyGreenEnergy_dev(g_now, faces[face_index].g_reference, faces[face_index].g_reference_inv, Ydev, nudev, hdev);

        ///Add energy to that face
        faces[face_index].energy += energy;
    }
}

void ComputeVertexCauchyGreenEnergy::compute_energy(void)
{

    ComputeVertexCauchyGreenEnergy_kernel(_system.Numfaces,
                                          &_system.faces[0], &_system.vertices[0], &Y[0], &nu[0], &h[0]);
}

// this is done in the cpu the user is responsible for calling get_host_mesh()
double ComputeVertexCauchyGreenEnergy::compute_edge_energy(int query_edge_index)
{
    //we need to loop the 2 faces that are connected to the edge_index
    host::vector<int> face_vec{_system.edges[query_edge_index].face_k, _system.edges[query_edge_index].face_l};
    //reset energy
    double edge_energy = 0.0;
    for (auto face_index : face_vec)
    {
        int type = _system.faces[face_index].type;
        int v1 = _system.faces[face_index].v1;
        int v2 = _system.faces[face_index].v2;
        int v3 = _system.faces[face_index].v3;

        double Ydev = Y[type];
        double nudev = nu[type];
        double hdev = h[type];

        // compute
        double g_now[3];
        host::compute_form_factor_triangle(g_now, _system.vertices[v1].r, _system.vertices[v2].r, _system.vertices[v3].r);
        edge_energy += ComputeVertexCauchyGreenEnergy_dev(g_now, _system.faces[face_index].g_reference, _system.faces[face_index].g_reference_inv, Ydev, nudev, hdev);
    }
    return edge_energy;
}

double ComputeVertexCauchyGreenEnergy::compute_vertex_energy(int query_vertex_index)
{
    double energy = 0.0;
    ///< get the triangle that this vertex is part of
    int he = _system.vertices[query_vertex_index]._hedge;
    int first = he;
    //std::cout<< "first " << first << "\n";
    int face_index, he_pair, he_pair_next;
    do
    {
        // DO SOMETHING WITH THAT FACE
        face_index = _system.halfedges[he].face;
        if (_system.faces[face_index].boundary == false) // Remember -1 is the virtual face outside of the mesh
        {
            int type = _system.faces[face_index].type;
            int v1 = _system.faces[face_index].v1;
            int v2 = _system.faces[face_index].v2;
            int v3 = _system.faces[face_index].v3;

            double Ydev = Y[type];
            double nudev = nu[type];
            double hdev = h[type];

            // compute
            double g_now[3];
            host::compute_form_factor_triangle(g_now, _system.vertices[v1].r, _system.vertices[v2].r, _system.vertices[v3].r);
            energy += ComputeVertexCauchyGreenEnergy_dev(g_now, _system.faces[face_index].g_reference, _system.faces[face_index].g_reference_inv, Ydev, nudev, hdev);
        }
        // MOVE TO THE NEXT FACE
        int he_prev = _system.halfedges[he].prev;
        he = _system.halfedges[he_prev].pair;
    } while ((he != first));
    return energy;
}

forceMatrix ComputeVertexCauchyGreenForce_dev(const real3 r1,
                                              const real3 r2,
                                              const real3 r3,
                                              const double *__restrict__ g_now,
                                              const double *__restrict__ g_reference,
                                              const double *__restrict__ g_reference_inv,
                                              const double Y,
                                              const double nu,
                                              const double h)
{

    double triangle_reference_area = host::compute_area_triangle_from_metric(g_reference);
    //     _                                         _     _                     _     _      _
    //    |   g_reference_inv[0]  g_reference_inv[1]  |   |   g_now[0]  g_now[1]  |   |  1  0  |
    //F = |                                           | x |                       | - |        |
    //    |_  g_reference_inv[1]  g_reference_inv[2] _|   |_  g_now[1]  g_now[2] _|   |_ 0  1 _|
    double F11 = g_reference_inv[0] * g_now[0] + g_reference_inv[1] * g_now[1] - 1.0;
    double F12 = g_reference_inv[0] * g_now[1] + g_reference_inv[1] * g_now[2];
    double F21 = g_reference_inv[1] * g_now[0] + g_reference_inv[2] * g_now[1];
    double F22 = g_reference_inv[1] * g_now[1] + g_reference_inv[2] * g_now[2] - 1.0;

    /*----------------------------------------------------------------------------------------------------------------*/
    /*-----------------------------------           FORCE MATRIX        ----------------------------------------------*/
    /*----------------------------------------------------------------------------------------------------------------*/
    double c1 = 0.25 * Y * h * triangle_reference_area / (1.0 - nu * nu);
    double c2 = F11 + nu * F22;
    double c3 = F22 + nu * F11;
    double c4 = (1.0 - nu);
    double c5 = F21;
    double c6 = F12;

    double alpha11 = g_reference_inv[0];
    double alpha12 = g_reference_inv[1];
    double alpha22 = g_reference_inv[2];

    real3 r12, r13;
    vsub(r12, r2, r1);
    vsub(r13, r3, r1);

    real3 forceM11, forceM12;
    forceM11.x = forceM11.y = forceM11.z = 0.0;
    forceM12.x = forceM12.y = forceM12.z = 0.0;

    //T1
    forceM11.x += c2 * (2.0 * alpha11 * r12.x + alpha12 * r13.x);
    forceM12.x += c2 * (alpha12 * r12.x);
    forceM11.y += c2 * (2.0 * alpha11 * r12.y + alpha12 * r13.y);
    forceM12.y += c2 * (alpha12 * r12.y);
    forceM11.z += c2 * (2.0 * alpha11 * r12.z + alpha12 * r13.z);
    forceM12.z += c2 * (alpha12 * r12.z);

    //T2
    forceM11.x += c3 * (alpha12 * r13.x);
    forceM12.x += c3 * (2.0 * alpha22 * r13.x + alpha12 * r12.x);
    forceM11.y += c3 * (alpha12 * r13.y);
    forceM12.y += c3 * (2.0 * alpha22 * r13.y + alpha12 * r12.y);
    forceM11.z += c3 * (alpha12 * r13.z);
    forceM12.z += c3 * (2.0 * alpha22 * r13.z + alpha12 * r12.z);

    //T3
    forceM11.x += c4 * (c5 * alpha11 * r13.x + c6 * 2.0 * alpha12 * r12.x + c6 * alpha22 * r13.x);
    forceM12.x += c4 * (c5 * 2.0 * alpha12 * r13.x + c5 * alpha11 * r12.x + c6 * alpha22 * r12.x);
    forceM11.y += c4 * (c5 * alpha11 * r13.y + c6 * 2.0 * alpha12 * r12.y + c6 * alpha22 * r13.y);
    forceM12.y += c4 * (c5 * 2.0 * alpha12 * r13.y + c5 * alpha11 * r12.y + c6 * alpha22 * r12.y);
    forceM11.z += c4 * (c5 * alpha11 * r13.z + c6 * 2.0 * alpha12 * r12.z + c6 * alpha22 * r13.z);
    forceM12.z += c4 * (c5 * 2.0 * alpha12 * r13.z + c5 * alpha11 * r12.z + c6 * alpha22 * r12.z);

    //T4
    forceM11.x *= c1;
    forceM12.x *= c1;
    forceM11.y *= c1;
    forceM12.y *= c1;
    forceM11.z *= c1;
    forceM12.z *= c1;

    forceM11.x *= -1.0;
    forceM12.x *= -1.0;
    forceM11.y *= -1.0;
    forceM12.y *= -1.0;
    forceM11.z *= -1.0;
    forceM12.z *= -1.0;

    forceMatrix result;

    result.forceM11 = forceM11;
    result.forceM12 = forceM12;

    return result;
}

void ComputeVertexCauchyGreenForce_kernel(const int Numfaces,
                                          HE_VertexProp *vertices,
                                          HE_FaceProp *faces,
                                          const double *__restrict__ _Y,
                                          const double *__restrict__ _nu,
                                          const double *__restrict__ _h)
{
    for (int face_index = 0; face_index < Numfaces; face_index++)
    {
        int type = faces[face_index].type;
        int v1 = faces[face_index].v1;
        int v2 = faces[face_index].v2;
        int v3 = faces[face_index].v3;

        double Ydev = _Y[type];
        double nudev = _nu[type];
        double hdev = _h[type];

        // compute
        double g_now[3];
        host::compute_form_factor_triangle(g_now, vertices[v1].r, vertices[v2].r, vertices[v3].r);

        forceMatrix fval = ComputeVertexCauchyGreenForce_dev(vertices[v1].r, vertices[v2].r, vertices[v3].r, g_now, faces[face_index].g_reference, faces[face_index].g_reference_inv, Ydev, nudev, hdev);

        /*----------------------------------------------------------------------------------------------------------------*/
        /*-----------------------------------           ACTUAL CALCULATION        ----------------------------------------*/
        /*----------------------------------------------------------------------------------------------------------------*/
        //v1
        vertices[v1].forceC.x += -1.0 * (fval.forceM11.x + fval.forceM12.x);
        vertices[v1].forceC.y += -1.0 * (fval.forceM11.y + fval.forceM12.y);
        vertices[v1].forceC.z += -1.0 * (fval.forceM11.z + fval.forceM12.z);

        //v2
        vertices[v2].forceC.x += fval.forceM11.x;
        vertices[v2].forceC.y += fval.forceM11.y;
        vertices[v2].forceC.z += fval.forceM11.z;

        //v3
        vertices[v3].forceC.x += fval.forceM12.x;
        vertices[v3].forceC.y += fval.forceM12.y;
        vertices[v3].forceC.z += fval.forceM12.z;
    }
}

void ComputeVertexCauchyGreenEnergy::compute(void)
{

    ComputeVertexCauchyGreenForce_kernel(_system.Numfaces,
                                         &_system.vertices[0],
                                         &_system.faces[0],
                                         &Y[0],
                                         &nu[0],
                                         &h[0]);
}