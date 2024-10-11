#include "potentialConstantAreaTriangle.hpp"

double ComputeVertexConstantAreaTriangleEnergy_dev(const double *__restrict__ g_now,
                                                   const double sigma,
                                                   const double target_area)
{
    double triangle_area = host::compute_area_triangle_from_metric(g_now);
    double energy = 0.5 * sigma * (triangle_area - target_area) * (triangle_area - target_area);
    return energy;
}

void ComputeVertexConstantAreaTriangleEnergy_kernel(const int Numfaces,
                                                    HE_FaceProp *faces,
                                                    const HE_VertexProp *vertices,
                                                    const double *__restrict__ _sigma,
                                                    const double *__restrict__ _target_area)
{
    for (int face_index = 0; face_index < Numfaces; face_index++)
    {
        int type = faces[face_index].type;
        int v1 = faces[face_index].v1;
        int v2 = faces[face_index].v2;
        int v3 = faces[face_index].v3;

        // compute
        double g_now[3];
        host::compute_form_factor_triangle(g_now, vertices[v1].r, vertices[v2].r, vertices[v3].r);
        double energy = ComputeVertexConstantAreaTriangleEnergy_dev(g_now, _sigma[type], _target_area[type]);

        ///Add energy to that face
        faces[face_index].energy += energy;
    }
}

void ComputeVertexConstantAreaTriangleEnergy::compute_energy(void)
{

    ComputeVertexConstantAreaTriangleEnergy_kernel(_system.Numfaces,
                                                   &_system.faces[0],
                                                   &_system.vertices[0],
                                                   &sigma[0],
                                                   &target_area[0]);
}

// this is done in the cpu the user is responsible for calling get_host_mesh()
double ComputeVertexConstantAreaTriangleEnergy::compute_edge_energy(int query_edge_index)
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

        // compute
        double g_now[3];
        host::compute_form_factor_triangle(g_now, _system.vertices[v1].r, _system.vertices[v2].r, _system.vertices[v3].r);
        edge_energy += ComputeVertexConstantAreaTriangleEnergy_dev(g_now, sigma[type], target_area[type]);
    }
    return edge_energy;
}

double ComputeVertexConstantAreaTriangleEnergy::compute_vertex_energy(int query_vertex_index)
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

            // compute
            double g_now[3];
            host::compute_form_factor_triangle(g_now, _system.vertices[v1].r, _system.vertices[v2].r, _system.vertices[v3].r);
            energy += ComputeVertexConstantAreaTriangleEnergy_dev(g_now, sigma[type], target_area[type]);
        }
        // MOVE TO THE NEXT FACE
        int he_prev = _system.halfedges[he].prev;
        he = _system.halfedges[he_prev].pair;
    } while ((he != first));
    return energy;
}

forceMatrix ComputeVertexConstantAreaTriangleForce_dev(const real3 r1,
                                                       const real3 r2,
                                                       const real3 r3,
                                                       const double *__restrict__ g_now,
                                                       const double sigma,
                                                       const double target_area)
{

    double triangle_area = host::compute_area_triangle_from_metric(g_now);

    /*----------------------------------------------------------------------------------------------------------------*/
    /*-----------------------------------           FORCE MATRIX        ----------------------------------------------*/
    /*----------------------------------------------------------------------------------------------------------------*/
    ///capture the face type
    real force_factor = sigma * (triangle_area - target_area);

    real3 r12, r13;
    vsub(r12, r2, r1);
    vsub(r13, r3, r1);

    real3 forceM11, forceM12;
    forceM11.x = forceM11.y = forceM11.z = 0.0;
    forceM12.x = forceM12.y = forceM12.z = 0.0;

    //T5
    forceM11.x += (0.25 * force_factor / triangle_area) * (g_now[2] * r12.x - g_now[1] * r13.x);
    forceM12.x += (0.25 * force_factor / triangle_area) * (g_now[0] * r13.x - g_now[1] * r12.x);
    forceM11.y += (0.25 * force_factor / triangle_area) * (g_now[2] * r12.y - g_now[1] * r13.y);
    forceM12.y += (0.25 * force_factor / triangle_area) * (g_now[0] * r13.y - g_now[1] * r12.y);
    forceM11.z += (0.25 * force_factor / triangle_area) * (g_now[2] * r12.z - g_now[1] * r13.z);
    forceM12.z += (0.25 * force_factor / triangle_area) * (g_now[0] * r13.z - g_now[1] * r12.z);

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

void ComputeVertexConstantAreaTriangleForce_kernel(const int Numfaces,
                                                   HE_VertexProp *vertices,
                                                   const HE_FaceProp *faces,
                                                   const double *__restrict__ _sigma,
                                                   const double *__restrict__ _target_area)
{
    for (int face_index = 0; face_index < Numfaces; face_index++)
    {
        int type = faces[face_index].type;
        int v1 = faces[face_index].v1;
        int v2 = faces[face_index].v2;
        int v3 = faces[face_index].v3;

        // compute
        double g_now[3];
        host::compute_form_factor_triangle(g_now, vertices[v1].r, vertices[v2].r, vertices[v3].r);

        forceMatrix fval = ComputeVertexConstantAreaTriangleForce_dev(vertices[v1].r,
                                                                      vertices[v2].r,
                                                                      vertices[v3].r,
                                                                      g_now,
                                                                      _sigma[type],
                                                                      _target_area[type]);

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

void ComputeVertexConstantAreaTriangleEnergy::compute(void)
{

    ComputeVertexConstantAreaTriangleForce_kernel(_system.Numfaces,
                                                  &_system.vertices[0],
                                                  &_system.faces[0],
                                                  &sigma[0],
                                                  &target_area[0]);
}