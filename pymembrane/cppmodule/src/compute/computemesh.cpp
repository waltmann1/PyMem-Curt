#include "computemesh.hpp"
#include "../system/systemclass.hpp"
#include "../evolver/evolverclass.hpp"
#include "../mesh/computegeometry.hpp"

void compute_edge_cotangents_kernel(const int Numedges,
                                    HE_EdgeProp *edges,
                                    HE_VertexProp *vertices,
                                    const HE_FaceProp *__restrict__ faces,
                                    const bool COMPUTEVERTEXNORMALS)

{
    for (int edge_index = 0;
         edge_index < Numedges;
         edge_index++)
    {
        real cot_alpha = 0.0;
        real cot_beta = 0.0;

        if (edges[edge_index].boundary == false)
        {
            int v0 = edges[edge_index].v0;
            int v1 = edges[edge_index].v1;
            int v2 = edges[edge_index].v2;
            int v3 = edges[edge_index].v3;

            /// vertices
            real3 r0 = vertices[v0].r;
            real3 r1 = vertices[v1].r;
            real3 r2 = vertices[v2].r;
            real3 r3 = vertices[v3].r;

            real3 r10, r12, r30, r32;
            vsub(r10, r0, r1);
            vsub(r12, r2, r1);
            vsub(r30, r0, r3);
            vsub(r32, r2, r3);

            real r10_dot_r12 = vdot(r10, r12);
            real3 vec_r10_cross_r12;
            vcross(vec_r10_cross_r12, r12, r10);
            real r10_cross_r12 = sqrt(vdot(vec_r10_cross_r12, vec_r10_cross_r12));
            //!< opposite cotangent angle of \vec{r0}-\vec{r1} and \vec{r2}-\vec{r1}
            cot_alpha = r10_dot_r12 / r10_cross_r12;

            real r30_dot_r32 = vdot(r30, r32);
            real3 vec_r30_cross_r32;
            vcross(vec_r30_cross_r32, r30, r32);
            real r30_cross_r32 = sqrt(vdot(vec_r30_cross_r32, vec_r30_cross_r32));
            //!< opposite cotangent angle of \vec{r0}-\vec{r3} and \vec{r2}-\vec{r3}
            cot_beta = r30_dot_r32 / r30_cross_r32;
        }
        else
        {
            int v0 = edges[edge_index].v0;
            int v1 = edges[edge_index].v1;
            int v2 = edges[edge_index].v2;
            int v3 = edges[edge_index].v3;

            int face_index = edges[edge_index].face_k;
            if (edges[edge_index].face_l > 0)
                face_index = edges[edge_index].face_l;

            int v;
            if (faces[face_index].v1 == v1)
                v = v1;
            else if (faces[face_index].v2 == v1)
                v = v1;
            else if (faces[face_index].v3 == v1)
                v = v1;
            else if (faces[face_index].v1 == v3)
                v = v3;
            else if (faces[face_index].v2 == v3)
                v = v3;
            else if (faces[face_index].v3 == v3)
                v = v3;

            /// vertices
            real3 r0 = vertices[v0].r;
            real3 r1 = vertices[v].r;
            real3 r2 = vertices[v2].r;

            real3 r10, r12;
            vsub(r10, r0, r1);
            vsub(r12, r2, r1);

            real r10_dot_r12 = vdot(r10, r12);
            real3 vec_r10_cross_r12;
            vcross(vec_r10_cross_r12, r12, r10);
            real r10_cross_r12 = sqrt(vdot(vec_r10_cross_r12, vec_r10_cross_r12));
            //!< opposite cotangent angle of \vec{r0}-\vec{r1} and \vec{r2}-\vec{r1}
            cot_alpha = r10_dot_r12 / r10_cross_r12;
        }
        edges[edge_index]._property.cot_alpha = cot_alpha;
        edges[edge_index]._property.cot_beta = cot_beta;
        //sum the the cotangent to i vertex
        if (COMPUTEVERTEXNORMALS)
        {
            real3 rij;
            vsub(rij, vertices[edges[edge_index].j].r, vertices[edges[edge_index].i].r);
            vertices[edges[edge_index].i].normal.x += (cot_alpha + cot_beta) * rij.x;
            vertices[edges[edge_index].i].normal.y += (cot_alpha + cot_beta) * rij.y;
            vertices[edges[edge_index].i].normal.z += (cot_alpha + cot_beta) * rij.z;
            vertices[edges[edge_index].j].normal.x += -(cot_alpha + cot_beta) * rij.x;
            vertices[edges[edge_index].j].normal.y += -(cot_alpha + cot_beta) * rij.y;
            vertices[edges[edge_index].j].normal.z += -(cot_alpha + cot_beta) * rij.z;
        }
    }
}

void compute_face_normals_kernel(const int Numfaces,
                                 HE_FaceProp *faces,
                                 HE_VertexProp *vertices,
                                 const bool COMPUTEVERTEXNORMALS,
                                 const bool AREAWeighted)
{
    for (int face_index = 0;
         face_index < Numfaces;
         face_index++)
    {
        int v1 = faces[face_index].v1;
        int v2 = faces[face_index].v2;
        int v3 = faces[face_index].v3;

        // compute
        real3 face_normal = host::compute_normal_triangle(vertices[v1].r, vertices[v2].r, vertices[v3].r);
        faces[face_index].normal = face_normal;
        if (COMPUTEVERTEXNORMALS)
        {
            real face_area_factor[3] = {1.0, 1.0, 1.0};
            if (AREAWeighted)
            {
                face_area_factor[0] = host::compute_area_triangle_from_vertex(vertices[v1].r, vertices[v2].r, vertices[v3].r) / 3.0;
                face_area_factor[1] = face_area_factor[0];
                face_area_factor[2] = face_area_factor[0];
            }
            else
            {
                face_area_factor[0] = host::compute_angle_vertex(vertices[v1].r, vertices[v2].r, vertices[v3].r);
                face_area_factor[1] = host::compute_angle_vertex(vertices[v2].r, vertices[v3].r, vertices[v1].r);
                face_area_factor[2] = host::compute_angle_vertex(vertices[v3].r, vertices[v1].r, vertices[v2].r);
            }
            //v1
            vertices[v1].normal.x += face_area_factor[0] * face_normal.x;
            vertices[v1].normal.y += face_area_factor[0] * face_normal.y;
            vertices[v1].normal.z += face_area_factor[0] * face_normal.z;
            //v2
            vertices[v2].normal.x += face_area_factor[1] * face_normal.x;
            vertices[v2].normal.y += face_area_factor[1] * face_normal.y;
            vertices[v2].normal.z += face_area_factor[1] * face_normal.z;
            //v3
            vertices[v3].normal.x += face_area_factor[2] * face_normal.x;
            vertices[v3].normal.y += face_area_factor[2] * face_normal.y;
            vertices[v3].normal.z += face_area_factor[2] * face_normal.z;
        }
    }
}

void ComputeMesh::compute_vertex_normals(bool vertex_normal_angle_weight)
{
    compute_face_normals_kernel(_system.Numfaces,
                                &_system.faces[0],
                                &_system.vertices[0],
                                true,
                                vertex_normal_angle_weight);
}


host::vector<real3> ComputeMesh::get_vertex_normals(void)
{
    host::vector<real3> normals(_system.get_vertices().size());
        for (int v = 0; v < _system.get_vertices().size(); v++)
        {
            normals[v].x = _system.get_vertices()[v].normal.x;
            normals[v].y = _system.get_vertices()[v].normal.y;
            normals[v].z = _system.get_vertices()[v].normal.z;
        }
    return normals;
}


host::vector<real> ComputeMesh::compute_vertex_area(void)
{
    host::vector<real> vertex_area(_system.Numvertices, 0.0);

    for (int face_index = 0;
         face_index < _system.Numfaces;
         face_index++)
    {
        int v1 = _system.faces[face_index].v1;
        int v2 = _system.faces[face_index].v2;
        int v3 = _system.faces[face_index].v3;

        // compute
        auto face_area = host::compute_area_triangle_from_vertex(_system.vertices[v1].r, _system.vertices[v2].r, _system.vertices[v3].r) / 3.0;
        //v1
        vertex_area[v1] += face_area;
        //v2
        vertex_area[v2] += face_area;
        //v3
        vertex_area[v3] += face_area;
    }
    return vertex_area;
}
void ComputeMesh::compute_face_normals(void)
{
    compute_face_normals_kernel(_system.Numfaces,
                                &_system.faces[0],
                                &_system.vertices[0],
                                false,
                                false);
}

void compute_edge_length_kernel(const int Numedges,
                                const HE_EdgeProp *__restrict__ edges,
                                const HE_VertexProp *__restrict__ vertices,
                                real *edge_lengths,
                                const BoxType _box)

{
    for (int edge_index = 0;
         edge_index < Numedges;
         edge_index++)
    {
        //real3 rij = host::minimum_image(vertices[edges[edge_index].i].r, vertices[edges[edge_index].j].r, _box);
        real3 rij;
        vsub(rij, vertices[edges[edge_index].i].r, vertices[edges[edge_index].j].r);
        edge_lengths[edge_index] = sqrt(vdot(rij, rij));
    }
}

host::vector<real> ComputeMesh::compute_edge_lengths(void)
{
    host::vector<real> edge_lengths(_system.Numedges, 0.0);
    compute_edge_length_kernel(_system.Numedges,
                               &_system.edges[0],
                               &_system.vertices[0],
                               &edge_lengths[0],
                               _system.get_box());
    return (host::copy(edge_lengths));
}

host::vector<real> ComputeMesh::gaussiancurvature(void)
{
    std::vector<double> curv;
    for (int vertex_index = 0; vertex_index < _system.Numvertices; vertex_index++)
    {
        double vertex_curv = 0.0;
        if (_system.vertices[vertex_index].boundary == false)
        {
            real3 r0 = _system.vertices[vertex_index].r;
            int first_he = _system.vertices[vertex_index]._hedge;
            int he = first_he;
            double phi = 0.0;
            do
            {
                int he_prev = _system.halfedges[he].prev;
                int v1_to = _system.halfedges[he].vert_to;
                real3 r1 = _system.vertices[v1_to].r;
                int v2_to = _system.halfedges[he_prev].vert_from;
                real3 r2 = _system.vertices[v2_to].r;

                real3 r10, r20;
                vsub(r10, r1, r0);
                vsub(r20, r2, r0);

                double r10_norm = sqrt(vdot(r10, r10));
                double r20_norm = sqrt(vdot(r20, r20));

                aXvec((1.0 / r10_norm), r10);
                aXvec((1.0 / r20_norm), r20);

                r10_norm = sqrt(vdot(r10, r10));
                r20_norm = sqrt(vdot(r20, r20));

                double r10_dot_r20 = vdot(r10, r20);
                if (r10_dot_r20 > 1.0)
                    r10_dot_r20 = 1.0;
                else if (r10_dot_r20 < -1.0)
                    r10_dot_r20 = -1.0;

                phi += acos(r10_dot_r20);
                he = _system.halfedges[he_prev].pair;
            } while (he != first_he);
            vertex_curv = 2.0 * M_PI - phi;
        }
        curv.push_back(vertex_curv);
    }
    return (curv);
}

host::vector<real> ComputeMesh::meancurvature(void)
{
    std::vector<double> curv;
    for (int vertex_index = 0; vertex_index < _system.Numvertices; vertex_index++)
    {
        real3 sigma, nv, nf;
        nv.x = nv.y = nv.z = 0.0;
        sigma.x = sigma.y = sigma.z = 0.0;
        double vertex_area = 0.0;
        int v0 = vertex_index, v1, v2;
        ///< get the triangle that this vertex is part of
        int he = _system.vertices[vertex_index]._hedge;
        int first = he;
        do
        {
            // DO SOMETHING WITH THAT EDGE
            int edge_index = _system.halfedges[he].edge;
            if (_system.edges[edge_index].boundary == false)
            {
                v1 = _system.halfedges[he].vert_to;
                //
                int he_next = _system.halfedges[he].next;
                //
                v2 = _system.halfedges[he_next].vert_to;
                nf = host::compute_normal_triangle(_system.vertices[v0].r, _system.vertices[v1].r, _system.vertices[v2].r);
                nv.x += nf.x;
                nv.y += nf.y;
                nv.z += nf.z;

                real3 r0 = _system.vertices[v0].r;
                real3 r1 = _system.vertices[v1].r;
                real3 r2 = _system.vertices[v2].r;
                real3 r01, r02, r10, r12, r20, r21;
                vsub(r01, r1, r0);
                vsub(r02, r2, r0);
                vsub(r10, r0, r1);
                vsub(r12, r2, r1);
                vsub(r20, r0, r2);
                vsub(r21, r1, r2);
                double r01_dot_r02 = vdot(r01, r02);
                double r10_dot_r12 = vdot(r10, r12);
                double r20_dot_r21 = vdot(r20, r21);
                real3 r01_cross_r02, r10_cross_r12, r20_cross_r21;
                vcross(r01_cross_r02, r01, r02);
                vcross(r10_cross_r12, r10, r12);
                vcross(r20_cross_r21, r21, r20);
                double r01_cross_r02n = sqrt(vdot(r01_cross_r02, r01_cross_r02));
                double r10_cross_r12n = sqrt(vdot(r10_cross_r12, r10_cross_r12));
                double r20_cross_r21n = sqrt(vdot(r20_cross_r21, r20_cross_r21));
                double cot_alpha = r10_dot_r12 / r10_cross_r12n;
                double cot_beta = r20_dot_r21 / r20_cross_r21n;
                double cot_weight = 0.5 * (cot_alpha + cot_beta);
                sigma.x += 0.5 * cot_alpha * r02.x + 0.5 * cot_beta * r01.x;
                sigma.y += 0.5 * cot_alpha * r02.y + 0.5 * cot_beta * r01.y;
                sigma.z += 0.5 * cot_alpha * r02.z + 0.5 * cot_beta * r01.z;
                if (r01_dot_r02 < 0 || r10_dot_r12 < 0 || r20_dot_r21 < 0)
                {
                    if (r01_dot_r02 < 0)
                        vertex_area += 0.5 * r01_cross_r02n;
                    else
                        vertex_area += 0.25 * r01_cross_r02n;
                }
                else
                    vertex_area += 0.125 * vdot(r02, r02) * cot_alpha + 0.125 * vdot(r01, r01) * cot_beta;
            }
            /*------------------------------------------------*/
            // MOVE TO THE NEXT EDGE
            int he_pair = _system.halfedges[he].pair;
            he = _system.halfedges[he_pair].next;
        } while ((he != first));
        double sign = (vdot(nv, sigma) > 0.) ? -1. : 1.;
        double H = sign * sqrt(vdot(sigma, sigma)) / vertex_area;
        // return 0.5*(sqrt(vdot(sigma, sigma)) / vertex_area)*(sqrt(vdot(sigma, sigma)));
        curv.push_back(H);
    }
    return (curv);
}

real ComputeMesh::compute_mesh_volume(void)
{
    real mesh_volume = 0.0;
    if (_system.close_surface == true)
    {
        for (int face_index = 0; face_index < _system.faces.size(); face_index++)
        {
            real3 r0 = _system.vertices[_system.faces[face_index].v1].r;
            real3 r1 = _system.vertices[_system.faces[face_index].v2].r;
            real3 r2 = _system.vertices[_system.faces[face_index].v3].r;
            real3 nt = host::compute_normal_triangle(r0, r1, r2);
            mesh_volume += vdot(r0, nt) / 6.0;
        }
    }
    else
    {
        std::cerr << "open surface doesn't have any volume\n";
    }
    return mesh_volume;
}

host::vector<real> ComputeMesh::compute_face_area(void)
{
    host::vector<real> face_area;
    for (int face_index = 0;
         face_index < _system.Numfaces;
         face_index++)
    {
        int v1 = _system.faces[face_index].v1;
        int v2 = _system.faces[face_index].v2;
        int v3 = _system.faces[face_index].v3;
        face_area.push_back(host::compute_area_triangle_from_vertex(_system.vertices[v1].r, _system.vertices[v2].r, _system.vertices[v3].r));
    }
    return face_area;
}

real ComputeMesh::compute_mesh_area(void)
{
    real mesh_area = 0.0;
    for (int face_index = 0;
         face_index < _system.Numfaces;
         face_index++)
    {
        int v1 = _system.faces[face_index].v1;
        int v2 = _system.faces[face_index].v2;
        int v3 = _system.faces[face_index].v3;
        mesh_area += host::compute_area_triangle_from_vertex(_system.vertices[v1].r, _system.vertices[v2].r, _system.vertices[v3].r);
    }
    return mesh_area;
}

host::vector<host::vector<real>> ComputeMesh::compute_face_metric(void)
{
    host::vector<host::vector<real>> metrics;
    for (int face_index = 0;
         face_index < _system.Numfaces;
         face_index++)
    {
        int v1 = _system.faces[face_index].v1;
        int v2 = _system.faces[face_index].v2;
        int v3 = _system.faces[face_index].v3;
        host::vector<real> g_now(3);
        host::compute_form_factor_triangle(&g_now[0], _system.vertices[v1].r, _system.vertices[v2].r, _system.vertices[v3].r);
        metrics.push_back(g_now);
    }
    return metrics;
}

std::map<std::string, real> ComputeMesh::compute_mesh_energy(EvolverClass &evolver)
{
    evolver.reset_mesh_energy();
    evolver.compute_mesh_energy();

    auto v_energy = 0.0;
    auto e_energy = 0.0;
    auto f_energy = 0.0;

    for(auto v : _system.vertices) v_energy+=v.energy;
    for(auto e : _system.edges) e_energy+=e.energy;
    for(auto f : _system.faces) f_energy+=f.energy;

    std::map<std::string, real> val;
    val["vertices"] = v_energy;
    val["edges"] = e_energy;
    val["faces"] = f_energy;
    return val;
}

host::vector<real3> ComputeMesh::compute_vertex_forces(EvolverClass &evolver)
{
    evolver.reset_mesh_forces();
    evolver.compute_mesh_forces();
    auto vertices = _system.get_vertices();
    host::vector<real3> forces(vertices.size());
    for (auto v : vertices)
        forces[v.id] = v.forceC;
    return forces;
}
realTensor ComputeMesh::compute_stresses(EvolverClass &evolver, const bool &include_velocities)
{
    evolver.reset_mesh_stresses();
    evolver.compute_mesh_stresses();
    return (this->get_stresses(evolver, include_velocities));
}

realTensor ComputeMesh::get_stresses(EvolverClass &evolver, const bool &include_velocities)
{
    realTensor zeroTensor;
    zeroTensor.xx = 0.0;
    zeroTensor.xy = 0.0;
    zeroTensor.xz = 0.0;
    zeroTensor.yx = 0.0;
    zeroTensor.yy = 0.0;
    zeroTensor.yz = 0.0;
    zeroTensor.zx = 0.0;
    zeroTensor.zy = 0.0;
    zeroTensor.zz = 0.0;

    realTensor total_stress = zeroTensor;

    if (evolver.has_vertex_forces == true)
    {
        auto stress_group = host::copy(_system.stress_group_vertices);
        for (auto vertex_index = 0; vertex_index < _system.Numvertices; vertex_index++)
        {
            total_stress.xx += stress_group[vertex_index].xx;
            total_stress.xy += stress_group[vertex_index].xy;
            total_stress.xz += stress_group[vertex_index].xz;
            total_stress.yx += stress_group[vertex_index].yx;
            total_stress.yy += stress_group[vertex_index].yy;
            total_stress.yz += stress_group[vertex_index].yz;
            total_stress.zx += stress_group[vertex_index].zx;
            total_stress.zy += stress_group[vertex_index].zy;
            total_stress.zz += stress_group[vertex_index].zz;
        }
    }
    if (evolver.has_edge_forces == true)
    {
        auto stress_group = host::copy(_system.stress_group_edges);
        for (auto edge_index = 0; edge_index < _system.Numedges; edge_index++)
        {
            total_stress.xx += stress_group[edge_index].xx;
            total_stress.xy += stress_group[edge_index].xy;
            total_stress.xz += stress_group[edge_index].xz;
            total_stress.yx += stress_group[edge_index].yx;
            total_stress.yy += stress_group[edge_index].yy;
            total_stress.yz += stress_group[edge_index].yz;
            total_stress.zx += stress_group[edge_index].zx;
            total_stress.zy += stress_group[edge_index].zy;
            total_stress.zz += stress_group[edge_index].zz;
        }
    }
    if (evolver.has_face_forces == true)
    {
        auto stress_group = host::copy(_system.stress_group_faces);
        for (auto face_index = 0; face_index < _system.Numfaces; face_index++)
        {
            total_stress.xx += stress_group[face_index].xx;
            total_stress.xy += stress_group[face_index].xy;
            total_stress.xz += stress_group[face_index].xz;
            total_stress.yx += stress_group[face_index].yx;
            total_stress.yy += stress_group[face_index].yy;
            total_stress.yz += stress_group[face_index].yz;
            total_stress.zx += stress_group[face_index].zx;
            total_stress.zy += stress_group[face_index].zy;
            total_stress.zz += stress_group[face_index].zz;
        }
    }

    if (include_velocities == true)
    {
        auto kinetic_energy_tensor = this->compute_kinetic_energy_tensor();
        total_stress.xx += kinetic_energy_tensor.xx;
        total_stress.xy += kinetic_energy_tensor.xy;
        total_stress.xz += kinetic_energy_tensor.xz;
        total_stress.yx += kinetic_energy_tensor.yx;
        total_stress.yy += kinetic_energy_tensor.yy;
        total_stress.yz += kinetic_energy_tensor.yz;
        total_stress.zx += kinetic_energy_tensor.zx;
        total_stress.zy += kinetic_energy_tensor.zy;
        total_stress.zz += kinetic_energy_tensor.zz;
    }

    total_stress.xx /= _system.Numvertices;
    total_stress.xy /= _system.Numvertices;
    total_stress.xz /= _system.Numvertices;
    total_stress.yx /= _system.Numvertices;
    total_stress.yy /= _system.Numvertices;
    total_stress.yz /= _system.Numvertices;
    total_stress.zx /= _system.Numvertices;
    total_stress.zy /= _system.Numvertices;
    total_stress.zz /= _system.Numvertices;

    return (total_stress);
}

realTensor ComputeMesh::compute_kinetic_energy_tensor(void)
{

    realTensor zeroTensor;
    zeroTensor.xx = 0.0;
    zeroTensor.xy = 0.0;
    zeroTensor.xz = 0.0;
    zeroTensor.yx = 0.0;
    zeroTensor.yy = 0.0;
    zeroTensor.yz = 0.0;
    zeroTensor.zx = 0.0;
    zeroTensor.zy = 0.0;
    zeroTensor.zz = 0.0;

    auto vertices = _system.get_vertices();
    realTensor kinetic_energy_tensor = zeroTensor;
    for (int vindex = 0;
         vindex < _system.Numvertices;
         vindex++)
    {
        auto vertex = vertices[vindex];
        realTensor kinetic_tensor;
        kinetic_energy_tensor.xx += vertex.mass * vertex.v.x * vertex.v.x;
        kinetic_energy_tensor.xy += vertex.mass * vertex.v.x * vertex.v.y;
        kinetic_energy_tensor.xz += vertex.mass * vertex.v.x * vertex.v.z;
        kinetic_energy_tensor.yx += vertex.mass * vertex.v.y * vertex.v.x;
        kinetic_energy_tensor.yy += vertex.mass * vertex.v.y * vertex.v.y;
        kinetic_energy_tensor.yz += vertex.mass * vertex.v.y * vertex.v.z;
        kinetic_energy_tensor.zx += vertex.mass * vertex.v.z * vertex.v.x;
        kinetic_energy_tensor.zy += vertex.mass * vertex.v.z * vertex.v.y;
        kinetic_energy_tensor.zz += vertex.mass * vertex.v.z * vertex.v.z;
    }
    return kinetic_energy_tensor;
}

real ComputeMesh::compute_kinetic_energy(void)
{
    auto vertices = _system.get_vertices();
    real KE = 0.0;
    for (auto vertex : vertices)
    {
        KE += 0.5 * vertex.mass * vdot(vertex.v, vertex.v);
    }
    return KE;
}

real ComputeMesh::compute_temperature(void)
{
    auto KE = this->compute_kinetic_energy();
    real T = (2.0 / 3.0) * KE / _system.Numvertices;
    //KE = ((dim/2.0) N kT)
    return T;
}

std::vector<real> ComputeMesh::compute_pressure(EvolverClass &evolver)
{
    auto stress = this->compute_stresses(evolver, true);

    real V = 1.0; //_system.get_box_volume();
    std::vector<real> P;
    P.push_back(_system.Numvertices * stress.xx / V);
    P.push_back(_system.Numvertices * stress.yy / V);
    P.push_back(_system.Numvertices * stress.zz / V);
    //auto P = (stress.xx + stress.yy + stress.zz) / (3.0 * _system.get_box_volume());
    return P;
}

realTensor ComputeMesh::compute_stresses_virial(EvolverClass &evolver, const bool &include_velocities)
{
    evolver.reset_mesh_atom_stresses();
    evolver.compute_mesh_atom_stresses();

    realTensor zeroTensor;
    zeroTensor.xx = 0.0;
    zeroTensor.xy = 0.0;
    zeroTensor.xz = 0.0;
    zeroTensor.yx = 0.0;
    zeroTensor.yy = 0.0;
    zeroTensor.yz = 0.0;
    zeroTensor.zx = 0.0;
    zeroTensor.zy = 0.0;
    zeroTensor.zz = 0.0;

    auto stress_virial_atom = host::copy(_system.stress_virial_atom);
    auto total_stress = zeroTensor;

    for (auto s : stress_virial_atom)
    {
        total_stress.xx += s.xx;
        total_stress.xy += s.xy;
        total_stress.xz += s.xz;
        total_stress.yx += s.yx;
        total_stress.yy += s.yy;
        total_stress.yz += s.yz;
        total_stress.zx += s.zx;
        total_stress.zy += s.zy;
        total_stress.zz += s.zz;
    }

    //total_stress = thrust::reduce(_system.stress_virial_atom.begin(), _system.stress_virial_atom.end(), zeroTensor, host::reduce_tensor<realTensor>());

    if (include_velocities == true)
    {
        auto kinetic_energy_tensor = this->compute_kinetic_energy_tensor();
        total_stress.xx += kinetic_energy_tensor.xx;
        total_stress.xy += kinetic_energy_tensor.xy;
        total_stress.xz += kinetic_energy_tensor.xz;
        total_stress.yx += kinetic_energy_tensor.yx;
        total_stress.yy += kinetic_energy_tensor.yy;
        total_stress.yz += kinetic_energy_tensor.yz;
        total_stress.zx += kinetic_energy_tensor.zx;
        total_stress.zy += kinetic_energy_tensor.zy;
        total_stress.zz += kinetic_energy_tensor.zz;
    }

    total_stress.xx /= _system.Numvertices;
    total_stress.xy /= _system.Numvertices;
    total_stress.xz /= _system.Numvertices;
    total_stress.yx /= _system.Numvertices;
    total_stress.yy /= _system.Numvertices;
    total_stress.yz /= _system.Numvertices;
    total_stress.zx /= _system.Numvertices;
    total_stress.zy /= _system.Numvertices;
    total_stress.zz /= _system.Numvertices;
    return (total_stress);
}

std::vector<realTensor> ComputeMesh::compute_stresses_atom(EvolverClass &evolver, const bool &include_velocities)
{
    evolver.reset_mesh_atom_stresses();
    evolver.compute_mesh_atom_stresses();

    std::vector<realTensor> stress_virial_atom = _system.stress_virial_atom;

    if (include_velocities == true)
    {
        for (auto i = 0; i < _system.Numvertices; i++)
        {
            stress_virial_atom[i].xx += _system.vertices[i].mass * _system.vertices[i].v.x * _system.vertices[i].v.x;;
            stress_virial_atom[i].xy += _system.vertices[i].mass * _system.vertices[i].v.x * _system.vertices[i].v.y;;
            stress_virial_atom[i].xz += _system.vertices[i].mass * _system.vertices[i].v.x * _system.vertices[i].v.z;;
            stress_virial_atom[i].yx += _system.vertices[i].mass * _system.vertices[i].v.y * _system.vertices[i].v.x;;
            stress_virial_atom[i].yy += _system.vertices[i].mass * _system.vertices[i].v.y * _system.vertices[i].v.y;;
            stress_virial_atom[i].yz += _system.vertices[i].mass * _system.vertices[i].v.y * _system.vertices[i].v.z;;
            stress_virial_atom[i].zx += _system.vertices[i].mass * _system.vertices[i].v.z * _system.vertices[i].v.x;;
            stress_virial_atom[i].zy += _system.vertices[i].mass * _system.vertices[i].v.z * _system.vertices[i].v.y;;
            stress_virial_atom[i].zz += _system.vertices[i].mass * _system.vertices[i].v.z * _system.vertices[i].v.z;;
        }
    }

    return (stress_virial_atom);
}




