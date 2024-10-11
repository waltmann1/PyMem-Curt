

#include "dumper.hpp"

#include <memory>
#include "../mesh/computegeometry.hpp"
#include "../mesh/halfedges.hpp"
#include "../types/globaltypes.hpp"
#include "../types/hostvector.hpp"
#include "../system/systemclass.hpp"

#include <vtkVersion.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolygon.h>
#include <vtkLine.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>

void DumperClass::user_vertex_data(std::string name, std::vector<real> &data)
{
    if (data.size() == _system.Numvertices)
        vertex_data_map[name] = data;
    else
        std::cerr << "data must have the same size that the number of vertices " << std::endl;
}

void DumperClass::user_face_data(std::string name, std::vector<real> &data)
{
    if (data.size() == _system.Numfaces)
        face_data_map[name] = data;
    else
        std::cerr << "data must have the same size that the number of faces " << std::endl;
}
void DumperClass::user_edge_data(std::string name, std::vector<real> &data)
{
    if (data.size() == _system.Numedges)
        edge_data_map[name] = data;
    else
        std::cerr << "data must have the same size that the number of edges " << std::endl;
}

//void DumperClass::user_vertex_data(std::string name, std::vector<real3> &data)
//{
//    if (data.size() == _system.Numvertices)
//        vertex_real3_data_map_data_map[name] = data;
//    else
//        std::cerr << "data must have the same size that the number of vertices " << std::endl;
//}

void DumperClass::user_vertex_data(std::string name, std::vector<realTensor> &data)
{
    if (data.size() == _system.Numvertices)
        vertex_tensor_data_map[name] = data;
    else
        std::cerr << "data must have the same size that the number of vertices " << std::endl;
}

void DumperClass::user_face_data(std::string name, std::vector<realTensor> &data)
{
    if (data.size() == _system.Numfaces)
        face_tensor_data_map[name] = data;
    else
        std::cerr << "data must have the same size that the number of faces " << std::endl;
}
void DumperClass::user_edge_data(std::string name, std::vector<realTensor> &data)
{
    if (data.size() == _system.Numedges)
        edge_tensor_data_map[name] = data;
    else
        std::cerr << "data must have the same size that the number of edges " << std::endl;
}
void DumperClass::mesh_vertices_faces(const std::string &file_name)
{
    auto vertices = _system.get_vertices();
    auto faces = _system.get_faces();

    std::string file_name2 = file_name + ".vtp";

    vtkSmartPointer<vtkPolyData> polydata_vtk = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> vertices_vtk = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> triangles_vtk = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkDoubleArray> normals = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> forces = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> velocities = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> coord = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> vnormal = vtkSmartPointer<vtkDoubleArray>::New();

    vtkSmartPointer<vtkPolygon> triangle_vtk = vtkSmartPointer<vtkPolygon>::New();
    vtkSmartPointer<vtkIntArray> ids_vtk = vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkDoubleArray> areas_vtk = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkIntArray> face_type_vtk = vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkIntArray> vertex_type_vtk = vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkDoubleArray> strain_xx = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> strain_xy = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> strain_yy = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> face_pressure = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> face_energy = vtkSmartPointer<vtkDoubleArray>::New();

    //Faces
    face_type_vtk->SetName("FaceType");
    face_type_vtk->SetNumberOfComponents(1);
    strain_xx->SetName("CauchyStrain_XX");
    strain_xx->SetNumberOfComponents(1);
    strain_xy->SetName("CauchyStrain_XY");
    strain_xy->SetNumberOfComponents(1);
    strain_yy->SetName("CauchyStrain_YY");
    strain_yy->SetNumberOfComponents(1);
    normals->SetName("FaceNormal");
    normals->SetNumberOfComponents(3);
    areas_vtk->SetName("FaceArea");
    areas_vtk->SetNumberOfComponents(1);
    face_pressure->SetName("face_pressure");
    face_pressure->SetNumberOfComponents(1);
    face_energy->SetName("face_energy");
    face_energy->SetNumberOfComponents(1);
    std::vector<vtkSmartPointer<vtkDoubleArray>> face_data_vtk;
    int facenumitems = 0;
    for (auto const &item : face_data_map)
    {
        face_data_vtk.resize(facenumitems + 1);
        face_data_vtk[facenumitems] = vtkSmartPointer<vtkDoubleArray>::New();
        const char *cstr = (item.first).c_str();
        face_data_vtk[facenumitems]->SetName(cstr);
        face_data_vtk[facenumitems]->SetNumberOfComponents(1);
        facenumitems++;
    }
    std::vector<vtkSmartPointer<vtkDoubleArray>> face_tensor_data_vtk;
    facenumitems = 0;
    for (auto const &item : face_tensor_data_map)
    {
        face_tensor_data_vtk.resize(facenumitems + 1);
        face_tensor_data_vtk[facenumitems] = vtkSmartPointer<vtkDoubleArray>::New();
        const char *cstr = (item.first).c_str();
        face_tensor_data_vtk[facenumitems]->SetName(cstr);
        face_tensor_data_vtk[facenumitems]->SetNumberOfComponents(9);
        facenumitems++;
    }
    //Vertices
    ids_vtk->SetName("Id");
    ids_vtk->SetNumberOfComponents(1);
    coord->SetName("VertexCoord");
    coord->SetNumberOfComponents(3);
    vertex_type_vtk->SetName("VertexType");
    vertex_type_vtk->SetNumberOfComponents(1);
    forces->SetName("VertexForce");
    forces->SetNumberOfComponents(3);
    velocities->SetName("VertexVel");
    velocities->SetNumberOfComponents(3);
    vnormal->SetName("VertexNormal");
    vnormal->SetNumberOfComponents(3);
    ///USER DATA
    std::vector<vtkSmartPointer<vtkDoubleArray>> vertex_data_vtk;
    int vertexnumitems = 0;
    for (auto const &item : vertex_data_map)
    {
        vertex_data_vtk.resize(vertexnumitems + 1);
        vertex_data_vtk[vertexnumitems] = vtkSmartPointer<vtkDoubleArray>::New();
        const char *cstr = (item.first).c_str();
        vertex_data_vtk[vertexnumitems]->SetName(cstr);
        vertex_data_vtk[vertexnumitems]->SetNumberOfComponents(1);
        vertexnumitems++;
    }
    std::vector<vtkSmartPointer<vtkDoubleArray>> vertex_tensor_data_vtk;
    vertexnumitems = 0;
    for (auto const &item : vertex_tensor_data_map)
    {
        vertex_tensor_data_vtk.resize(vertexnumitems + 1);
        vertex_tensor_data_vtk[vertexnumitems] = vtkSmartPointer<vtkDoubleArray>::New();
        const char *cstr = (item.first).c_str();
        vertex_tensor_data_vtk[vertexnumitems]->SetName(cstr);
        vertex_tensor_data_vtk[vertexnumitems]->SetNumberOfComponents(9);
        vertexnumitems++;
    }

    for (unsigned int i = 0; i < faces.size(); i++)
    {
        int v1 = faces[i].v1;
        int v2 = faces[i].v2;
        int v3 = faces[i].v3;
        auto r1 = vertices[v1].r;
        auto r2 = vertices[v2].r;
        auto r3 = vertices[v3].r;

        triangle_vtk->GetPointIds()->SetNumberOfIds(3);
        triangle_vtk->GetPointIds()->SetId(0, v1);
        triangle_vtk->GetPointIds()->SetId(1, v2);
        triangle_vtk->GetPointIds()->SetId(2, v3);
        triangles_vtk->InsertNextCell(triangle_vtk);
        face_type_vtk->InsertNextValue(faces[i].type);

        real3 normal = host::compute_normal_triangle(vertices[v1].r, vertices[v2].r, vertices[v3].r);
        real normal_norm = sqrt(vdot(normal, normal));
        real face_area = (0.5 * normal_norm);
        areas_vtk->InsertNextValue(face_area);

        normal.x /= normal_norm;
        normal.y /= normal_norm;
        normal.z /= normal_norm;
        real normal_tuple[3] = {normal.x, normal.y, normal.z};
        normals->InsertNextTuple(normal_tuple);

        //< strains
        real metric_now[3];
        host::compute_form_factor_triangle(metric_now, vertices[v1].r, vertices[v2].r, vertices[v3].r);

        strain_xx->InsertNextValue((metric_now[0] - faces[i].g_reference[0]));
        strain_xy->InsertNextValue((metric_now[1] - faces[i].g_reference[1]));
        strain_yy->InsertNextValue((metric_now[2] - faces[i].g_reference[2]));

        face_energy->InsertNextValue(faces[i].energy);
        face_pressure->InsertNextValue(((metric_now[0] - faces[i].g_reference[0]) + (metric_now[2] - faces[i].g_reference[2])));
        ///USER DATA
        int item_index = 0;
        for (auto const &item : face_data_map)
        {
            face_data_vtk[item_index]->InsertNextValue((item.second)[i]);
            item_index++;
        }
        item_index = 0;
        for (auto const &item : face_tensor_data_map)
        {
            real tensor_data[9] = {((item.second)[i]).xx, (item.second)[i].xy, (item.second)[i].xz,
                                   ((item.second)[i]).yx, (item.second)[i].yy, (item.second)[i].yz,
                                   ((item.second)[i]).zx, (item.second)[i].zy, (item.second)[i].zz};
            face_tensor_data_vtk[item_index]->InsertNextTuple(tensor_data);
            item_index++;
        }
    }

    polydata_vtk->SetPolys(triangles_vtk);
    //polydata_vtk->GetCellData()->AddArray(ids);
    polydata_vtk->GetCellData()->AddArray(areas_vtk);
    polydata_vtk->GetCellData()->SetNormals(normals);
    polydata_vtk->GetCellData()->AddArray(face_type_vtk);
    polydata_vtk->GetCellData()->AddArray(strain_xx);
    polydata_vtk->GetCellData()->AddArray(strain_xy);
    polydata_vtk->GetCellData()->AddArray(strain_yy);
    polydata_vtk->GetCellData()->AddArray(face_pressure);
    polydata_vtk->GetCellData()->AddArray(face_energy);
    ///USER DATA
    for (auto vtk_user_data : face_data_vtk)
        polydata_vtk->GetCellData()->AddArray(vtk_user_data);
    for (auto vtk_user_data : face_tensor_data_vtk)
        polydata_vtk->GetCellData()->AddArray(vtk_user_data);

    for (unsigned int i = 0; i < vertices.size(); i++)
    {
        vertices_vtk->InsertNextPoint(vertices[i].r.x, vertices[i].r.y, vertices[i].r.z);
        ids_vtk->InsertNextValue(vertices[i].id);
        real pos[3] = {vertices[i].r.x, vertices[i].r.y, vertices[i].r.z};
        coord->InsertNextTuple(pos);
        vertex_type_vtk->InsertNextValue(vertices[i].type);

        real f[3] = {vertices[i].forceC.x, vertices[i].forceC.y, vertices[i].forceC.z};
        forces->InsertNextTuple(f);

        real vel[3] = {vertices[i].v.x, vertices[i].v.y, vertices[i].v.z};
        velocities->InsertNextTuple(vel);

        real norm[3] = {vertices[i].normal.x, vertices[i].normal.y, vertices[i].normal.z};
        real norm_norm = sqrt(norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2]);
        if (norm_norm > 0.0)
        {
            norm[0] /= norm_norm;
            norm[1] /= norm_norm;
            norm[2] /= norm_norm;
        }
        vnormal->InsertNextTuple(norm);
        ///USER DATA
        int item_index = 0;
        for (auto const &item : vertex_data_map)
        {
            vertex_data_vtk[item_index]->InsertNextValue((item.second)[i]);
            item_index++;
        }
        item_index = 0;
        for (auto const &item : vertex_tensor_data_map)
        {
            real tensor_data[9] = {((item.second)[i]).xx, (item.second)[i].xy, (item.second)[i].xz,
                                   ((item.second)[i]).yx, (item.second)[i].yy, (item.second)[i].yz,
                                   ((item.second)[i]).zx, (item.second)[i].zy, (item.second)[i].zz};
            vertex_tensor_data_vtk[item_index]->InsertNextTuple(tensor_data);
            item_index++;
        }
    }
    polydata_vtk->SetPoints(vertices_vtk);
    polydata_vtk->GetPointData()->AddArray(ids_vtk);
    polydata_vtk->GetPointData()->AddArray(coord);
    polydata_vtk->GetPointData()->AddArray(vertex_type_vtk);
    polydata_vtk->GetPointData()->AddArray(forces);
    polydata_vtk->GetPointData()->AddArray(velocities);
    polydata_vtk->GetPointData()->AddArray(vertex_type_vtk);
    polydata_vtk->GetPointData()->AddArray(vnormal);
    ///USER DATA
    for (auto vtk_user_data : vertex_data_vtk)
        polydata_vtk->GetPointData()->AddArray(vtk_user_data);

    for (auto vtk_user_data : vertex_tensor_data_vtk)
        polydata_vtk->GetPointData()->AddArray(vtk_user_data);
    // Write the file
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(file_name2.c_str());
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(polydata_vtk);
#else
    writer->SetInputData(polydata_vtk);
#endif
    if (vtkLegacy == true)
        writer->SetDataModeToAscii();
    else
        writer->SetDataModeToBinary();
    writer->Write();
}

void DumperClass::mesh_edge_vtk(const std::string &file_name)
{
    auto vertices = _system.get_vertices();
    auto edges = _system.get_edges();

    std::string file_name2 = file_name + "_edges.vtp";

    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

    vtkSmartPointer<vtkLine> edge = vtkSmartPointer<vtkLine>::New();
    vtkSmartPointer<vtkIntArray> ids = vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkDoubleArray> lens = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkIntArray> type = vtkSmartPointer<vtkIntArray>::New();

    ids->SetName("Id");
    ids->SetNumberOfComponents(1);
    lens->SetName("Length");
    lens->SetNumberOfComponents(1);
    type->SetName("edge_type");
    type->SetNumberOfComponents(1);

    ///USER DATA
    std::vector<vtkSmartPointer<vtkDoubleArray>> edge_data_vtk;
    int edgenumitems = 0;
    for (auto const &item : edge_data_map)
    {
        edge_data_vtk.resize(edgenumitems + 1);
        edge_data_vtk[edgenumitems] = vtkSmartPointer<vtkDoubleArray>::New();
        const char *cstr = (item.first).c_str();
        edge_data_vtk[edgenumitems]->SetName(cstr);
        edge_data_vtk[edgenumitems]->SetNumberOfComponents(1);
        edgenumitems++;
    }

    for (unsigned int i = 0; i < vertices.size(); i++)
    {
        points->InsertNextPoint(vertices[i].r.x, vertices[i].r.y, vertices[i].r.z);
        ids->InsertNextValue(vertices[i].id);
    }

    polydata->SetPoints(points);
    polydata->GetPointData()->AddArray(ids);

    for (unsigned int i = 0; i < edges.size(); i++)
    {
        edge->GetPointIds()->SetId(0, edges[i].i);
        edge->GetPointIds()->SetId(1, edges[i].j);
        lines->InsertNextCell(edge);
        real3 rij;
        vsub(rij, vertices[edges[i].j].r, vertices[edges[i].i].r);
        real length = sqrt(vdot(rij, rij));
        lens->InsertNextValue(length);
        type->InsertNextValue(edges[i].type);
        ///USER DATA
        int item_index = 0;
        for (auto const &item : edge_data_map)
        {
            edge_data_vtk[item_index]->InsertNextValue((item.second)[i]);
            item_index++;
        }
    }
    polydata->SetLines(lines);
    polydata->GetCellData()->AddArray(lens);
    polydata->GetCellData()->AddArray(type);
    ///USER DATA
    for (auto vtk_user_data : edge_data_vtk)
        polydata->GetPointData()->AddArray(vtk_user_data);
    // Write the file
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(file_name2.c_str());
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(polydata);
#else
    writer->SetInputData(polydata);
#endif
    if (vtkLegacy == true)
        writer->SetDataModeToAscii();
    else
        writer->SetDataModeToBinary();
    writer->Write();
}

void DumperClass::mesh_ply(const std::string &file_name)
{
    auto vertices = _system.get_vertices();
    auto edges = _system.get_edges();
    auto faces = _system.get_faces();
    auto halfedges = _system.get_halfedges();

    int Numvertices = vertices.size();
    int Numedges = edges.size();
    int Numfaces = faces.size();

    std::string ply_file2 = file_name + ".ply";
    std::ofstream myfile(ply_file2);
    if (myfile.is_open())
    {
        /*
            PLY file
            ply
            format ascii 1.0
            comment created by platoply
            element vertex 8
            property float32 x
            property float32 y
            property float32 z
            element face 6
            property list uint8 int32 vertex_indices
            end_header
            -1 -1 -1 
            1 -1 -1 
            1 1 -1 
            -1 1 -1 
            -1 -1 1 
            1 -1 1 
            1 1 1 
            -1 1 1 
            4 0 1 2 3 
            4 5 4 7 6 
            4 6 2 1 5 
            4 3 7 4 0 
            4 7 3 2 6 
            4 5 1 0 4 
            */
        myfile << "ply\n";
        myfile << "format ascii 1.0\n";
        myfile << "comment created by D. A. Matoz-Fernandez\n";
        myfile << "element vertex " << Numvertices << "\n";
        myfile << "property float32 x\n";
        myfile << "property float32 y\n";
        myfile << "property float32 z\n";
        myfile << "element face " << Numfaces << "\n";
        myfile << "property list uint8 int32 vertex_indices\n";
        myfile << "property float32 nx\n";
        myfile << "property float32 ny\n";
        myfile << "property float32 nz\n";
        //myfile << "property uchar red\n";
        //myfile << "property uchar green\n";
        //myfile << "property uchar blue\n";
        myfile << "element edge " << Numedges << "\n";
        myfile << "property int vertex1\n";
        myfile << "property int vertex2\n";
        myfile << "end_header\n";
        for (int vert_index = 0; vert_index < Numvertices; vert_index++)
        {
            myfile << vertices[vert_index].r.x << " " << vertices[vert_index].r.y << " " << vertices[vert_index].r.z << "\n";
        }
        for (int face_index = 0; face_index < Numfaces; face_index++)
        {
            int v1 = faces[face_index].v1; ///< vertices index
            int v2 = faces[face_index].v2; ///< vertices index
            int v3 = faces[face_index].v3; ///< vertices index

            //int red = 255;
            //int blue = 0;
            //int green = 255;
            //myfile << "3" << " " << v1 << " " << v2 << " " << v3 << " " << red << " " << green << " " << blue << "\n";

            real3 normal, v12, v13;
            vsub(v12, vertices[v2].r, vertices[v1].r);
            vsub(v13, vertices[v3].r, vertices[v1].r);
            vcross(normal, v12, v13);
            real normal_norm = sqrt(vdot(normal, normal));
            normal.x /= normal_norm;
            normal.y /= normal_norm;
            normal.z /= normal_norm;

            myfile << "3"
                   << " " << v1 << " " << v2 << " " << v3 << " ";
            myfile << normal.x << " " << normal.y << " " << normal.z << std::endl;
        }

        for (int id = 0; id < Numedges; id++)
        {
            int v1 = edges[id].i;
            int v2 = edges[id].j;
            myfile << v1 << " " << v2 << "\n";
        }
        myfile.close();
    }
    else
        std::cerr << "Unable to open dump mesh file\n";
}

/**
    * Json dumper
    **/
void DumperClass::mesh_json(const std::string &file_name)
{

    auto vertices = _system.get_vertices();
    auto edges = _system.get_edges();
    auto faces = _system.get_faces();
    auto halfedges = _system.get_halfedges();

    int Numvertices = vertices.size();
    int Numedges = edges.size();
    int Numfaces = faces.size();
    int Numhalfedges = halfedges.size();

    using json = nlohmann::json;
    std::map<std::string, json> j_vertices;
    std::map<std::string, json> j_edges;
    std::map<std::string, json> j_faces;
    std::map<std::string, json> j_halfedges;

    json j_mesh;
    using namespace std::chrono;
    system_clock::time_point today = system_clock::now();
    time_t tt = system_clock::to_time_t(today);
    j_mesh["@@"] = "Created with membrane code. All right reserved to Dr. D.A Matoz-Fernandez, for more information contact fdamatoz@gmail.com";
    j_mesh["@@Date"] = ctime(&tt);
    j_mesh["@Numvertices"] = Numvertices;
    j_mesh["@Numedges"] = Numedges;
    j_mesh["@Numfaces"] = Numfaces;
    j_mesh["@Numhalfedges"] = Numhalfedges;

    /*
################### HE_Vertex   ###################
####################################################
real3 r;              //!< position of the vertex in a 3D space
int id;               //!< unique id
bool boundary;        //!< if true, vertex is on boundary
int  coordination;    //!< number of neighbours this vertex has
//handle
int _hedge;           //!< HANDLE INDEX OF: one of the half-edges emanating from the vertex
PropertyVertices _property;
################### properties   ###################
####################################################
int type;                //!< Vertex/Edge/Face Type
real3 v;                 //!< Vertex/Edge/Face Velocity
real age;                //!< Vertex/Edge/Face Age
real3 forceC;            //!< Vertex/Edge/Face conservative force
real3 forceD;            //!< Vertex/Edge/Face dissipative force
real energy;             //!< Vertex/Edge/Face energy
real mass;               //!< Vertex mass
*/

    for (unsigned int i = 0; i < Numvertices; i++)
    {
        json properties;
        properties["type"] = vertices[i].type;
        properties["v"] = {vertices[i].v.x, vertices[i].v.y, vertices[i].v.z};
        //properties["age"] = vertices[i].age;
        properties["forceC"] = {vertices[i].forceC.x, vertices[i].forceC.y, vertices[i].forceC.z};
        properties["forceD"] = {vertices[i].forceD.x, vertices[i].forceD.y, vertices[i].forceD.z};
        properties["energy"] = vertices[i].energy;
        properties["mass"] = vertices[i].mass;

        json jvertex;
        jvertex["r"] = {vertices[i].r.x, vertices[i].r.y, vertices[i].r.z};
        jvertex["id"] = vertices[i].id;
        jvertex["boundary"] = vertices[i].boundary;
        jvertex["coordination"] = vertices[i].coordination;
        jvertex["hedge"] = vertices[i]._hedge;
        jvertex["Property"] = properties;
        j_vertices[std::to_string(i)] = jvertex;
    }
    j_mesh["Vertices"] = j_vertices;

    /*
################### HE_Edge   ###################
####################################################
int i, j;               //!< indices of two vertices that connect the edge
int id;                 //!< unique id
bool boundary;          //!< if true, edge is a boundary edge
int face_k;             //!< index to one of the faces shared by this edge -1 if is outer face
int face_l;             //!< i01_disclination   ###################
####################################################
int type;                //!< Vertex/Edge/Face Type
real age;                //!< Vertex/Edge/Face Age
real energy;             //!< Vertex/Edge/Face energy
*/
    for (unsigned int i = 0; i < Numedges; i++)
    {
        json properties;
        properties["type"] = edges[i].type;
        //properties["age"] = edges[i].age;
        properties["energy"] = edges[i].energy;

        json jedge;
        jedge["i"] = edges[i].i;
        jedge["j"] = edges[i].j;
        jedge["id"] = edges[i].id;
        jedge["boundary"] = edges[i].boundary;
        jedge["face_k"] = edges[i].face_k;
        jedge["face_l"] = edges[i].face_l;
        jedge["v0"] = edges[i].v0;
        jedge["v1"] = edges[i].v1;
        jedge["v2"] = edges[i].v2;
        jedge["v3"] = edges[i].v3;
        jedge["hedge"] = edges[i]._hedge;
        jedge["Property"] = properties;
        j_edges[std::to_string(i)] = jedge;
    }
    j_mesh["Edges"] = j_edges;
    /*
################### HE_Face   ###################
####################################################
//properties
int id;                 //!< face id
bool outer;             //!< if true, face is a ghost outer face
int nsides;             //!< number of sides face has
real3 normal;           //!< normal to that face
real area;              //!< area of the triangle
real g_reference[3];    //!< reference metric tensor [0] |(v2-v1)|^2 [1] (v2-v1)*(v3-v1) [2]|(v3-v1)|^2
real g_reference_inv[3];//!< reference inverse metric tensor
real3 normal_reference; //!< reference normal
int v1,v2,v3;           //!< the 3 vertices that belong to that face
bool boundary;          //!< if true, the face is a boundary
//handle
int _hedge;             //!< HANDLE INDEX OF: one of the half-edges bordering the face
PropertyFaces _property;
################### properties   ###################
####################################################
int type;                //!< Vertex/Edge/Face Type
real age;                //!< Vertex/Edge/Face Age
real energy;             //!< Vertex/Edge/Face energy
real g_reference0_growth[3];//!< reference metric tensor [0] |(v2-v1)|^2 [1] (v2-v1)*(v3-v1) [3]|(v3-v1)|^2 used for growth
*/
    for (unsigned int i = 0; i < Numfaces; i++)
    {
        json properties;
        properties["type"] = faces[i].type;
        //properties["age"] = faces[i].age;
        properties["energy"] = faces[i].energy;
        //properties["g_reference0_growth"] = {faces[i].g_reference0_growth[0], faces[i].g_reference0_growth[1], faces[i].g_reference0_growth[2]};

        json jface;
        jface["id"] = faces[i].id;
        jface["outer"] = faces[i].outer;
        jface["nsides"] = faces[i].nsides;
        jface["boundary"] = faces[i].boundary;
        ///
        int v1 = faces[i].v1; ///< vertices index
        int v2 = faces[i].v2; ///< vertices index
        int v3 = faces[i].v3; ///< vertices index

        real3 normal, v12, v13;
        vsub(v12, vertices[v2].r, vertices[v1].r);
        vsub(v13, vertices[v3].r, vertices[v1].r);
        vcross(normal, v12, v13);
        real normal_norm = sqrt(vdot(normal, normal));
        normal.x /= normal_norm;
        normal.y /= normal_norm;
        normal.z /= normal_norm;

        real face_area = (0.5 * normal_norm);
        ///
        jface["normal"] = {normal.x, normal.y, normal.z};
        jface["area"] = face_area;
        jface["g_reference"] = {faces[i].g_reference[0], faces[i].g_reference[1], faces[i].g_reference[2]};
        jface["g_reference_inv"] = {faces[i].g_reference_inv[0], faces[i].g_reference_inv[1], faces[i].g_reference_inv[2]};
        jface["normal_reference"] = {faces[i].normal_reference.x, faces[i].normal_reference.y, faces[i].normal_reference.z};
        jface["v1"] = faces[i].v1;
        jface["v2"] = faces[i].v2;
        jface["v3"] = faces[i].v3;
        jface["hedge"] = faces[i]._hedge;
        jface["Property"] = properties;
        j_faces[std::to_string(i)] = jface;
    }
    j_mesh["Faces"] = j_faces;
    /*
################### HE_Edge   ###################
####################################################
//handle
int index;           //!< Index, used for half-edge sorting
int vert_from;       //!< HANDLE INDEX OF: vertex at the beginning of the half-edge
int vert_to;         //!< HANDLE INDEX OF: vertex at the end of the half-edge
int edge;            //!< HANDLE INDEX OF: edge this he is part of
int face;            //!< HANDLE INDEX OF: face the half-edge borders
int pair;            //!< HANDLE INDEX OF: oppositely oriented adjacent half-edge
int next;            //!< HANDLE INDEX OF: next half-edge around the face
int prev;            //!< HANDLE INDEX OF: previous half-edge around the face
bool boundary;       //!< if true, the halfedge is a boundary
PropertyEdges _property;
*/
    for (unsigned int i = 0; i < Numhalfedges; i++)
    {
        json jhedge;
        jhedge["vert_from"] = halfedges[i].vert_from;
        jhedge["vert_to"] = halfedges[i].vert_to;
        jhedge["edge"] = halfedges[i].edge;
        jhedge["face"] = halfedges[i].face;
        jhedge["pair"] = halfedges[i].pair;
        jhedge["next"] = halfedges[i].next;
        jhedge["prev"] = halfedges[i].prev;
        jhedge["boundary"] = halfedges[i].boundary;
        j_halfedges[std::to_string(i)] = jhedge;
    }
    j_mesh["Halfedges"] = j_halfedges;

    // write prettified JSON to another file
    std::string file_name2 = file_name + ".json";
    std::ofstream output_file(file_name2);
    output_file << std::setw(4) << j_mesh << std::endl;
}