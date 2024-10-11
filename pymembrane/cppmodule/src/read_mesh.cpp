#include "read_mesh.hpp"
#include "mesh/computegeometry.hpp"
#include "../external/json_lib/json.hpp"
//operators
#include "mesh/meshoperators.hpp"

#include <pybind11/pybind11.h>
namespace py = pybind11;

void ReadMesh::__add_vertex(int id, real x, real y, real z, int type)
{
    //properties
    HE_VertexProp vertex;
    vertex.id = vertices.size(); //!< unique id
    vertex.r.x = x;              //!< position of the vertex in a 3D space
    vertex.r.y = y;              //!< position of the vertex in a 3D space
    vertex.r.z = z;              //!< position of the vertex in a 3D space
    vertex.boundary = true;      //!< by default all the vertices are not boundary
    vertex.coordination = -1;    //!< number of neighbours
    vertex._hedge = -1;          //!< no halfedge is assigned
    vertex.type = type;
    vertex.ip.x = 0;
    vertex.ip.y = 0;
    vertex.ip.z = 0;
    vertex.normal.x = 0.;
    vertex.normal.y = 0.;
    vertex.normal.z = 0.;
    vertices.push_back(vertex);
}
void ReadMesh::__add_face(int id, std::vector<int> verts, int orientation, int type)
{
    HE_FaceProp face; //Create new face
    face.id = id;
    face.outer = false;
    face.type = type;
    face.normal_reference.x = 0.0;
    face.normal_reference.y = 0.0;
    face.normal_reference.z = 1.0;

    int v1 = verts[0];
    int v2 = verts[1];
    int v3 = verts[2];

    if(orientation==-1)
    {
        int v_dummy = v2;
        v2 = v3;
        v3 = v_dummy;
    }
    verts[0] = v1;
    verts[1] = v2;
    verts[2] = v3;
    face.v1 = v1;
    face.v2 = v2;
    face.v3 = v3;
    //std::cout << v1 << " " << v2 << " " << v3 << std::endl;

    std::vector<HE_HalfEdgeProp> he_list;                         // Local storage for prev/next bookkeeping

    int N = verts.size();
    for (int i = 0; i < N; i++) // Loop over all vertices in the face
    {
        auto vi = verts[i];                                        //  First vertex
        auto vj = verts[(i + 1) % N];                              //  and its follower, that can be wrapped around
        //std::cout<< vi << " " << vj << std::endl;
        bool new_pair = false;                                     //  flag that checks if we created a new pair of half-edges

        HE_HalfEdgeProp he;     he.boundary=false;
        HE_HalfEdgeProp hep;    hep.boundary=false;
        std::pair<int,int> vij(vi,vj);
        std::pair<int,int> vji(vj,vi);
        if (hep_map.find(vij) == hep_map.end()) // if the half-edge does not exist
        {
            he.index = halfedges.size();        //create it
            hep.index = he.index + 1;           //and create its neighbour
            he.pair = hep.index;                //set them as pairs of each other
            hep.pair = he.index;
            hep.boundary = true;                //per construction, pair is boundary 
            hep_map[vji] = hep.index;
            new_pair = true;                    //set new_pair flag to True for post-processing*/
        }
        else                                    //if the half-edge exists retrieve it
        {
            he = halfedges[hep_map[vij]];
            he.boundary = false;                //per construction, it cannat be boundary 
            hep_map.erase(vij);
        }
        //Do the bookkeeping of the connectivity 
        he.vert_from = vi;                  
        he.vert_to = vj;                    
        he.face = face.id;                  
        he_list.push_back(he);
        if(i > 0)
        {
            if(he_list[i-1].index>=halfedges.size()) std::cout << "error here if(i > 0)" << std::endl;
            he.prev = he_list[i-1].index;
            halfedges[he_list[i-1].index].next = he.index;
        }        
        if(i == N-1)
        {
            if(he_list[0].index>=halfedges.size()) std::cout << "error here if(i == N-1)" << std::endl;
            he.next = he_list[0].index;
            halfedges[he_list[0].index].prev = he.index;
            face._hedge = he.index;
        }
        vertices[vi]._hedge = he.index;
        //Add new pair to the list of all half-edges
        if(new_pair)
        {
            halfedges.push_back(he);
            halfedges.push_back(hep);
        }
        else //Update the halfedge
        {
            halfedges[he.index] = he;
        }
        //construct the face and the vertices at the boundary 
        if(he.boundary == false)
            face.outer = false;
    }
    faces.push_back(face);
}

int ReadMesh::__add_edge(const HE_HalfEdgeProp& he, const HE_HalfEdgeProp& he_pair)
{
    HE_EdgeProp edge;
    edge.id = edges.size();
    edge.i = he.vert_from;
    edge.j = he.vert_to;
    edge.boundary = false;
    if((he.boundary==true) || (he_pair.boundary==true))  edge.boundary = true; 
    edge.face_k = he.face;             //!< index to one of the faces shared by this edge -1 if is outer face
    edge.face_l = he_pair.face;        //!< index to one of the faces shared by this edge -1 if is outer face
    edge._hedge = he.index;  
    edge.type = 0;
    edge.v0 = he.vert_to;
    edge.v2 = he.vert_from;
    int he_next = he.next;
    edge.v1 = halfedges[he_next].vert_to;
    int he_index_pair_next = he_pair.next;
    edge.v3 = halfedges[he_index_pair_next].vert_to;
    edges.push_back(edge);
    return(edge.id);
}
void ReadMesh::__build_boundary(void)
{
    // After all inner faces have been adde, what is left in the hep_map dictionary
    // are boundary edges. 
    // We need to connect them to each other in the proper order
    close_surface = true;
    for(auto item : hep_map)
    {
        py::print("something in the hep map");
        close_surface = false;
        auto vij = item.first;
        auto he = halfedges[hep_map[vij]];
        int i = vij.first;
        int j = vij.second;
        he.vert_from = i;
        he.vert_to = j; 
        he.face = -1;
        auto v = vertices[j];
        auto hev = halfedges[v._hedge];
        while(!hev.boundary)
        {
            hev = halfedges[halfedges[hev.prev].pair];
        }
        he.next = hev.index; 
        hev.prev = he.index;
        halfedges[he.index] = he;
        halfedges[hev.index] = hev;
    }
    Numhalfedges = halfedges.size();
}

void ReadMesh::__build_edges(void)
{
    std::map<std::pair<int,int>, int> edges_list;
    for(int index=0;index<halfedges.size();index++)
    {
        auto he = halfedges[index];
        auto he_pair = halfedges[he.pair];
        std::pair<int, int> vij(he.vert_from, he.vert_to);
        std::pair<int, int> vji(he.vert_to, he.vert_from);
        
        if (edges_list.find(vij) == edges_list.end())   // if the edge does not exist
        {
            ///update the edge reference in the halfedge            
            int edge_id = __add_edge(he, he_pair);
            edges_list[vij] = edge_id;
            edges_list[vji] = edge_id;
            halfedges[index].edge = edge_id;
        }
        else                                            //if the edge exists retrieve it
        {
            halfedges[index].edge = edges_list[vij];
        }
    }

}

void ReadMesh::__build_mesh(std::string &faces_file, std::string &vertices_file)
{
    /*! read vertex */
    std::ifstream file_vertex(vertices_file, std::ios::in);
    if (file_vertex.good())
    {
        std::string str;
        int num_triangles = 0;
        while (getline(file_vertex, str))
        {
            std::istringstream ss(str);
            //std::cout<< str << "\n";
            real num;
            std::vector<real> aux;
            while (ss >> num)
            {
                aux.push_back(num);
                //std::cout<< num << "\t";
            }
            this->__add_vertex(aux[0], aux[1], aux[2], aux[3], aux[4]);
            ///< sanity check
            //std::cout<< "vertex[" << vert.size()-1 << "] = < " << vert[vert.size()-1].x << "," << vert[vert.size()-1].y << "," << vert[vert.size()-1].z << " >\n";
        }
    }
    else
    {
        py::print("Error no vertex file is provided\n");
        exit(1);
    }

    py::print("finisehd vertices\n");

    /*! read triangles */
    std::ifstream file_triangles(faces_file, std::ios::in);
    if (file_triangles.good())
    {
        std::string str;
        int num_triangles = 0;
        while (getline(file_triangles, str))
        {
            std::istringstream ss(str);
            //std::cout<< str << "\n";
            int num;
            std::vector<int> aux;
            while (ss >> num)
            {
                aux.push_back(num);
            }
            std::vector<int> triangle = {aux[1], aux[2], aux[3]};
            __add_face(aux[0], triangle, aux[4], aux[5]);
        }
    }
    else
    {
        py::print("Error no triangles file is provided\n");
        exit(1);
    }

    py::print("finisehd triangles\n");
    this->__build_boundary();
    py::print("built boundary\n");
    this->__build_edges();
    py::print("built edges\n");
    std::transform(vertices.begin(), vertices.end(), vertices.begin(), host::reset_vertex_forces());
    std::transform(vertices.begin(), vertices.end(), vertices.begin(), host::reset_vertex_energy());
    std::transform(vertices.begin(), vertices.end(), vertices.begin(), host::reset_vertex_velocities());
    std::transform(edges.begin(), edges.end(), edges.begin(), host::reset_edge_energy());
    std::transform(faces.begin(), faces.end(), faces.begin(), host::reset_face_energy());
    py::print("finisehd initial read_mesh\n");
}


void ReadMesh::__build_mesh_common(std::vector<triangle_type>& tri,
                                   std::vector<vertices_type>& vert)
{
    //build vertices
    for(int vertex_index=0; vertex_index<vert.size(); vertex_index++)
    {
        __add_vertex(vert[vertex_index].id, vert[vertex_index].r.x, vert[vertex_index].r.y, vert[vertex_index].r.z, vert[vertex_index].type);
    }
    Numvertices = vertices.size();
    //build faces and halfedges
    for(int face_index=0; face_index<tri.size(); face_index++)
    {
        std::vector<int> verts;
        verts = {tri[face_index].v1, tri[face_index].v2, tri[face_index].v3};
        __add_face(tri[face_index].id, verts, 1, tri[face_index].type);
    }
    Numfaces = faces.size();                                   //!<  Num of faces
    //build boundaries
    //# After all inner faces have been added, what is left in the hep_map dictionary
    //# are boundary edges.
    //# We need to connect them to each other in the proper order
    close_surface = true;
    for(auto idx_map : hep_map)
    {
        close_surface = false;
        int he_index = idx_map.second;
        std::pair<int,int> vij = idx_map.first;
        int vi = vij.first;
        int vj = vij.second;

        //update the boundary_vertices
        vertices[vi].boundary = true;      //!< the vertex is not at the boundary

        //printf("vij < %i, %i> -> he_index = %i\n", vij.first, vij.second, he_index);
        halfedges[he_index].vert_from = vi;
        halfedges[he_index].vert_to = vj;
        halfedges[he_index].face = -1; //# This is the fictitious outer face.

        HE_HalfEdgeProp hev = halfedges[vertices[vj]._hedge];
        while(hev.boundary != true)
        {
            //printf("hev.boundary[%i] = %s\n", hev.index, hev.boundary ? "true" : "false");
            hev = halfedges[halfedges[hev.prev].pair];
        }
        //printf("End loop --- hev.boundary[%i] = %s\n", hev.index,  hev.boundary ? "true" : "false");
        halfedges[he_index].next = hev.index;
        halfedges[hev.index].prev = he_index;
    }
    Numhalfedges = halfedges.size();
    //build edges
    Numedges = 0;
    for(int he_index=0; he_index<Numhalfedges; he_index+=2)
    {
        int he_index_pair = halfedges[he_index].pair;
        __add_edge(halfedges[he_index], halfedges[he_index_pair]);
        //printf("index =%i \t pair = %i\n",halfedges[he_index].index, halfedges[he_index].pair);
    }
    Numedges = edges.size();

}

void ReadMesh::__json_mesh(std::string& json_file)
{
    using json = nlohmann::json;
    bool vertices_flag = false; //< check if everything was loaded
    bool faces_flag = false;
    // read a JSON file
    std::ifstream file_json(json_file);
    if (file_json.good())
    //if(true)
    {
        std::vector<triangle_type> tri;
        std::vector<vertices_type> vert;

        json json_mesh;
        file_json >> json_mesh;
        std::cout<< "--------------------------------------------------------" << std::endl;
        std::cout<< "Reading Mesh from JSON file " << std::endl;
        std::cout<< json_mesh["@@"] << std::endl;
        std::cout<< json_mesh["@@Date"] << std::endl;
        //std::cout<< mesh.dump() << std::endl;
        if (json_mesh.find("Vertices") != json_mesh.end())
        {
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

            json j_vertices = json_mesh["Vertices"];
            //json j_vertices = json_mesh["Vertices"];

            //std::cout<< " Numvertices " << Numvertices << std::endl;
            for(auto jvertex:j_vertices)
            {
                vert.resize(vert.size()+1);
                vert[vert.size()-1].id = jvertex["id"];
                vert[vert.size()-1].r.x = jvertex["r"][0];
                vert[vert.size()-1].r.y = jvertex["r"][1];
                vert[vert.size()-1].r.z = jvertex["r"][2];
                //vertices[i].id = jvertex["id"];
                //vertices[i].boundary = false;//jvertex["boundary"];
                //vertices[i].coordination = -1;//jvertex["coordination"];
                //vertices[i]._hedge = -1;//jvertex["hedge"];

                json properties = jvertex["Property"];
                //std::cout<< " vertex type " << properties["type"] << " vertex id " << jvertex["id"] <<std::endl;
                vert[vert.size()-1].type = properties["type"];
                vert[vert.size()-1].mass = properties["mass"];
                vert[vert.size()-1].energy = properties["energy"];
                //vert[vert.size()-1].v.x = 0.0;//properties["v"][0];
                //vert[vert.size()-1].v.y = 0.0;//properties["v"][1];
                //vert[vert.size()-1].v.z = 0.0;//properties["v"][2];
                //vert[vert.size()-1].age = properties["age"];
                vert[vert.size()-1].forceC.x = 0.0;//properties["forceC"][0];
                vert[vert.size()-1].forceC.y = 0.0;//properties["forceC"][1];
                vert[vert.size()-1].forceC.z = 0.0;//properties["forceC"][2];
                vert[vert.size()-1].forceD.x = 0.0;//properties["forceD"][0];
                vert[vert.size()-1].forceD.y = 0.0;//properties["forceD"][1];
                vert[vert.size()-1].forceD.z = 0.0;//properties["forceD"][2];
            }
            vertices_flag = true;
        }
        else { std::cerr<< "vertices not found"; exit(1); }
        if (json_mesh.find("Faces") != json_mesh.end() && vertices_flag==true)
        {
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



            json j_faces = json_mesh["Faces"];
            for(auto jface:j_faces)
            {
                tri.resize(tri.size()+1);
                tri[tri.size()-1].id = jface["id"];
                std::cout<< " face id " << jface["id"] <<std::endl;
                tri[tri.size()-1].v1 = jface["v1"];
                tri[tri.size()-1].v2 = jface["v2"];
                tri[tri.size()-1].v3 = jface["v3"];
                tri[tri.size()-1].orientation = 1;
                tri[tri.size()-1].is_reference_normal = true;
                tri[tri.size()-1].normal_reference.x = jface["normal_reference"][0];
                tri[tri.size()-1].normal_reference.y = jface["normal_reference"][1];
                tri[tri.size()-1].normal_reference.z = jface["normal_reference"][2];
                tri[tri.size()-1].is_reference_g = true;
                tri[tri.size()-1].g_reference[0] = jface["g_reference"][0];
                tri[tri.size()-1].g_reference[1] = jface["g_reference"][1];
                tri[tri.size()-1].g_reference[2] = jface["g_reference"][2];


                json properties = jface["Property"];
                PropertyFaces _tri_prop;
                //tri[tri.size()-1].type = properties["type"];
                tri[tri.size()-1].type = 0;
                //tri[tri.size()-1].age = properties["age"];
                tri[tri.size()-1].energy = properties["energy"];
                //std::cout<< " face energy " <<  properties["energy"] << " face id " << jface["id"] <<std::endl;
                //tri[tri.size()-1].prop.g_reference0_growth[0] = properties["g_reference0_growth"][0];
                //tri[tri.size()-1].prop.g_reference0_growth[1] = properties["g_reference0_growth"][1];
                //tri[tri.size()-1].prop.g_reference0_growth[2] = properties["g_reference0_growth"][2];
            }
            faces_flag = true;
        }
        else { std::cerr<< "faces not found"; exit(1); }
        if(faces_flag == true && vertices_flag == true)
        {
            std::sort(vert.begin(), vert.end());
            std::sort(tri.begin(), tri.end());
            for(auto tt:tri) std::cout<< tt.id << "\n";
            __build_mesh_common(tri, vert);
        }
    }
    else
    {
        std::cerr << json_file << std::endl;
        std::cerr<<"no json mesh provided"<<std::endl;
        exit(1);
    }
}
