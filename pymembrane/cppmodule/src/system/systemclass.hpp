/************************************************************************************
* MIT License                                                                       *
*                                                                                   *
* Copyright (c) 2020 Dr. Daniel Alejandro Matoz Fernandez                           *
*               fdamatoz@gmail.com                                                  *
* Permission is hereby granted, free of charge, to any person obtaining a copy      *
* of this software and associated documentation files (the "Software"), to deal     *
* in the Software without restriction, including without limitation the rights      *
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell         *
* copies of the Software, and to permit persons to whom the Software is             *
* furnished to do so, subject to the following conditions:                          *
*                                                                                   *
* The above copyright notice and this permission notice shall be included in all    *
* copies or substantial portions of the Software.                                   *
*                                                                                   *
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR        *
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,          *
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE       *
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER            *
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,     *
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE     *
* SOFTWARE.                                                                         *
*************************************************************************************/
#ifndef __systemclass_hpp__
#define __systemclass_hpp__

#include <algorithm>

#include "../types/hostvector.hpp"
#include "../box/box.hpp"
#include "../box/pbc.hpp"

//meshes
#include "../mesh/halfedges.hpp"
#include "../mesh/computegeometry.hpp"
#include "../mesh/meshoperators.hpp"
#include "read_mesh.hpp"

//compute
#include "../compute/computemesh.hpp"

//dumper
#include "../dumper/dumper.hpp"
#include <pybind11/pybind11.h>
namespace py = pybind11;

class SystemClass
{
public:
    SystemClass() : compute_mesh(*this), dumper(*this)
    {
        Numhalfedges = 0; //!< Number of halfedges
        Numvertices = 0;  //!< Number of vertices
        Numedges = 0;     //!< Number of edges
        Numfaces = 0;     //!< Number of faces
    }
    /**
     * Create a particles from given vectors
    */
    SystemClass(const BoxType &box) : compute_mesh(*this), dumper(*this), _box(box)
    {
        Numhalfedges = 0; //!< Number of halfedges
        Numvertices = 0;  //!< Number of vertices
        Numedges = 0;     //!< Number of edges
        Numfaces = 0;     //!< Number of faces
    }
    ~SystemClass() {}

    std::string get_mesh_info(void)
    {
        std::string info = "Mesh\n";
        info += " Numvertices = " + std::to_string(vertices.size()) + "\n";
        info += " NumFaces = " + std::to_string(faces.size()) + "\n";
        info += " NumEdges = " + std::to_string(edges.size()) + "\n";
        info += " NumHEdges = " + std::to_string(halfedges.size()) + "\n";
        return info;
    }

    /**
     * @brief get the particles loaded in the system
     * @param void 
     * @return const BoxType&
     */
    const BoxType &get_box(void) { return _box; }
    /**
     * @brief read the mesh from a file
     * @param std::map<std::string, std::string>& files 
     */
    void read_mesh_from_files(std::map<std::string, std::string> &files)
    {
        ReadMesh reader(files);
        host::vector<HE_Face<PropertyFaces>> _faces = reader.get_faces();             //!< Faces list
        host::vector<HE_Vertex<PropertyVertices>> _vertices = reader.get_vertices();  //!< Vertex list
        host::vector<HE_Edge<PropertyEdges>> _edges = reader.get_edges();             //!< Edges list
        host::vector<HE_HalfEdge<PropertyEdges>> _halfedges = reader.get_halfedges(); //!< Half-edge list
        close_surface = reader.is_close_surface();
        this->set(_faces, _vertices, _edges, _halfedges);
        //the mesh is no copy to host
        copy_in_host = false;
        compute_mesh.compute_vertex_normals();
        compute_mesh.compute_face_normals();
    }
    void read_mesh_from_json(std::string &json_file)
    {
        ReadMesh reader(json_file);
        host::vector<HE_Face<PropertyFaces>> _faces = reader.get_faces();             //!< Faces list
        host::vector<HE_Vertex<PropertyVertices>> _vertices = reader.get_vertices();  //!< Vertex list
        host::vector<HE_Edge<PropertyEdges>> _edges = reader.get_edges();             //!< Edges list
        host::vector<HE_HalfEdge<PropertyEdges>> _halfedges = reader.get_halfedges(); //!< Half-edge list
        close_surface = reader.is_close_surface();
        this->set(_faces, _vertices, _edges, _halfedges);
        //the mesh is no copy to host
        copy_in_host = false;
        compute_mesh.compute_vertex_normals();
        compute_mesh.compute_face_normals();
    }
    /**
     * @brief get the halfedges loaded in the system
     * @param void 
     * @return host::vector<HE_HalfEdge<PropertyEdges>> 
     */
    host::vector<HE_HalfEdge<PropertyEdges>> get_halfedges(void) { return (host::copy(halfedges)); }
    /**
     * @brief get the vertices loaded in the system
     * @param void 
     * @return host::vector<HE_Vertex<PropertyVertices>> 
     */
    host::vector<HE_Vertex<PropertyVertices>> get_vertices(void) { return (host::copy(vertices)); }
    /**
     * @brief get the edges loaded in the system
     * @param void 
     * @return host::vector<HE_Edge<PropertyEdges>> 
     */
    host::vector<HE_Edge<PropertyEdges>> get_edges(void) { return (host::copy(edges)); }
    /**
     * @brief get the faces loaded in the system
     * @param void 
     * @return host::vector<HE_Face<PropertyFaces>>
     */
    host::vector<HE_Face<PropertyFaces>> get_faces(void) { return (host::copy(faces)); }
    /**
     * @brief set the halfedges loaded in the system
     * @param host::vector<HE_HalfEdge<PropertyEdges>>& _halfedges
     * @return void
     */
    void set_halfedges(host::vector<HE_HalfEdge<PropertyEdges>> &_halfedges)
    {
        halfedges = _halfedges;
        Numhalfedges = _halfedges.size();
    }
    /**
     * @brief get the vertices loaded in the system
     * @param host::vector<HE_Vertex<PropertyVertices>>& _vertices 
     * @return  
     */
    void set_vertices(host::vector<HE_Vertex<PropertyVertices>> &_vertices)
    {
        vertices = _vertices;
        Numvertices = _vertices.size();
        compute_mesh.compute_vertex_normals();
        compute_mesh.compute_face_normals();
    }
    /**
     * @brief get the edges loaded in the system
     * @param host::vector<HE_Edge<PropertyEdges>>& _edges 
     * @return  
     */
    void set_edges(host::vector<HE_Edge<PropertyEdges>> &_edges)
    {
        edges = _edges;
        Numedges = _edges.size();
        //compute_mesh.compute_vertex_normals();
        //compute_mesh.compute_face_normals();
    }
    /**
     * @brief get the faces loaded in the system
     * @param void 
     * @return host::vector<HE_Face<PropertyFaces>>
     */
    void set_faces(host::vector<HE_Face<PropertyFaces>> &_faces)
    {
        faces = _faces;
        Numfaces = _faces.size();
        compute_mesh.compute_vertex_normals();
        compute_mesh.compute_face_normals();
    }
    void set(host::vector<HE_Face<PropertyFaces>> &_faces,         //!< Faces list
             host::vector<HE_Vertex<PropertyVertices>> &_vertices, //!< Vertex list
             host::vector<HE_Edge<PropertyEdges>> &_edges,         //!< Edges list
             host::vector<HE_HalfEdge<PropertyEdges>> &_halfedges  //!< Half-edge list
    )
    {
        Numhalfedges = _halfedges.size(); //!< Number of halfedges
        Numvertices = _vertices.size();   //!< Number of vertices
        Numedges = _edges.size();         //!< Number of edges
        Numfaces = _faces.size();         //!< Number of faces

        halfedges = _halfedges;
        vertices = _vertices;
        edges = _edges;
        faces = _faces;

        py::print("Mesh");
        py::print(" Numvertices ", vertices.size());
        py::print(" NumFaces ", faces.size());
        py::print(" NumEdges ", edges.size());
        py::print(" NumHEdges ", halfedges.size());
        
        //init the stresses
        this->init_stresses();
    }
    //other useful functions
    host::vector<int> get_edge_neighbours_host(int edge_index)
    {
        host::vector<int> edge_index_vec(5);
        int he0 = edges[edge_index]._hedge;
        int he0_next = halfedges[he0].next;
        int he0_prev = halfedges[he0].prev;
        int he0_pair = halfedges[he0].pair;
        int he0_pair_next = halfedges[he0_pair].next;
        int he0_pair_prev = halfedges[he0_pair].prev;
        edge_index_vec[0] = edge_index;
        edge_index_vec[1] = halfedges[he0_next].edge;
        edge_index_vec[2] = halfedges[he0_prev].edge;
        edge_index_vec[3] = halfedges[he0_pair_next].edge;
        edge_index_vec[4] = halfedges[he0_pair_prev].edge;
        return edge_index_vec;
    }
    //return compute mesh
    const ComputeMesh &get_compute_mesh(void) { return compute_mesh; }
    const DumperClass &get_dumper(void) { return dumper; }

    void init_stresses(void)
    {
        stress_group_faces.resize(Numfaces);
        stress_group_vertices.resize(Numvertices);
        stress_group_edges.resize(Numedges);
        stress_virial_atom.resize(Numvertices);
        stress_kinetic_atom.resize(Numvertices);
        std::transform(stress_group_faces.begin(), stress_group_faces.end(), stress_group_faces.begin(), host::reset_tensor<realTensor>());
        std::transform(stress_group_vertices.begin(), stress_group_vertices.end(), stress_group_vertices.begin(), host::reset_tensor<realTensor>());
        std::transform(stress_group_edges.begin(), stress_group_edges.end(), stress_group_edges.begin(), host::reset_tensor<realTensor>());
        std::transform(stress_virial_atom.begin(), stress_virial_atom.end(), stress_virial_atom.begin(), host::reset_tensor<realTensor>());
        std::transform(stress_kinetic_atom.begin(), stress_kinetic_atom.end(), stress_kinetic_atom.begin(), host::reset_tensor<realTensor>());
    }
    /**
     * @brief get the STRESSES in the vertices loaded in the system
     * @param void 
     * @return host::vector<HE_Edge<PropertyEdges>> 
     */
    host::vector<realTensor> get_stress_vertices(void) { return (host::copy(stress_group_vertices)); }
    /**
     * @brief get the STRESSES in the edges loaded in the system
     * @param void 
     * @return host::vector<HE_Edge<PropertyEdges>> 
     */
    host::vector<realTensor> get_stress_edges(void) { return (host::copy(stress_group_edges)); }
    /**
     * @brief get the faces loaded in the system
     * @param void 
     * @return host::vector<HE_Face<PropertyFaces>>
     */
    host::vector<realTensor> get_stress_faces(void) { return (host::copy(stress_group_faces)); }
    /**
     * @brief get the faces loaded in the system
     * @param void 
     * @return host::vector<HE_Face<PropertyFaces>>
     */
    host::vector<realTensor> get_stress_virial(void) { return (host::copy(stress_virial_atom)); }
    /**
     * @brief get the faces loaded in the system
     * @param void 
     * @return host::vector<HE_Face<PropertyFaces>>
     */
    host::vector<realTensor> get_stress_kinetic(void) { return (host::copy(stress_kinetic_atom)); }
    //private:
    //< mesh
    host::vector<HE_Face<PropertyFaces>> faces;         //!< Faces list
    host::vector<HE_Vertex<PropertyVertices>> vertices; //!< Vertex list
    host::vector<HE_HalfEdge<PropertyEdges>> halfedges; //!< Half-edge list
    host::vector<HE_Edge<PropertyEdges>> edges;         //!< Edges list

    //<group of stresses
    host::vector<realTensor> stress_group_faces;
    host::vector<realTensor> stress_group_vertices;
    host::vector<realTensor> stress_group_edges;
    host::vector<realTensor> stress_virial_atom;
    host::vector<realTensor> stress_kinetic_atom;

    int Numvertices;    //!< Number of vertices
    int Numedges;       //!< Number of edges
    int Numfaces;       //!< Number of faces
    int Numhalfedges;   //!< Number of halfedges
    bool close_surface; //!< Mesh is a close surface
    bool copy_in_host;
    //compute
    ComputeMesh compute_mesh;
    //dumper
    DumperClass dumper;

private:
    BoxType _box;
};

#endif