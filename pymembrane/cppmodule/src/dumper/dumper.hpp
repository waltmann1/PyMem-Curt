

#ifndef __dumper_hpp__
#define __dumper_hpp__

#include <string>
#include <memory>
#include <vector>
#include <map>
#include <iostream>

#include "../types/globaltypes.hpp"

#include "json_lib/json.hpp"

#include <chrono>
#include <fstream>

class SystemClass;

class DumperClass
{
public:
    DumperClass(SystemClass& system) : _system(system), vtkLegacy(true)
    {
    }
    void mesh_vertices_faces(const std::string &file_name);
    void mesh_edge_vtk(const std::string &file_name);
    void mesh_vtk(const std::string &file_name)
    {
        this->mesh_vertices_faces(file_name);
        this->mesh_edge_vtk(file_name);
    }
    void mesh_json(const std::string &file_name);
    void mesh_ply(const std::string &file_name);
    void setvtkLegacyFormat(const bool& flag)
    {
        vtkLegacy = flag;
    }
   
    void user_vertex_data(std::string name, std::vector<real>&data);
    void user_face_data(std::string name, std::vector<real>&data);
    void user_edge_data(std::string name, std::vector<real>&data);

    //void user_vertex_data(std::string name, std::vector<real3>&data);

    void user_vertex_data(std::string name, std::vector<realTensor>&data);
    void user_face_data(std::string name, std::vector<realTensor>&data);
    void user_edge_data(std::string name, std::vector<realTensor>&data);
    
private:
    SystemClass& _system;
    bool vtkLegacy;
    
    ///USER VTK DATA DUMP
    std::map<std::string,std::vector<realTensor>> vertex_tensor_data_map;
    std::map<std::string,std::vector<realTensor>> face_tensor_data_map;
    std::map<std::string,std::vector<realTensor>> edge_tensor_data_map;

    //std::map<std::string,std::vector<real3>> vertex_real3_data_map;

    std::map<std::string,std::vector<real>> vertex_data_map;
    std::map<std::string,std::vector<real>> face_data_map;
    std::map<std::string,std::vector<real>> edge_data_map;
    
};


#endif