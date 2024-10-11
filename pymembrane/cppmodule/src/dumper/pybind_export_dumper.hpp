#ifndef __pybind_export_dumper_hpp__
#define __pybind_export_dumper_hpp__

#include "dumper.hpp"

void export_DumperClass(py::module &m)
{
    py::class_<DumperClass>(m, "Dumper")
        .def(py::init<SystemClass &>())
        .def("vtk", &DumperClass::mesh_vtk)
        .def("edge_vtk", &DumperClass::mesh_edge_vtk)
        .def("json", &DumperClass::mesh_json)
        .def("ply", &DumperClass::mesh_ply)
        .def("setvtkLegacyFormat", &DumperClass::setvtkLegacyFormat)
        
        .def("user_data_vertex", (void (DumperClass::*)(std::string, std::vector<real>&)) & DumperClass::user_vertex_data)
        .def("user_data_vertex", (void (DumperClass::*)(std::string, std::vector<realTensor>&)) & DumperClass::user_vertex_data)

        .def("user_data_face", (void (DumperClass::*)(std::string, std::vector<real>&)) & DumperClass::user_face_data)
        .def("user_data_face", (void (DumperClass::*)(std::string, std::vector<realTensor>&)) & DumperClass::user_face_data)
        
        .def("user_data_edge", (void (DumperClass::*)(std::string, std::vector<real>&)) & DumperClass::user_edge_data)
        .def("user_data_edge", (void (DumperClass::*)(std::string, std::vector<realTensor>&)) & DumperClass::user_edge_data)

        //.def("user_data_vertex", (void (DumperClass::*)(std::string, std::vector<real3>&)) & DumperClass::user_vertex_data)
        //.def("user_data_face", &DumperClass::user_face_data)
        //.def("user_data_edge", &DumperClass::user_edge_data);
        ;
}

#endif