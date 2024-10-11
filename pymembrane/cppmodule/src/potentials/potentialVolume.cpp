#include "potentialVolume.hpp"


double ComputeVertexVolumeEnergy::compute_volume(void)
{
    double vol = 0.0;
    for (int face_index = 0; face_index < _system.Numfaces; face_index++)
    {
        int v1 = _system.faces[face_index].v1;
        int v2 = _system.faces[face_index].v2;
        int v3 = _system.faces[face_index].v3;

        real3 normal = host::compute_normal_triangle(_system.vertices[v1].r, _system.vertices[v2].r, _system.vertices[v3].r);

        vol += vdot(normal, _system.vertices[v1].r)/6.;
    }
    return vol;
}
void ComputeVertexVolumeEnergy::compute_energy(void)
{    
    double vol = compute_volume();
    double totenergy = lambda * (1.-vol/vol0) * (1.-vol/vol0);
    //for (int face_index = 0; face_index < _system.Numfaces; face_index++)
      //  _system.faces[face_index]._property.energy += totenergy/(_system.Numfaces);
}

double ComputeVertexVolumeEnergy::compute_vertex_energy(int vertex_index)
{
    double vol = compute_volume();
    return lambda * (1.-vol/vol0) * (1.-vol/vol0);
}


void ComputeVertexVolumeEnergy::compute_vertex_forces(void)
{ 
    double vol = compute_volume();
    for(int face_index=0;face_index<_system.Numfaces; face_index++)
    {
        int v1 = _system.faces[face_index].v1;
        int v2 = _system.faces[face_index].v2;
        int v3 = _system.faces[face_index].v3;

        real3 normal = host::compute_normal_triangle(_system.vertices[v1].r, _system.vertices[v2].r, _system.vertices[v3].r);
        double force = lambda*(1.-vol/vol0)/3.0/vol0;
        double fx = force*normal.x;
        double fy = force*normal.y;
        double fz = force*normal.z;
        //v1
        //_system.vertices[v1]._property.forceC.x+=fx;
        //_system.vertices[v1]._property.forceC.y+=fy;
        //_system.vertices[v1]._property.forceC.z+=fz;

        //v2
        //_system.vertices[v2]._property.forceC.x+=fx;
        //_system.vertices[v2]._property.forceC.y+=fy;
        //_system.vertices[v2]._property.forceC.z+=fz;

        //v3
        //_system.vertices[v3]._property.forceC.x+=fx;
        //_system.vertices[v3]._property.forceC.y+=fy;
        //_system.vertices[v3]._property.forceC.z+=fz;
    }
}
