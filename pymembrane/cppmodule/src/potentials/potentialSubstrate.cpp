#include "potentialSubstrate.hpp"

void ComputeVertexLimitForce_kernel(int Numvertices,
                                    HE_VertexProp *vertices,
                                    real *_kz1,
                                    real *_kz2)
{
    for (int vertex_index = 0; vertex_index < Numvertices; vertex_index++)
    {
        int type = vertices[vertex_index].type;
        auto r = vertices[vertex_index].r;
        
        auto fval = 0.0;
        
        if(r.z>=0.0)
            fval = -_kz1[type]*r.z;
        else
            fval = -_kz2[type]*r.z;

        vertices[vertex_index].forceC.z+= fval;
    }
}
void ComputeVertexSubstrateEnergy::compute(void)
{

    ComputeVertexLimitForce_kernel(_system.Numedges,
                                   &_system.vertices[0],
                                   &kz1[0],
                                   &kz2[0]);
}
