#ifndef __computegeometry_hpp__
#define __computegeometry_hpp__

#pragma once

#include "../types/globaltypes.hpp"
#include "../box/pbc.hpp"

namespace host
{
    /** @defgroup ComputeGPUfn Geometry
*  @brief ComputeGeometry computes properties in Mesh or the space
*  @{
*/
    inline
    real3 unit_vector(const real3 &v)
    {
        real norm = sqrt(vdot(v,v));
        real3 v_unit = v;
        v_unit.x/=norm;
        v_unit.y/=norm;
        v_unit.z/=norm;
        return(v_unit);
    }
  
    inline 
    real3 vector_cross(const real3 &v1, const real3 &v2)
    {
        real3 v;
        vcross(v, v1, v2);
        return v;
    }

    inline 
    real3 vector_sum(const real3 &v1, const real3 &v2)
    {
        real3 v;
        v.x = v1.x + v2.x;
        v.y = v1.y + v2.y;
        v.z = v1.z + v2.z;
        return v;
    }

    inline 
    real3 vector_subtract(const real3 &v1, const real3 &v2)
    {
        real3 v;
        v.x = v1.x - v2.x;
        v.y = v1.y - v2.y;
        v.z = v1.z - v2.z;
        return v;
    }

    inline 
    real3 vector_subtract(const real3 &v1, const real3 &v2, const BoxType &box)
    {
        return (host::minimum_image(v2, v1, box));
        /*real3 v;
        v.x = v1.x - v2.x;
        v.y = v1.y - v2.y;
        v.z = v1.z - v2.z;
        return v;*/
    }

    inline 
    void compute_form_factor_triangle(real *_metric, const real3 &r1, const real3 &r2, const real3 &r3, const BoxType &box)
    {
        real3 v12, v13;
        v12 = host::minimum_image(r1, r2, box);
        v13 = host::minimum_image(r1, r3, box);
        _metric[0] = vdot(v12, v12);
        _metric[1] = vdot(v12, v13);
        _metric[2] = vdot(v13, v13);
    }

    inline 
    void compute_form_factor_triangle(real *_metric, const real3 &r1, const real3 &r2, const real3 &r3)
    {
        real3 v12, v13;
        vsub(v12, r2, r1);
        vsub(v13, r3, r1);
        _metric[0] = vdot(v12, v12);
        _metric[1] = vdot(v12, v13);
        _metric[2] = vdot(v13, v13);
    }

    inline 
    real3 compute_normal_triangle(const real3 &r1, const real3 &r2, const real3 &r3, const BoxType &box)
    {
        real3 normal, v12, v13;
        v12 = host::minimum_image(r1, r2, box);
        v13 = host::minimum_image(r1, r3, box);
        vcross(normal, v12, v13);
        return normal;
    }

    inline 
    real3 compute_normal_triangle(const real3 &r1, const real3 &r2, const real3 &r3)
    {
        real3 normal, v12, v13;
        vsub(v12, r2, r1);
        vsub(v13, r3, r1);
        vcross(normal, v12, v13);
        return normal;
    }

    inline 
    real3 compute_normal_triangle_unit(const real3 &r1, const real3 &r2, const real3 &r3, const BoxType &box)
    {
        return (unit_vector(compute_normal_triangle(r1, r2, r3, box)));
    }

    inline 
    real3 compute_normal_triangle_unit(const real3 &r1, const real3 &r2, const real3 &r3)
    {
        return (unit_vector(compute_normal_triangle(r1, r2, r3)));
    }

    inline 
    real compute_area_triangle_from_vertex(const real3 &r1, const real3 &r2, const real3 &r3, const BoxType &box)
    {
        real3 normal = compute_normal_triangle(r1, r2, r3, box);
        return (0.5 * sqrt(vdot(normal, normal)));
    }

    inline 
    real compute_area_triangle_from_vertex(const real3 &r1, const real3 &r2, const real3 &r3)
    {
        real3 normal = compute_normal_triangle(r1, r2, r3);
        return (0.5 * sqrt(vdot(normal, normal)));
    }

    inline 
    real compute_area_triangle_from_metric(const real *_metric)
    {
        return (0.5 * sqrt(_metric[0] * _metric[2] - _metric[1] * _metric[1]));
    }

    inline 
    real compute_angle_vertex(const real3 &r1, const real3 &r2, const real3 &r3, const BoxType &box)
    {
        real3 v12, v13;
        v12 = host::minimum_image(r1, r2, box);
        v13 = host::minimum_image(r1, r3, box);
        real angle = acos(vdot(v12, v13) / sqrt(vdot(v12, v12) * vdot(v13, v13)));
        return angle;
    }

    inline 
    real compute_angle_vertex(const real3 &r1, const real3 &r2, const real3 &r3)
    {
        real3 v12, v13;
        vsub(v12, r2, r1);
        vsub(v13, r3, r1);
        real angle = acos(vdot(v12, v13) / sqrt(vdot(v12, v12) * vdot(v13, v13)));
        return angle;
    }

    inline 
    void compute_matrix_F(real *F, const real *g_reference_inv, const real *g_now)
    {

        /*
         _                                         _     _                     _     _      _
        |   g_reference_inv[0]  g_reference_inv[1]  |   |   g_now[0]  g_now[1]  |   |  1  0  |
    F = |                                           | x |                       | - |        |
        |_  g_reference_inv[1]  g_reference_inv[2] _|   |_  g_now[1]  g_now[2] _|   |_ 0  1 _|
    */
        F[0] = g_reference_inv[0] * g_now[0] + g_reference_inv[1] * g_now[1] - 1.0; //F11
        F[1] = g_reference_inv[0] * g_now[1] + g_reference_inv[1] * g_now[2];       //F12
        F[2] = g_reference_inv[1] * g_now[0] + g_reference_inv[2] * g_now[1];       //F21
        F[3] = g_reference_inv[1] * g_now[1] + g_reference_inv[2] * g_now[2] - 1.0; //F22
    }

    inline 
    real3 cmassT(const real3 &r1, const real3 &r2, const real3 &r3)
    {
        real3 vcm;
        vcm.x = (r1.x + r2.x + r3.x) / 3.0;
        vcm.y = (r1.y + r2.y + r3.y) / 3.0;
        vcm.z = (r1.z + r2.z + r3.z) / 3.0;
        return vcm;
    }

    inline 
    void RefMatrixFromCartesian(const real theta, const real *__restrict__ grefCart, real *grefCylin)
    {
        real a = sin(theta);
        real b = cos(theta);
        real grr = b * (b * grefCart[0] + a * grefCart[1]) + a * (b * grefCart[1] + a * grefCart[2]); //grr
        real grt = b * (b * grefCart[1] + a * grefCart[2]) - a * (b * grefCart[0] + a * grefCart[1]); //grt
        real gtt = b * (b * grefCart[2] - a * grefCart[1]) - a * (b * grefCart[1] - a * grefCart[0]); //gtt
        grefCylin[0] = grr;
        grefCylin[1] = grt;
        grefCylin[2] = gtt;
    }

    inline 
    void RefMatrixFromCylindrical(const real theta, const real *__restrict__ grefCylin, real *grefCart)
    {
        real a = sin(theta);
        real b = cos(theta);
        real g11 = b * (b * grefCylin[0] - a * grefCylin[1]) - a * (b * grefCylin[1] - a * grefCylin[2]);
        real g12 = a * (b * grefCylin[0] - a * grefCylin[1]) + b * (b * grefCylin[1] - a * grefCylin[2]);
        //real g21 = b*(a*grefCylin[0] + b*grefCylin[1]) - a*(a*grefCylin[1] + b*grefCylin[2]);
        real g22 = a * (a * grefCylin[0] + b * grefCylin[1]) + b * (a * grefCylin[1] + b * grefCylin[2]);
        grefCart[0] = g11;
        grefCart[1] = g12;
        grefCart[2] = g22;
    }
    
    /*! @} */
} // namespace host
#endif
/*! @} */
