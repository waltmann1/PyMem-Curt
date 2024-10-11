/************************************************************************************
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
#ifndef __globaltypes_hpp__
#define __globaltypes_hpp__

#include <iostream>
#include <sstream>
#include <math.h>

/** @} */
#define BIG_ENERGY_LIMIT 1e15 //!< Effectively acts as infinity
#define defPI 3.141592653589793
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

/** @brief real type **/
using real = double;

/** @brief 2D xy type **/
template <typename T>
struct xyType
{
    T x, y, z, w;
};
//using real3 = double4;
using real3 = xyType<real>;
using inth3 = xyType<int>;
using bool3 = xyType<bool>;

template <typename T>
struct TensorType
{
    T xx, xy, xz, yx, yy, yz, zx, zy, zz;
};
using realTensor = TensorType<real>;

/** @brief Vector addition. */
#define vset(v, val) \
    (v.x = val),     \
        (v.y = val), \
        (v.z = val)

/** @brief Vector addition. */
#define vsum(v, v1, v2)      \
    (v.x = v1.x + v2.x),     \
        (v.y = v1.y + v2.y), \
        (v.z = v1.z + v2.z)

/** @brief Vector subtraction. */
#define vsub(v, v1, v2)      \
    (v.x = v1.x - v2.x),     \
        (v.y = v1.y - v2.y), \
        (v.z = v1.z - v2.z)

/** @brief Vector dot product. */
#define vdot(v1, v2) (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z)

/** @brief Vector cross product. */
#define vcross(v, v1, v2)                  \
    (v.x = v1.y * v2.z - v1.z * v2.y),     \
        (v.y = v1.z * v2.x - v1.x * v2.z), \
        (v.z = v1.x * v2.y - v1.y * v2.x)

/** brief Constant times a vector **/
#define aXvec(a, v)      \
    (v.x = (a)*v.x),     \
        (v.y = (a)*v.y), \
        (v.z = (a)*v.z)


/**
 * @brief Force matrix type used for edge to vertex force
 */
typedef struct
{
    real3 forceM11, forceM12, forceM13;
} forceMatrix;

#endif
