/*
    FAmath.cpp
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2021/08/17
    Description: Useful math constants and functions.
*/

#include <cmath>

#include "FAdefines.h"
#include "FAmath.h"

real norm2(real x, real y, real z)
{
    return x * x + y * y + z * z;
}

real norm(real x, real y, real z)
{
    return std::sqrt(x * x + y * y + z * z);
}

real vecnorm(rvec3& v)
{
    return norm(v[XX], v[YY], v[ZZ]);
}

real vecnorm(rvec4& v)
{
    v[XYZ] = norm(v[XX], v[YY], v[ZZ]);
    return v[XYZ];
}

real vecdot(const rvec3& vi, const rvec3& vj)
{
    return vi[XX] * vj[XX] + vi[YY] * vj[YY] + vi[ZZ] * vj[ZZ];
}

real vecdot(const rvec3& vi, const rvec4& vj)
{
    return vi[XX] * vj[XX] + vi[YY] * vj[YY] + vi[ZZ] * vj[ZZ];
}

real vecdot(const rvec4& vi, const rvec3& vj)
{
    return vi[XX] * vj[XX] + vi[YY] * vj[YY] + vi[ZZ] * vj[ZZ];
}

real vecdot(const rvec4& vi, const rvec4& vj)
{
    return vi[XX] * vj[XX] + vi[YY] * vj[YY] + vi[ZZ] * vj[ZZ];
}

rvec3 veccross(const rvec3& vi, const rvec3& vj)
{
    rvec3 v;
    v[XX] = vi[YY] * vj[ZZ] - vi[ZZ] * vj[YY];
    v[YY] = vi[ZZ] * vj[XX] - vi[XX] * vj[ZZ];
    v[ZZ] = vi[XX] * vj[YY] - vi[YY] * vj[XX];
    return v;
}

rvec3 veccross(const rvec4& vi, const rvec4& vj)
{
    rvec3 v;
    v[XX] = vi[YY] * vj[ZZ] - vi[ZZ] * vj[YY];
    v[YY] = vi[ZZ] * vj[XX] - vi[XX] * vj[ZZ];
    v[ZZ] = vi[XX] * vj[YY] - vi[YY] * vj[XX];
    return v;
}

real vecdist(const rvec3& vi, const rvec3& vj)
{
    return norm(vi[XX] - vj[XX], vi[YY] - vj[YY], vi[ZZ] - vj[ZZ]);
}

real vecdist(const rvec4& vi, const rvec4& vj)
{
    return norm(vi[XX] - vj[XX], vi[YY] - vj[YY], vi[ZZ] - vj[ZZ]);
}

rmat3 matmul(const rvec3& vi, const rvec3& vj)
{
    rmat3 m;
    for (int i = 0; i < DIM; ++i)
        for (int j = 0; j < DIM; ++j)
            m[i][j] = vi[i] * vj[j];
    return m;
}

rmat3 matmul(const rvec3& vi, const rvec4& vj)
{
    rmat3 m;
    for (int i = 0; i < DIM; ++i)
        for (int j = 0; j < DIM; ++j)
            m[i][j] = vi[i] * vj[j];
    return m;
}

rmat3 matmul(const rvec4 &vi, const rvec3 &vj)
{
    rmat3 m;
    for (int i = 0; i < DIM; ++i)
        for (int j = 0; j < DIM; ++j)
            m[i][j] = vi[i] * vj[j];
    return m;
}

rmat3 matmul(const rvec4& vi, const rvec4& vj)
{
    rmat3 m;
    for (int i = 0; i < DIM; ++i)
        for (int j = 0; j < DIM; ++j)
            m[i][j] = vi[i] * vj[j];
    return m;
}

rvec3 matmul(const rmat3& m, const rvec3& v)
{
    rvec3 r;
    for (int i = 0; i < DIM; ++i)
        for (int j = 0; j < DIM; ++j)
            r[i] += m[i][j] * v[j];
    return r;
}

rvec3 matdiag(const rmat3& m)
{
    rvec3 diag;
    for (int i = 0; i < 3; ++i)
        diag[i] = m[i][i];
    return diag;
}

rvec4 matdiag(const rmat4& m)
{
    rvec4 diag;
    for (int i = 0; i < 4; ++i)
        diag[i] = m[i][i];
    return diag;
}

real trace(const rmat3 &m)
{
    return m[XX][XX] + m[YY][YY] + m[ZZ][ZZ];
}

real trace(const rmat4 &m)
{
    return m[XX][XX] + m[YY][YY] + m[ZZ][ZZ] + m[XYZ][XYZ];
}

void initvec(VecVec3 &vec, ulong len)
{
    vec.resize(len, rvec3());
}

void initvec(VecVec4 &vec, ulong len)
{
    vec.resize(len, rvec4());
}

void initvec(VecMat3 &vec, ulong len)
{
    vec.resize(len, rmat3());
}

void initvec(VecReal &vec, ulong len, real c)
{
    vec.resize(len, c);
}

void initmat(DVecVec3 &mat, ulong nrow, ulong ncol)
{
    mat.resize(nrow, VecVec3(ncol, rvec3()));
}

void initmat(DVecVec4 &mat, ulong nrow, ulong ncol)
{
    mat.resize(nrow, VecVec4(ncol, rvec4()));
}

void initmat(DVecMat3 &mat, ulong nrow, ulong ncol)
{
    mat.resize(nrow, VecMat3(ncol, rmat3()));
}

void initmat(DVecReal &mat, ulong nrow, ulong ncol, real c)
{
    mat.resize(nrow, VecReal(ncol, c));
}

void clearvec(VecVec3& vec)
{
    for (ulong i = 0; i < vec.size(); ++i)
        vec[i] = rvec3();
}

void clearmat(DVecVec3& mat)
{
    for (ulong i = 0; i < mat.size(); ++i)
        clearvec(mat[i]);
}

double* rvec2d(const rvec3& v)
{
    double* dv = new double[3];
    for (int i = 0; i < 3; ++i)
        dv[i] = v[i];
    return dv;
}

double* rvec2d(const rvec4& v)
{
    double* dv = new double[4];
    for (int i = 0; i < 4; ++i)
        dv[i] = v[i];
    return dv;
}

double *rmat2d(const rmat3 &m)
{
    double* dm = new double[3 * 3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            dm[i * 3 + j] = m[i][j];
    return dm;
}

double* rmat2d(const rmat4& m)
{
    double* dm = new double[4 * 4];
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            dm[i * 4 + j] = m[i][j];
    return dm;
}

void d2rvec(double* dv, rvec3& v)
{
    for (int i = 0; i < 3; ++i)
        v[i] = dv[i];
    delete[] dv;
}

void d2rvec(double* dv, rvec4& v)
{
    for (int i = 0; i < 4; ++i)
        v[i] = dv[i];
    delete[] dv;
}

void d2rmat(double *dm, rmat3 &m)
{
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            m[i][j] = dm[i * 3 + j];
    delete[] dm;
}

void d2rmat(double* dm, rmat4& m)
{
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            m[i][j] = dm[i * 4 + j];
    delete[] dm;
}

void rotate(const rmat3& U, const rvec3& r1center, const rvec3& r2center, const VecVec3& r1, VecVec3& rr)
{
    ulong i = 0, rN = r1.size();
    rr.resize(rN, rvec3());
    for (; i < rN; ++i)
        rr[i] = r2center + matmul(U, r1[i] - r1center);
}

void rotate(const rmat3& U, const rvec3& r1center, const rvec3& r2center, const rvec4& r1, rvec4& rr)
{
    rvec3 r1xyz, rrxyz;
    r1xyz = {r1[XX], r1[YY], r1[ZZ]};
    rrxyz = r2center + matmul(U, r1xyz - r1center);
    rr = {rrxyz[XX], rrxyz[YY], rrxyz[ZZ], r1[XYZ]};
}

void rotate(const rmat3 &U, const rvec4 &f0, rvec4 &f)
{
    f[XX] = U[XX][XX] * f0[XX] + U[XX][YY] * f0[YY] + U[XX][ZZ] * f0[ZZ];
    f[YY] = U[YY][XX] * f0[XX] + U[YY][YY] * f0[YY] + U[YY][ZZ] * f0[ZZ];
    f[ZZ] = U[ZZ][XX] * f0[XX] + U[ZZ][YY] * f0[YY] + U[ZZ][ZZ] * f0[ZZ];
}

void rotate(const rmat3 &U, rvec4 &f)
{
    real fx, fy, fz;
    fx = U[XX][XX] * f[XX] + U[XX][YY] * f[YY] + U[XX][ZZ] * f[ZZ];
    fy = U[YY][XX] * f[XX] + U[YY][YY] * f[YY] + U[YY][ZZ] * f[ZZ];
    fz = U[ZZ][XX] * f[XX] + U[ZZ][YY] * f[YY] + U[ZZ][ZZ] * f[ZZ];
    f[XX] = fx;
    f[YY] = fy;
    f[ZZ] = fz;
}


void rot_force(const rmat3 &U, VecVec4 &Fi)
{
    for (auto fi = Fi.begin(); fi != Fi.end(); ++fi)
        rotate(U, *fi);
}

void rot_force(const rmat3 &U, DMapVec4 &Fij)
{
    for (auto fi = Fij.begin(); fi != Fij.end(); ++fi)
        for (auto fij = fi->second.begin(); fij != fi->second.end(); ++fij)
            rotate(U, Fij[fi->first][fij->first]);
}

real fitcos(real cosx)
{
    return std::max<real>(std::min<real>(cosx, 1.0), -1.0);
}

real rad2deg(real radian)
{
    return radian / pi * 180;
}

real deg2rad(real degree)
{
    return degree * pi / 180;
}
