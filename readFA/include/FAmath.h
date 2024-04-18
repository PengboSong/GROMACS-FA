/*
    FAmath.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2021/08/17
    Description: Useful math constants and functions.
*/

#ifndef SRC_READFA_FAMATH_H_
#define SRC_READFA_FAMATH_H_

#include "real.h"
#include "rvec.h"
#include "rmat.h"

real norm2(real fx, real fy, real fz);

real norm(real fx, real fy, real fz);

real vecnorm(rvec3& v);
real vecnorm(rvec4& v);

// Calculate dot product of two vectors
real vecdot(const rvec3& vi, const rvec3& vj);
real vecdot(const rvec3& vi, const rvec4& vj);
real vecdot(const rvec4& vi, const rvec3& vj);
real vecdot(const rvec4& vi, const rvec4& vj);

// Calculate cross product of two vectors
rvec3 veccross(const rvec3& vi, const rvec3& vj);
rvec3 veccross(const rvec4& vi, const rvec4& vj);

// Calculate distance between two vectors
real vecdist(const rvec3& vi, const rvec3& vj);
real vecdist(const rvec4& vi, const rvec4& vj);

// Calculate matrix multiply result of two vectors
rmat3 matmul(const rvec3& vi, const rvec3& vj);
rmat3 matmul(const rvec3& vi, const rvec4& vj);
rmat3 matmul(const rvec4& vi, const rvec3& vj);
rmat3 matmul(const rvec4& vi, const rvec4& vj);

rvec3 matmul(const rmat3& m, const rvec3& v);

// 
rvec3 matdiag(const rmat3& m);
rvec4 matdiag(const rmat4& m);

// 
real trace(const rmat3& m);
real trace(const rmat4& m);

void initvec(VecVec3& vec, ulong len);
void initvec(VecVec4& vec, ulong len);
void initvec(VecMat3& vec, ulong len);
void initvec(VecReal& vec, ulong len, real c = 0.0f);
void initmat(DVecVec3& mat, ulong nrow, ulong ncol);
void initmat(DVecVec4& mat, ulong nrow, ulong ncol);
void initmat(DVecMat3& mat, ulong nrow, ulong ncol);
void initmat(DVecReal& mat, ulong nrow, ulong ncol, real c = 0.0f);

void clearvec(VecVec3& vec);
void clearmat(DVecVec3& mat);

double* rvec2d(const rvec3& v);
double* rvec2d(const rvec4& v);
double* rmat2d(const rmat3& m);
double* rmat2d(const rmat4& m);
void d2rvec(double* dv, rvec3& v);
void d2rvec(double* dv, rvec4& v);
void d2rmat(double* dm, rmat3& m);
void d2rmat(double* dm, rmat4& m);

// Rotate r2 to fit r1: r2 = U(r1 - r1_center) + r2_center
void rotate(const rmat3& U, const rvec3& r1center, const rvec3& r2center, const VecVec3& r1, VecVec3& rr);
void rotate(const rmat3& U, const rvec3& r1center, const rvec3& r2center, const rvec4& r1, rvec4& rr);

// Force rotation without displacement: f2 = U f1
void rotate(const rmat3& U, const rvec4& f0, rvec4& f);
void rotate(const rmat3& U, rvec4& f);

void rot_force(const rmat3& U, VecVec4& Fi);
void rot_force(const rmat3& U, DMapVec4& Fij);

// Fit cosine value to range [-1.0, 1.0]
real fitcos(real cosx);

real rad2deg(real radian);
real deg2rad(real degree);

#endif /* SRC_READFA_FAMATH_H_ */
