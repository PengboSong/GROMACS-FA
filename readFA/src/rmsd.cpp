/*
    rmsd.cpp
    Author: Zhirong Liu [LiuZhiRong@pku.edu.cn]
            Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2014/06/12
    Description: Functions to calculate RMSD. Rewritten in C++.
*/

#include <cmath>

#include "rmsd.h"
#include "FAmath.h"

int EigenRvec3(rmat3 &a, rvec3 &d, rmat3 &z)
{
    rvec3 e;
    double *da = rmat2d(a), *dd = rvec2d(d), *de = rvec2d(e), *dz = rmat2d(z);
    int v = myEigen(da, 3, dd, de, dz);
    d2rmat(da, a);
    d2rvec(dd, d);
    d2rvec(de, e);
    d2rmat(dz, z);
    return v;
}

int EigenRvec4(rmat4 &a, rvec4 &d, rvec4 &e, rmat4 &z)
{
    double *da = rmat2d(a), *dd = rvec2d(d), *de = rvec2d(e), *dz = rmat2d(z);
    int v = myEigen(da, 4, dd, de, dz);
    d2rmat(da, a);
    d2rvec(dd, d);
    d2rvec(de, e);
    d2rmat(dz, z);
    return v;
}

void GetFfR(rmat3 &RR, rmat4 &FF)
{
    FF[0][0] = RR[0][0] + RR[1][1] + RR[2][2];
    FF[0][1] = RR[1][2] - RR[2][1];
    FF[0][2] = RR[2][0] - RR[0][2];
    FF[0][3] = RR[0][1] - RR[1][0];

    FF[1][0] = FF[0][1];
    FF[1][1] = RR[0][0] - RR[1][1] - RR[2][2];
    FF[1][2] = RR[0][1] + RR[1][0];
    FF[1][3] = RR[0][2] + RR[2][0];

    FF[2][0] = FF[0][2];
    FF[2][1] = FF[1][2];
    FF[2][2] = -RR[0][0] + RR[1][1] - RR[2][2];
    FF[2][3] = RR[1][2] + RR[2][1];

    FF[3][0] = FF[0][3];
    FF[3][1] = FF[1][3];
    FF[3][2] = FF[2][3];
    FF[3][3] = -RR[0][0] - RR[1][1] + RR[2][2];
}

void GetUfq(real q0, real q1, real q2, real q3, rmat3 &U)
{
    U[0][0] = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3;
    U[0][1] = 2 * (q1 * q2 - q0 * q3);
    U[0][2] = 2 * (q1 * q3 + q0 * q2);

    U[1][0] = 2 * (q1 * q2 + q0 * q3);
    U[1][1] = q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3;
    U[1][2] = 2 * (q2 * q3 - q0 * q1);

    U[2][0] = 2 * (q1 * q3 - q0 * q2);
    U[2][1] = 2 * (q2 * q3 + q0 * q1);
    U[2][2] = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;
}

int MoveRMSD(VecVec3 &r1, VecVec3 &r2, VecVec3 &rm)
{
    real rmsd;
    rvec3 r1center, r2center, dR, rdR;
    rvec4 lambda;
    rmat3 U;
    rmat4 q;
    int atomn = r1.size();
    VecVec3 g(atomn, rvec3());

    if (CalRMSD(rmsd, r1, r2, r1center, r2center, lambda, q, true, U, false, g) == -1)
    {
        rm = r1;
        return -1;
    }

    rm.resize(atomn, rvec3());
    for (int ai = 0; ai < atomn; ++ai)
    {
        dR = r1[ai] - r1center;
        for (int i = 0; i < 3; ++i)
            rdR[i] = vecdot(U[i], dR);
        rm[ai] = r2center + rdR;
    }
    return 0;
}

int CalRMSD(real &rmsd, VecVec3 &r1, VecVec3 &r2,
            rvec3 &r1center, rvec3 &r2center,
            rvec4 &lambda, rmat4 &q,
            bool calc_U, rmat3 &U,
            bool calc_g, VecVec3 &g)
{
    real r1norm, r2norm, sqrmsd;
    real q0, q1, q2, q3;
    real revn, revrmsd;
    rmat3 RR;
    rmat4 FF;
    rvec4 ee;
    rvec3 r1c, r2c, Ur2c;
    int nL, atomn;

    atomn = r1.size();
    revn = 1.0 / atomn;
    r1center = {0., 0., 0.};
    r2center = {0., 0., 0.};
    for (int ai = 0; ai < atomn; ++ai)
    {
        r1center += r1[ai];
        r2center += r2[ai];
    }
    r1center *= revn;
    r2center *= revn;

    // norm and R = r1::r2 (centric)
    r1norm = 0.;
    r2norm = 0.;
    for (int ai = 0; ai < atomn; ++ai)
    {
        r1c = r1[ai] - r1center;
        r2c = r2[ai] - r2center;
        r1norm += r1c.norm2();
        r2norm += r2c.norm2();
        RR += matmul(r1c, r2c);
    }

    // F
    GetFfR(RR, FF);

    // Slove the eigen problem: Fq=lambda.q
    if (EigenRvec4(FF, lambda, ee, q) != 1) // fail
    {
        printf("Failed to solve the eigen problem.\n");
        return -1;
    }

    // rmsd
    nL = 0; // No. of lambda for rmsd
    for (int i = 1; i < 4; ++i)
        if (lambda[i] > lambda[nL])
            nL = i;
    sqrmsd = (r1norm + r2norm - 2. * lambda[nL]) * revn;
    rmsd = (sqrmsd > 0.) ? std::sqrt(sqrmsd) : 0.;
    revrmsd = (rmsd > 0.) ? 1. / rmsd : 0.;

    // U
    if (calc_U || calc_g)
    {
        q0 = q[0][nL];
        q1 = q[1][nL];
        q2 = q[2][nL];
        q3 = q[3][nL];

        GetUfq(q0, q1, q2, q3, U);
    }

    // g =(r1c-U'.r2c)/(n.rmsd)
    if (calc_g)
        for (int ai = 0; ai < atomn; ++ai)
        {
            r1c = r1[ai] - r1center;
            r2c = r2[ai] - r2center;
            if (rmsd > 0.)
            {
                // Ur2c[i] = Σ_j U[j][i] * r2c[j]
                Ur2c = {0., 0., 0.};
                for (int i = 0; i < 3; ++i)
                    Ur2c += U[i] * r2c[i];
                g[ai] = (r1c - Ur2c) * revn * revrmsd;
            }
            else
                g[ai] = {0., 0., 0.};
        }

    return nL;
}