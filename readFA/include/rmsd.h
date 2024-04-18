/*
    rmsd.h
    Author: Zhirong Liu [LiuZhiRong@pku.edu.cn]
            Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2014/06/12
    Description: Functions to calculate RMSD. Rewritten in C++.
*/

#ifndef SRC_CLIB_RMSD_H_
#define SRC_CLIB_RMSD_H_

#include "real.h"
#include "rvec.h"
#include "rmat.h"
#include "FAdefines.h"

extern "C" {
#include "NRecipes-UT.h"
}

int EigenRvec3(rmat3& a, rvec3& d, rmat3& z);

int EigenRvec4(rmat4 &a, rvec4 &d, rvec4 &e, rmat4 &z);

// Get matrix F from R
void GetFfR(rmat3 &RR, rmat4 &FF);

// Get matrix U from q
void GetUfq(real q0, real q1, real q2, real q3, rmat3 &U);

/* Move r1 to r2 by rmsd, save at rm
return -1 for error and rm = r1
*/
int MoveRMSD(VecVec3 &r1, VecVec3 &r2, VecVec3 &rm);

/* Calculate the rmsd between r1 and r2
   input: r1, r2, 
          calc_U (whether to calculate the rotation matrix U)
          calc_g (whether to calculate the gradient, g, of rmsd to r1
   output: r1center, r2center (barycenter)
           rmsd,  =(1/n)(|U.r1-r2|^2)^1/2
           lambda[4], q[4][4] (all eingen values and vectors of the quaternion)
                            note that the i_th eigenvector is q[][i]
           U (rotation matrix), g (rmsd gradient to r1)
   return: No. of eigenvalue, lambda, to calculate the rmsd
           -1 if error
*/
/* reference:
[1] Evangelos A. Coutsias, Chaok Seok, and Ken A. Dill.
  Using quaternions to calculate RMSD.
  J. Comput. Chem. {\bf 25}, 1849-1857 (2004).
[2] Simon K. Kearsley.
  On the orthogonal transformation used for structural comparisons.
  Acta Cryst. A {\bf 45}, 208-210 (1989).
[3] my own worm-tech.tex and worm-tech.ps
*/
int CalRMSD(real &rmsd, VecVec3 &r1, VecVec3 &r2,
            rvec3 &r1center, rvec3 &r2center,
            rvec4 &lambda, rmat4 &q,
            bool calc_U, rmat3 &U,
            bool calc_g, VecVec3 &g);

#endif /* SRC_CLIB_RMSD_H_ */