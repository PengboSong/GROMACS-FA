/*
    NRecipes-UT.h
    Author: Zhirong Liu [LiuZhiRong@pku.edu.cn]
            Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2007/06/27
    Description: Functions to calculate RMSD. Rewritten in C++.
*/

#ifndef SRC_CLIB_NRECIPESUT_H_
#define SRC_CLIB_NRECIPESUT_H_

#include <math.h>
#include <malloc.h>

int myEigen(double *a, int n, double d[], double e[], double *z); // Eigen problems

double pythag(double a, double b);

int tqli(double d[], double e[], int n, double *z);

void tred2(double *a, int n, double d[], double e[]);

int ludcmp(double *a, int n, int *indx, double *d);

void lubksb(double *a, int n, int *indx, double *b, int nb);

int mySolLinear(double *a, int n, double *b, int nb); // solve linear equation AX=b with LU decomposition

double myDet(double *a, int n);                       // calculate the determinant of the matrix

#endif /* SRC_CLIB_NRECIPESUT_H_ */