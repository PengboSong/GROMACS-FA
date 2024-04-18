/*
    real.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2021/08/19
    Description: Definition of real.
*/

#ifndef HAS_REAL
#define HAS_REAL

#ifdef DOUBLE_PRECISION
using real = double;
#else
using real = float;
#endif

#endif /* HAS_REAL */
