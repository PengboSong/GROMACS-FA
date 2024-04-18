/*
    rvecTest.cpp
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2021/08/23
    Description: Test functionalities of rvec and rmat.
*/

#include <iostream>
#include <iomanip>

#include "real.h"
#include "rvec.h"
#include "rmat.h"
#include "rmsd.h"
#include "FAmath.h"
#include "FAdefines.h"
#include "FAutility.h"

constexpr bool real_equal(const real r1, const real r2)
{
    return std::abs(r1 - r2) < 1e-3;
}

bool check_rvec3(const rvec3 v, const real x, const real y, const real z)
{
    return real_equal(v[XX], x) && real_equal(v[YY], y) && real_equal(v[ZZ], z);
}

bool check_rvec4(const rvec4 v, const real x, const real y, const real z, const real xyz)
{
    return real_equal(v[XX], x) && real_equal(v[YY], y) && real_equal(v[ZZ], z) && real_equal(v[XYZ], xyz);
}

bool check_rmat3(const rmat3 m,
                 const real xx, const real xy, const real xz,
                 const real yx, const real yy, const real yz,
                 const real zx, const real zy, const real zz)
{
    return real_equal(m[XX][XX], xx) && real_equal(m[XX][YY], xy) && real_equal(m[XX][ZZ], xz) && \
           real_equal(m[YY][XX], yx) && real_equal(m[YY][YY], yy) && real_equal(m[YY][ZZ], yz) && \
           real_equal(m[ZZ][XX], zx) && real_equal(m[ZZ][YY], zy) && real_equal(m[ZZ][ZZ], zz);
}

void perform_test(bool test, const int id, const std::string desc)
{
    std::string status = test ? "Success" : "Failed";
    std::cout << "Test " << id << " (" << desc << "): " << status << std::endl;
}

int main(int argc, char **argv)
{
    int tid = 0;

    rvec<9> nums = {9.537, 5.701, 6.564, 9.857, 3.061, 9.184, 9.462, 0.577, 5.098};

    real a = 1.648, b = 8.138, c = 2.411;
    rvec3 va, vb, vc, vd;
    rmat3 m;

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            m[i][j] = nums[i * 3 + j];

    for (int i = 0; i < 3; ++i)
    {
        va[i] = nums[i];
        vb[i] = nums[i + 3];
        vc[i] = nums[i + 6];
    }

    perform_test(check_rvec3(va,           9.537,  5.701,  6.564), ++tid, "[rvec3] value assignment va");
    perform_test(check_rvec3(vb,           9.857,  3.061,  9.184), ++tid, "[rvec3] value assignment vb");
    perform_test(check_rvec3(vc,           9.462,  0.577,  5.098), ++tid, "[rvec3] value assignment vc");
    perform_test(check_rvec3(va + vb,     19.394,  8.762, 15.748), ++tid, "[rvec3] addition va + vb");
    perform_test(check_rvec3(va - vb,     -0.320,  2.640, -2.620), ++tid, "[rvec3] subtraction va - vb");
    perform_test(check_rvec3(vc * c,      22.812,  1.391, 12.291), ++tid, "[rvec3] scalar-multiplication vc * c");
    perform_test(check_rvec3(va * vc,     90.239,  3.289, 33.463), ++tid, "[rvec3] vector-multiplication va * vc");
    perform_test(check_rvec3(vb * vc,     93.267,  1.766, 46.820), ++tid, "[rvec3] vector-multiplication vb * vc");

    vd = {4.255, 0.943, 3.793};
    perform_test(check_rvec3(vd,           4.255,  0.943,  3.793), ++tid, "[rvec3] value assignment vd");
    vd += va;
    perform_test(check_rvec3(vd,          13.792,  6.644, 10.357), ++tid, "[rvec3] increment vd += va");
    vd -= vb;
    perform_test(check_rvec3(vd,           3.935,  3.583,  1.173), ++tid, "[rvec3] decrement vd -= vb");
    vd *= c;
    perform_test(check_rvec3(vd,           9.487,  8.639,  2.828), ++tid, "[rvec3] scalar-multiplication vd *= c");
    perform_test(check_rvec3(vd.square(), 90.008, 74.626,  7.998), ++tid, "[rvec3] vector square vd ** 2");
    perform_test(check_rvec3(vd.sqrt(),    3.080,  2.939,  1.682), ++tid, "[rvec3] vector sqrt vd ** .5");

    perform_test(real_equal(vd.norm2(), 172.632), ++tid, "[rvec3] vector norm2 ||vd||2");

    perform_test(check_rmat3(m,      9.537,  5.701,  6.564,  9.857,  3.061,  9.184,  9.462,  0.577,  5.098), ++tid, "[rmat3] assign values to m");
    m += a;
    perform_test(check_rmat3(m,     11.185,  7.349,  8.212, 11.505,  4.709, 10.832, 11.110,  2.225,  6.746), ++tid, "[rmat3] increment m += a");
    m -= b;
    perform_test(check_rmat3(m,      3.047, -0.789,  0.074,  3.367, -3.429,  2.694,  2.972, -5.913, -1.392), ++tid, "[rmat3] decrement m -= b");
    m *= c;
    perform_test(check_rmat3(m,      7.346, -1.902,  0.178,  8.118, -8.267,  6.495,  7.165,-14.256, -3.356), ++tid, "[rmat3] scalar-multiplication m *= c");
    perform_test(check_rmat3(m * c, 17.712, -4.586,  0.430, 19.572,-19.932, 15.660, 17.276,-34.372, -8.092), ++tid, "[rmat3] scalar-multiplication m * c");
    
    rvec4 ve = {0.321, 0.640, 0.521};
    perform_test(check_rvec4(ve,           0.321,  0.640,  0.521,  0.000), ++tid, "[rvec4] value assignment ve");
    ve += {0.154, 0.649, 0.405};
    perform_test(check_rvec4(ve,           0.475,  1.289,  0.926,  0.000), ++tid, "[rvec4] increment ve");
    ve -= {0.737, 0.401, 0.886};
    perform_test(check_rvec4(ve,          -0.262,  0.888,  0.040,  0.000), ++tid, "[rvec4] decrement ve");

    const rvec4 vf = {0.752, 0.743, 0.707, 1.272};
    rvec4 vg = vf;
    perform_test(check_rvec4(vg,           0.752,  0.743,  0.707,  1.272), ++tid, "[rvec4] value assignment vg");

    perform_test(real_equal(vecdot(va, vb), 171.741), ++tid, "[rvec3] vector dot product va . vb");
    perform_test(real_equal(vecdot(ve, vc),  -1.763), ++tid, "[rvec3] vector dot product ve . vc");
    perform_test(real_equal(vecdot(ve, vf),   0.491), ++tid, "[rvec3] vector dot product ve . vf");
    perform_test(check_rvec3(veccross(va, vb), 32.266,-22.886,-27.002), ++tid, "[rvec3] value cross product va x vb");
    perform_test(check_rvec3(veccross(ve, vf),  0.598,  0.215, -0.862), ++tid, "[rvec3] value cross product ve x vf");
    perform_test(real_equal(vecdist(va, vb),   3.733), ++tid, "[rvec3] vector distance |va - vb|");
    perform_test(real_equal(vecdist(ve, vf),   1.222), ++tid, "[rvec3] vector distance |ve - vf|");
    perform_test(check_rmat3(matmul(va, vb),   94.006, 29.193, 87.588, 56.195, 17.451, 52.358, 64.701, 20.092, 60.284), ++tid, "[rmat3] matrix multiplication va * vb");
    perform_test(check_rmat3(matmul(va, ve),   -2.499,  8.469,  0.381, -1.494,  5.062,  0.228, -1.720,  5.829,  0.263), ++tid, "[rmat3] matrix multiplication va * ve");
    perform_test(check_rmat3(matmul(ve, vf),   -0.197, -0.195, -0.185,  0.668,  0.660,  0.628,  0.030,  0.030,  0.028), ++tid, "[rmat3] matrix multiplication ve * vf");
    perform_test(check_rvec3(matmul(m, va), 60.388, 72.922,-34.967), ++tid, "[rvec3] matrix multiplication m * va");
    perform_test(check_rvec3(matdiag(m),     7.346, -8.267, -3.356), ++tid, "[rvec3] matrix diagonal diag(m)");

    double *da = rvec2d(va), *de = rvec2d(ve);
    rvec3 dva;
    rvec4 dve, rve;
    d2rvec(da, dva);
    d2rvec(de, dve);
    perform_test(check_rvec3(dva,  9.537,  5.701,  6.564), ++tid, "[rvec3] c array conversion va");
    perform_test(check_rvec4(dve, -0.262,  0.888,  0.040,  0.000), ++tid, "[rvec4] c array conversion ve");

    rmat4 AA, zz, dAA;
    rvec4 dd, ee;
    AA[0] = {8.696, 4.104, 2.739, 4.115};
    AA[1] = {4.104, 7.145, 5.172, 4.985};
    AA[2] = {2.739, 5.172, 3.276, 1.540};
    AA[3] = {4.115, 4.985, 1.540, 5.186};
    EigenRvec4(AA, dd, ee, zz);
    std::cout << "Input 4x4 matrix = [" << AA << "]" << std::endl;
    std::cout << "Eigenvalues & Eigenvectors = [" << dd << zz << "]" << std::endl;

    m[0] = { 9.528, -16.518, -23.355};
    m[1] = {-1.216, -18.435, -59.63 };
    m[2] = { 2.337,  -7.133, -17.566};
    std::cout << "Input 3x3 matrix = [" << m << "]" << std::endl;
    std::cout << "Eigenvalues & Eigenvectors = [";
    outmatpca(std::cout, m);
    std::cout << "]" << std::endl;
    m[0] = {12.822, 4.187, 0.543};
    m[1] = { 4.187, 1.367, 0.177};
    m[2] = { 0.543, 0.177, 0.023};
    std::cout << "Input 3x3 matrix = [" << m << "]" << std::endl;
    std::cout << "Eigenvalues & Eigenvectors = [";
    outmatpca(std::cout, m);
    std::cout << "]" << std::endl;

    return 0;
}
