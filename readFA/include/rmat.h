/*
    rmat.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2021/08/19
    Description: Real matrix containers.
*/

#ifndef SRC_READFA_RMAT_H_
#define SRC_READFA_RMAT_H_

#include <cstdint>
#include <iostream>
#include <map>
#include <vector>

#include "real.h"
#include "rvec.h"

template <ulong nrow, ulong ncol>
class rmat
{

    using row = rvec<ncol>;

public:
    rmat()
    {
        m.resize(nrow, row());
    }

    const row &operator[](const ulong i) const
    {
        if (i < nrow)
            return m[i];
        else
            throw std::runtime_error("rmat index out of range.");
    }

    row &operator[](const ulong i)
    {
        if (i < nrow)
            return m[i];
        else
            throw std::runtime_error("rmat index out of range.");
    }

    void operator+=(const real c)
    {
        for (ulong i = 0; i < nrow; ++i)
            m[i] += c;
    }

    void operator+=(const rmat<nrow, ncol> mat)
    {
        for (ulong i = 0; i < nrow; ++i)
            m[i] += mat[i];
    }

    void operator-=(const real c)
    {
        for (ulong i = 0; i < nrow; ++i)
            m[i] -= c;
    }

    void operator-=(const rmat<nrow, ncol> mat)
    {
        for (ulong i = 0; i < nrow; ++i)
            m[i] -= mat[i];
    }

    void operator*=(const real c)
    {
        for (ulong i = 0; i < nrow; ++i)
            m[i] *= c;
    }

    const rmat<nrow, ncol> operator+(const real c) const
    {
        rmat<nrow, ncol> nmat;
        for (ulong i = 0; i < nrow; ++i)
            nmat[i] = m[i] + c;
        return nmat;
    }

    const rmat<nrow, ncol> operator+(const rmat<nrow, ncol> mat) const
    {
        rmat<nrow, ncol> nmat;
        for (ulong i = 0; i < nrow; ++i)
            nmat[i] = m[i] + mat[i];
        return nmat;
    }

    const rmat<nrow, ncol> operator-(const real c) const
    {
        rmat<nrow, ncol> nmat;
        for (ulong i = 0; i < nrow; ++i)
            nmat[i] = m[i] - c;
        return nmat;
    }

    const rmat<nrow, ncol> operator-(const rmat<nrow, ncol> mat) const
    {
        rmat<nrow, ncol> nmat;
        for (ulong i = 0; i < nrow; ++i)
            nmat[i] = m[i] - mat[i];
        return nmat;
    }

    const rmat<nrow, ncol> operator*(const real c) const
    {
        rmat<nrow, ncol> nmat;
        for (ulong i = 0; i < nrow; ++i)
            nmat[i] = m[i] * c;
        return nmat;
    }

    const rmat<nrow, ncol> operator/(const rmat<nrow, ncol> mat) const
    {
        rmat<nrow, ncol> nmat;
        for (ulong i = 0; i < nrow; ++i)
            nmat[i] = m[i] / mat[i];
        return nmat;
    }

    friend std::ostream &operator<<(std::ostream &os, rmat<nrow, ncol> m)
    {
        for (ulong i = 0; i < nrow; ++i)
            for (ulong j = 0; j < ncol; ++j)
                os << ',' << m[i][j];
        return os;
    }

private:
    std::vector<row> m;
};

using rmat2 = rmat<2, 2>;                    // rmat2[2][2]
using rmat3 = rmat<3, 3>;                    // rmat3[3][3]
using rmat4 = rmat<4, 4>;                    // rmat4[4][4]
using VecMat3 = std::vector<rmat3>;          // VecMat3[i][3][3]
using DVecMat3 = std::vector<VecMat3>;       // DVecMat3[i][j][3][3]
using MapMat3 = std::map<int32_t, rmat3>;    // MapMat3[i][3][3]
using DMapMat3 = std::map<int32_t, MapMat3>; // DMapMat3[i][j][3][3]

using VecReal = std::vector<real>;           // VecReal[i]
using DVecReal = std::vector<VecReal>;       // DVecReal[i][j]
using MapReal = std::map<int32_t, real>;     // MapReal[i]
using DMapReal = std::map<int32_t, MapReal>; // DMapReal[i][j]

#endif /* SRC_READFA_RMAT_H_ */
