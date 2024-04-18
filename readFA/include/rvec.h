/*
    rvec.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2021/08/19
    Description: Real vector containers.
*/

#ifndef SRC_READFA_RVEC_H_
#define SRC_READFA_RVEC_H_

#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <iostream>
#include <map>
#include <vector>

#include "real.h"

template <ulong len>
class rvec
{
public:
    rvec()
    {
        v.resize(len, 0.);
    }

    rvec(const std::initializer_list<real> lst)
    {
        v.resize(len, 0.);
        const real *it = lst.begin();
        for (ulong i = 0; i < std::max<ulong>(len, lst.size());)
            v[i++] = *it++;
    }

    const real operator[](const ulong i) const
    {
        if (i < len)
            return v[i];
        else
            throw std::runtime_error("rvec index out of range.");
    }

    real &operator[](const ulong i)
    {
        if (i < len)
            return v[i];
        else
            throw std::runtime_error("rvec index out of range.");
    }

    void operator+=(const real c)
    {
        for (ulong i = 0; i < len; ++i)
            v[i] += c;
    }

    void operator+=(const rvec<len> &f)
    {
        for (ulong i = 0; i < len; ++i)
            v[i] += f[i];
    }

    void operator+=(const std::initializer_list<real> lst)
    {
        const real *it = lst.begin();
        for (ulong i = 0; i < std::max<ulong>(len, lst.size());)
            v[i++] += *it++;
    }

    void operator-=(const real c)
    {
        for (ulong i = 0; i < len; ++i)
            v[i] -= c;
    }

    void operator-=(const rvec<len> &f)
    {
        for (ulong i = 0; i < len; ++i)
            v[i] -= f[i];
    }

    void operator-=(const std::initializer_list<real> lst)
    {
        const real *it = lst.begin();
        for (ulong i = 0; i < std::max<ulong>(len, lst.size());)
            v[i++] -= *it++;
    }

    void operator*=(const real c)
    {
        for (ulong i = 0; i < len; ++i)
            v[i] *= c;
    }

    const rvec<len> operator+(const rvec<len> &f) const
    {
        rvec<len> nvec;
        for (ulong i = 0; i < len; ++i)
            nvec[i] = v[i] + f[i];
        return nvec;
    }

    const rvec<len> operator-(const rvec<len> &f) const
    {
        rvec<len> nvec;
        for (ulong i = 0; i < len; ++i)
            nvec[i] = v[i] - f[i];
        return nvec;
    }

    const rvec<len> operator*(const real c) const
    {
        rvec<len> nvec;
        for (ulong i = 0; i < len; ++i)
            nvec[i] = c * v[i];
        return nvec;
    }

    const rvec<len> operator*(const rvec<len> &f) const
    {
        rvec<len> nvec;
        for (ulong i = 0; i < len; ++i)
            nvec[i] = v[i] * f[i];
        return nvec;
    }

    const rvec<len> operator/(const real c) const
    {
        rvec<len> nvec;
        for (ulong i = 0; i < len; ++i)
            nvec[i] = c / v[i];
        return nvec;
    }

    const rvec<len> operator/(const rvec<len> &f) const
    {
        rvec<len> nvec;
        for (ulong i = 0; i < len; ++i)
            nvec[i] = v[i] / f[i];
        return nvec;
    }

    const rvec<len> square() const
    {
        rvec<len> nvec;
        for (ulong i = 0; i < len; ++i)
            nvec[i] = v[i] * v[i];
        return nvec;
    }

    const rvec<len> sqrt() const
    {
        rvec<len> nvec;
        for (ulong i = 0; i < len; ++i)
            nvec[i] = std::sqrt(std::abs(v[i]));
        return nvec;
    }

    const real norm2() const
    {
        real value = 0.;
        for (ulong i = 0; i < len; ++i)
            value += v[i] * v[i];
        return value;
    }

    real *data()
    {
        return v.data();
    }

    friend std::ostream &operator<<(std::ostream &os, rvec<len> v)
    {
        for (ulong i = 0; i < len; ++i)
            os << ',' << v[i];
        return os;
    }

private:
    std::vector<real> v;
};

using rvec3 = rvec<3>;                         // rvec3[3]
using VecVec3 = std::vector<rvec3>;            // VecVec3[i][3]
using DVecVec3 = std::vector<VecVec3>;         // DVecVec3[i][j][3]
using MapVecVec3 = std::map<int32_t, VecVec3>; // MapVecVec3[i][j][3]
using MapVec3 = std::map<int32_t, rvec3>;      // MapVec3[i][3]
using DMapVec3 = std::map<int32_t, MapVec3>;   // DMapVec3[i][j][3]

using rvec4 = rvec<4>;                         // rvec4[4]
using VecVec4 = std::vector<rvec4>;            // VecVec4[i][4]
using DVecVec4 = std::vector<VecVec4>;         // DVecVec4[i][j][4]
using MapVec4 = std::map<int32_t, rvec4>;      // MapVec4[i][4]
using DMapVec4 = std::map<int32_t, MapVec4>;   // DMapVec3[i][j][4]

// Total count of interaction types, borrowed from GROMACS-FA ForceAnal module
static constexpr uint64_t Interact_COUNT = 7;
static constexpr uint64_t Interact_FORCEVEC_LEN = 4 * (Interact_COUNT + 1);
static constexpr uint64_t Interact_ANAL = 10;

using ifvec = rvec<Interact_FORCEVEC_LEN>;     // ifvec[32]
using MapIFVec = std::map<int32_t, ifvec>;     // MapIFVec[i][32]
using DMapIFVec = std::map<int32_t, MapIFVec>; // DMapIFVec[i][j][32]


#endif /* SRC_READFA_RVEC_H_ */
