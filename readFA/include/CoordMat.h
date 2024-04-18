/*
    CoordMat.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2023/03/22
    Description: Container for coordinate matrix and index map.
*/

#ifndef SRC_READFA_COORDMAT_H_
#define SRC_READFA_COORDMAT_H_

#include "FAdefines.h"
#include "rvec.h"

template <typename T>
class IndexMap
{
    using IndexMapContainer = std::map<int32_t, uint32_t>;
    using DataContainer = std::vector<T>;

public:
    IndexMap() : n(0), i(0)
    {
    }

    IndexMap(uint32_t N) : n(N), i(0)
    {
        m.resize(N, T());
    }

    bool empty() const
    {
        return n == 0;
    }

    uint32_t size() const
    {
        return n;
    }

    DataContainer& mat()
    {
        return m;
    }

    const DataContainer mat() const
    {
        return m;
    }

    T& operator[](const int32_t ri)
    {
        return m[reindex_map.at(ri)];
    }

    const T operator[](const int32_t ri) const
    {
        return m[reindex_map.at(ri)];
    }

    void copy_map(IndexMapContainer& other_map)
    {
        for (IndexMapContainer::iterator it = reindex_map.begin(); it != reindex_map.end(); ++it)
            reindex_map[it->first] = it->second;
    }

    void update_map(const int32_t ri)
    {
        reindex_map[ri] = i++;
    }

    bool check_map() const
    {
        if (i != n) return false;
        uint32_t bits = 0;
        for (uint32_t j = 0; j < n; ++j)
            bits ^= (j + 1);
        for (IndexMapContainer::iterator it = reindex_map.begin(); it != reindex_map.end(); ++it)
            bits ^= (it->second + 1);
        if (bits != 0) return false;
        return true;
    }

private:
    uint32_t n, i;

    IndexMapContainer reindex_map;

    DataContainer m;
};

using IMapVec3 = IndexMap<rvec3>;
using DIMapVec3 = IndexMap<IMapVec3>;
using IMapReal = IndexMap<real>;
using DIMapReal = IndexMap<IMapReal>;


#endif
