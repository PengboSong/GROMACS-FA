/*
    AtomForce.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2021/07/16
    Description: Container for atom force.
*/

#ifndef SRC_READFA_ATOMFORCE_H_
#define SRC_READFA_ATOMFORCE_H_

#include <cmath>
#include <vector>
#include <map>

#include "FAdefines.h"

class AtomForce
{
public:
    AtomForce() : atomn(0)
    {
    }

    ~AtomForce()
    {
        std::vector<InteractionType>().swap(itypes);
        std::vector<real>().swap(forces);
    }

    const uint32_t get_atomn() { return atomn; }

    void clear()
    {
        atomn = 0;
        idmap.clear();
        reverse_idmap.clear();
        itypes.clear();
        forces.clear();
    }

    int32_t find_ai(int32_t ai)
    {
        if (idmap.find(ai) == idmap.end())
        {
            idmap.insert(std::pair<int32_t, int32_t>(ai, atomn));
            reverse_idmap.insert(std::pair<int32_t, int32_t>(atomn, ai));
            ++atomn;
            itypes.resize(atomn, 0);
            forces.resize(4 * atomn, 0.);
        }
        return idmap.at(ai);
    }

    void add(int32_t ai, InteractionType itype, real fx, real fy, real fz)
    {
        int32_t vecid = find_ai(ai);
        itypes[vecid] |= itype;

        int32_t idx = 4 * vecid;
        forces[idx + XX] += fx;
        forces[idx + YY] += fy;
        forces[idx + ZZ] += fz;
    }

    void get(int32_t ai, InteractionType &itype, real &fx, real &fy, real &fz, bool *ok)
    {
        if (idmap.find(ai) != idmap.end())
        {
            int32_t vecid = idmap.at(ai);
            itype = itypes[vecid];

            int32_t idx = 4 * vecid;
            fx = forces[idx + XX];
            fy = forces[idx + YY];
            fz = forces[idx + ZZ];
            *ok = true;
        }
        else
        {
            itype = 0;
            fx = 0.;
            fy = 0.;
            fz = 0.;
            *ok = false;
        }
    }

    void norm()
    {
        real fx, fy, fz;
        for (uint32_t i = 0; i < 4 * atomn;)
        {
            fx = forces[i++];
            fy = forces[i++];
            fz = forces[i++];
            forces[i++] = std::sqrt(fx * fx + fy * fy + fz * fz);
        }
    }

    void sum(real &fx, real &fy, real &fz, real &f)
    {
        fx = 0.;
        fy = 0.;
        fz = 0.;
        for (uint32_t i = 0; i < 4 * atomn; ++i)
        {
            fx += forces[i++];
            fy += forces[i++];
            fz += forces[i++];
        }
        f = std::sqrt(fx * fx + fy * fy + fz * fz);
    }

    void operator*=(real factor)
    {
        for (uint32_t i = 0; i < 4 * atomn; ++i)
            forces[i] *= factor;
    }

protected:
    uint32_t atomn;

    std::map<int32_t, int32_t> idmap; // Mapping from real atom/residue index to vector index

    std::map<int32_t, int32_t> reverse_idmap; // Mapping from vector index to real atom/residue index

    std::vector<InteractionType> itypes;

    std::vector<real> forces;
};

using AtomForceMap = std::map<int32_t, AtomForce>;

#endif /* SRC_READFA_ATOMFORCE_H_ */