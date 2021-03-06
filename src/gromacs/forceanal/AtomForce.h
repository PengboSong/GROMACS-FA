/*
 * AtomForce.h
 *
 *  Created on: June 16, 2021
 *      Author: Pengbo Song
 */

#ifndef SRC_GROMACS_FORCEANAL_ATOMFORCE_H_
#define SRC_GROMACS_FORCEANAL_ATOMFORCE_H_

#include <vector>

#include "gromacs/math/extended_vec.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"
#include "ForceAnalConst.h"
#include "InteractionForce.h"
#include "InteractionType.h"

namespace ForceAnal {
    template<typename T>
    class OffsetVector
    {
        public:

        OffsetVector() : offset(0), length(0)
        {}

        ~OffsetVector()
        {
            std::vector<T>().swap(container);
        }

        void init()
        {
            /*
            container.reserve(length);
            for (uint64_t i = 0; i < length; ++i)
                container[i] = T();
            */
            container.resize(length, T());
        }

        void reset(int64_t sloc, int64_t eloc)
        {
            offset = -sloc;
            length = eloc - sloc + 1;
            init();
        }

        T& operator[](const int64_t idx)
        {
            return container[idx + offset];
        }

        int64_t offset;

        uint64_t length;

        std::vector<T> container;
    };
    
    using DetailedAtomForce = OffsetVector<InteractionForce>;
    using DetailedForce = OffsetVector<DetailedAtomForce>;

    class AtomForce
    {
        public:

        AtomForce() : offset(0), atomn(0), forcelen(0)
        {}

        ~AtomForce()
        {
            std::vector<InteractionType>().swap(itypes);
            std::vector<real>().swap(forces);
        }

        void init()
        {
            /*
            itypes.reserve(atomn);
            for (uint64_t i = 0; i < atomn; ++i)
                itypes[i] = 0;
            forces.reserve(forcelen);
            for (uint64_t i = 0; i < forcelen; ++i)
                forces[i] = 0.;
            */
            itypes.resize(atomn, 0);
            forces.resize(forcelen, 0.);
        }

        void reset(int64_t sloc, int64_t eloc)
        {
            offset = -sloc;
            atomn = eloc - sloc + 1;
            forcelen = 4 * atomn;
            init();
        }

        void add(int64_t ai, InteractionType itype, rvec force)
        {
            int64_t idx = 4 * (ai + offset);
            forces[idx + XX] += force[XX];
            forces[idx + YY] += force[YY];
            forces[idx + ZZ] += force[ZZ];
            itypes[ai + offset] &= itype;
        }

        void norm2()
        {
            real fx, fy, fz;
            for (uint64_t i = 0; i < forcelen;)
            {
                fx = forces[i++];
                fy = forces[i++];
                fz = forces[i++];
                forces[i++] = fx * fx + fy * fy + fz * fz;
            }
        }

        void operator*=(double factor)
        {
            for (uint64_t i = 0; i < forcelen; ++i)
                forces[i] *= factor;
        }

        int64_t offset;

        uint64_t atomn;

        uint64_t forcelen;

        std::vector<InteractionType> itypes;

        std::vector<real> forces;
    };

    using SummedForce = OffsetVector<AtomForce>;
}

#endif /* SRC_GROMACS_FORCEANAL_ATOMFORCE_H_ */
