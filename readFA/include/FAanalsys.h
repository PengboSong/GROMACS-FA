/*
    FAanalsys.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2023/07/11
    Description: Collections of force and coordinate analysis module.
*/

#ifndef SRC_READFA_FAANALSYS_H_
#define SRC_READFA_FAANALSYS_H_

#include "rvec.h"
#include "rmat.h"
#include "FAdefines.h"
#include "FAmath.h"
#include "FAotypes.h"
#include "FAsettings.h"
#include "FAafanal.h"
#include "FAcrossanal.h"
#include "FApfanal.h"
#include "FAgraphanal.h"
#include "FAspfanal.h"

class AtomForceAnal;
class PairwiseForceAnal;
class CrossAnal;
class SummedPairwiseForceAnal;
class GraphAnal;

class AnalSystem
{
public:
    AnalSystem(const FAargs& inargs, std::string tag = "");

    inline void reset_tag(std::string tag) { prefix = tag; }

    void add_detailed_force(int32_t i, int32_t j, rvec4& force);

    void avg_forces(const rmat3& U);

    void digest(VecVec4* avgRi, VecVec4* avgRj, DVecVec4* avgRij);

    void update(const rmat3& U, const VecVec4* Ri, const VecVec4* Rj);

    void clear();

    void write_anal();

protected:
    const FAargs& args;

    std::string prefix;

    VecVec4 Fi;
    VecVec4 Fj;
    DMapVec4 Fij;

    VecVec4 avgFi;
    VecVec4 avgFj;
    DMapVec4 avgFij;

    AtomForceAnal FiAnal;
    AtomForceAnal FjAnal;
    PairwiseForceAnal FijAnal;
    CrossAnal FiFjAnal;
    SummedPairwiseForceAnal SjFijAnal;
    GraphAnal GAnal;
};

#endif