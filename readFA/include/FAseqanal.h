/*
    FAseqanal.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2023/07/04
    Description: Time sequence force and coordinate analysis module.
*/

#ifndef SRC_READFA_FASEQANAL_H_
#define SRC_READFA_FASEQANAL_H_

#include <map>
#include <utility>
#include <vector>
#include <unordered_map>

#include "rvec.h"
#include "rmat.h"
#include "FAdefines.h"
#include "FAmath.h"
#include "FAotypes.h"
#include "FAsettings.h"

using PINT = std::pair<int32_t, int32_t>;
using TVec3 = VecVec3;
using TVec4 = VecVec4;
using MapTVec4 = std::map<int32_t, TVec4>;
using PMapTVec4 = std::map<PINT, TVec4>;

class SequenceAnal
{
public:
    SequenceAnal(const FAsettings &args);

protected:
    /* Containers for force data
    */
    MapTVec4 Fi;
    PMapTVec4 Fij;

    /* Containers for conformation analysis
    */
    MapTVec4 Ri;
    PMapTVec4 Rij;
};

#endif