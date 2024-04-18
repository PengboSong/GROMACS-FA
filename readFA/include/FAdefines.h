/*
    FAdefines.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2021/07/21
    Description: Frequently used aliases and defines.
*/

#ifndef SRC_READFA_FADEFINES_H_
#define SRC_READFA_FADEFINES_H_

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <map>
#include <tuple>
#include <utility>

#include "real.h"

static const int XX = 0;
static const int YY = 1;
static const int ZZ = 2;
static const int XYZ = 3;
static const int DIM = 3;

static const real pi = 3.14159265f;
static const real pN = 1.66053878f;

using InteractionType = uint8_t;

using BlockLoc = std::tuple<int, uint32_t, uint64_t>;
using ForceDetail = std::tuple<uint8_t, real, real, real, real>;
using AtomResTable = std::map<int32_t, int32_t>;

static const char* WORDXYZ = "xyz";

#endif /* SRC_READFA_FADEFINES_H_ */