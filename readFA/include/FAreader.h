/*
    FAreader.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2021/08/18
    Description: Core module of readFA.
*/

#ifndef SRC_READFA_FAREADER_H_
#define SRC_READFA_FAREADER_H_

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "real.h"
#include "rvec.h"
#include "rmat.h"
#include "rmsd.h"
#include "CoordMat.h"
#include "FAdefines.h"
#include "FAmath.h"
#include "FAotypes.h"
#include "FAsettings.h"
#include "FAafanal.h"
#include "FApfanal.h"
#include "FAcrossanal.h"
#include "FAgraphanal.h"
#include "FAanalsys.h"

class AnalSystem;

class FAReader : public FAsettings
{
    friend class AnalSystem;

public:
    FAReader(FAargs& inargs);

protected:
    void read_listed_forces();

    void read_detailed_forces();

    void read_atom_forces();

    void read_frame(uint32_t k, VecVec4* Ri, VecVec4* Rj);

    void fit_frames();

    /*
    U - 
    */
    VecVec4 avgRi;
    VecVec4 avgRj;
    DVecVec4 avgRij;
    MapVecVec3 mergeRi;
    MapVecVec3 mergeRj;
};

#endif /* SRC_READFA_FAREADER_H_ */