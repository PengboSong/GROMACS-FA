/*
    FAspfanal.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2023/09/08
    Description: Summed pairwise force analysis module.
*/

#ifndef SRC_READFA_FASPFANAL_H_
#define SRC_READFA_FASPFANAL_H_

#include "rvec.h"
#include "rmat.h"
#include "FAdefines.h"
#include "FAmath.h"
#include "FAotypes.h"
#include "FApfanal.h"

class PairwiseForceAnal;

class SummedPairwiseForceAnal
{
public:
    SummedPairwiseForceAnal(const PairwiseForceAnal* pf);

    void write_header(std::ostream& out);

    void write_data(std::ostream& out);

    void tofile(std::string fnm);

protected:
    const PairwiseForceAnal *Fij;
};

#endif