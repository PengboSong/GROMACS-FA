/*
    FAgraphanal.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2023/07/10
    Description: Force graph analysis module.
*/

#ifndef SRC_READFA_FAGRAPHANAL_H_
#define SRC_READFA_FAGRAPHANAL_H_

#include "real.h"
#include "rvec.h"
#include "rmat.h"
#include "FAdefines.h"
#include "FAmath.h"
#include "FAotypes.h"
#include "FApfanal.h"

class PairwiseForceAnal;

class GraphAnal
{
public:
    GraphAnal(const PairwiseForceAnal* pf, const real lowpf);

    void write_header(std::ostream& out);

    void write_data(std::ostream& out);

    void tofile(std::string fnm);

protected:
    const PairwiseForceAnal *Fij;

    // Lower bound of pairwise forces as edges of force graph
    real g_lowpf;
};

#endif