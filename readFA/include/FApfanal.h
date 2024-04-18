/*
    FApfanal.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2023/06/28
    Description: Pairwise force analysis module.
*/

#ifndef SRC_READFA_FAPFANAL_H_
#define SRC_READFA_FAPFANAL_H_

#include "rvec.h"
#include "rmat.h"
#include "FAdefines.h"
#include "FAmath.h"
#include "FAotypes.h"
#include "FAsettings.h"
#include "FAgraphanal.h"
#include "FAspfanal.h"

class SummedPairwiseForceAnal;

class GraphAnal;

class PairwiseForceAnal
{
    friend class SummedPairwiseForceAnal;

    friend class GraphAnal;

public:
    PairwiseForceAnal(const FAargs& args);

    void update(const DMapVec4* F, const VecVec4* Ri, const VecVec4* Rj);

    void clear();

    void digest(DMapVec4* avgF, DVecVec4* avgR);

    void write_header(std::ofstream& out);

    void write_data(std::ostream& out);

    void tofile(std::string fnm);

protected:
    // This identifier is used to determine how results are calculated and displayed
    OutputType otype;

    // This flag determines whether standard deviation calculation should be disabled
    bool onesample;

    // This flag determines whether averages of forces and coordinates should be calculated first
    bool avgreq;

    // Container mapping atom id to compact vector index
    IdxMap M2Vi, M2Vj;

    // Reversed map converting vector index to atom id
    RevIdxMap V2Mi, V2Mj;

    // Particles number count for group
    uint32_t Ni, Nj;

    // Block length
    uint32_t K;
    
    // Coefficient to average forces (=1/K)
    real avgc;

    // Coefficient to calculate sample std. dev. of forces (=1/sqrt(K-1))
    real stdc;

    // Output data precision
    uint32_t ndigits;

    /* Containers for force data
    */
    DMapVec4 Fij;
    DMapVec4 SFij;
    DMapVec4 SFij2;

    DMapVec4 *avgFij;
    DMapVec4 DFij;
    DMapReal SDFij;
    DMapMat3 SDFijDFij;

    /* Containers for conformation analysis
    */
    DVecVec4 Rij;
    DVecVec4 SRij;
    DVecVec4 SRij2;

    DVecVec4 *avgRij;
    DVecVec4 DRij;
    DVecReal SDRij;
    DVecMat3 SDRijDRij;

    /* Containers for correlation analysis
    */
    DMapReal SFijProj;
    DMapReal SFijProj2;

    DMapMat3 SFijRij;
    DMapReal SCOSFijRij;
    DMapReal SDEGFijRij;

    DMapReal SCOSFijDRij;
    DMapReal SDEGFijDRij;
    DMapReal SCOSDFijDRij;
    DMapReal SDEGDFijDRij;
};

#endif