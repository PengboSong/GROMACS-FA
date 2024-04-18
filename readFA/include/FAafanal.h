/*
    FAafanal.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2023/06/28
    Description: Atom force analysis module.
*/

#ifndef SRC_READFA_FAAFANAL_H_
#define SRC_READFA_FAAFANAL_H_

#include "rvec.h"
#include "rmat.h"
#include "FAdefines.h"
#include "FAmath.h"
#include "FAotypes.h"
#include "FAsettings.h"
#include "FAcrossanal.h"

class CrossAnal;

class AtomForceAnal
{
    friend class CrossAnal;

public:
    AtomForceAnal(const FAargs &args, bool primary);

    void update(const VecVec4* F, const VecVec4* R);

    void clear();

    void digest(VecVec4* avgF, VecVec4* avgR);

    void write_header(std::ostream& out);

    void write_data(std::ostream& out);

    void tofile(std::string fnm);

protected:
    // This identifier is used to determine how results are calculated and displayed
    OutputType otype;

    // This flag determines whether standard deviation calculation should be disabled
    bool onesample;

    // This flag determines whether averages of forces and coordinates should be calculated first
    bool avgreq;

    // This flag determines whether memory-cost cross analysis should be performed
    bool crossanal;

    // Container mapping atom id to compact vector index
    IdxMap atom2vec;

    // Reversed map converting vector index to atom id
    RevIdxMap vec2atom;

    // Particles number count for group
    uint32_t N;

    // Block length
    uint32_t K;
    
    // Coefficient to average forces (=1/K)
    real avgc;

    // Coefficient to calculate sample std. dev. of forces (=1/sqrt(K-1))
    real stdc;

    // Output data precision
    uint32_t ndigits;

    /* Containers for force analysis
    Fi   - 
    SFi  - 
    SFi2 - 
    */
    VecVec4 Fi;
    VecVec4 stdFi;
    VecVec4 SFi;
    VecVec4 SFi2;

    VecVec4 *avgFi;
    VecVec4 DFi;
    VecReal SDFi;
    VecMat3 SDFiDFi;

    /* Containers for conformation analysis
    DRi   - 
    DRj   - 
    SDRi  - 
    SDRi2 - 
    SDR   - 
    SDR2  - 
    */
    VecVec4 Ri;
    VecVec4 stdRi;
    VecVec4 SRi;
    VecVec4 SRi2;

    VecVec4 *avgRi;
    VecVec4 DRi;
    VecReal SDRi;
    VecMat3 SDRiDRi;

    /* Containers for force-conformation correlation analysis
    
    */
    VecMat3 SFiRi;
    VecMat3 SFiDRi;
    VecReal SPFiDRi;
    VecReal SCOSFiDRi;
    VecReal SDEGFiDRi;
    VecReal SCOSDFiDRi;
    VecReal SDEGDFiDRi;
};

#endif