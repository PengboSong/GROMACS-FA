/*
    FAcrossanal.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2023/06/28
    Description: Cross force and coordinate analysis module.
*/

#ifndef SRC_READFA_FACROSSANAL_H_
#define SRC_READFA_FACROSSANAL_H_

#include "rvec.h"
#include "rmat.h"
#include "rmsd.h"
#include "FAdefines.h"
#include "FAmath.h"
#include "FAotypes.h"
#include "FAafanal.h"

class AtomForceAnal;

class CrossAnal
{
public:
    CrossAnal(const AtomForceAnal *af1, const AtomForceAnal *af2);
    
    void update();

    void write_header(std::ostream& out);

    void write_data(std::ostream& out);

    void tofile(std::string fnm);

protected:
    const AtomForceAnal *fi, *fj;

    /* Containers for correlation analysis
    */
    DVecMat3 SFiFj;
    DVecMat3 SFiRj;
    DVecMat3 SRiFj;
    DVecMat3 SRiRj;
    DVecReal SDij;
    
    /* Containers for angle analysis
    */
    DVecReal SCOSFiFj;
    DVecReal SCOSFiDRj;
    DVecReal SCOSDRiFj;
    DVecReal SCOSDFiDFj;
    DVecReal SCOSDFiDRj;
    DVecReal SCOSDRiDFj;
    DVecReal SCOSDRiDRj;

    DVecReal SDEGFiFj;
    DVecReal SDEGFiDRj;
    DVecReal SDEGDRiFj;
    DVecReal SDEGDFiDFj;
    DVecReal SDEGDFiDRj;
    DVecReal SDEGDRiDFj;
    DVecReal SDEGDRiDRj;
};

#endif