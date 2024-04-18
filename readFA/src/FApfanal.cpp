/*
    FApfanal.cpp
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2023/06/28
    Description: Pairwise force analysis module.
*/

#include <iomanip>

#include "regression.h"
#include "rmsd.h"
#include "FApfanal.h"
#include "FAutility.h"

PairwiseForceAnal::PairwiseForceAnal(const FAargs& args)
 : otype(args.otype),
   onesample(args.K <= 1),
   avgreq(args.avgreq),
   M2Vi(args.M2V.first), M2Vj(args.M2V.second),
   V2Mi(args.V2M.first), V2Mj(args.V2M.second),
   Ni(args.N.first), Nj(args.N.second),
   K(args.K), avgc(args.avgc), stdc(args.stdc),
   ndigits(args.ndigits),
   avgFij(nullptr), avgRij(nullptr)
{
    if (otype.raw()) return;

    if (otype.coord())
    {
        initmat(Rij, Ni, Nj);
        initmat(SRij, Ni, Nj);
        initmat(SRij2, Ni, Nj);
    }
}

void PairwiseForceAnal::update(const DMapVec4 *F, const VecVec4 *Ri, const VecVec4 *Rj)
{
    if (otype.raw()) return;

    if (((otype.pf() || otype.fcoord()) && F == nullptr) || ((otype.coord() || otype.fcoord()) && (Ri == nullptr || Rj == nullptr)))
        throw std::runtime_error("Pairwise force analysis module failed for empty data.");
    if ((otype.coord() || otype.fcoord()) && (Ri->size() != Ni || Rj->size() != Nj))
        throw std::runtime_error("Pairwise force analysis module failed for mismatched data length.");

    int32_t mi, mj;
    uint32_t i, j;
    if (otype.coord())
    {
        for (uint32_t i = 0; i < Ni; ++i)
            for (uint32_t j = 0; j < Nj; ++j)
            {
                Rij[i][j] = (*Rj)[j] - (*Ri)[i];
                vecnorm(Rij[i][j]);
                SRij[i][j] += Rij[i][j];
                SRij2[i][j] += Rij[i][j].square();
                if (avgreq)
                {
                    DRij[i][j] = Rij[i][j] - (*avgRij)[i][j];
                    SDRij[i][j] += vecnorm(DRij[i][j]);
                    SDRijDRij[i][j] += matmul(DRij[i][j], DRij[i][j]);
                }
            }
    }

    for (auto fi = F->cbegin(); fi != F->cend(); ++fi)
    {
        mi = fi->first;
        i = M2Vi[mi];
        for (auto fij = fi->second.cbegin(); fij != fi->second.cend(); ++fij)
        {
            mj = fij->first;
            j = M2Vj[mj];
            if (otype.pf())
            {
                Fij[mi][mj] = fij->second;
                SFij[mi][mj] += Fij[mi][mj];
                SFij2[mi][mj] += Fij[mi][mj].square();
                if (avgreq)
                {
                    DFij[mi][mj] = fij->second - (*avgFij)[mi][mj];
                    SDFij[mi][mj] += vecnorm(DFij[mi][mj]);
                    SDFijDFij[mi][mj] += matmul(DFij[mi][mj], DFij[mi][mj]);
                }
            }
            if (otype.fcoord())
            {
                rvec4 r = otype.coord() ? Rij[i][j] : (*Rj)[j] - (*Ri)[i];
                real rnorm = otype.coord() ? Rij[i][j][XYZ] : vecnorm(r);
                real fijproj = vecdot(fij->second, r) / rnorm;
                SFijProj[mi][mj] += fijproj;
                SFijProj2[mi][mj] += fijproj * fijproj;
                SFijRij[mi][mj] += matmul(fij->second, r);
                real COSFijRij = fitcos(fijproj / fij->second[XYZ]);
                SCOSFijRij[mi][mj] += COSFijRij;
                SDEGFijRij[mi][mj] += rad2deg(std::acos(COSFijRij));
                if (avgreq)
                {
                    rvec4 df = otype.pf() ? DFij[mi][mj] : fij->second - (*avgFij)[mi][mj];
                    rvec4 dr = otype.coord() ? DRij[i][j] : r - (*avgRij)[i][j];
                    real dfnorm = otype.pf() ? DFij[mi][mj][XYZ] : vecnorm(df);
                    real drnorm = otype.coord() ? DRij[i][j][XYZ] : vecnorm(dr);
                    real COSFijDRij = fitcos(vecdot(fij->second, dr) / fij->second[XYZ] / drnorm);
                    SCOSFijDRij[mi][mj] += COSFijDRij;
                    SDEGFijDRij[mi][mj] += rad2deg(std::acos(COSFijDRij));
                    real COSDFijDRij = fitcos(vecdot(df, dr) / dfnorm / drnorm);
                    SCOSDFijDRij[mi][mj] += COSDFijDRij;
                    SDEGDFijDRij[mi][mj] += rad2deg(std::acos(COSDFijDRij));
                }
            }
        }
    }
}

void PairwiseForceAnal::clear()
{
    if (otype.raw()) return;
    
    if (otype.pf())
    {
        DMapVec4 emptyFij;
        Fij.swap(emptyFij);
        if (avgreq)
        {
            DMapVec4 emptyDFij;
            DFij.swap(emptyDFij);
        }
    }
}

void PairwiseForceAnal::digest(DMapVec4 *avgF, DVecVec4 *avgR)
{
    if (otype.raw()) return;
    
    if (((otype.pf() || otype.fcoord()) && avgF == nullptr) || ((otype.coord() || otype.fcoord()) && avgR == nullptr)) avgreq = false;
    if ((otype.coord() || otype.fcoord()) && (avgR->size() != Ni || (*avgR)[0].size() != Nj)) avgreq = false;
    if (!avgreq) return;

    if (otype.coord())
    {
        initmat(DRij, Ni, Nj);
        initmat(SDRij, Ni, Nj);
        initmat(SDRijDRij, Ni, Nj);
    }
    if (otype.pf() || otype.fcoord()) avgFij = avgF;
    if (otype.coord() || otype.fcoord()) avgRij = avgR;
}

void PairwiseForceAnal::write_header(std::ofstream &out)
{
    out << "i,j";
    if (otype.pf())
    {
        out << ",<Fijx>,<Fijy>,<Fijz>,<|Fij|>,<|Fij|^2>^1/2,|<Fij>|";
        if (!onesample)
            out << ",var(Fijx),var(Fijy),var(Fijz),var(|Fij|),var(Fij),std(Fijx),std(Fijy),std(Fijz),std(|Fij|),|std(Fij)|";
        if (avgreq)
        {
            out << ",<DFij>";
            if (!onesample) out << ",var(DFij),std(DFij)";
            out << mtxheader("<DFijDFij>") << ",trace(<DFijDFij>)" << pcaheader("<DFijDFij>");
        }
    }
    if (otype.coord())
    {
        out << ",<Rijx>,<Rijy>,<Rijz>,<|Rij|>,<|Rij|^2>^1/2,|<Rij>|";
        if (!onesample)
            out << ",var(Rijx),var(Rijy),var(Rijz),var(|Rij|),var(Rij),std(Rijx),std(Rijy),std(Rijz),std(|Rij|),|std(Rij)|";
        if (avgreq)
        {
            out << ",<DRij>";
            if (!onesample) out << ",var(DRij),std(DRij)";
            out << mtxheader("<DRijDRij>") << ",trace(<DRijDRij>)" << pcaheader("<DRijDRij>");
        }
    }
    if (otype.fcoord())
    {
        out << ",<FijProj>";
        if (!onesample)
        {
            out << ",var(FijProj),std(FijProj)";
            out << ",var(FijProj)/var(|Rij|),|Rij|-FijProj/b,|Rij|-FijProj/k,|Rij|-FijProj/var(E(|Rij|)),|Rij|-FijProj/var(E(FijProj)),|Rij|-FijProj/var(Eerr(|Rij|)),|Rij|-FijProj/var(Eerr(FijProj)),|Rij|-FijProj/lb(CI95(b)),|Rij|-FijProj/ub(CI95(b)),|Rij|-FijProj/lb(CI95(k)),|Rij|-FijProj/ub(CI95(k))";
        }
        out << mtxheader("<FijRij>") << ",trace(<FijRij>)" << pcaheader("<FijRij>");
        out << mtxheader("<DFijDRij>") << ",trace(<DFijDRij>)" << pcaheader("<DFijDRij>");
        if (!onesample) out << mtxheader("cov(Fij;Rij)") << mtxheader("p(Fij;Rij)");
        out << ",cos(Fij;Rij),deg(Fij;Rij)";
        if (avgreq) out << ",cos(Fij;DRij),deg(Fij;DRij),cos(DFij;DRij),deg(DFij;DRij)";
    }
    out << std::endl;
}

void PairwiseForceAnal::write_data(std::ostream &out)
{
    int32_t mi, mj;
    uint32_t i, j;
    rvec4 avgF, varF, stdF, avgR, varR, stdR;
    rmat3 DFijDFij, DRijDRij, FijRij, DFijDRij;
    real varc = stdc * stdc;   // =1/K-1
    for (auto fi = SFij.cbegin(); fi != SFij.cend(); ++fi)
    {
        mi = fi->first;
        i = M2Vi[mi];
        for (auto fij = fi->second.cbegin(); fij != fi->second.cend(); ++fij)
        {
            mj = fij->first;
            j = M2Vj[mj];
            out << mi << ',' << mj;
            if (otype.pf())
            {
                avgF = fij->second * avgc;
                out << avgF << ',' << std::sqrt(SFij2[mi][mj][XYZ] * avgc) << ',' << norm(avgF[XX], avgF[YY], avgF[ZZ]);   // <Fijx>,<Fijy>,<Fijz>,<|Fij|>,<|Fij|^2>^1/2,|<Fij>|
                if (!onesample)
                {
                    varF = (SFij2[mi][mj] - avgF.square() * K) * varc;
                    out << varF << ',' << varF[XX] + varF[YY] + varF[ZZ];   // var(Fijx),var(Fijy),var(Fijz),var(|Fij|),var(Fij)
                    stdF = varF.sqrt();
                    out << stdF << ',' << norm(stdF[XX], stdF[YY], stdF[ZZ]);   // std(Fijx),std(Fijy),std(Fijz),std(|Fij|),|std(Fij)|
                }
                if (avgreq)
                {
                    out << ',' << SDFij[mi][mj] * avgc;   // <DFij>
                    if (!onesample)
                    {
                        real varDFij = (trace(SDFijDFij[mi][mj]) - SDFij[mi][mj] * SDFij[mi][mj] * avgc) * varc; 
                        out << ',' << varDFij << ',' << std::sqrt(varDFij);   // var(DFij),std(DFij)
                    }
                    DFijDFij = SDFijDFij[mi][mj] * avgc;
                    out << DFijDFij;   // <DFijDFij>
                    out << ',' << trace(DFijDFij);   // trace(<DFijDFij>)
                    if (outmatpca(out, DFijDFij) != 1)   // pca(<DFijDFij>)
                        printf("Warning: Failed to solve eigen problems on DFijDFij at i=%d, j=%d.", mi, mj);
                }
            }
            if (otype.coord())
            {
                avgR = SRij[i][j] * avgc;
                out << avgR << ',' << std::sqrt(SRij2[i][j][XYZ] * avgc) << ',' << norm(avgR[XX], avgR[YY], avgR[ZZ]);   // <Rijx>,<Rijy>,<Rijz>,<|Rij|>,<|Rij|^2>^1/2,|<Rij>|
                if (!onesample)
                {
                    varR = (SRij2[i][j] - avgR.square() * K) * varc;
                    out << varR << ',' << varR[XX] + varR[YY] + varR[ZZ];   // var(Rijx),var(Rijy),var(Rijz),var(|Rij|),var(Rij)
                    stdR = varR.sqrt();
                    out << stdR << ',' << norm(stdR[XX], stdR[YY], stdR[ZZ]);   // std(Rijx),std(Rijy),std(Rijz),std(|Rij|),|std(Rij)|
                }
                if (avgreq)
                {
                    out << ',' << SDRij[i][j] * avgc;   // <DRij>
                    if (!onesample)
                    {
                        real varDRij = (trace(SDRijDRij[i][j]) - SDRij[i][j] * SDRij[i][j] * avgc) * varc;
                        out << ',' << varDRij << ',' << std::sqrt(varDRij);   // var(DRij),std(DRij)
                    }
                    DRijDRij = SDRijDRij[i][j] * avgc;
                    out << DRijDRij;   // <DRijDRij>
                    out << ',' << trace(DRijDRij);   // trace(<DRijDRij>)
                    if (outmatpca(out, DRijDRij) != 1)   // pca(<DRijDRij>)
                        printf("Warning: Failed to solve eigen problems on DRijDRij at i=%d, j=%d.", mi, mj);
                }
            }
            if (otype.fcoord())
            {
                real avgFijProj = SFijProj[mi][mj] * avgc;
                out << ',' << avgFijProj;   // <FijProj>
                if (!onesample)
                {
                    real varFijProj = (SFijProj2[mi][mj] - avgFijProj * avgFijProj * K) * varc;
                    out << ',' << varFijProj << ',' << std::sqrt(varFijProj);   // var(FijProj), std(FijProj)
                    if (otype.coord())
                    {
                        FitResult fit_Rij_FijProj;   // Deming fitting result for <|Rij|> & FijProj
                        DemingRegression::fit(K, avgR[XYZ], avgFijProj, SRij2[i][j][XYZ] * avgc, SFijProj2[mi][mj] * avgc, trace(SFijRij[mi][mj]), varR[XYZ], varFijProj, fit_Rij_FijProj);
                        out << fit_Rij_FijProj;   // var(FijProj)/var(|Rij|),|Rij|-FijProj/b,|Rij|-FijProj/k,|Rij|-FijProj/var(E(|Rij|)),|Rij|-FijProj/var(E(FijProj)),|Rij|-FijProj/var(Eerr(|Rij|)),|Rij|-FijProj/var(Eerr(FijProj)),|Rij|-FijProj/lb(CI95(b)),|Rij|-FijProj/ub(CI95(b)),|Rij|-FijProj/lb(CI95(k)),|Rij|-FijProj/ub(CI95(k))
                    }
                }
                FijRij = SFijRij[mi][mj] * avgc;
                out << FijRij;   // <FijRij>
                out << ',' << trace(FijRij);   // trace(<FijRij>)
                if (outmatpca(out, FijRij) != 1)   // pca(<FijRij>)
                    printf("Warning: Failed to solve eigen problems on FijRij at i=%d, j=%d.", mi, mj);
                DFijDRij = FijRij - matmul(avgF, avgR);
                out << DFijDRij;   // <DFijDRij>
                out << ',' << trace(DFijDRij);   // trace(<DFijDRij>)
                if (outmatpca(out, DFijDRij) != 1)   // pca(<DFijDRij>)
                    printf("Warning: Failed to solve eigen problems on DFijDRij at i=%d, j=%d.", mi, mj);
                if (!onesample)
                {
                    out << DFijDRij * (1.0 + varc);   // cov(Fij;Rij)
                    out << DFijDRij / matmul(stdF, stdR) * (1.0 + varc);   // p(Fij;Rij)
                }
                out << ',' << SCOSFijRij[mi][mj] * avgc << ',' << SDEGFijRij[mi][mj] * avgc;   // cos(Fij;Rij),deg(Fij;Rij)
                if (avgreq)
                    out << ',' << SCOSFijDRij[mi][mj] * avgc << ',' << SDEGFijDRij[mi][mj] * avgc << ',' << SCOSDFijDRij[mi][mj] * avgc << ',' << SDEGDFijDRij[mi][mj] * avgc;   // cos(Fij;DRij),deg(Fij;DRij),cos(DFij;DRij),deg(DFij;DRij)
            }
            out << std::endl;
        }
    }
}

void PairwiseForceAnal::tofile(std::string fnm)
{
    if (otype.raw()) return;
    
    std::ofstream out;
    outfile(out, fnm, "pairwise analysis results");
    out << std::setiosflags(std::ios::fixed) << std::setprecision(ndigits);
    write_header(out);
    write_data(out);
    out.close();
}
