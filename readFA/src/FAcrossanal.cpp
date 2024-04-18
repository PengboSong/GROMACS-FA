/*
    FAcrossanal.cpp
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2023/06/28
    Description: Cross force and coordinate analysis module.
*/

#include "FAcrossanal.h"
#include "FAutility.h"

CrossAnal::CrossAnal(const AtomForceAnal *af1, const AtomForceAnal *af2)
 : fi(af1), fj(af2)
{
    if (!af1->crossanal || !af2->crossanal) return;
    
    if (af1->otype.atomf() && af2->otype.atomf())
    {
        initmat(SFiFj, af1->N, af2->N);
        initmat(SCOSFiFj, af1->N, af2->N);
        initmat(SDEGFiFj, af1->N, af2->N);
        if (af1->avgreq && af2->avgreq)
        {
            initmat(SCOSDFiDFj, af1->N, af2->N);
            initmat(SDEGDFiDFj, af1->N, af2->N);
        }
    }
    if (af1->otype.coord() && af2->otype.coord())
    {
        initmat(SDij, af1->N, af2->N);
        initmat(SRiRj, af1->N, af2->N);
        if (af1->avgreq && af2->avgreq)
        {
            initmat(SCOSDRiDRj, af1->N, af2->N);
            initmat(SDEGDRiDRj, af1->N, af2->N);
        }
    }
    if (af1->otype.fcoord() && af2->otype.fcoord())
    {
        if (af1->otype.coord())
        {
            initmat(SRiFj, af1->N, af2->N);
            if (af1->avgreq)
            {
                initmat(SCOSDRiFj, af1->N, af2->N);
                initmat(SDEGDRiFj, af1->N, af2->N);
            }
            if (af1->avgreq && af2->avgreq)
            {
                initmat(SCOSDRiDFj, af1->N, af2->N);
                initmat(SDEGDRiDFj, af1->N, af2->N);
            }
        }
        if (af2->otype.coord())
        {
            initmat(SFiRj, af1->N, af2->N);
            if (af2->avgreq)
            {
                initmat(SCOSFiDRj, af1->N, af2->N);
                initmat(SDEGFiDRj, af1->N, af2->N);
            }
            if (af1->avgreq && af2->avgreq)
            {
                initmat(SCOSDFiDRj, af1->N, af2->N);
                initmat(SDEGDFiDRj, af1->N, af2->N);
            }
        }
    }
}

void CrossAnal::update()
{
    if (!fi->crossanal || !fj->crossanal) return;

    for (uint32_t i = 0; i < fi->N; ++i)
        for (uint32_t j = 0; j < fj->N; ++j)
        {
            // Force-Force correlation
            if (fi->otype.atomf() && fj->otype.atomf())
            {
                SFiFj[i][j] += matmul(fi->Fi[i], fj->Fi[j]);
                real COSFiFj = fitcos(vecdot(fi->Fi[i], fj->Fi[j]) / fi->Fi[i][XYZ] / fj->Fi[j][XYZ]);
                SCOSFiFj[i][j] += COSFiFj;
                SDEGFiFj[i][j] += rad2deg(std::acos(COSFiFj));
                if (fi->avgreq && fj->avgreq)
                {
                    real COSDFiDFj = fitcos(vecdot(fi->DFi[i], fj->DFi[j]) / fi->DFi[i][XYZ] / fj->DFi[j][XYZ]);
                    SCOSDFiDFj[i][j] += COSDFiDFj;
                    SDEGDFiDFj[i][j] += rad2deg(std::acos(COSDFiDFj));
                }
            }
            // Coord-Coord correlation
            if (fi->otype.coord() && fj->otype.coord())
            {
                rvec4 Rij = fj->Ri[j] - fi->Ri[i];
                SDij[i][j] += norm(Rij[XX], Rij[YY], Rij[ZZ]);
                SRiRj[i][j] += matmul(fi->Ri[i], fj->Ri[j]);
                if (fi->avgreq && fj->avgreq)
                {
                    real COSDRiDRj = fitcos(vecdot(fi->DRi[i], fj->DRi[j]) / fi->DRi[i][XYZ] / fj->DRi[j][XYZ]);
                    SCOSDRiDRj[i][j] += COSDRiDRj;
                    SDEGDRiDRj[i][j] += rad2deg(std::acos(COSDRiDRj));
                }
            }
            // Force-Coord correlation
            if (fi->otype.fcoord() && fj->otype.fcoord())
            {                
                if (fi->otype.coord())
                {
                    SRiFj[i][j] += matmul(fi->Ri[i], fj->Fi[j]);
                    if (fi->avgreq)
                    {
                        real COSDRiFj = fitcos(vecdot(fi->DRi[i], fj->Fi[j]) / fi->DRi[i][XYZ] / fj->Fi[j][XYZ]);
                        SCOSDRiFj[i][j] += COSDRiFj;
                        SDEGDRiFj[i][j] += rad2deg(std::acos(COSDRiFj));
                    }
                    if (fi->avgreq && fj->avgreq)
                    {
                        real COSDRiDFj = fitcos(vecdot(fi->DRi[i], fj->DFi[j]) / fi->DRi[i][XYZ] / fj->DFi[j][XYZ]);
                        SCOSDRiDFj[i][j] += COSDRiDFj;
                        SDEGDRiDFj[i][j] += rad2deg(std::acos(COSDRiDFj));
                    }
                }
                if (fj->otype.coord())
                {
                    SFiRj[i][j] += matmul(fi->Fi[i], fj->Ri[j]);
                    if (fj->avgreq)
                    {
                        real COSFiDRj = fitcos(vecdot(fi->Fi[i], fj->DRi[j]) / fi->Fi[i][XYZ] / fj->DRi[j][XYZ]);
                        SCOSFiDRj[i][j] += COSFiDRj;
                        SDEGFiDRj[i][j] += rad2deg(std::acos(COSFiDRj));
                    }
                    if (fi->avgreq && fj->avgreq)
                    {
                        real COSDFiDRj = fitcos(vecdot(fi->DFi[i], fj->DRi[j]) / fi->DFi[i][XYZ] / fj->DRi[j][XYZ]);
                        SCOSDFiDRj[i][j] += COSDFiDRj;                    
                        SDEGDFiDRj[i][j] += rad2deg(std::acos(COSDFiDRj));
                    }
                }
            }
        }
}

void CrossAnal::write_header(std::ostream &out)
{
    out << "i,j";
    // Force-Force correlation
    if (fi->otype.atomf() && fj->otype.atomf())
    {
        out << mtxheader("<FiFj>") << ",trace(<FiFj>)" << pcaheader("<FiFj>");
        out << mtxheader("<DFiDFj>") << ",trace(<DFiDFj>)" << pcaheader("<DFiDFj>");
        if (!fi->onesample && !fj->onesample) out << mtxheader("cov(Fi;Fj)") << mtxheader("p(Fi;Fj)");
        out << ",cos(Fi;Fj),deg(Fi;Fj)";
        if (fi->avgreq && fj->avgreq) out << ",cos(DFi;DFj),deg(DFi;DFj)";
    }
    // Coord-Coord correlation
    if (fi->otype.coord() && fj->otype.coord())
    {
        out << ",<Dij>" << mtxheader("<DRiDRj>") << ",trace(<DRiDRj>)" << pcaheader("<DRiDRj>");
        if (!fi->onesample && !fj->onesample) out << mtxheader("cov(Ri;Rj)") << mtxheader("p(Ri;Rj)");
        if (fi->avgreq && fj->avgreq) out << ",cos(DRi;DRj),deg(DRi;DRj)";
    }
    // Force-Coord correlation
    if (fi->otype.fcoord())
    {
        if (fi->otype.coord())
        {
            out << mtxheader("<DRiFj>") << ",trace(<DRiFj>)" << pcaheader("<DRiFj>");
            if (!fi->onesample && !fj->onesample) out << mtxheader("cov(Ri;Fj)") << mtxheader("p(Ri;Fj)");
            if (fi->avgreq) out << ",cos(DRi;Fj),deg(DRi;Fj)";
            if (fi->avgreq && fj->avgreq) out << ",cos(DRi;DFj),deg(DRi;DFj)";
        }
        if (fj->otype.coord())
        {
            out << mtxheader("<FiDRj>") << ",trace(<FiDRj>)" << pcaheader("<FiDRj>");
            if (!fi->onesample && !fj->onesample) out << mtxheader("cov(Fi;Rj)") << mtxheader("p(Fi;Rj)");
            if (fj->avgreq) out << ",cos(Fi;DRj),deg(Fi;DRj)";
            if (fi->avgreq && fj->avgreq) out << ",cos(DFi;DRj),deg(DFi;DRj)";
        }
    }
    out << std::endl;
}

void CrossAnal::write_data(std::ostream &out)
{
    rmat3 FiFj, DFiDFj, DRiFj, FiDRj, DRiDRj;
    real avgc2 = fi->avgc * fi->avgc, cstdc2 = 1.0F + fi->stdc * fi->stdc;   // avgc=1/K^2, cstdc2=K/K-1
    for (int i = 0; i < fi->N; ++i)
        for (int j = 0; j < fj->N; ++j)
        {
            out << fi->vec2atom[i] << ',' << fj->vec2atom[j];
            // Force-Force correlation
            if (fi->otype.atomf() && fj->otype.atomf())
            {
                FiFj = SFiFj[i][j] * fi->avgc;
                DFiDFj = FiFj - matmul(fi->SFi[i], fj->SFi[j]) * avgc2;
                out << FiFj;   // <FiFj>
                out << ',' << trace(FiFj);   // trace(<FiFj>
                if (outmatpca(out, FiFj) != 1)   // pca(<FiFj>)
                    printf("Warning: Failed to solve eigen problems on FiFj at i=%d, j=%d.", fi->vec2atom[i], fj->vec2atom[j]);
                out << DFiDFj;   // <DFiDFj>
                out << ',' << trace(DFiDFj);   // trace(<DFiDFj>)
                if (outmatpca(out, DFiDFj) != 1)   // pca(<DFiDFj>)
                    printf("Warning: Failed to solve eigen problems on DFiDFj at i=%d, j=%d.", fi->vec2atom[i], fj->vec2atom[j]);
                if (!fi->onesample && !fj->onesample)
                {
                    out << DFiDFj * cstdc2;   // cov(Fi;Fj)
                    out << DFiDFj / matmul(fi->stdFi[i], fj->stdFi[j]) * cstdc2;   // p(Fi;Fj)
                }
                out << ',' << SCOSFiFj[i][j] * fi->avgc << ',' << SDEGFiFj[i][j] * fi->avgc;   // cos(Fi;Fj),deg(Fi;Fj)
                if (fi->avgreq && fj->avgreq)
                    out << ',' << SCOSDFiDFj[i][j] * fi->avgc << ',' << SDEGDFiDFj[i][j] * fi->avgc;   // cos(DFi;DFj),deg(DFi;DFj)
            }
            // Coord-Coord correlation
            if (fi->otype.coord() && fj->otype.coord())
            {
                DRiDRj = SRiRj[i][j] * fi->avgc - matmul(fi->SRi[i], fj->SRi[j]) * avgc2;
                out << ',' << SDij[i][j] * fi->avgc;   // <Dij>
                out << DRiDRj;   // <DRiDRj>
                out << ',' << trace(DRiDRj);   // trace(<DRiDRj>)
                if (outmatpca(out, DRiDRj) != 1)   // pca(<DRiDRj>)
                    printf("Warning: Failed to solve eigen problems on DRiDRj at i=%d, j=%d.", fi->vec2atom[i], fj->vec2atom[j]);
                if (!fi->onesample && !fj->onesample)
                {
                    out << DRiDRj * cstdc2;   // cov(Ri;Rj)
                    out << DRiDRj / matmul(fi->stdRi[i], fj->stdRi[j]) * cstdc2;   // p(Ri;Rj)
                }
                if (fi->avgreq && fj->avgreq)
                    out << ',' << SCOSDRiDRj[i][j] * fi->avgc << ',' << SDEGDRiDRj[i][j] * fi->avgc;   // cos(DRi;DRj),deg(DRi;DRj)
            }            
            // Force-Coord correlation
            if (fi->otype.fcoord())
            {
                if (fi->otype.coord())
                {
                    DRiFj = SRiFj[i][j] * fi->avgc - matmul(fi->SRi[i], fj->SFi[j]) * avgc2;
                    out << DRiFj;   // <DRiFj>
                    out << ',' << trace(DRiFj);   // trace(<DRiFj>)
                    if (outmatpca(out, DRiFj) != 1)   // pca(<DRiFj>)
                        printf("Warning: Failed to solve eigen problems on DRiFj at i=%d, j=%d.", fi->vec2atom[i], fj->vec2atom[j]);
                    if (!fi->onesample && !fj->onesample)
                    {
                        out << DRiFj * cstdc2;   // cov(Ri;Fj)
                        out << DRiFj / matmul(fi->stdRi[i], fj->stdFi[j]) * cstdc2;   // p(Ri;Fj)
                    }
                    if (fi->avgreq)
                        out << ',' << SCOSDRiFj[i][j] * fi->avgc << ',' << SDEGDRiFj[i][j] * fi->avgc;   // cos(DRi;Fj),deg(DRi;Fj)
                    if (fi->avgreq && fj->avgreq)
                        out << ',' << SCOSDRiDFj[i][j] * fi->avgc << ',' << SDEGDRiDFj[i][j] * fi->avgc;   // cos(DRi;DFj),deg(DRi;DFj)
                }
                if (fj->otype.coord())
                {
                    FiDRj = SFiRj[i][j] * fi->avgc - matmul(fi->SFi[i], fj->SRi[j]) * avgc2;
                    out << FiDRj;   // <FiDRj>
                    out << ',' << trace(FiDRj);   // trace(<FiDRj>)
                    if (outmatpca(out, FiDRj) != 1)   // pca(<FiDRj>)
                        printf("Warning: Failed to solve eigen problems on FiDRj at i=%d, j=%d.", fi->vec2atom[i], fj->vec2atom[j]);
                    if (!fi->onesample && !fj->onesample)
                    {
                        out << FiDRj * cstdc2;   // cov(Fi;Rj)
                        out << FiDRj / matmul(fi->stdFi[i], fj->stdRi[j]) * cstdc2;   // p(Fi;Rj)
                    }
                    if (fj->avgreq)
                        out << ',' << SCOSFiDRj[i][j] * fi->avgc << ',' << SDEGFiDRj[i][j] * fi->avgc;   // cos(Fi;DRj),deg(Fi;DRj)
                    if (fi->avgreq && fj->avgreq)
                        out << ',' << SCOSDFiDRj[i][j] * fi->avgc << ',' << SDEGDFiDRj[i][j] * fi->avgc;   // cos(DFi;DRj),deg(DFi;DRj)
                }
            }
            out << std::endl;
        }
}

void CrossAnal::tofile(std::string fnm)
{
    if (!fi->crossanal || !fj->crossanal) return;

    std::ofstream out;
    outfile(out, fnm, "correlation analysis results");
    out << std::setiosflags(std::ios::fixed) << std::setprecision(fi->ndigits);
    write_header(out);
    write_data(out);
    out.close();
}
