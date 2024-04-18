/*
    FAafanal.cpp
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2023/06/28
    Description: Atom force analysis module.
*/

#include "FAafanal.h"
#include "FAutility.h"

AtomForceAnal::AtomForceAnal(const FAargs &args, bool primary)
 : otype(args.otype),
   onesample(args.K <= 1),
   avgreq(args.avgreq), crossanal(args.crossanal),
   atom2vec(primary ? args.M2V.first : args.M2V.second),
   vec2atom(primary ? args.V2M.first : args.V2M.second),
   N(primary ? args.N.first : args.N.second),
   K(args.K), avgc(args.avgc), stdc(args.stdc),
   ndigits(args.ndigits),
   avgFi(nullptr), avgRi(nullptr)
{
    if (otype.raw()) return;

    if (otype.atomf())
    {
        initvec(Fi, N);
        initvec(stdFi, N);
        initvec(SFi, N);
        initvec(SFi2, N);
    }
    if (otype.coord())
    {
        initvec(Ri, N);
        initvec(stdRi, N);
        initvec(SRi, N);
        initvec(SRi2, N);
    }
    if (otype.fcoord())
    {
        initvec(SFiRi, N);
        initvec(SFiDRi, N);
        initvec(SPFiDRi, N);
    }
}

void AtomForceAnal::update(const VecVec4 *F, const VecVec4 *R)
{
    if (otype.raw()) return;

    if (((otype.atomf() || otype.fcoord()) && F == nullptr) || ((otype.coord() || otype.fcoord()) && R == nullptr))
        throw std::runtime_error("Atom force analysis module failed for empty data.");
    if (((otype.atomf() || otype.fcoord()) && F->size() != N) || ((otype.coord() || otype.fcoord()) && R->size() != N))
        throw std::runtime_error("Atom force analysis module failed for mismatched data length.");
    
    for (int i = 0; i < N; ++i)
    {
        if (otype.atomf())
        {
            Fi[i] = (*F)[i];
            SFi[i] += Fi[i];
            SFi2[i] += Fi[i].square();
            if (avgreq)
            {
                DFi[i] = Fi[i] - (*avgFi)[i];
                SDFi[i] += vecnorm(DFi[i]);
                SDFiDFi[i] += matmul(DFi[i], DFi[i]);
            }
        }
        if (otype.coord())
        {
            Ri[i] = (*R)[i];
            SRi[i] += Ri[i];
            SRi2[i] += Ri[i].square();
            if (avgreq)
            {
                DRi[i] = Ri[i] - (*avgRi)[i];
                SDRi[i] += vecnorm(DRi[i]);
                SDRiDRi[i] += matmul(DRi[i], DRi[i]);
            }
        }
        if (otype.fcoord())
        {
            SFiRi[i] += matmul((*F)[i], (*R)[i]);
            if (avgreq)
            {
                rvec4 df = otype.atomf() ? DFi[i] : ((*F)[i] - (*avgFi)[i]);
                rvec4 dr = otype.coord() ? DRi[i] : ((*R)[i] - (*avgRi)[i]);
                real dfnorm = otype.atomf() ? DFi[i][XYZ] : vecnorm(df);
                real drnorm = otype.coord() ? DRi[i][XYZ] : vecnorm(dr);
                SFiDRi[i] += matmul(((*F)[i]), dr);
                SPFiDRi[i] += vecdot((*F)[i], dr) / (*F)[i][XYZ] / drnorm;
                real COSFiDRi = fitcos(vecdot((*F)[i], dr) / (*F)[i][XYZ] / drnorm);
                SCOSFiDRi[i] += COSFiDRi;
                SDEGFiDRi[i] += rad2deg(std::acos(COSFiDRi));
                real COSDFiDRi = fitcos(vecdot(df, dr) / dfnorm / drnorm);
                SCOSDFiDRi[i] += COSDFiDRi;
                SDEGDFiDRi[i] += rad2deg(std::acos(COSDFiDRi));
            }
        }
    }
}

void AtomForceAnal::clear()
{
    // No cleaning
}

void AtomForceAnal::digest(VecVec4 *avgF, VecVec4 *avgR)
{
    if (otype.raw()) return;

    if (((otype.atomf() || otype.fcoord()) && avgF == nullptr) || ((otype.coord() || otype.fcoord()) && avgR == nullptr)) avgreq = false;
    if (((otype.atomf() || otype.fcoord()) && avgF->size() != N) || ((otype.coord() || otype.fcoord()) && avgR->size() != N)) avgreq = false;
    if (!avgreq) return;

    if (otype.atomf())
    {
        initvec(DFi, N);
        initvec(SDFi, N);
        initvec(SDFiDFi, N);
    }
    if (otype.coord())
    {
        initvec(DRi, N);
        initvec(SDRi, N);
        initvec(SDRiDRi, N);
    }
    if (otype.fcoord())
    {
        initvec(SCOSFiDRi, N);
        initvec(SDEGFiDRi, N);
        initvec(SCOSDFiDRi, N);
        initvec(SDEGDFiDRi, N);
    }
    if (otype.atomf() || otype.fcoord()) avgFi = avgF;
    if (otype.coord() || otype.fcoord()) avgRi = avgR;
}

void AtomForceAnal::write_header(std::ostream &out)
{
    out << "i";
    if (otype.atomf())
    {
        out << ",<Fix>,<Fiy>,<Fiz>,<|Fi|>,<|Fi|^2>^1/2,|<Fi>|";
        if (!onesample)
            out << ",var(Fix),var(Fiy),var(Fiz),var(|Fi|),var(Fi),std(Fix),std(Fiy),std(Fiz),std(|Fi|),|std(Fi)|";
        if (avgreq)
        {
            out << ",<DFi>";
            if (!onesample) out << ",var(DFi),std(DFi)";
            out << mtxheader("<DFiDFi>") << ",trace(<DFiDFi>)" << pcaheader("<DFiDFi>");
        }
    }
    if (otype.coord())
    {
        out << ",<Rix>,<Riy>,<Riz>,<|Ri|>,<|Ri|^2>^1/2,|<Ri>|";
        if (!onesample)
            out << ",var(Rix),var(Riy),var(Riz),var(|Ri|),var(Ri),std(Rix),std(Riy),std(Riz),std(|Ri|),|std(Ri)|";
        if (avgreq)
        {
            out << ",<DRi>";
            if (!onesample) out << ",var(DRi),std(DRi)";
            out << mtxheader("<DRiDRi>") << ",trace(<DRiDRi>)" << pcaheader("<DRiDRi>");
        }
    }
    if (otype.fcoord())
    {
        out << mtxheader("<FiDRi>") << ",trace(<FiDRi>),p(Fi;Ri)" << pcaheader("<FiDRi>");
        if (!onesample) out << mtxheader("cov(Fi;Ri)") << mtxheader("p(Fi;Ri)");
        if (avgreq) out << ",cos(Fi;DRi),deg(Fi;DRi),cos(DFi;DRi),deg(DFi;DRi)";
    }
    out << std::endl;
}

void AtomForceAnal::write_data(std::ostream &out)
{
    rvec4 avgF, varF, stdF, avgR, varR, stdR;
    rmat3 DFiDFi, DRiDRi, FiDRi;
    real varc = stdc * stdc;   // =1/K-1
    for (int i = 0; i < N; ++i)
    {
        out << vec2atom[i];
        if (otype.atomf())
        {
            avgF = SFi[i] * avgc;
            out << avgF << ',' << std::sqrt(SFi2[i][XYZ] * avgc) << ',' << norm(avgF[XX], avgF[YY], avgF[ZZ]);   // <Fix>,<Fiy>,<Fiz>,<|Fi|>,<|Fi|^2>^1/2,|<Fi>|
            if (!onesample)
            {
                varF = (SFi2[i] - avgF.square() * K) * varc;
                out << varF << ',' << varF[XX] + varF[YY] + varF[ZZ];   // var(Fix),var(Fiy),var(Fiz),var(|Fi|),var(Fi)
                stdF = varF.sqrt();
                stdFi[i] = stdF;
                out << stdF << ',' << std::sqrt(varF[XX] + varF[YY] + varF[ZZ]);   // std(Fix),std(Fiy),std(Fiz),std(|Fi|),|std(Fi)|
            }
            if (avgreq)
            {
                out << ',' << SDFi[i] * avgc;   // <DFi>
                if (!onesample)
                {
                    real varDFi = (trace(SDFiDFi[i]) - SDFi[i] * SDFi[i] * avgc) * varc;
                    out << ',' << varDFi << ',' << std::sqrt(varDFi);   // var(DFi), std(DFi)
                }
                DFiDFi = SDFiDFi[i] * avgc;
                out << DFiDFi;   // <DFiDFi>
                out << ',' << trace(DFiDFi);   // trace(<DFiDFi>)
                if (outmatpca(out, DFiDFi) != 1)   // pca(<DFiDFi>)
                    printf("Warning: Failed to solve eigen problems on DFiDFi at i=%d.", vec2atom[i]);
            }
        }
        if (otype.coord())
        {
            avgR = SRi[i] * avgc;
            out << avgR << ',' << std::sqrt(SRi2[i][XYZ] * avgc) << ',' << norm(avgR[XX], avgR[YY], avgR[ZZ]);   // <Rix>,<Riy>,<Riz>,<|Ri|>,<|Ri|^2>^1/2,|<Ri>|
            if (!onesample)
            {
                varR = (SRi2[i] - avgR.square() * K) * varc;
                out << varR << ',' << varR[XX] + varR[YY] + varR[ZZ];   // var(Rix),var(Riy),var(Riz),var(|Ri|),var(Ri)
                stdR = varR.sqrt();
                stdRi[i] = stdR;
                out << stdR << ',' << std::sqrt(varR[XX] + varR[YY] + varR[ZZ]);   // std(Rix),std(Riy),std(Riz),|std(Ri)|
            }
            if (avgreq)
            {
                out << ',' << SDRi[i] * avgc;   // <DRi>
                if (!onesample)
                {
                    real varDRi = (trace(SDRiDRi[i]) - SDRi[i] * SDRi[i] * avgc) * varc;
                    out << ',' << varDRi << ',' << std::sqrt(varDRi);   // var(DRi),std(DRi)
                }
                DRiDRi = SDRiDRi[i] * avgc;
                out << DRiDRi;   // <DRiDRi>
                out << ',' << trace(DRiDRi);   // trace(<DRiDRi>)
                if (outmatpca(out, DRiDRi) != 1)   // pca(<DRiDRi>)
                    printf("Warning: Failed to solve eigen problems on DRiDRi at i=%d.", vec2atom[i]);
            }
        }
        if (otype.fcoord())
        {
            // FiDRi = SFiRi[i] * avgc - matmul(avgF, avgR);
            FiDRi = SFiDRi[i] * avgc;
            out << FiDRi;   // <FiDRi>
            out << ',' << trace(FiDRi) << ',' << SPFiDRi[i] * avgc;   // trace(<FiDRi>)
            if (outmatpca(out, FiDRi) != 1)   // pca(<FiDRi>)
                printf("Warning: Failed to solve eigen problems on FiDRi at i=%d.", vec2atom[i]);
            if (!onesample)
            {
                out << FiDRi * (1.0 + varc);   // cov(Fi,Ri)
                out << FiDRi / matmul(stdF, stdR) * (1.0 + varc);   // p(Fi,Ri)
            }
            if (avgreq)
            {
                out << ',' << SCOSFiDRi[i] * avgc << ',' << SDEGFiDRi[i] * avgc;   // cos(Fi,DRi),deg(Fi,DRi)
                out << ',' << SCOSDFiDRi[i] * avgc << ',' << SDEGDFiDRi[i] * avgc;  // cos(DFi,DRi),deg(DFi,DRi)
            }
        }
        out << std::endl;
    }
}

void AtomForceAnal::tofile(std::string fnm)
{
    if (otype.raw()) return;

    std::ofstream out;
    outfile(out, fnm, "atom/residue-level analysis results");
    out << std::setiosflags(std::ios::fixed) << std::setprecision(ndigits);
    write_header(out);
    write_data(out);
    out.close();
}
