/*
    FAspfanal.cpp
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2023/09/08
    Description: Summed pairwise force analysis module.
*/

#include "FAspfanal.h"
#include "FAutility.h"

SummedPairwiseForceAnal::SummedPairwiseForceAnal(const PairwiseForceAnal *pf)
 : Fij(pf)
{
}

void SummedPairwiseForceAnal::write_header(std::ostream &out)
{
    out << "i";
    if (Fij->otype.pf())
        out << ",sum_j(|<Fij>|),sum_j(<|Fij|^2>^1/2)";
    if (Fij->otype.fcoord())
        out << ",sum_j(<FijProj>)";
    out << std::endl;
}

void SummedPairwiseForceAnal::write_data(std::ostream &out)
{
    int32_t i, j;
    real ravgc = 1.0F / Fij->avgc;
    for (auto fi = Fij->SFij.cbegin(); fi != Fij->SFij.cend(); ++fi)
    {
        i = fi->first;
        real sumj_stress = 0.0F, sumj_fnorm = 0.0F, sumj_fproj = 0.0F;
        for (auto fij = fi->second.cbegin(); fij != fi->second.cend(); ++fij)
        {
            j = fij->first;
            if (i == j) continue;
            if (Fij->otype.pf())
            {
                sumj_stress += norm(fij->second[XX], fij->second[YY], fij->second[ZZ]);
                sumj_fnorm += std::sqrt(Fij->SFij2.at(i).at(j)[XYZ] * Fij->avgc);
            }
            if (Fij->otype.fcoord())
                sumj_fproj += Fij->SFijProj.at(i).at(j);
        }
        out << i;
        if (Fij->otype.pf())
        {
            out << ',' << sumj_stress * Fij->avgc;   // sum_j(|<Fij>|)
            out << ',' << sumj_fnorm;   // sum_j(<|Fij|^2>^1/2)
        }
        if (Fij->otype.fcoord())
            out << ',' << sumj_fproj * Fij->avgc;   // sum_j(<FijProj>)
        out << std::endl;
    }
}

void SummedPairwiseForceAnal::tofile(std::string fnm)
{
    if (!Fij->otype.pf() && !Fij->otype.fcoord()) return;

    std::ofstream out;
    outfile(out, fnm, "summed pairwise analysis results");
    out << std::setiosflags(std::ios::fixed) << std::setprecision(Fij->ndigits);
    write_header(out);
    write_data(out);
    out.close();
}
