/*
    FAutility.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2023/07/12
    Description: Collections of useful functions.
*/

#ifndef SRC_READFA_FAUTILITY_H_
#define SRC_READFA_FAUTILITY_H_

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "FAdefines.h"
#include "rmsd.h"

static inline std::string mtxheader(const char* s)
{
    std::stringstream out;
    for (int i = 0; i < DIM; ++i)
        for (int j = 0; j < DIM; ++j)
            out << ',' << s << WORDXYZ[i] << WORDXYZ[j];
    return out.str();
}

static inline std::string pcaheader(const char* s)
{
    std::vector<const char*> WORDABC = {"α", "β", "γ"};
    std::stringstream out;
    for (int i = 0; i < DIM; ++i) out << ',' << s << "λ" << WORDABC[i];
    for (int i = 0; i < DIM; ++i)
        for (int j = 0; j < DIM; ++j)
            out << ',' << s << 'U' << WORDABC[i] << WORDXYZ[j];
    return out.str();
}

static inline std::string modfnm(std::string fnm, std::string prefix, std::string suffix)
{
    return prefix + (suffix.empty() ? fnm : (fnm.find_first_of('.') == std::string::npos ? fnm + suffix : fnm.substr(0, fnm.find_first_of('.')) + suffix + fnm.substr(fnm.find_first_of('.'))));
}

static inline void outfile(std::ofstream& out, std::string fnm, std::string hint)
{
    if (hint.empty()) hint = "force analysis data";
    out.open(fnm);
    if (!out.is_open())
        throw std::runtime_error("Can not write " + hint + " to file " + fnm);
}

static inline int outmatpca(std::ostream& out, const rmat3& mat)
{
    rvec3 lam;
    rmat3 U, M = mat;
    int v = EigenRvec3(M, lam, U);
    int s = 0, m = 1, l = 2;
    if (lam[0] > lam[1])
    {
        std::swap(lam[0], lam[1]);
        std::swap(s, m);
    }
    if (lam[0] > lam[2])
    {
        std::swap(lam[0], lam[2]);
        std::swap(s, l);
    }
    if (lam[1] > lam[2])
    {
        std::swap(lam[1], lam[2]);
        std::swap(m, l);
    }
    rvec3 eigvecA({U[XX][s], U[YY][s], U[ZZ][s]});
    rvec3 eigvecB({U[XX][m], U[YY][m], U[ZZ][m]});
    rvec3 eigvecR({U[XX][l], U[YY][l], U[ZZ][l]});
    out << lam << eigvecA << eigvecB << eigvecR;
    return v;
}

#endif /* SRC_READFA_FAUTILITY_H_ */