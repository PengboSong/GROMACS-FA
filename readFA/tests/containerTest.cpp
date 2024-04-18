/*
    containerTest.cpp
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2021/08/26
    Description: Test functionalities of containers.
*/

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <utility>

#include "real.h"
#include "rvec.h"
#include "rmat.h"
#include "FAmath.h"

void deltaforce(const rvec4 f, real deltaf, rvec4& df)
{
    for (int i = 0; i < 3; ++i)
        df[i] = f[i] + deltaf;
}

int main(int argc, char **argv)
{
    const rvec4 f1 = {9.537, 5.701, 6.564};
    const rvec4 f2 = {9.857, 3.061, 9.184};
    const rvec4 f3 = {9.462, 0.577, 5.098};

    int nframe = 10;

    MapVec4 Fi;
    MapVec4 SFi;
    MapVec4 SFi2;
    DMapVec4 Fij;
    DMapVec4 SFij;
    DMapVec4 SFij2;
    DMapMat3 SFiFj;
    
    // df1 - f01, df2 - f02, df3 - f12
    int ai, aj;
    real deltaf;
    rvec4 df1, df2, df3, ftmp;
    for (int k = 0; k < nframe; ++k)
    {
        deltaf = .15 * k;
        deltaforce(f1, deltaf, df1);
        deltaforce(f2, deltaf, df2);
        deltaforce(f3, deltaf, df3);

        Fi[0] += df1;
        Fi[1] -= df1;
        Fij[0][1] += df1;
        Fij[1][0] -= df1;

        Fi[0] += df2;
        Fi[2] -= df2;
        Fij[0][2] += df2;
        Fij[2][0] -= df2;

        Fi[1] += df3;
        Fi[2] -= df3;
        Fij[1][2] += df3;
        Fij[2][1] -= df3;

        for (MapVec4::iterator fi = Fi.begin(); fi != Fi.end(); ++fi)
        {
            ai = fi->first;
            ftmp = fi->second;
            ftmp[3] = norm(ftmp[0], ftmp[1], ftmp[2]);
            SFi[ai] += ftmp;
            SFi2[ai] += ftmp.square();
            std::cout << k << ' ' << ai << ' ' << ftmp << ' ' << ftmp.square() << std::endl;
        }
        for (DMapVec4::iterator fi = Fij.begin(); fi != Fij.end(); ++fi)
        {
            ai = fi->first;
            for (MapVec4::iterator fij = fi->second.begin(); fij != fi->second.end(); ++fij)
            {
                aj = fij->first;
                ftmp = fij->second;
                ftmp[3] = norm(ftmp[0], ftmp[1], ftmp[2]);
                SFij[ai][aj] += ftmp;
                SFij2[ai][aj] += ftmp.square();
                SFiFj[ai][aj] += matmul(Fi[ai], Fi[aj]);
            }
            fi->second.clear();
        }
        Fij.clear();
        Fi.clear();
    }
    
    real avgc = 1.0 / nframe;
    real stdc = std::sqrt(1.0 / (nframe - 1));
    rvec4 avgf, stdf;
    for (const std::pair<int32_t, rvec4> &sfi : SFi)
    {
        ai = sfi.first;
        avgf = sfi.second * avgc;
        stdf = (SFi2[ai] - avgf.square() * nframe).sqrt() * stdc;
        std::cout << ai << avgf << stdf << std::endl;
    }
    
    for (const std::pair<int32_t, MapVec4> &sfi : SFij)
    {
        ai = sfi.first;
        for (const std::pair<int32_t, rvec4> &sfij : sfi.second)
        {
            aj = sfij.first;
            avgf = sfij.second * avgc;
            stdf = (SFij2[ai][aj] - avgf.square() * nframe).sqrt() * stdc;
            std::cout << ai << ',' << aj << avgf << stdf << SFiFj[ai][aj] * avgc << std::endl;
        }
    }
    
    return 0;
}
