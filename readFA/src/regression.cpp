#include "regression.h"
/*
    regression.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2023/09/09
    Description: Functions for analysis of linear regression.
*/

#include "regression.h"

rmat2 DemingRegression::cov(real avgXt, real avgYt, real avgXtsq, real avgYtsq, real avgXtYt)
{
    rmat2 mcov;
    // Population covariance
    mcov[XX][XX] = avgXtsq - avgXt * avgXt;
    mcov[YY][YY] = avgYtsq - avgYt * avgYt;
    mcov[XX][YY] = mcov[YY][XX] = avgXtYt - avgXt * avgYt;
    return mcov;
}

rmat2 DemingRegression::estcov(int n, real avgXt, real delta, real b1, real varXt, real varXterr)
{
    rmat2 estmcov;
    real S = (delta + b1 * b1) * varXterr * (n - 1) / (n - 2);
    estmcov[YY][YY] = (S / std::sqrt(varXt) + S * varXterr / varXt - b1 * b1 * varXterr * varXterr / varXt) / (n - 1);
    estmcov[XX][XX] = S / n + avgXt * avgXt * estmcov[YY][YY];
    estmcov[XX][YY] = estmcov[YY][XX] = -avgXt * estmcov[YY][YY];
    return estmcov;
}

int DemingRegression::fit(int n, real avgXt, real avgYt, real avgXtsq, real avgYtsq, real avgXtYt, real varXt, real varYt, FitResult &res)
{
    int flag = 0;
    res.delta = varYt / varXt;
    rmat2 mcov = DemingRegression::cov(avgXt, avgYt, avgXtsq, avgYtsq, avgXtYt);
    // Sample covariance
    mcov *= 1.0 * n / (n - 1);
    if (mcov[XX][YY] == 0)
    {
        if (mcov[YY][YY] > res.delta * mcov[XX][XX]) return -1;
        else 
        {
            flag = 1;
            res.b1 = 0;
            res.b0 = avgYt;
            res.varYt = mcov[XX][XX] * res.delta - mcov[YY][YY];
            res.varXt = res.varYt / res.delta;
            res.varYterr = mcov[YY][YY];
            res.varXterr = res.varYterr / res.delta;
        }
    }
    else
    {
        real mcovyypdxx = mcov[YY][YY] + mcov[XX][XX] * res.delta;   // COVyy + δ COVxx
        real mcovyymdxx = mcov[YY][YY] - mcov[XX][XX] * res.delta;   // COVyy - δ COVxx
        real DD = std::sqrt(mcovyymdxx * mcovyymdxx + 4 * res.delta * mcov[XX][YY] * mcov[XX][YY]);   // ((COVyy - δ COVxx)^2 + 4δ COVxy^2)^1/2
        res.b1 = (mcovyymdxx + DD) / 2 / mcov[XX][YY];
        res.b0 = avgYt - res.b1 * avgXt;
        res.varXt = -(mcovyymdxx - DD) / 2 / res.delta;
        res.varYt = res.varXt * res.delta;
        res.varXterr = (mcovyypdxx - DD) / 2 / res.delta;
        res.varYterr = res.varXterr * res.delta;
    }

    rmat2 estmcov = estcov(n, avgXt, res.delta, res.b1, res.varXt, res.varXterr);
    res.b0CI95lb = res.b0 - Z95 * estmcov[XX][XX];
    res.b0CI95ub = res.b0 + Z95 * estmcov[XX][XX];
    res.b1CI95lb = res.b1 - Z95 * estmcov[YY][YY];
    res.b1CI95ub = res.b1 + Z95 * estmcov[YY][YY];
    return flag;
}
