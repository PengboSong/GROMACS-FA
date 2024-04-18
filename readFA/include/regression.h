/*
    regression.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2023/09/09
    Description: Functions for analysis of regression.
*/

#ifndef SRC_READFA_REGRESSION_H_
#define SRC_READFA_REGRESSION_H_

#include "real.h"
#include "rmat.h"
#include "FAdefines.h"

typedef struct 
{
    real delta;      // Ratio of variance of Yt to variance of Xt
    real b0;         // Estimated intercept
    real b1;         // Estimated slope
    real varXt;      // Estimated variance of Xt
    real varYt;      // Estimated variance of Yt
    real varXterr;   // Estimated variance of Xt error
    real varYterr;   // Estimated variance of Yt error
    real b0CI95lb;   // Lower bound of 95% confidence interval of intercept
    real b0CI95ub;   // Upper bound of 95% confidence interval of intercept
    real b1CI95lb;   // Lower bound of 95% confidence interval of slope
    real b1CI95ub;   // Upper bound of 95% confidence interval of slope
} FitResult;

static inline std::ostream& operator<<(std::ostream& out, const FitResult& res)
{
    out << ',' << res.delta << ',' << res.b0 << ',' << res.b1;
    out << ',' << res.varXt << ',' << res.varYt << ',' << res.varXterr << ',' << res.varYterr;
    out << ',' << res.b0CI95lb << ',' << res.b0CI95ub << ',' << res.b1CI95lb << ',' << res.b1CI95ub;
    return out;
}

static const real Z90 = 1.645F;
static const real Z95 = 1.96F;
static const real Z99 = 2.58F;

class DemingRegression
{
public:

    static rmat2 cov(real avgXt, real avgYt, real avgXtsq, real avgYtsq, real avgXtYt);

    static rmat2 estcov(int n, real avgXt, real delta, real b1, real varXt, real varXterr);

    static int fit(int n, real avgXt, real avgYt, real avgXtsq, real avgYtsq, real avgXtYt, real varXt, real varYt, FitResult& res);
};

#endif