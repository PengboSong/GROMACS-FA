/*
    FAgraphanal.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2023/07/10
    Description: Time sequence force and coordinate analysis module.
*/

#include <cmath>
#include <tuple>
#include <utility>

#include "FAgraphanal.h"
#include "FAutility.h"

GraphAnal::GraphAnal(const PairwiseForceAnal *pf, const real lowpf)
 : Fij(pf), g_lowpf(lowpf)
{
}

void GraphAnal::write_header(std::ostream &out)
{
    printf("Assign edge weights according to <|Fij|>\n");
    out << "i,j,Wij" << std::endl;
}

void GraphAnal::write_data(std::ostream &out)
{
    real eij, wij;
    double logprod = 0.0;
    uint32_t ecnt = 0;
    int32_t i, j;
    std::vector<std::tuple<int32_t, int32_t, real>> edges;
    for (auto fi = Fij->SFij.cbegin(); fi != Fij->SFij.cend(); ++fi)
    {
        i = fi->first;
        for (auto fij = fi->second.cbegin(); fij != fi->second.cend(); ++fij)
        {                
            j = fij->first;
            eij = fij->second[XYZ] * Fij->avgc;
            if (eij < g_lowpf) continue;
            logprod += std::log(static_cast<double>(eij));
            edges.emplace_back(std::make_tuple(i, j, eij));
            ++ecnt;
        }
    }
    double prod = std::exp(logprod / static_cast<double>(ecnt));
    for (auto& edge : edges)
    {
        std::tie(i, j, eij) = edge;
        wij = prod / eij;
        out << i << ',' << j << ',' << wij << std::endl;
    }
}

void GraphAnal::tofile(std::string fnm)
{
    if ((!Fij->otype.pf() && !Fij->otype.fcoord()) || !Fij->otype.fgraph()) return;

    std::ofstream out;
    outfile(out, fnm, "graph edges data");
    out << std::setiosflags(std::ios::fixed) << std::setprecision(Fij->ndigits);
    write_header(out);
    write_data(out);
    out.close();
}
