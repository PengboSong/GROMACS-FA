/*
    FAotypes.cpp
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2021/08/18
    Description: Control output types and data processing functions.
*/

#include "FAotypes.h"

OutputType::OutputType() : c(OutSummedForce)
{
}

OutputType::OutputType(const codetype otype) : c(otype)
{
}

void OutputType::operator=(const codetype otype)
{
    c = otype;
}

void OutputType::operator+=(const codetype otype)
{
    c &= otype;
}

void OutputType::operator-=(const codetype otype)
{
    c -= (c & otype);
}

bool OutputType::raw() const
{
    return c == OutRawData;
}

bool OutputType::atomf() const
{
    return (c & OutSummedForce) != 0;
}

bool OutputType::pf() const
{
    return (c & OutPairwiseForce) != 0;
}

bool OutputType::coord() const
{
    return (c & OutCoord) != 0;
}

bool OutputType::fcoord() const
{
    return (c & OutForceCoord) != 0;
}

bool OutputType::fgraph() const
{
    return (c & OutForceGraph) != 0;
}
