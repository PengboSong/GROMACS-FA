/*
    FAotypes.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2021/08/18
    Description: Control output types and data processing functions.
*/

#ifndef SRC_READFA_FAOTYPES_H_
#define SRC_READFA_FAOTYPES_H_

#include <cstdint>

class OutputType
{

using codetype = uint8_t;

public:
    OutputType();

    OutputType(const codetype otype);

    void operator=(const codetype otype);

    void operator+=(const codetype otype);

    void operator-=(const codetype otype);

    // This flag is used to determine whether data processing is disabled
    bool raw() const;

    bool atomf() const;

    // This flag is used to determine whether memory cost pairwise force module is activated
    bool pf() const;

    // This flag is used to determine whether coordinate analysis module is activated
    bool coord() const;

    // This flag is used to determine whether force-coordinate analysis module is activated
    bool fcoord() const;

    // This flag is used to determine whether force-graph module is activated
    bool fgraph() const;

    static const codetype OutRawData         = 0;
    static const codetype OutSummedForce     = 1 << 0;
    static const codetype OutPairwiseForce   = 1 << 1;
    static const codetype OutCoord           = 1 << 2;
    static const codetype OutForceCoord      = 1 << 3;
    static const codetype OutForceGraph      = 1 << 4;
    static const codetype SummedForceOnly    = OutSummedForce;
    static const codetype PairwiseForceOnly  = OutPairwiseForce;
    static const codetype CoordOnly          = OutCoord;
    static const codetype ForceCoordOnly     = OutForceCoord;
    static const codetype ForceGraphOnly     = OutForceGraph;
    static const codetype OutForce           = OutSummedForce | OutPairwiseForce;
    static const codetype NoPairwiseForce    = OutSummedForce | OutCoord | OutForceCoord;
    static const codetype NoForceCoord       = OutSummedForce | OutPairwiseForce | OutCoord;
    static const codetype OutEverything      = OutForce | OutCoord | OutForceCoord | OutForceGraph;

private:
    codetype c;
};


#endif /* SRC_READFA_FAOTYPES_H_ */