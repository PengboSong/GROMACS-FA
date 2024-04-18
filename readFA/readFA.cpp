/*
    readFA.cpp
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2020/09/30
    Description: Core module of readFA.
*/

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <utility>

#include "FAsettings.h"
#include "FAotypes.h"
#include "FAreader.h"

int main(int argc, char **argv)
{
    FAargs args;
    args.otype = OutputType::SummedForceOnly;
    args.g_lowpf = 0.1;
    args.sblock = 0;
    args.eblock = -1;
    args.ndigits = 3;
    args.forcerot = true;
    args.avgreq = true;
    args.crossanal = true;
    for (int i = 0; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg == "-f" && i < argc - 1)
            args.infnm = argv[++i];
        else if (arg == "-o" && i < argc - 1)
            args.outfnm = argv[++i];
        else if (arg == "--ref" && i < argc - 1)
            args.refnm = argv[++i];
        else if (arg == "--map" && i < argc - 2)
        {
            std::string mapfnmA(argv[++i]);
            std::string mapfnmB(argv[++i]);
            args.mapfnm = std::make_pair(mapfnmA, mapfnmB);
        }
        else if (arg == "--coord" && i < argc - 2)
        {
            std::string coordfnmA(argv[++i]);
            std::string coordfnmB(argv[++i]);
            args.coordfnm = std::make_pair(coordfnmA, coordfnmB);
        }
        else if (arg == "--block" && i < argc - 2)
        {
            args.sblock = std::stoi(argv[++i]);
            args.eblock = std::stoi(argv[++i]);
        }
        else if (arg == "--out" && i < argc - 1)
        {
            std::string otypestr = argv[++i];
            if (otypestr == "raw")
                args.otype = OutputType::OutRawData;
            else if (otypestr == "summed")
                args.otype = OutputType::SummedForceOnly;
            else if (otypestr == "pairwise")
                args.otype = OutputType::PairwiseForceOnly;
            else if (otypestr == "coord")
                args.otype = OutputType::CoordOnly;
            else if (otypestr == "forcecoord")
                args.otype = OutputType::ForceCoordOnly;
            else if (otypestr == "force-graph")
                args.otype = OutputType::ForceGraphOnly;
            else if (otypestr == "force")
                args.otype = OutputType::OutForce;
            else if (otypestr == "no-pairwise")
                args.otype = OutputType::NoPairwiseForce;
            else if (otypestr == "no-forcecoord")
                args.otype = OutputType::NoForceCoord;
            else if (otypestr == "all")
                args.otype = OutputType::OutEverything;
        }
        else if (arg == "--force-graph-edge" && i < argc - 1)
            args.g_lowpf = std::stof(argv[++i]);
        else if (arg == "--out-digits" && i < argc - 1)
            args.ndigits = std::stoi(argv[++i]);
        else if (arg == "--no-forcerot")
            args.forcerot = false;
        else if (arg == "--no-prerun")
            args.avgreq = false;
        else if (arg == "--no-crossanal")
            args.crossanal = false;
    }

    FAReader readFA(args);

    return 0;
}
