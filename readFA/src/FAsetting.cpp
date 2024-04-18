/*
    FAsettings.cpp
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2021/08/18
    Description: Settings base class for FAreader.
*/

#include <algorithm>
#include <iomanip>

#include "FAsettings.h"

FAsettings::FAsettings(FAargs& inargs)
 : args(inargs), coordsecflag(true)
{
    if (args.otype.raw())
    {
        args.forcerot = false;
        args.avgreq = false;
        args.crossanal = false;
    }
    if (args.otype.fgraph() && !args.otype.pf())
    {
        args.otype -= OutputType::OutForceGraph;
        printf("Warning: Pairwise force analysis is required for force graph analysis. Force graph analysis is auto off.\n");
    }

    // Reading mapping data
    if (!read_map(args.mapfnm.first, args.M2V.first, args.V2M.first, args.N.first) || !read_map(args.mapfnm.second, args.M2V.second, args.V2M.second, args.N.second))
        throw std::runtime_error("Mapping data is required for further analysis.");

    // Reading coordinates data
    if (args.coordfnm.first == args.coordfnm.second && args.mapfnm.first == args.mapfnm.second) coordsecflag = false;
    coordflag = false;
    if (args.otype.raw() || args.otype.coord() || args.otype.fcoord())
    {
        // FLags for reading coordinates (true for success, false for failed or empty)
        coordflag = read_coord(args.coordfnm.first, coordpri, args.N.first);
        if (coordsecflag)
            coordflag = coordflag && read_coord(args.coordfnm.second, coordsec, args.N.second);
        if (!coordflag)
        {
            // Disable coord and force-coord results output
            args.otype -= OutputType::OutCoord;
            args.otype -= OutputType::OutForceCoord;
            printf("Warning: No coordinates data. Coordinate analysis is auto off.\n");
        }
    }

    // Prepare input file stream
    indata.open(args.infnm, std::ios::binary);
    if (!indata.is_open())
        throw std::runtime_error("Can not read data from file " + args.infnm);
    
    // Read the first byte from file as filecode
    indata.read((char *)&args.filecode, sizeof(uint8_t));

    // Get total filesize
    int64_t startpos = indata.tellg();
    indata.seekg(0, std::ios::end);
    int64_t endpos = indata.tellg();
    indata.seekg(startpos);
    uint64_t blocklen = 0;
    if (endpos != -1)
        blocklen = endpos - startpos;
    else
        throw std::runtime_error("Failed to read filesize of file " + args.infnm);
    
    // Determine block number and start address of each block
    uint32_t frameid = 0;
    uint32_t forces = 0;
    uint64_t current_blockloc = global_headersize;
    uint64_t headersize = local_headersize;
    uint64_t pairsize = 0;
    switch (args.filecode)
    {
    case 1:
        pairsize = list_pairsize;
        break;
    case 2:
        pairsize = detailed_pairsize;
        break;
    case 3:
        pairsize = list_pairsize;
        break;
    case 4:
        pairsize = atom_pairsize;
        args.otype -= OutputType::OutPairwiseForce;
        args.otype -= OutputType::OutForceGraph;
        break;
    default:
        // TODO
        break;
    }
    while (indata.good())
    {
        indata.read((char *)&frameid, sizeof(uint32_t));
        indata.read((char *)&forces, sizeof(uint32_t));
        current_blockloc += headersize;
        blocklocs.push_back(std::make_tuple(frameid, forces, current_blockloc));
        uint64_t blocksize = forces * pairsize;
        if (blocksize > blocklen)
            throw std::runtime_error("Got too large number of forces " + std::to_string(forces) + " when reading frame " + std::to_string(frameid) + ". The file may be corrupted.");
        current_blockloc += blocksize;
        if (current_blockloc >= blocklen)
            break;
        indata.seekg(current_blockloc);
    }
    args.blockn = blocklocs.size();
    printf("Got %d blocks in force data file.\n", args.blockn);

    // Check force block length and coordinate matrix length
    if (coordflag && (args.blockn != coordpri.size()))
        throw std::runtime_error("Got unequal length of force blocks and coordinate matrices for primary group.");
    if (coordsecflag && coordflag && (args.blockn != coordsec.size()))
        throw std::runtime_error("Got unequal length of force blocks and coordinate matrices for secondary group.");

    // Reset block index range
    int32_t lastframe = args.blockn - 1;
    if (args.sblock < 0)
        args.sblock += args.blockn;
    if (args.eblock < 0)
        args.eblock += args.blockn;
    args.sblock = (lastframe > args.sblock) ? args.sblock : lastframe;
    args.eblock = (lastframe > args.eblock) ? args.eblock : lastframe;
    if (args.sblock > args.eblock)
        std::swap<uint32_t>(args.sblock, args.eblock);
    ++args.eblock;
    args.K = args.eblock - args.sblock;
    args.avgc = 1.0 / args.K;
    args.stdc = std::sqrt(1.0 / (args.K - 1));
    printf("Read frames in range from %d (included) to %d (not included)\n", args.sblock, args.eblock);

    // Fitting reference structures
    if (args.forcerot && !read_refcoord(args.refnm)) args.forcerot = false;
    if (!args.forcerot)
    {
        rmat3 unimat;
        unimat[XX][XX] = unimat[YY][YY] = unimat[ZZ][ZZ] = 1.0;
        for (int32_t k = args.sblock; k < args.eblock; ++k) U[k] = unimat;
    }
}

FAsettings::~FAsettings()
{
    if (indata.is_open()) indata.close();
}

bool FAsettings::read_map(std::string mapfnm, IdxMap &idxmap, RevIdxMap &revmap, uint32_t& N)
{
    // Read atom map
    if (mapfnm.empty()) return false;
    std::ifstream mapstream(mapfnm, std::ios::binary);
    if (!mapstream.is_open()) return false;
    
    uint32_t atomn, ai;
    mapstream.read((char *)&atomn, sizeof(uint32_t));
    printf("Got %d nodes from mapping file %s.\n", atomn, mapfnm.c_str());
    N = atomn;
    revmap.resize(atomn, 0);
    mapstream.read((char *)revmap.data(), sizeof(int32_t) * atomn);
    for (ai = 0; ai < atomn; ++ai)
        idxmap[revmap[ai]] = ai;
    mapstream.close();
    printf("Successfully reading mapping data from file %s.\n", mapfnm.c_str());
    return true;
}

bool FAsettings::read_coord(std::string coordfnm, DVecVec3 &coordmat, const uint32_t N)
{
    // Read coordinate data
    if (coordfnm.empty()) return false;
    std::ifstream coordstream(coordfnm, std::ios::binary);
    if (!coordstream.is_open()) return false;

    uint32_t framen, fid, atomn, ai;
    coordstream.read((char*)&framen, sizeof(uint32_t));
    coordstream.read((char*)&atomn, sizeof(uint32_t));
    printf("Got %d frames and %d nodes from coordinate file %s.\n", framen, atomn, coordfnm.c_str());
    if (N != 0 && atomn != N)
        throw std::runtime_error("Got unequal atom count from atom map and coordinates.");
    initmat(coordmat, framen, atomn);
    for (fid = 0; fid < framen; ++fid)
        for (ai = 0; ai < atomn; ++ai)
            coordstream.read((char*)coordmat[fid][ai].data(), sizeof(float)*DIM);
    coordstream.close();
    printf("Successfully reading coordinate data from file %s.\n", coordfnm.c_str());
    return true;
}

bool FAsettings::read_refcoord(std::string refnm)
{
    if (!args.forcerot) return false;

    DVecVec3 refmat;   // Container for reference structures
    if (!read_coord(args.refnm, refmat, 0))
    {
        printf("Warning: Reference structures are required for force rotation. Force rotation is auto off.\n");
        return false;
    }

    /*
    RMSD - Calculated RMSD
    g    - Calculated gradient
    lam  - Calculated lambda
    q    - Calculated q
    */
    real RMSD = 0., avgRMSD = 0.;
    rvec3 r0ctr, rctr;
    rvec4 lam;
    rmat4 q;
    VecVec3 rg;
    initvec(rg, refmat[args.sblock].size());
    for (int32_t k = args.sblock; k < args.eblock; ++k)
        if (CalRMSD(RMSD, refmat[k], refmat[args.sblock], rctr, r0ctr, lam, q, true, U[k], false, rg) == -1)
        {
            printf("Warning: Fitting reference structures failed at frame %d. Force rotation is auto off.\n", k);
            return false;
        }
        else avgRMSD += RMSD;
    avgRMSD *= args.avgc;
    printf("Average RMSD of fitted reference structures = %.3f A.\n", avgRMSD);
    return true;
}

bool FAsettings::idx_found(IdxMap map, int32_t idx, uint32_t &mapidx)
{
    IdxMap::iterator it = map.find(idx);
    if (it != map.end())
    {
        mapidx = it->second;
        return true;
    }
    else
        return false;
}
