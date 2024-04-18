#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <unordered_set>

using GroupIdx = std::unordered_set<uint32_t>;

void convert_listed_forces(std::istream &is, std::ostream &os, GroupIdx grp)
{
    uint32_t frameid = 0;
    uint32_t forces = 0;
    uint32_t outforces = 0;
    int32_t i, j;
    std::streamsize forcebytesn = 17;
    std::vector<char> forcebytes(forcebytesn, 0);

    while (is.good())
    {
        outforces = 0;
        is.read((char *)&frameid, sizeof(uint32_t));
        if (!is.good()) break;
        is.read((char *)&forces, sizeof(uint32_t));
        os.write((char *)&frameid, sizeof(uint32_t));
        std::streampos rewriteloc = os.tellp();
        os.write((char *)&outforces, sizeof(uint32_t));
        for (uint32_t fi = 0; fi < forces; ++fi)
        {
            is.read((char *)&i, sizeof(int32_t));
            is.read((char *)&j, sizeof(int32_t));
            is.read(forcebytes.data(), forcebytesn);
            
            if (grp.count(i) && grp.count(j))
            {
                os.write((char *)&i, sizeof(int32_t));
                os.write((char *)&j, sizeof(int32_t));
                os.write(forcebytes.data(), forcebytesn);
                ++outforces;
            }
        }
        std::streampos endloc = os.tellp();
        os.seekp(rewriteloc);
        os.write((char *)&outforces, sizeof(uint32_t));
        os.seekp(endloc);
    }
}

void convert_atom_forces(std::istream &is, std::ostream &os, GroupIdx grp)
{
    uint32_t frameid = 0;
    uint32_t forces = 0;
    uint32_t outforces = grp.size();
    std::streamsize forcebytesn = 16;
    std::vector<char> forcebytes(forcebytesn, 0);

    while (is.good())
    {
        is.read((char *)&frameid, sizeof(uint32_t));
        if (!is.good()) break;
        is.read((char *)&forces, sizeof(uint32_t));
        os.write((char *)&frameid, sizeof(uint32_t));
        os.write((char *)&outforces, sizeof(uint32_t));
        for (uint32_t ai = 0; ai < forces; ++ai)
        {
            is.read(forcebytes.data(), forcebytesn);            
            if (grp.count(ai))
                os.write(forcebytes.data(), forcebytesn);
        }
    }
}

void convert_force_data(const std::string infnm, const std::string outfnm, const std::string grpfnm)
{
    GroupIdx grp;
    uint8_t filecode = 0;

    // Read group information
    std::ifstream ingrp(grpfnm, std::ios::binary);
    if (ingrp.is_open())
    {
        uint32_t atomn = 0;
        ingrp.read((char *)&atomn, sizeof(uint32_t));
        std::vector<uint32_t> grplst(atomn, 0);
        ingrp.read((char *)grplst.data(), sizeof(uint32_t) * atomn);
        grp = GroupIdx(grplst.begin(), grplst.end());
    }

    // Prepare input file stream
    std::ifstream indata(infnm, std::ios::binary);
    if (!indata.is_open())
        throw std::runtime_error("Can not read data from file " + infnm);

    std::ofstream outdata(outfnm, std::ios::binary);
    if (!outdata.is_open())
        throw std::runtime_error("Can not write data to file " + outfnm);

    // Read the first byte from file as filecode
    indata.read((char *)&filecode, sizeof(uint8_t));
    outdata.write((char *)&filecode, sizeof(uint8_t));
    switch (filecode)
    {
        case 1:
            convert_listed_forces(indata, outdata, grp);
            break;
        case 2:
            // DO NOTHING
            break;
        case 3:
            convert_listed_forces(indata, outdata, grp);
            break;
        case 4:
            convert_atom_forces(indata, outdata, grp);
            break;
        default:
            // DO NOTHING
            break;
    }
}

int main(int argc, char **argv)
{
    std::string infnm, outfnm, grpfnm;
    bool atomlist = false;
    for (int i = 0; i < argc; ++i)
    {
        std::string arg(argv[i]);
        if (arg == "-f" && i < argc - 1)
            infnm = argv[++i];
        else if (arg == "-o" && i < argc - 1)
            outfnm = argv[++i];
        else if (arg == "-g" && i < argc - 1)
            grpfnm = argv[++i];
    }
    convert_force_data(infnm, outfnm, grpfnm);

    return 0;
}
