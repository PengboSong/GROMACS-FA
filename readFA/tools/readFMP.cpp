#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <utility>

using real = float;

void read_atom_map(std::string fnm, bool atomlist)
{
    uint8_t maptype = 0;
    uint32_t atomn = 0;

    std::ifstream indata(fnm, std::ios::binary);
    if (!indata.is_open())
        throw std::runtime_error("Failed to open file.");

    indata.read((char*)&maptype, sizeof(uint8_t));
    indata.read((char*)&atomn, sizeof(uint32_t));
    int32_t* buffer = new int32_t[atomn];
    indata.read((char*)buffer, sizeof(int32_t) * atomn);
    indata.close();
    std::vector<int32_t> container(buffer, buffer + atomn);

    std::string identifier = (maptype == 2) ? "resi." : "mol";
    if (atomlist)
    {
        std::cout << "atomi," + identifier + 'i' << std::endl;
        for (uint32_t ai = 0; ai < container.size(); ++ai)
            std::cout << ai + 1 << ',' << container[ai] << std::endl;
    }
    else
    {
        std::map<int32_t, std::vector<int32_t>> atom_map;
        int32_t ri = 0, ai = 0;
        for (std::size_t i = 0; i < container.size(); ++i)
        {
            ai = static_cast<int32_t>(i);
            ri = container.at(ai);
            if (atom_map.find(ri) != atom_map.end())
                atom_map.at(ri).push_back(ai);
            else
                atom_map[ri] = {ai};
        }
        
        for (const std::pair<int32_t, std::vector<int32_t>>& pair : atom_map)
        {
            std::cout << identifier << ' ' << pair.first << std::endl;
            for (const int32_t& ai : pair.second)
                std::cout << ai << ' ';
            std::cout << std::endl;
        }
    }
    
}

int main(int argc, char **argv)
{
    std::string datafnm;
    bool atomlist = false;
    for (int i = 0; i < argc; ++i)
    {
        std::string arg(argv[i]);
        if (arg == "-f" && i < argc - 1)
            datafnm = argv[++i];
        else if (arg == "--atom-list")
            atomlist = true;
    }
    read_atom_map(datafnm, atomlist);

    return 0;
}
