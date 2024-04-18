/*
    FAsettings.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2021/08/18
    Description: Settings base class for FAreader.
*/

#ifndef SRC_READFA_FASETTINGS_H_
#define SRC_READFA_FASETTINGS_H_

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <utility>

#include "real.h"
#include "rmsd.h"
#include "rvec.h"
#include "CoordMat.h"
#include "FAdefines.h"
#include "FAmath.h"
#include "FAotypes.h"

using IdxMap = std::map<int32_t, uint32_t>;
using RevIdxMap = std::vector<int32_t>;

static constexpr uint64_t global_headersize = sizeof(uint8_t);
static constexpr uint64_t local_headersize = 2 * sizeof(uint32_t);
static constexpr uint64_t list_pairsize = 2 * sizeof(int32_t) + sizeof(InteractionType) + 4 * sizeof(real);
static constexpr uint64_t detailed_pairsize = 2 * sizeof(int32_t) + Interact_FORCEVEC_LEN * sizeof(real);
static constexpr uint64_t atom_pairsize = 4 * sizeof(real);

typedef struct {
    // Input force data filename
    std::string infnm;

    // Basic output force analysis results filename
    std::string outfnm;

    // Input packed reference structure filename
    std::string refnm;

    // Input atom/residue map filename
    std::pair<std::string, std::string> mapfnm;
    
    // Input atom/residue coordinate data filename
    std::pair<std::string, std::string> coordfnm;
    
    // Force data format identifier
    // 0 - None
    // 1 - Summed Mode
    // 2 - Detailed Mode
    // 3 - Listed Mode
    // 4 - Atom Force Mode
    uint8_t filecode;

    // This identifier is used to determine how results are calculated and displayed
    OutputType otype;

    // This flag determines whether forces are rotated based on fitting
    bool forcerot;

    // This flag determines whether averages of forces and coordinates should be calculated first
    bool avgreq;

    // This flag determines whether memory-cost cross analysis should be performed
    bool crossanal;
    
    // Container mapping atom/residue id to compact vector index
    // Maps for group 1 and group 2 are required
    std::pair<IdxMap, IdxMap> M2V;

    // Reversed map converting vector index to atom id
    // Reversed maps for group 1 and group 2 are required
    std::pair<RevIdxMap, RevIdxMap> V2M;

    // Particles (usually called atom in program) number count for group 1 and group 2
    std::pair<uint32_t, uint32_t> N;

    // Total block length
    uint32_t blockn;

    // Start/End block index
    uint32_t sblock, eblock;

    // Block length in range [sblock, eblock)
    uint32_t K;

    // Coefficient to average forces (=1/K)
    real avgc;

    // Coefficient to calculate sample std. dev. of forces (=1/sqrt(K-1))
    real stdc;

    // Output data precision
    uint32_t ndigits;

    // Lower bound of pairwise forces as edges of force graph
    real g_lowpf;

} FAargs;

class AtomForceAnal;
class PairwiseForceAnal;
class CrossAnal;
class SequenceAnal;

class FAsettings
{
    friend class AtomForceAnal;
    friend class PairwiseForceAnal;
    friend class CrossAnal;
    friend class SequenceAnal;

public:
    FAsettings(FAargs& inargs);

    ~FAsettings();

    bool read_map(std::string mapfnm, IdxMap& idxmap, RevIdxMap& revmap, uint32_t& N);

    bool read_coord(std::string coordfnm, DVecVec3& coordmat, const uint32_t N);

    bool read_refcoord(std::string refnm);

    bool idx_found(IdxMap map, int32_t idx, uint32_t& mapidx);

protected:
    // Collections of arguments
    FAargs args;

    // Input file handler to read force data
    std::ifstream indata;

    // Container for block indexes, block lengthes and block addresses
    std::vector<BlockLoc> blocklocs;

    // FLags for reading coordinates (true for success, false for failed or empty)
    bool coordflag;

    // Container for rotation matrices
    MapMat3 U;

    // Container for primary coordinate matrix
    DVecVec3 coordpri;

    // Container for secondary coordinate matrix
    DVecVec3 coordsec;

    // Flags for whether the secondary coordinate matrix is activated
    bool coordsecflag;

};

#endif /* SRC_READFA_FASETTINGS_H_ */
