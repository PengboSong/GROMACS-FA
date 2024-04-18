/*
    FAanalsys.cpp
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2023/07/11
    Description: Collections of force and coordinate analysis module.
*/

#include "FAanalsys.h"
#include "FAutility.h"

AnalSystem::AnalSystem(const FAargs& inargs, std::string tag)
 : args(inargs), prefix(tag),
   FiAnal(AtomForceAnal(args, true)), FjAnal(args, false), FijAnal(args),
   FiFjAnal(&FiAnal, &FjAnal), SjFijAnal(&FijAnal), GAnal(&FijAnal, args.g_lowpf)
{
    if (!args.otype.raw())
    {
        initvec(Fi, args.N.first);
        initvec(Fj, args.N.second);
        if (args.avgreq)
        {
            initvec(avgFi, args.N.first);
            initvec(avgFj, args.N.second);
        }
    }
}

void AnalSystem::add_detailed_force(int32_t i, int32_t j, rvec4 &force)
{
    if (!args.M2V.first.count(i))
        throw std::runtime_error("Node " + std::to_string(i) + " can not be found in group 1.");
    
    if (!args.M2V.second.count(j))
        throw std::runtime_error("Node " + std::to_string(j) + " can not be found in group 2.");

    Fi[args.M2V.first.at(i)] += force;
    Fj[args.M2V.second.at(j)] -= force;
    if (args.otype.pf()) Fij[i][j] = force;
}

void AnalSystem::avg_forces(const rmat3 &U)
{
    if (args.forcerot)
    {
        rot_force(U, Fi);
        rot_force(U, Fj);
        rot_force(U, Fij);
    }

    for (uint32_t i = 0; i < args.N.first; ++i)
    {
        vecnorm(Fi[i]);
        avgFi[i] += Fi[i];
    }
    for (uint32_t j = 0; j < args.N.second; ++j)
    {
        vecnorm(Fj[j]);
        avgFj[j] += Fj[j];
    }
    for (auto& fi : Fij)
        for (auto& fij : fi.second)
        {
            vecnorm(fij.second);
            avgFij[fi.first][fij.first] += fij.second;
        }
}

void AnalSystem::digest(VecVec4 *avgRi, VecVec4 *avgRj, DVecVec4 *avgRij)
{
    for (uint32_t i = 0; i < args.N.first; ++i) avgFi[i] *= args.avgc;
    for (uint32_t j = 0; j < args.N.second; ++j) avgFj[j] *= args.avgc;
    for (auto& avgfi : avgFij)
        for (auto& avgfij : avgfi.second)
            avgfij.second *= args.avgc;

    FiAnal.digest(&avgFi, avgRi);
    FjAnal.digest(&avgFj, avgRj);
    FijAnal.digest(&avgFij, avgRij);
}

void AnalSystem::update(const rmat3& U, const VecVec4* Ri, const VecVec4* Rj)
{
    if (args.forcerot)
    {
        rot_force(U, Fi);
        rot_force(U, Fj);
        rot_force(U, Fij);
    }

    for (uint32_t i = 0; i < args.N.first; ++i) vecnorm(Fi[i]);
    for (uint32_t j = 0; j < args.N.second; ++j) vecnorm(Fj[j]);
    for (auto& fi : Fij)
        for (auto& fij : fi.second)
            vecnorm(fij.second);
    
    FiAnal.update(&Fi, Ri);
    FjAnal.update(&Fj, Rj);
    FijAnal.update(&Fij, Ri, Rj);
    FiFjAnal.update();
}

void AnalSystem::clear()
{
    // Clear reused force data containers
    std::fill(Fi.begin(), Fi.end(), rvec4());
    std::fill(Fj.begin(), Fj.end(), rvec4());
    if (args.otype.pf())
    {
        DMapVec4 emptyFij;
        Fij.swap(emptyFij);
    }
}

void AnalSystem::write_anal()
{
    printf("Writing %s series of force analysis results to file.\n", prefix.c_str());
    FiAnal.tofile(modfnm(args.outfnm, prefix, "_group1"));
    FjAnal.tofile(modfnm(args.outfnm, prefix, "_group2"));
    FijAnal.tofile(modfnm(args.outfnm, prefix, "_pairwise"));
    FiFjAnal.tofile(modfnm(args.outfnm, prefix, "_corr"));
    SjFijAnal.tofile(modfnm(args.outfnm, prefix, "_summed-pairwise"));
    GAnal.tofile(modfnm(args.outfnm, prefix, "_graph"));
}
