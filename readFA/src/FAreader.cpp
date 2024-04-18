/*
    FAreader.cpp
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2021/08/18
    Description: Core module of readFA.
*/

#include "FAreader.h"
#include "FAanalsys.h"
#include "FAutility.h"

FAReader::FAReader(FAargs& inargs)
 : FAsettings(inargs)
{
    if (args.otype.coord() || args.otype.fcoord()) fit_frames();

    /*
    std::ofstream outcoord;
    outfile(outcoord, "debug_coord.csv", "rotation & coordinate matrix");
    outcoord << std::setiosflags(std::ios::fixed) << std::setprecision(6);
    outcoord << mtxheader("U");
    for (uint32_t i = 0; i < args.N.first; ++i) outcoord << ",x" << i << ",y" << i << ",z" << i;
    outcoord << std::endl;
    for (uint32_t k = args.sblock; k < args.eblock; ++k)
    {
        outcoord << k << U[k];
        for (uint32_t i = 0; i < args.N.first; ++i) outcoord << mergeRi[k][i];
        outcoord << std::endl;
    }
    outcoord.close();
    return;
    */
    
    switch (args.filecode)
    {
    case 1:        
        read_listed_forces();
        break;
    case 2:
        read_detailed_forces();
        break;
    case 3:
        read_listed_forces();
        break;
    case 4:
        read_atom_forces();
        break;
    default:
        // Do Nothing
        break;
    }
}

void FAReader::read_listed_forces()
{
    uint32_t frameid = 0;
    uint32_t forces = 0;
    uint64_t current_blockloc = 0;

    int32_t i, j;
    InteractionType itype;
    real fx, fy, fz, f;
    
    AnalSystem Anal(args);
    VecVec4 Ri, Rj;
    if (!args.otype.raw() && (args.otype.coord() || args.otype.fcoord()))
    {
        initvec(Ri, args.N.first);
        if (coordsecflag) initvec(Rj, args.N.second);
    }

    std::ofstream outdata;
    if (args.otype.raw())
    {
        outfile(outdata, modfnm(args.outfnm, "", "_raw"), "raw pairwise data");
        outdata << "frame,i,j,fx,fy,fz,f,type";
        if (coordflag) outdata << ",xi,yi,zi,xj,yj,zj";
        outdata << std::endl;
        outdata << std::setiosflags(std::ios::fixed) << std::setprecision(args.ndigits);
    }

    if (args.avgreq)
    {
        printf("Prerun stage for average forces.\n");
        for (uint32_t k = args.sblock; k < args.eblock; ++k)
        {
            // Locate block start address of block k
            std::tie(frameid, forces, current_blockloc) = blocklocs[k];
            printf("Reading block %5d, with force number = %9d\n", k, forces);
            indata.seekg(current_blockloc);
            for (uint32_t fnum = 0; fnum < forces; ++fnum)
            {
                indata.read((char *)&i, sizeof(int32_t));
                indata.read((char *)&j, sizeof(int32_t));
                indata.read((char *)&itype, sizeof(InteractionType));
                indata.read((char *)&fx, sizeof(real));
                indata.read((char *)&fy, sizeof(real));
                indata.read((char *)&fz, sizeof(real));
                indata.read((char *)&f, sizeof(real));

                rvec4 force({fx, fy, fz, f});
                force *= pN;   // Convert to unit pN
                
                Anal.add_detailed_force(i, j, force);
            }
            Anal.avg_forces(U[k]);
            Anal.clear();
        }
        if (args.otype.coord() || args.otype.fcoord()) Anal.digest(&avgRi, &(coordsecflag ? avgRj : avgRi), &avgRij);
        else Anal.digest(nullptr, nullptr, nullptr);
    }

    printf("\nReading forces for detailed analysis.\n");
    for (uint32_t k = args.sblock; k < args.eblock; ++k)
    {
        // Locate block start address of block k
        std::tie(frameid, forces, current_blockloc) = blocklocs[k];
        printf("Reading block %5d, with force number = %9d\n", k, forces);
        indata.seekg(current_blockloc);
        for (uint32_t fnum = 0; fnum < forces; ++fnum)
        {
            indata.read((char *)&i, sizeof(int32_t));
            indata.read((char *)&j, sizeof(int32_t));
            indata.read((char *)&itype, sizeof(InteractionType));
            indata.read((char *)&fx, sizeof(real));
            indata.read((char *)&fy, sizeof(real));
            indata.read((char *)&fz, sizeof(real));
            indata.read((char *)&f, sizeof(real));

            rvec4 force({fx, fy, fz, f});
            force *= pN;   // Convert to unit pN

            if (args.otype.raw())
            {
                outdata << k << ',' << i << ',' << j << force << ',' << static_cast<int>(itype);
                if (coordflag)
                    outdata << coordpri[k][i] << (coordsecflag ? coordsec : coordpri)[k][j];
                outdata << std::endl;
                continue;
            }
            
            Anal.add_detailed_force(i, j, force);
        }
        if (args.otype.raw()) continue;
        if (args.otype.coord() || args.otype.fcoord())
        {
            read_frame(k, &Ri, &Rj);
            Anal.update(U[k], &Ri, &(coordsecflag ? Rj : Ri));
        }
        else Anal.update(U[k], nullptr, nullptr);
        Anal.clear();
    }
    if (args.otype.raw())
    {
        outdata.close();
        return;
    }
    Anal.write_anal();
}

void FAReader::read_detailed_forces()
{
    uint32_t frameid = 0, forces = 0;
    uint64_t current_blockloc = 0;

    int32_t i, j;
    rvec4 ftmp;
    ifvec fvec;
    
    uint64_t itype;
    std::vector<AnalSystem> Collections(Interact_ANAL, args);
    std::vector<std::string> tags = {"bond_", "polar_", "angle_", "dihedral_", "1-4_", "coulomb_", "vdw_", "bonded_", "nonbonded_", "total_"};
    for (itype = 0; itype < Interact_ANAL; ++itype)
        Collections[itype].reset_tag(tags[itype]);

    VecVec4 Ri, Rj;
    if (!args.otype.raw() && args.otype.coord())
    {
        initvec(Ri, args.N.first);
        if (coordsecflag) initvec(Rj, args.N.second);
    }

    std::ofstream outdata;
    if (args.otype.raw())
    {
        outfile(outdata, modfnm(args.outfnm, "", "_raw"), "raw pairwise data");
        outdata << "frame,i,j,bond:fx,bond:fy,bond:fz,bond:f,polar:fx,polar:fy,polar:fz,polar:f,angle:fx,angle:fy,angle:fz,angle:f,dihedral:fx,dihedral:fy,dihedral:fz,dihedral:f,1-4:fx,1-4:fy,1-4:fz,1-4:f,coulomb:fx,coulomb:fy,coulomb:fz,coulomb:f,vdw:fx,vdw:fy,vdw:fz,vdw:f,fx,fy,fz,f";
        if (coordflag) outdata << ",xi,yi,zi,xj,yj,zj";
        outdata << std::endl;
        outdata << std::setiosflags(std::ios::fixed) << std::setprecision(args.ndigits);
    }

    if (args.avgreq)
    {
        printf("Prerun stage for average forces.\n");
        for (uint32_t k = args.sblock; k < args.eblock; ++k)
        {
            // Locate block start address of block k
            std::tie(frameid, forces, current_blockloc) = blocklocs[k];
            printf("Reading block %5d, with force number = %9d\n", k, forces);
            indata.seekg(current_blockloc);
            for (uint32_t fnum = 0; fnum < forces; ++fnum)
            {
                indata.read((char *)&i, sizeof(int32_t));
                indata.read((char *)&j, sizeof(int32_t));
                indata.read((char *)fvec.data(), sizeof(real) * Interact_FORCEVEC_LEN);
                fvec *= pN;

                rvec4 bonded, nonbonded, total;
                for (itype = 0; itype < Interact_COUNT; ++itype)
                {
                    rvec4 force = {fvec[4 * itype + XX], fvec[4 * itype + YY], fvec[4 * itype + ZZ], fvec[4 * itype + XYZ]};
                    (itype < 4 ? bonded : nonbonded) += force;
                    total += force;
                    Collections[itype].add_detailed_force(i, j, force);
                }
                Collections[itype++].add_detailed_force(i, j, bonded);
                Collections[itype++].add_detailed_force(i, j, nonbonded);
                Collections[itype++].add_detailed_force(i, j, total);
            }
            for (auto Anal = Collections.begin(); Anal != Collections.end(); ++Anal)
            {
                Anal->avg_forces(U[k]);
                Anal->clear();
            }
        }
        for (auto Anal = Collections.begin(); Anal != Collections.end(); ++Anal)
        {
            if (args.otype.coord() || args.otype.fcoord()) Anal->digest(&avgRi, &(coordsecflag ? avgRj : avgRi), &avgRij);
            else Anal->digest(nullptr, nullptr, nullptr);
        }
    }
    
    printf("\nReading forces for detailed analysis.\n");
    for (uint32_t k = args.sblock; k < args.eblock; ++k)
    {
        // Locate block start address of block k
        std::tie(frameid, forces, current_blockloc) = blocklocs[k];
        printf("Reading block %5d, with force number = %9d\n", k, forces);
        indata.seekg(current_blockloc);
        for (uint32_t fnum = 0; fnum < forces; ++fnum)
        {
            indata.read((char *)&i, sizeof(int32_t));
            indata.read((char *)&j, sizeof(int32_t));
            indata.read((char *)fvec.data(), sizeof(real) * Interact_FORCEVEC_LEN);
            fvec *= pN;

            if (args.otype.raw())
            {
                outdata << k << ',' << i << ',' << j << fvec;
                if (coordflag)
                    outdata << coordpri[k][i] << (coordsecflag ? coordsec : coordpri)[k][j];
                outdata << std::endl;
                continue;
            }

            rvec4 bonded, nonbonded, total;
            for (itype = 0; itype < Interact_COUNT; ++itype)
            {
                rvec4 force = {fvec[4 * itype + XX], fvec[4 * itype + YY], fvec[4 * itype + ZZ], fvec[4 * itype + XYZ]};
                (itype < 4 ? bonded : nonbonded) += force;
                total += force;
                Collections[itype].add_detailed_force(i, j, force);
            }
            Collections[itype++].add_detailed_force(i, j, bonded);
            Collections[itype++].add_detailed_force(i, j, nonbonded);
            Collections[itype++].add_detailed_force(i, j, total);
        }
        if (args.otype.raw()) continue;
        if (args.otype.coord() || args.otype.fcoord()) read_frame(k, &Ri, &Rj);
        for (auto Anal = Collections.begin(); Anal != Collections.end(); ++Anal)
        {
            if (args.otype.coord() || args.otype.fcoord()) Anal->update(U[k], &Ri, &(coordsecflag ? Rj : Ri));
            else Anal->update(U[k], nullptr, nullptr);
            Anal->clear();
        }
    }
    if (args.otype.raw())
    {
        outdata.close();
        return;
    }
    for (auto Anal = Collections.begin(); Anal != Collections.end(); ++Anal)
        Anal->write_anal();
}

void FAReader::read_atom_forces()
{
    uint32_t frameid = 0;
    uint32_t atomn = 0;
    uint64_t current_blockloc = 0;

    char *line;

    rvec4 force;

    VecVec4 avgFi, Fi, Ri;

    AtomForceAnal FiAnal(args, true);
    CrossAnal FiFjAnal(&FiAnal, &FiAnal);

    std::ofstream outdata;
    if (args.otype.raw())
    {
        outfile(outdata, modfnm(args.outfnm, "", "_raw"), "raw atom data");
        outdata << "frame,i,fx,fy,fz,f";
        if (coordflag) outdata << ",x,y,z";
        outdata << std::endl;
        outdata << std::setiosflags(std::ios::fixed) << std::setprecision(args.ndigits);
    }
    else
    {
        initvec(Fi, args.N.first);
        if (args.avgreq) initvec(avgFi, args.N.first);
        if (args.otype.coord()) initvec(Ri, args.N.first);
    }

    if (args.avgreq)
    {
        printf("Prerun stage for average forces.\n");
        for (uint32_t k = args.sblock; k < args.eblock; ++k)
        {
            std::tie(frameid, atomn, current_blockloc) = blocklocs[k];
            printf("Reading block %5d, with force number = %9d\n", k, atomn);
            if (atomn > args.N.first)
                throw std::runtime_error("Got mismatched nodes number from mapping data and force file.");
            indata.seekg(current_blockloc);
            for (uint32_t a = 0; a < atomn; ++a)
            {
                indata.read((char *)force.data(), sizeof(real) * 4);
                force *= pN;   // Convert to unit pN
                Fi[a] = force;   // Forces assigned to vectors instead of accumulating, NO clear
            }
            if (args.forcerot) rot_force(U[k], Fi);
            for (uint32_t a = 0; a < atomn; ++a) avgFi[a] += Fi[a] * args.avgc;
        }
        FiAnal.digest(&avgFi, (args.otype.coord() ? &avgRi : nullptr));
    }

    printf("\nReading forces for detailed analysis.\n");
    for (uint32_t k = args.sblock; k < args.eblock; ++k)
    {
        std::tie(frameid, atomn, current_blockloc) = blocklocs[k];
        printf("Reading block %5d, with force number = %9d\n", k, atomn);
        if (atomn > args.N.first)
            throw std::runtime_error("Got mismatched nodes number from mapping data and force file.");
        indata.seekg(current_blockloc);
        for (uint32_t a = 0; a < atomn; ++a)
        {
            indata.read((char *)force.data(), sizeof(real) * 4);
            force *= pN;   // Convert to unit pN

            if (args.otype.raw())
            {
                outdata << k << ',' << a << force;
                if (coordflag) outdata << coordpri[k][a];
                outdata << std::endl;
                continue;
            }

            Fi[a] = force;   // Forces assigned to vectors instead of accumulating, NO clear
        }
        
        if (args.otype.raw()) continue;

        if (args.otype.coord()) read_frame(k, &Ri, nullptr);
        if (args.forcerot) rot_force(U[k], Fi);
        FiAnal.update(&Fi, &Ri);
        FiFjAnal.update();
    }
    if (args.otype.raw())
    {
        outdata.close();
        return;
    }
    FiAnal.tofile(args.outfnm);
    FiFjAnal.tofile(modfnm(args.outfnm, "", "_corr"));
}

void FAReader::read_frame(uint32_t k, VecVec4* Ri, VecVec4* Rj)
{
    if (!args.otype.coord()) return;

    if (Ri != nullptr)
    {
        VecVec3& r  = mergeRi.at(k);
        for (uint32_t i = 0; i < args.N.first; ++i)
        {
            (*Ri)[i][XX]  = r[i][XX];
            (*Ri)[i][YY]  = r[i][YY];
            (*Ri)[i][ZZ]  = r[i][ZZ];
            (*Ri)[i][XYZ] = vecnorm(r[i]);
        }
    }
    if (coordsecflag && Rj != nullptr)
    {
        VecVec3& rr = mergeRj.at(k);
        for (uint32_t j = 0; j < args.N.second; ++j)
        {
            (*Rj)[j][XX]  = rr[j][XX];
            (*Rj)[j][YY]  = rr[j][YY];
            (*Rj)[j][ZZ]  = rr[j][ZZ];
            (*Rj)[j][XYZ] = vecnorm(rr[j]);
        }
    }
}

void FAReader::fit_frames()
{
    /*
    rctr   - Coordinate center of state (Group 1)
    rrctr  - Coordinate center of state (Group 2)
    r0ctr  - Coordinate center of referenced state (Group 1)
    rr0ctr - Coordinate center of referenced state (Group 2)
    r      - Coordinate array of state (Group 1)
    rr     - Coordinate array of state (Group 2)
    r0     - Coordinate array of referenced state (Group 1)
    rr0    - Coordinate array of referenced state (Group 2)
    */
    if (!args.otype.coord() && !args.otype.fcoord()) return;
    
    if (args.avgreq)
    {
        initvec(avgRi, args.N.first);
        if (coordsecflag) initvec(avgRj, args.N.second);
        if (args.otype.pf()) initmat(avgRij, args.N.first, args.N.second);
    }

    rvec3 r0ctr;
    VecVec3& r0 = coordpri.at(args.sblock);

    // Calculate center of coordinate matrix
    for (uint32_t i = 0; i < args.N.first; ++i) r0ctr += r0[i];
    r0ctr *= 1.0 / args.N.first;

    for (int32_t k = args.sblock; k < args.eblock; ++k)
    {
        rvec3 rctr;
        VecVec3& r = coordpri.at(k);
        if (!args.forcerot) mergeRi[k] = r;
        else
        {
            for (uint32_t i = 0; i < args.N.first; ++i) rctr += r[i];
            rctr *= 1.0 / args.N.first;
            rotate(U[k], rctr, r0ctr, r, mergeRi[k]);
        }


        if (!args.avgreq) continue;
        for (uint32_t i = 0; i < args.N.first; ++i)
        {
            avgRi[i][XX] += mergeRi[k][i][XX];
            avgRi[i][YY] += mergeRi[k][i][YY];
            avgRi[i][ZZ] += mergeRi[k][i][ZZ];
        }
    }
    if (args.avgreq)
    {
        for (uint32_t i = 0; i < args.N.first; ++i)
        {
            avgRi[i] *= args.avgc;
            vecnorm(avgRi[i]);
        }
    }

    if (coordsecflag)
    {
        rvec3 rr0ctr;
        VecVec3& rr0 = coordsec.at(args.sblock);

        // Calculate center of coordinate matrix
        for (uint32_t j = 0; j < args.N.second; ++j) rr0ctr += rr0[j];
        rr0ctr *= 1.0 / args.N.second;

        for (int32_t k = args.sblock; k < args.eblock; ++k)
        {
            rvec3 rrctr;
            VecVec3& rr = coordsec.at(k);

            if (!args.forcerot) mergeRj[k] = rr;
            else
            {
                for (uint32_t j = 0; j < args.N.second; ++j) rrctr += rr[j];
                rrctr *= 1.0 / args.N.first;
                rotate(U[k], rrctr, rr0ctr, rr, mergeRj[k]);
            }

            if (!args.avgreq) continue;
            for (uint32_t j = 0; j < args.N.second; ++j)
            {
                avgRj[j][XX] += mergeRj[k][j][XX];
                avgRj[j][YY] += mergeRj[k][j][YY];
                avgRj[j][ZZ] += mergeRj[k][j][ZZ];
            }
        }
        if (args.avgreq)
        {
            for (uint32_t j = 0; j < args.N.second; ++j)
            {
                avgRj[j] *= args.avgc;
                vecnorm(avgRj[j]);
            }
        }
    }

    if (!args.avgreq || !args.otype.pf()) return;
    for (int32_t k = args.sblock; k < args.eblock; ++k)
        for (uint32_t i = 0; i < args.N.first; ++i)
            for (uint32_t j = 0; j < args.N.second; ++j)
            {
                avgRij[i][j][XX] += (coordsecflag ? mergeRj : mergeRi)[k][j][XX] - mergeRi[k][i][XX];
                avgRij[i][j][YY] += (coordsecflag ? mergeRj : mergeRi)[k][j][YY] - mergeRi[k][i][YY];
                avgRij[i][j][ZZ] += (coordsecflag ? mergeRj : mergeRi)[k][j][ZZ] - mergeRi[k][i][ZZ];
            }
    for (uint32_t i = 0; i < args.N.first; ++i)
        for (uint32_t j = 0; j < args.N.second; ++j)
        {
            avgRij[i][j] *= args.avgc;
            vecnorm(avgRij[i][j]);
        }
}
