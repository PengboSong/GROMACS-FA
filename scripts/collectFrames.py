#!/usr/bin/env python3
# coding = utf-8


import argparse
import os
import re
import struct


class ElementData:
    """Get mass data from element name.

    Element number, symbol and atom weights from H(1) to U(92) are taken from 
    IUPAC ATOMIC WEIGHTS OF THE ELEMENTS 2021 at https://iupac.qmul.ac.uk/AtWt/.\

    Three decimal places are reserved. If significant digits are insufficient
    for three decimal places, zeros are appended for completion.

    Another term of LP is added with 0 element number and 0.000 atom weight
    for it is usually treated as a virtual particle.
    This term is also preserved for unrecognized elements.
    """

    def __init__(self, checkInitLetterFirst = True):
        self.elements = {}
        self.elementIndexes = {}
        self.checkInitLetterFirst = checkInitLetterFirst

        # Elements detailed information are added to built-in tables
        elementslst = [(0, 'LP', 0.000), (1, 'H', 1.008), (2, 'He', 4.003), (3, 'Li', 6.940), (4, 'Be', 9.012), (5, 'B', 10.810), (6, 'C', 12.011), (7, 'N', 14.007), (8, 'O', 15.999), (9, 'F', 18.998), (10, 'Ne', 20.180), (11, 'Na', 22.990), (12, 'Mg', 24.305), (13, 'Al', 26.982), (14, 'Si', 28.085), (15, 'P', 30.974), (16, 'S', 32.060), (17, 'Cl', 35.450), (18, 'Ar', 39.950), (19, 'K', 39.098), (20, 'Ca', 40.078), (21, 'Sc', 44.956), (22, 'Ti', 47.867), (23, 'V', 50.942), (24, 'Cr', 51.996), (25, 'Mn', 54.938), (26, 'Fe', 55.845), (27, 'Co', 58.933), (28, 'Ni', 58.693), (29, 'Cu', 63.546), (30, 'Zn', 65.382), (31, 'Ga', 69.723), (32, 'Ge', 72.631), (33, 'As', 74.922), (34, 'Se', 78.972), (35, 'Br', 79.904), (36, 'Kr', 83.798), (37, 'Rb', 85.468), (38, 'Sr', 87.621), (39, 'Y', 88.906), (40, 'Zr', 91.224), (41, 'Nb', 92.906), (42, 'Mo', 95.951), (43, 'Tc', 97.000), (44, 'Ru', 101.072), (45, 'Rh', 102.905), (46, 'Pd', 106.420), (47, 'Ag', 107.868), (48, 'Cd', 112.414), (49, 'In', 114.818), (50, 'Sn', 118.711), (51, 'Sb', 121.760), (52, 'Te', 127.603), (53, 'I', 126.904), (54, 'Xe', 131.294), (55, 'Cs', 132.905), (56, 'Ba', 137.328), (57, 'La', 138.905), (58, 'Ce', 140.116), (59, 'Pr', 140.908), (60, 'Nd', 144.242), (61, 'Pm', 145.000), (62, 'Sm', 150.362), (63, 'Eu', 151.964), (64, 'Gd', 157.253), (65, 'Tb', 158.925), (66, 'Dy', 162.500), (67, 'Ho', 164.930), (68, 'Er', 167.259), (69, 'Tm', 168.934), (70, 'Yb', 173.045), (71, 'Lu', 174.967), (72, 'Hf', 178.487), (73, 'Ta', 180.948), (74, 'W', 183.841), (75, 'Re', 186.207), (76, 'Os', 190.233), (77, 'Ir', 192.217), (78, 'Pt', 195.085), (79, 'Au', 196.967), (80, 'Hg', 200.592), (81, 'Tl', 204.380), (82, 'Pb', 207.210), (83, 'Bi', 208.980), (84, 'Po', 209.000), (85, 'At', 210.000), (86, 'Rn', 222.000), (87, 'Fr', 223.000), (88, 'Ra', 226.000), (89, 'Ac', 227.000), (90, 'Th', 232.038), (91, 'Pa', 231.036), (92, 'U', 238.029)]
        self.addElements(elementslst)
    
    def addElement(self, index, element, mass):
        # Uppercase symbols are stored for consistency
        element = element.upper()
        self.elements[element] = dict(index=index, mass=mass)
        self.elementIndexes[index] = dict(element=element, mass=mass)
    
    def addElements(self, infolst):
        for elementinfo in infolst:
            self.addElement(*elementinfo)
    
    def guessElement(self, atomnm):
        # Strip digits in the head or in the tail
        atomnm = re.sub(r'^\d+|\d+$', '', atomnm).upper()
        if atomnm in self.elements:
            return atomnm
        elif self.checkInitLetterFirst:
            if atomnm[0] in self.elements:
                return atomnm[0]
            elif atomnm[:2] in self.elements:
                return atomnm[:2]
        else:
            if atomnm[:2] in self.elements:
                return atomnm[:2]
            elif atomnm[0] in self.elements:
                return atomnm[0]
        # Failed to guess element from given atom name
        print("[WARNING] Failed to guess element from atom name %s." % atomnm)
        # A reserved term with zero atom weight is returned
        return "LP"
    
    def getMass(self, atomnm):
        return self.elements[self.guessElement(atomnm)]["mass"]


GetElementMass = ElementData()


def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--workdir',
        type=str,
        default='.',
        help="Path to the workdir where holds PDB files extracted from GROMACS trajectory"
    )
    parser.add_argument(
        '--filename',
        type=str,
        default='frames',
        help="Filename of output binary packed coordinate file"
    )
    parser.add_argument(
        '--mode',
        type=str,
        default='all',
        help="Frame coordinates extraction method. Should be one of all (extract"
             "coordinates of all atoms), com (calculate COM of each residue/molecules)"
    )

    FLAGS, unparsed = parser.parse_known_args()

    return FLAGS


def extractCoordFromLine(line):
    # Coordinate X at column 30-37
    # Coordinate Y at column 38-45
    # Coordinate Z at column 46-53
    # Coordinate unit conversion: A (PDB) -> nm (GROMACS)
    return [.1 * float(line[30:38]), .1 * float(line[38:46]), .1 * float(line[46:54])]


def readXYZFromPDB(filenm):
    xyzdata = []
    with open(filenm, 'r') as f:
        for ln, line in enumerate(f.readlines()):
            if line[0:6] not in ("ATOM  ", "HETATM"): continue
            xyzdata.append(extractCoordFromLine(line))
    return xyzdata


def calcCOM(resdata, resatoms):
    assert len(resdata) == len(resatoms)
    wx = wy = wz = sw = 0.0
    for (x, y, z), atomnm in zip(resdata, resatoms):
        w = GetElementMass.getMass(atomnm)
        wx += w * x
        wy += w * y
        wz += w * z
        sw += w
    return [wx / sw, wy / sw, wz / sw]


def readCOMFromPDB(filenm):
    comdata = []
    with open(filenm, 'r') as f:
        reskey = (0, '')
        resdata = []
        resatoms = []
        for line in f.readlines():
            if line[0:6] not in ("ATOM  ", "HETATM"): continue
            atomnm = line[12:16].strip()
            resid = int(line[22:26])
            resnm = line[17:20].strip()
            if (reskey[0] != resid) or (reskey[1] != resnm):
                # Start a new residue
                if len(resdata) != 0:
                    comdata.append(calcCOM(resdata, resatoms))
                    resdata = []
                    resatoms = []
            resdata.append(extractCoordFromLine(line))
            resatoms.append(atomnm)
            reskey = (resid, resnm)
        if len(resdata) != 0:
            comdata.append(calcCOM(resdata, resatoms))        
    return comdata  


def packXYZ(filenm, mode = "all"):
    xyzbytes = b''
    N = 0
    # Default `mode` keyword to `all`
    readfunc = readCOMFromPDB if mode == "com" else readXYZFromPDB
    for xyz in readfunc(filenm):
        xyzbytes += struct.pack('fff', *xyz)
        N += 1
    return N, xyzbytes


def main():
    FLAGS = parseArgs()
    assert os.path.isdir(FLAGS.workdir), "The specified workdir is not an available directory."

    frames = []
    framefnms = [f for f in os.listdir(FLAGS.workdir) if os.path.splitext(f)[1] == '.pdb']
    framepattern = re.compile(r'[A-Za-z]*([0-9]+)\.pdb')
    for fnm in framefnms:
        matchres = re.match(framepattern, fnm)
        if matchres:
            frameid = int(matchres.groups()[0])
            frames.append((frameid, os.path.join(FLAGS.workdir, fnm)))
    frames.sort(key=lambda obj: obj[0])
    framen = len(frames)
    globalN = 0
    frameid, framefnm = frames[0]
    globalN, xyzbytes = packXYZ(filenm=framefnm, mode=FLAGS.mode)
    with open(os.path.join(FLAGS.workdir, FLAGS.filename + ".xyz"), 'wb') as fxyz:
        fxyz.write(struct.pack('II', framen, globalN))
        print("Got %d frames and %d atoms/residues/molecules in total" % (framen, globalN))
        print("Writing frame %d/%d from %s ..." % (frameid + 1, framen, framefnm))
        fxyz.write(xyzbytes)
        for frameid, framefnm in frames[1:]:
            N, xyzbytes = packXYZ(filenm=framefnm, mode=FLAGS.mode)
            if (N != globalN):
                raise ValueError("Got unequal atom/residue/molecule number (Got %d, expected %d)" % (N, globalN))
            print("Writing frame %d/%d from %s ..." % (frameid + 1, framen, framefnm))
            fxyz.write(xyzbytes)
    

if __name__ == '__main__':
    main()