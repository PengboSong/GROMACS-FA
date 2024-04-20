#!/usr/bin/env python3

from __future__ import print_function
import os
import math
from typing import List

import numpy as np
import pandas as pd
from pymol import cgo, cmd


def check_sele(sele):
    xyz = cmd.get_coords(sele, 1)
    if xyz is None:
        print("[WARNING] Selection " + sele + " has no atoms. Skipped.")
        return None
    else:
        if xyz.shape[0] > 1:
            print("[WARNING] Selection " + sele +
                  " has multiple atoms. Automatically choose the first atom in selection.")
        return xyz[0]


def draw_arrow(xyz1: List[float], xyz2: List[float], radius: float=0.2, cap: float=0.3, gap: float=0.1, rgbcolor: List[float]=[1.0, 0.0, 0.0]):
    vec = xyz2 - xyz1
    dist = math.sqrt(vec.dot(vec))
    gap = min(gap, .05 * dist)
    gapvec = vec * (gap / dist)
    cap = min(cap, .1 * dist)
    capvec = vec * (cap / dist)

    tail = xyz1 + gapvec
    middle = xyz2 - gapvec - capvec
    head = xyz2 - gapvec

    objs = []
    # CYLINDER: SAUSAGE, XYZ1, XYZ2, RADIUS, RGBCOLOR1, RGBCOLOR2
    objs.extend([cgo.SAUSAGE, *tail, *middle, radius, *rgbcolor, *rgbcolor])
    # CAP: CONE, XYZ1, XYZ2, RADIUS1, RADIUS2, RGBCOLOR1, RGBCOLOR2, CAPS1, CAPS2
    objs.extend([cgo.CONE, *middle, *head, radius * 1.5, 0., *rgbcolor, *rgbcolor, 1., 0.])
    return objs


def draw_ellipsoid(xyz: List[float], scale: float, lam: List[float], pc1: List[float], pc2: List[float], pc3: List[float], rgbcolor: List[float]):
    coeff1 = lam[0] / lam[0] / np.linalg.norm(pc1)
    coeff2 = lam[1] / lam[0] / np.linalg.norm(pc2)
    coeff3 = lam[2] / lam[0] / np.linalg.norm(pc3)
    size = lam[0] * scale
    for i in range(3):
        pc1[i] *= coeff1
        pc2[i] *= coeff2
        pc3[i] *= coeff3

    # COLOR: COLOR, R, G, B
    # ELLIPSOID: COLOR, R, G, B, ELLIPSOID, X, Y, Z, SIZE, PC1X, PC1Y, PC1Z, PC2X, PC2Y, PC2Z, PC3X, PC3Y, PC3Z
    return [cgo.COLOR, *rgbcolor, cgo.ELLIPSOID, *xyz, size, *pc1, *pc2, *pc3]


def draw_cylinder(xyz1: List[float], xyz2: List[float], radius: float=0.1, gap: float=0.1, rgbcolor: List[float]=[1.0, 0.0, 0.0]):
    vec = xyz2 - xyz1
    dist = math.sqrt(vec.dot(vec))
    gap = min(gap, .05 * dist)
    gapvec = vec * (gap / dist)

    tail = xyz1 + gapvec
    head = xyz2 - gapvec
    # CYLINDER: SAUSAGE, XYZ1, XYZ2, RADIUS, RGBCOLOR1, RGBCOLOR2
    return [cgo.SAUSAGE, *tail, *head, radius, *rgbcolor, *rgbcolor]
    

def draw_pairwise_force(forcedata: str, grpAunit: str = 'residue', grpBunit: str = 'residue', grpAoffset: int = 1, grpBoffset: int = 1, minforce: float = 100.0, grppairs: int = 200, mindseq: int = 1, width: float = 1.0):
    grpAoffset = int(grpAoffset)
    grpBoffset = int(grpBoffset)
    minforce = float(minforce)
    grppairs = int(grppairs)
    mindseq = int(mindseq)
    width = float(width)

    if not os.path.isfile(forcedata):
        print("Can not open force data file %s." % forcedata)

    df = pd.read_csv(forcedata, sep=',', header=0)
    df["dij"] = df.apply(lambda row: abs(row.i - row.j), axis=1)
    df["fijprojval"] = df.apply(lambda row: abs(row["<FijProj>"]), axis=1)
    # Pairwise forces have fij = -fji, only one of (i,j) and (j,i) is required
    df = df.loc[df.i < df.j].loc[df.fijprojval > minforce].loc[df.dij >= mindseq]
    # Sort pairwise forces by projected force in descending order
    df = df.sort_values(by=["fijprojval"], ascending=False)
    maxf = df.fijprojval.max()
    lstA = list(set(map(int, df.i.values)))
    lstB = list(set(map(int, df.j.values)))
    
    npairs = len(df)   # Total number of force vector pairs to plot
    print(f"Number of pairwise forces > {minforce:.1f} pN: {npairs}")

    # Map force data residue ids to pymol selections
    if grpAunit == "residue":
        mapA = dict(zip(lstA, ['resi %d and name CA' % (ri + grpAoffset) for ri in lstA]))
    elif grpAunit == "atom":
        mapA = dict(zip(lstA, ['index %d' % (ai + grpAoffset) for ai in lstA]))

    if grpBunit == "residue":
        mapB = dict(zip(lstB, ['resi %d and name CA' % (ri + grpBoffset) for ri in lstB]))
    elif grpBunit == "atom":
        mapB = dict(zip(lstB, ['index %d' % (ai + grpBoffset) for ai in lstB]))

    grpn = npairs // grppairs
    if npairs % grppairs != 0:
        grpn += 1
    print(f"Number of force groups: {grpn}")
    cnt, gi = 0, 0
    cylobjs = []
    for _, row in df.iterrows():
        i = mapA[int(row.i)]
        j = mapB[int(row.j)]
        f = row["<FijProj>"]
        r = width * abs(f) / maxf
        # Plot pairwise force fij from res i to j
        xyzi = check_sele(i)
        xyzj = check_sele(j)
        if xyzi is not None and xyzj is not None:
            cylcgo = draw_cylinder(xyzi, xyzj, radius=r, rgbcolor=[1., 0., 0.] if f < 0 else [0., 0., 1.])
            if cylcgo:
                cylobjs.extend(cylcgo)
        cnt += 1
        if cnt >= grppairs:
            cnt = 0
            gi += 1
            cmd.load_cgo(cylobjs, f"Force_Group_{gi}")
            cylobjs.clear()
    if cnt != 0:
        gi += 1
        cmd.load_cgo(cylobjs, f"Force_Group_{gi}")


cmd.extend("draw_pairwise_force", draw_pairwise_force)
