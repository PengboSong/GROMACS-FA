#!/usr/bin/env python3
# coding=utf-8


import argparse
import os

import pandas as pd


def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--data',
        type=str,
        help="Filename of CSV data table with correlation matrix."
    )
    parser.add_argument(
        '--pdb',
        type=str,
        help="Filename of PDB format file to insert ANISOU domains."
    )
    parser.add_argument(
        '--type',
        type=str,
        default="force",
        help="Matrix type name for which to insert into ANISOU domains (force, coord supported)."
    )
    parser.add_argument(
        '--unit',
        type=str,
        default="atom",
        help="Unit type name for which to insert into ANISOU domains (atom, residue supported)."
    )

    FLAGS, unparsed = parser.parse_known_args()

    return FLAGS


def main():
    FLAGS = parseArgs()
    assert os.path.isfile(FLAGS.data), "Input force data file does not exist."
    assert os.path.isfile(FLAGS.pdb), "Input PDB file does not exist."
    assert FLAGS.type in ["force", "coord"]
    assert FLAGS.unit in ["atom", "residue"]
    df = pd.read_csv(FLAGS.data, sep=',', header=0)
    assert (FLAGS.type == "force" and "<DFiDFi>xx" in df.columns) or (FLAGS.type == "coord" and "<DRiDRi>xx" in df.columns)
    c = 1E-2 if FLAGS.type == "force" else 1E4
    mtxheader = "<DFiDFi>" if FLAGS.type == "force" else "<DRiDRi>"
    content = []
    i = 0
    with open(FLAGS.pdb, 'r') as f:
        for line in f.readlines():
            if line.startswith("ATOM") or line.startswith("HETATM"):
                content.append(line)
                if FLAGS.unit == "residue" and line[12:16].strip() != 'CA': continue
                anisou_domain = f"{round(c*df.loc[i, mtxheader + 'xx']):7d}{round(c*df.loc[i, mtxheader + 'yy']):7d}{round(c*df.loc[i, mtxheader + 'zz']):7d}{round(c*df.loc[i, mtxheader + 'xy']):7d}{round(c*df.loc[i, mtxheader + 'xz']):7d}{round(c*df.loc[i, mtxheader + 'yz']):7d}"
                i += 1
                content.append("ANISOU" + line[6:28] + anisou_domain + line[70:])
            else:
                content.append(line)
    pdbdir, pdbfnm = os.path.split()
    with open(os.path.join(pdbdir, "anisou_" + pdbfnm), 'w') as fw:
        fw.writelines(content)


if __name__ == '__main__':
    main()        
