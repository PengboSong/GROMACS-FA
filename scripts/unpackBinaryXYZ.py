#!/usr/bin/env python3
# coding = utf-8


import argparse
import os
import struct


def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f', '--file',
        type=str,
        default='.',
        help="Filename of output binary packed coordinate file"
    )

    FLAGS, unparsed = parser.parse_known_args()

    return FLAGS


def main():
    FLAGS = parseArgs()
    assert os.path.isfile(FLAGS.file), "The specified file not exists."

    with open(FLAGS.file, 'rb') as f:
        framen, N = struct.unpack("II", f.read(struct.calcsize("II")))
        print("Got %s frames and %d atoms/residues/molecules in total" % (framen, N))
        blockfmt = "fff" * N
        for framei in range(framen):
            print("Packed XYZ data at frame %d" % framei)
            print("i     X          Y          Z          ")
            data = struct.unpack(blockfmt, f.read(struct.calcsize(blockfmt)))
            for i in range(N):
                print("%5d %10.3f %10.3f %10.3f" % (i, data[3*i], data[3*i+1], data[3*i+2]))
    

if __name__ == '__main__':
    main()
