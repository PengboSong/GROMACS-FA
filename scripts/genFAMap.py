#!/usr/bin/env python3
# coding = utf-8

import argparse
import struct


def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-s',
        type=int,
        default=0,
        help="Start integer of the required map (included)."
    )
    parser.add_argument(
        '-e',
        type=int,
        help="End integer of the required map (excluded)."
    )
    parser.add_argument(
        '--excl',
        type=list,
        default=[],
        help="List of integers to be excluded from "
    )
    parser.add_argument(
        '-o',
        type=str,
        default="atoms.map",
        help="Filename of atom map file generated."
    )

    FLAGS, unparsed = parser.parse_known_args()

    return FLAGS


def main():
    FLAGS = parseArgs()
    s = max(0, FLAGS.s)
    e = max(1 + s, FLAGS.e)
    exclLst = [int(x) for x in FLAGS.excl]
    # Create a list of unsigned integers
    unsignedInts = [x for x in range(s, e) if x not in exclLst]
    # 
    unsignedInts.insert(0, len(unsignedInts))
    # Create a binary data stream using struct.pack()
    buff = struct.pack('I' * len(unsignedInts), *unsignedInts)
    # Write the binary data to file
    with open(FLAGS.o, 'wb') as f:
        f.write(buff)


if __name__ == '__main__':
    main()
