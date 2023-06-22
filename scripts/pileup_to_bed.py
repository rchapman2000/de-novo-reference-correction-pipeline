import sys
import os
import argparse as ap

def main():

    parser = ap.ArgumentParser()

    parser.add_argument('-i', '--input', required = True, type=str, \
        help='[Required] - path to pileup file', \
        action = 'store', dest = 'infile')
    parser.add_argument('-o', '--output', required = True, type=str, \
        help='[Required] - Output bed file', action='store', dest='outfile')
    parser.add_argument('--minCov', required=False, type=int, \
        help='[Required] The minimum coverage required for a position to be output.', \
        action='store', dest='minCov')

    args = parser.parse_args()

    minCov = 0
    if (args.minCov):
        minCov = args.minCov

    f = open(args.infile, "r")
    o = open(args.outfile, "w+")

    for l in f:
        split = l.strip().split("\t")
        if int(split[3]) >= minCov:
            chrom = split[0]
            pos = int(split[1])
            o.write("{0}\t{1}\t{2}\n".format(chrom, pos - 1, pos))

    f.close()
    o.close()

if __name__ == "__main__":
    main()