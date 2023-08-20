import sys
import os
import re
import argparse as ap


def main():
    # Creates an argument parser to define/handle user supplied arguments.
    parser = ap.ArgumentParser()

    parser.add_argument('-i', '--input', required = True, type=str, \
        help='[Required] - Vcf file', \
        action = 'store', dest = 'infile')
    parser.add_argument('-o', '--output', required = True, type=str, \
        help='[Required] - Output vcf file name', action='store', dest='outfile')
    
    args = parser.parse_args()

    # Handles the input file parameter
    infile = ""
    # Checks whether the file supplied by the user exists.
    if not os.path.exists(args.infile):
        # If not, notify the user and exit.
        sys.exit("ERROR: Input file {0} does not exist. Please supply an existing file.".format(args.infile))
    else:
        # If the file exists, create a file stream that opens it.
        infile = open(args.infile, "r")


    # Create an output file stream that writes the output file supplied by the user.
    # Use the "w+" directive to create the file if it does not exist.
    outfile = open(args.outfile, "w+")

    # Loop over each line in the input VCF file.
    for line in infile:
        # The purpose of this script is to remove any genotyping information from the VCF
        # file. Thus, we need to both remove those columns from each entry line, and modify the header to remove
        # any of the FORMAT definition lines.

        # To identify whether a line is part of the header, use regex
        # to check that the line begins with two "##" characters,
        if re.match("##", line):
            # If so, check whether it is a FORMAT definition line by using
            # another regex.
            if not re.match("##FORMAT", line):
                # If the line is not a FORMAT definition line, we can
                # write it in the output file.
                outfile.write(line)
        else:
            # If the line is not a header line, then it will be
            # tab-delimited. Remove the newline character and
            # split the line at the tab characters.
            split = line.strip("\n").split("\t")

            # Remove the fields related to the genotype information,
            # FORMAT and sample, which are columns 9 and 10
            genotypeRemoved = split[:8]

            # Write the modified line to the output file.
            outfile.write("\t".join(genotypeRemoved) + "\n")
    
    # Close the input and output files.
    infile.close()
    outfile.close()


    


if __name__ == "__main__":
    main()