##!/usr/bin/env python3

import numpy as np
import pandas as pd

def arguments():
    import argparse, sys
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", help="Input vcf file, default: stdin.",
    	nargs='?', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument("-o", "--output_file", help="Output vcf file, default: stdout",
    	nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    ArgP = parser.parse_args()
    return ArgP

def main():
    # get input arguments
    args = arguments()

    # read in the input file
    df = pd.read_csv(args.input_file, 
                     sep="\t", 
                     header=None, 
                     names=["scaffold", "start", "stop"])
    
    # create new columns:
    # source, this script
    df["source"] = "bed2gff3"

    # feature type, set to peak
    df["feature"] = "peak"

    # score, strand and frame. none of there are known so we set them as "."
    df["score"] = "."
    df["strand"] = "."
    df["frame"] = "."

    # attribute, the ID of the peak, combine scaffold, start and stop
    df["attribute"] = "ID=" + df["scaffold"] + ":" + df["start"].astype(str) + "-" + df["stop"].astype(str)

    # order the columns
    df = df[["scaffold", "source", "feature", "start", "stop", "score", "strand", "frame", "attribute"]]

    # convert datafram to a tab separated string
    df_str = df.to_csv(index=False, header=False, sep='\t')

    # write the output file
    with args.output_file as outfile:
        outfile.write("##gff-version 3\n")
        outfile.write(df_str)

if __name__ == '__main__':
    main()