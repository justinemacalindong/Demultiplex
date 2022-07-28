#!/usr/bin/env python

import argparse
import numpy as np
import bioinfo
import gzip
import matplotlib.pyplot as plt

#2: Setting global variables
def get_args():
    parser = argparse.ArgumentParser(description="Determine the quality score of each nucleotide.")
    parser.add_argument("-f", "--filename", help="Specify the filename")
    parser.add_argument("-r", "--readlength", help="Specify the length of the read", type=int)
    parser.add_argument("-pn", "--plotname", help="Specify the name of your plot file")
    parser.add_argument("-pt", "--plottitle", help="Specify the file name in the title of your plot")
    return parser.parse_args()

args = get_args()

the_list = [0]*args.readlength
        
def populate_list(file: str) -> tuple[list, int]:
    with gzip.open(file,"r") as fh:
        line_count = 0 
        for line in fh:
            line = line.decode("ascii")
            line = line.strip('\n')
            line_count += 1
            if line_count % 4 == 0:
                for count, letter in enumerate(line):
                    the_list[count] += bioinfo.convert_phred(letter) 
    return (the_list, line_count)

my_list, num_lines = populate_list(args.filename)

for index, total in enumerate(my_list):
    my_list[index] = total/(num_lines/4)
    print(index,'    ',my_list[index])

plt.bar(range(0,args.readlength), my_list)
plt.title("Mean Quality Score of vs # Base Pair " + args.plottitle)
plt.ylabel("Mean Quality Score")
plt.xlabel("# Base Pair")
plt.savefig(args.plotname)