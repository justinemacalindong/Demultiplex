#!/usr/bin/env python
import argparse
import numpy as np

#2: Setting global variables
def get_args():
    parser = argparse.ArgumentParser(description="Determine the quality score of each nucleotide.")
    parser.add_argument("-f", "--filename", help="Specify the filename")
    parser.add_argument("-r", "--readlength", help="Specify the length of the read")

args = get_args()

num_lines = 0 
with open(args.filename, "r") as fh:
    for lines in fh:
        num_lines += 1

all_qscores = np.zeros((args.readlength,int(num_lines//4)))
mean = np.zeros((args.readlength), dtype=float)


with open(args.filename,"r") as fh:
    line_count = 0 
    for line in fh:
        line = line.strip('\n')
        line_count += 1
        if line_count % 4 == 0:
            counter_array = line_count//4
            for count, letter in enumerate(line):
                all_qscores[count,counter_array-1] = (bioinfo.convert_phred(letter))

for idx, value in enumerate(all_qscores):
    mean[idx] = np.mean(value)

print(mean)