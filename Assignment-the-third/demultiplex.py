#!/usr/bin/env python

import argparse
import bioinfo
import gzip
import numpy as np

def get_args():
    parser = argparse.ArgumentParser(description="Determine the quality score of each nucleotide.")
    parser.add_argument("-r1", "--read1", help="Specify the read1 file")
    parser.add_argument("-r2", "--read2", help="Specify read2 file")
    parser.add_argument("-i1", "--index1", help="Specify the index1 file")
    parser.add_argument("-i2", "--index2", help="Specify the index2 file")
    parser.add_argument("-k", "--known", help="Specify known index file")
    return parser.parse_args()

args = get_args()

read1 = gzip.open(args.read1,"rt") # Reading and opening read1 file
read2 = gzip.open(args.read2,"rt") # Reading and opening read2 file
index1 = gzip.open(args.index1,"rt") # Reading and opening index1 file 
index2 = gzip.open(args.index2,"rt") # Reading and opening index2 file
known_indexes =  open(args.known, "r") # Reading and opening the known indexes file

r1_unknown = open("read1_unknown.fq", "w") # Opening file to write read1 UNKNOWN records to
r2_unknown = open("read2_unknown.fq", "w") # Opening file to write read2 UNKNOWN records to
r1_hopped = open("read1_hopped.fq", "w") # Opening file to write read1 INDEX HOPPED records to
r2_hopped = open("read2_hopped.fq", "w") # Opening file to write read 2 INDEX HOPPED records to
results = open("results.tsv", "w") # Opening file to write counts to

known = set() # Creating a set to store known indexes into 
known_indexes.readline() # Skipping header line
for line in known_indexes: 
        line = line.split("\t") 
        known.add(line[4].strip()) # Adding each index in the file to the set of known indexes

# Opening dictionaries to store indexes and opening statements
r1_open = {}
r2_open = {}

# Setting known index as key and opening a file for both reads with its sequence and storing it in values
for idx in known:
    r1_open[idx] = open("read1_" + idx + ".fq", "w")
    r2_open[idx] = open("read2_" + idx + ".fq", "w")

# Creating a list of the open files to write to them later
r1_files = list(r1_open.values())
r2_files = list(r2_open.values())

# Initializing counters
index_dict = {}
hopped_dict = {}
unknown_count = 0
matched_count = 0
hopped_count = 0 

while True:
    # Reading header lines first
    read1_header = read1.readline().strip() 
    read2_header = read2.readline().strip()
    index1_header = index1.readline().strip()
    index2_header = index2.readline().strip()
    if read1_header.startswith("@"): # Making sure the first line read is the header line
        # Storing the four lines of each record to an array
        r1_array = np.array([read1_header,read1.readline().strip(),read1.readline().strip(),read1.readline().strip()])
        r2_array = np.array([read2_header,read2.readline().strip(),read2.readline().strip(),read2.readline().strip()])
        i1_array = np.array([index1_header,index1.readline().strip(),index1.readline().strip(),index1.readline().strip()])
        i2_array = np.array([index2_header,index2.readline().strip(),index2.readline().strip(),index2.readline().strip()])
        rev_comp = bioinfo.reverse_complement(i2_array[1]) # Storing the reverse complement of the sequence from index2 in a variable

        if i1_array[1] in known: # If the sequence from index1 is in the known set
            i1_count = 0 # Initializing running sum for converted index1 phred scores
            i2_count = 0 # Initializing running sum for converted index2 phred scores

            # Calculating phred score mean for index1 record
            for i1_letter in i1_array[3]:
                i1_count += bioinfo.convert_phred(i1_letter)
            i1_phred_avg = i1_count/8
            if i1_phred_avg >= 30: # If the index1 record phred score mean is greater than or equal to 30

                # Calculating phred score mean for index2 record
                for i2_letter in i2_array[3]:
                    i2_count += bioinfo.convert_phred(i2_letter)
                i2_phred_avg = i2_count/8
                if i2_phred_avg >= 30: # If the index2 record phred score mean is greater than or equal to 30

                    if i1_array[1] == rev_comp: # If the sequence from index1 matches the reverse complement of the sequence from index2

                        dict_count = 0 # Initializing counter for dictionary lines
                        for i in r1_open.keys(): 
                            dict_count += 1
                            if i1_array[1] == i: # If sequence matches the key in dictionary, write that record to file
                                bioinfo.outputting(r1_files[dict_count-1], r2_files[dict_count-1], r1_array, r2_array, i1_array, rev_comp)
                                matched_count += 1
                                if str(i1_array[1] + ':' + i1_array[1]) in index_dict: 
                                   index_dict[str(i1_array[1] + ':' + i1_array[1])] += 1
                                else:
                                   index_dict[str(i1_array[1] + ':' + i1_array[1])] = 1

                    else: # If the sequence from index1 DOES NOT match the reverse complement of the sequence from index2

                        if rev_comp in known: # If the reverse comp is in known, write to hopped file and add to hopped count
                            bioinfo.outputting(r1_hopped, r2_hopped, r1_array, r2_array, i1_array, rev_comp)
                            hopped_count += 1

                            # Counting the number of records for hopped index
                            if str(i1_array[1] + ':' + rev_comp) in hopped_dict: 
                                hopped_dict[str(i1_array[1] + ':' + rev_comp)] += 1
                            else:
                                hopped_dict[str(i1_array[1] + ':' + rev_comp)] = 1
                        
                        else: # If the reverse comp is NOT in known, write to unknown file and add to unknown count
                            bioinfo.outputting(r1_unknown, r2_unknown, r1_array, r2_array, i1_array, rev_comp)
                            unknown_count += 1

                else: # If the index2 record phred score is less than 30, write to unknown files and add to unknown count
                    bioinfo.outputting(r1_unknown, r2_unknown, r1_array, r2_array, i1_array, rev_comp)
                    unknown_count += 1

            elif i1_phred_avg < 30: # If the index1 record phred score is less than 30, write to unknown files and add to unknown count
                bioinfo.outputting(r1_unknown, r2_unknown, r1_array, r2_array, i1_array, rev_comp)
                unknown_count += 1

        else: # If the sequence from index1 is NOT in the known set, write to unkown file and add to unknown count
            bioinfo.outputting(r1_unknown, r2_unknown, r1_array, r2_array, i1_array, rev_comp)
            unknown_count += 1

    else: # Exiting the loop
        break 

# Calculating total number of records
total_count = unknown_count + hopped_count + matched_count
        
# Printing and writing the counts of unknown, hopped, and matched records
results.write("Total number of unknown index pairs:" + '\t' + str(unknown_count) + '\n')
results.write("Total number of hopped index pairs:" + '\t' + str(hopped_count) + '\n')
results.write("Total number of matched index pairs:" + '\t' + str(matched_count) + '\n')
results.write("Total number of records:" + '\t' + str(total_count) + '\n' + '\n')

sorted_idict = sorted(index_dict.items())
sorted_hdict = sorted(hopped_dict.items())

results.write(str("Matched indexes" + '\n'))
results.write(str("Index" + "\t" + "Count" + "\t" + "Percentage" + "\n"))
for item1 in sorted_idict: 
    percent1 = (int(item1[1])/363246735) * 100
    results.write(str(item1[0]) + '\t' + str(item1[1]) + '\t' + str(percent1) + '\n')

results.write(str("\n" + "Hopped indexes" + '\n'))
results.write(str("Index1:Index2" + "\t" + "Count" + "\n"))
for item2 in sorted_hdict: 
    results.write(str(item2[0]) + '\t' + str(item2[1]) + '\n')

print("All done!")

# Closing all files
read1.close()
read2.close()
index1.close()
index2.close()
r1_unknown.close()
r2_unknown.close()
r1_hopped.close()
r2_hopped.close()
results.close()

for file1 in r1_files:
    file1.close()

for file2 in r2_files:
    file2.close()


