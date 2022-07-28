#!/usr/bin/env python

# Author: Justine Macalindong justinem@uoregon.edu


'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''
__version__ = "0.3"         
                            

DNA_bases = "ATGCN"
RNA_bases = "AUGCN"

def convert_phred(letter: str) -> int:
    """Converts a single character into a phred score"""
    return ord(letter) - 33

def qual_score(phred_score: str) -> float:
    """Calculates the average quality score of an entire string"""
    sum = 0 
    length = 0 
    for letter in phred_score:
        length += 1
        sum += convert_phred(letter)
    return(sum/length)

def validate_base_seq(seq,RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(seq)<=(RNAbases if RNAflag else DNAbases)

def gc_content():
    '''This function calculates the GC content of a DNA sequence.'''
    DNA = DNA.upper()
    return (DNA.count("G") + DNA.count("C")) / len(DNA)

def oneline_fasta(fasta, new_file):
    '''This function takes all sequence lines of a FASTA file and puts them all on one line.'''
    sequence = ''
    output = open(new_file, "w")
    with open(fasta, "r") as file:
        for line in file:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    if sequence != "":
                        # output.write(str(line) + '\n')
                        output.write(str(sequence) + '\n')
                    output.write(str(line) + '\n')
                else:
                    sequence+=line
        if sequence:
            output.write(str(sequence) + '\n')

if __name__ == "__main__":
    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"
    print("Passed DNA and RNA tests")

    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")

    assert qual_score(phred_score) == 37.62105263157895, "wrong average phred score"
    print("You calcluated the correct average phred score")

    assert gc_content("GCGCGC") == 1
    assert gc_content("AATTATA") == 0
    assert gc_content("GCATGCAT") == 0.5
    print("correctly calculated GC content")