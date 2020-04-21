#!/bin/python

# Author: Wyatt Eng
# Date 4/7/20
# Goal: Identify sequences that have been index swapped based on the output of an initial demultiplexing step.
#       Track the sequence combinations that result in an index swap and output the resulting data to a file.
#       Input Files: 2 Undetermined files from R1 and R2, 1 known index file
#       Output Files: 2 Identified index swapped files (R1 and R2), 2 Identified error-corrected index swapped files (R1 and R2), \
#       and 2 Unknown files (R1 and R2)

import argparse
import re
import gzip
from Bio import SeqIO
from mimetypes import guess_type
from functools import partial

# Set input variables
def get_args():
    parser = argparse.ArgumentParser(
        description="Identify the records that contain an index swap or are undidentifiable. \
                     Output files: R1_IndexSwapped.fastq, R2_IndexSwapped.fastq, Unknown_R1.fastq, Unknown_R2.fastq")
    parser.add_argument("-U1", "--undetermined_R1", help="path to a fastq file containing the undetermined reads from R1", type=str)
    parser.add_argument("-U2", "--undetermined_R2", help="path to a fastq file containing the undetermined reads from R2", type=str)
    parser.add_argument("-I", "--known_indexes", help="path to a file containing a list of the known indexes", type=str)
    return parser.parse_args()

args = get_args()
# Setting user inputted variables
U1_File = args.undetermined_R1
U2_File = args.undetermined_R2
Index_File = args.known_indexes

##### !!!!!DELETE LATER!!!!! #####
U1_File = "Undetermined_R1_First100.fastq.gz"
U2_File = "Undetermined_R2_First100.fastq.gz"
Index_File = "index_seq.txt"

#Define global variables
input_files = [U1_File, U2_File]
unknown_files = ["Unknown_R1.fastq", "Unknown_R2.fastq"]
Known_indexes = set() # Create a set to store all of the known indexes
RC_Known_indexes=set() # A set to store the reverse complement of all of the known indexes

# handle both regular text and gzipped files / From: https://stackoverflow.com/questions/42757283/seqio-parse-on-a-fasta-gz
encoding = guess_type(U1_File)[1]  # uses file extension
_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

def hamming_distance(s1, s2):
    '''Return the Hamming distance between equal-length sequences /
    From:https://biology.stackexchange.com/questions/23523/hamming-distance-between-two-dna-strings'''
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

#Used in reverse complement function
def reverse(string):
    '''Reverse a string. Used in reverse_complement'''
    str = ""
    for char in string:
        str = char + str
    return str

# Dictionary used in reverse_complement function, defined dictionary outside
# of function to maximize efficiency.
complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

def reverse_complement(sequence):
    '''Finds the reverse complement of a given string'''
    temp_string = str()
    for character in sequence:
        temp_string += complement_dict[character]
    return reverse(temp_string) #Returning the reverse of the complement

def load_known_indexes(file):
    '''Load the known indexes into a set and the reverse complement of the known indexes to a set which can later be compared to.\
    Returns true if the current index is one nucleotied away from a known index; otherwise false'''
    with open(file) as ind:
        for line in ind:
            line = line.strip()
            Known_indexes.add(line)
            rc_index=reverse_complement(line)
            RC_Known_indexes.add(rc_index)

def One_Away(current_index, int):
    '''Checks to see if its 1 away from the known indexs, int 1 = forward indexes, 2 = reverse complement index'''
    counter_list = [] # Counter list stores the number of nucleotide mismatches per known index
    one_away_bool = False # bool variable will be returned to signify if the current index is one away from a known index

    # Set which index dictionary to use, either known indexes or the reverse complement of each known index
    if int == 1:
        index_dict = Known_indexes
    elif int == 2:
        index_dict = RC_Known_indexes

    # loop through every nucleotide of each index and count the number of mismatches between the current index and known indexes
    for index in index_dict:
        counter = 0  # Counter holds the number of mismatches bwetween current index and the known index
        for i in range(len(current_index)):
            if current_index[i] != index[i]:
                counter = counter + 1
        counter_list.append(counter) # storing the number of mismatches in a list to be referenced later

    # Checking to see if current index was one nucleotide away from a known index
    if 1 in counter_list or 0 in counter_list:
        one_away_bool = True

    # return true if the current index is one nucleotied away from a known index; otherwise false
    return one_away_bool


def main():
    # Load indexes (including reverse complemented indexes)
    load_known_indexes(Index_File)

    # file_number stores the file number in order to name the correct output files
    file_number = 1
    for file in input_files:
        #Name output files according to R1 or R2
        IndexSwap_OutFile = "R" + str(file_number) + "_IndexSwapped.fastq"
        Unknown_OutFile = "Unknown_R" + str(file_number) + ".fastq"
        EC_IH_OutFile = "R" + str(file_number) + "_ErrorCorrected_IndexSwapped.fastq"

        #Open write files for R1 or R2
        Index_output_handle = open(IndexSwap_OutFile, 'w')
        Unknown_output_handle = open(Unknown_OutFile, 'w')
        EC_IH_output_handle = open(EC_IH_OutFile, 'w')

        # Sort the Undetermined file into index swapped or unknown records
        with _open(file) as handle:
            for record in SeqIO.parse(handle, "fastq"):
                indexes = re.search('(\w{8})[+](\w{8}$)', record.description) # Extracting the 2 indexes from the header line of record
                index1 = indexes.group(1)
                index2 = indexes.group(2)

                index1_one_away = One_Away(index1, 1) # Stores true if index1 is one away from a known index
                index2_one_away = One_Away(index2, 2) # Stores true if index2 is one away from a known index

                if index1 in Known_indexes and index2 in RC_Known_indexes: # checking for index hopped if index 1 and index 2 have no errors
                    Index_output_handle.write(record.format("fastq"))
                elif index1_one_away and index2_one_away:   # Check to see if both records are 1 away from a known index
                    EC_IH_output_handle.write(record.format("fastq"))
                else:
                    Unknown_output_handle.write(record.format("fastq"))
        handle.close()

        #Increment the Read file number to 2
        file_number = file_number + 1

        #Close write files
        Index_output_handle.close()
        Unknown_output_handle.close()
        EC_IH_output_handle.close()

main()




















