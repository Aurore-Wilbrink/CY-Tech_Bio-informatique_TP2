#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Aurore WILBRINK"
__copyright__ = "CY Tech"
__credits__ = ["Aurore WILBRINK"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Aurore Wilbrink"
__email__ = "wilbrinkau@eisti.eu"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):
    """ Generate Sequences from fasta file
    :param amplicon_file: file fasta.gz
    :param minseqlen: minimal length of sequences (int)
    :return: Sequence (len(Sequence) >= minseqlen)
    """
    with gzip.open(amplicon_file, 'rt') as myfile:
        sub_sequence = []
        sequences = myfile.readlines()
        for i, sequence in enumerate(sequences):
            sequence = sequence[:-1]
            if len(sequence) == 0: continue  # ligne vide
            elif (
                    (
                        sequence[0] == "T"
                        or sequence[0] == "G"
                        or sequence[0] == "A"
                        or sequence[0] == "C"
                    )
                    and len(sequence) > 0
            ):
                sub_sequence.append(sequence)
            elif (
                    sequence[0] != "T"
                    or sequence[0] != "G"
                    or sequence[0] != "A"
                    or sequence[0] != "C"
            ) or i == len(sequences) - 1:
                sequence = "".join(sub_sequence)
                sub_sequence = []
                if len(sequence) >= minseqlen:
                    yield sequence

def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """ Generate unique sequences that appear at least mincount and return them (descending order)
    :param amplicon_file: fasta file
    :param minseqlen: minimal length of a sequence
    :param mincount: minimal number of occurancies
    :return: [sequence, number of occurancies]
    """
    sequences = []
    for sequence in read_fasta(amplicon_file, minseqlen):
        if len(sequences) == 0: sequences.append([sequence, 1])
        else:
            for i in range(len(sequences)):
                if sequence == sequences[i][0]:
                    sequences[i][1] = sequences[i][1]+1
                    break
                else: count = False
            if count == False :
                sequences.append([sequence, 1])
    sequences.sort(key=lambda x: x[1], reverse=True)

    for sequence in sequences :
        if sequence[1] >= mincount :
            yield [sequence[0], sequence[1]]


def get_chunks(sequence, chunk_size):
    """
    Give a list of sub-sequences (at leat 4 segments / sequence)
    :param sequence: sequence
    :param chunk_size: length of a segment
    :return: list of sub-sequences (len = chunk_size) non-overlapping
    """
    #on split la sequence en au moins 4 segments (ou plus) de taille chunk-size
    seq = list(sequence) # "AAAAA" -> A,A,A,A,A
    start = 0
    sub_sequence = []
    while (start+chunk_size < len(sequence)) and (len(sub_sequence) < 4):
        sub_seq = []
        for i in range(start, start+chunk_size) :
            sub_seq.append(seq[i])
        sub_sequence.append("".join(sub_seq)) # A,A,A -> "AAA"
        start = start+chunk_size # on update start
    return sub_sequence

def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2): 
    return list(set(lst1) & set(lst2))

def cut_kmer(sequence, kmer_size):
    """ From a sequence, generate kmers (length k)
    :param sequence: sequence
    :param kmer_size: length of a kmer
    :return: kmers
    """
    seq = list(sequence)
    for i in range(len(sequence)-kmer_size+1) :
        kmer=[]
        for j in range(i, i+kmer_size) :
            kmer.append(seq[j])
        yield "".join(kmer)

def get_identity(alignment_list):
    """ Give the percentage of identity between 2 sequences
    :param alignment_list: alignment (sequence1, sequence2)
    :return: percentage of identity between the two sequences
    """
    alignment = nw.global_align(alignment_list[0], alignment_list[1], gap_open=-1, gap_extend=-1,
                    matrix=os.path.abspath(os.path.join(os.path.dirname(__file__), "MATCH")))
    sequence1 = list(alignment[0])
    sequence2 = list(alignment[1])
    identity = 0
    for i in range(len(sequence1)):
        if sequence1[i] == sequence2[i]:
            identity = identity + 1/len(sequence1)
    return identity

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """ Fait appel au générateur fourni par dereplication_fulllength et retourne un générateur des séquences non chimérique au format: yield [sequence, count]
    :param amplicon_file:
    :param minseqlen:
    :param mincount:
    :param chunk_size:
    :param kmer_size:
    :return:
    """
    segments = []
    for sequence in dereplication_fulllength(amplicon_file, minseqlen, mincount):
        segments.append(get_chunks(sequence[0], chunk_size))


           # tu prends le segment 1 test, tu prends les kmers, et tu testes identity avec les segments 1 parents (n-1 avant le test)
    for num_sequence in range(2,len(segments)) : # pour toutes les séquences de 3 à la fin
        for i in range(len(segments[0])): # on teste les segments 1, puis 2, puis 3, puis 4
            for kmer in cut_kmer(segments[num_sequence][i], kmer_size) :
                if get_identity(kmer, segments[0]) :
                    True






def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    pass

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    pass
#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()


if __name__ == '__main__':
    #idres = get_identity(("TGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAG", "TGGGGAATA--GCACAATGGGCGCAAGCCTCTAGCAG"))
    #print(idres)
    #print(idres == 86.5)
    #amplicon_file = "/Users/aurorewilbrink/Documents/EISTI3A/agc-tp/tests/test_sequences.fasta.gz"
    #minseqlen = 200
    #mincount = 3
    #chunk_size = 50
    #kmer_size = 8
    #chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size)

