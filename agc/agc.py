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
from collections import Counter
import numpy


# import nwalign3 as nw

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
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default=400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default=10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default=100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default=8,
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
        for sequence in sequences:
            if (
                    len(sequence) == 0
                    or not sequence[0] in "TGAC"
            ):
                sequence = "".join(sub_sequence).replace("\n", "")
                sub_sequence = []
                if len(sequence) >= minseqlen:
                    yield sequence
            else:
                sub_sequence.append(sequence)


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """ Generate unique sequences that appear at least mincount and return them (descending order)
    :param amplicon_file: fasta file
    :param minseqlen: minimal length of a sequence
    :param mincount: minimal number of occurancies
    :return: [sequence, number of occurancies]
    """
    sequences = {}
    for sequence in read_fasta(amplicon_file, minseqlen):
        if not sequence in sequences:
            sequences[sequence] = 1
        else:
            sequences[sequence] += 1
    sequences = {k: v for k, v in sorted(sequences.items(),
                                         key=lambda item: item[1],
                                         reverse=True)
                 if v >= mincount}
    for sequence in sequences:
        yield [sequence, sequences[sequence]]


def get_chunks(sequence, chunk_size):
    """
    Give a list of sub-sequences (at least 4 segments / sequence)
    :param sequence: sequence
    :param chunk_size: length of a segment
    :return: list of sub-sequences (len = chunk_size) non-overlapping
    """
    # on split la sequence en au moins 4 segments (ou plus) de taille chunk-size
    print("sequence: " + str(sequence))
    seq = list(sequence)  # "AAAAA" -> A,A,A,A,A
    start = 0
    sub_sequence = []
    while (start + chunk_size < len(sequence)) and (len(sub_sequence) < 4):
        sub_seq = []
        for i in range(start, start + chunk_size):
            sub_seq.append(seq[i])
        sub_sequence.append("".join(sub_seq))  # A,A,A -> "AAA"
        start = start + chunk_size  # on update start
    return sub_sequence


def get_unique(ids):
    """ Function predefined by the teacher
    :param ids:
    :return:
    """
    return {}.fromkeys(ids).keys()


def common(lst1, lst2):
    """
    Function predefined by the teacher
    :param lst1:
    :param lst2:
    :return:
    """
    return list(set(lst1) & set(lst2))


def cut_kmer(sequence, kmer_size):
    """ From a sequence, generate kmers (length k)
    :param sequence: sequence
    :param kmer_size: length of a kmer
    :return: kmers
    """
    seq = list(sequence)
    for i in range(len(sequence) - kmer_size + 1):
        kmer = []
        for j in range(i, i + kmer_size):
            kmer.append(seq[j])
        yield "".join(kmer)


def get_identity(alignment_list):
    """ Give the percentage of identity between 2 sequences
    :param alignment_list: alignment (sequence1, sequence2)
    :return: percentage of identity between the two sequences
    """
    alignment = alignment_list
    sequence1 = list(alignment[0])
    sequence2 = list(alignment[1])
    identity = 0
    for i in range(len(sequence1)):
        if (
                sequence1[i] == sequence2[i]
                and sequence1[i] != "-"
        ):
            identity = identity + 1 / len(sequence1)
    return identity * 100


def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    """ Select only unique kmers in a dictionary of kmers
    :param kmer_dict: dictionary of kmers
    :param sequence: sequence
    :param id_seq: integer
    :param kmer_size: integer
    :return: dictionary of kmers
    """
    for kmer in cut_kmer(sequence, kmer_size):
        if kmer in kmer_dict:
            if not id_seq in kmer_dict[kmer]:
                kmer_dict[kmer].append(id_seq)
        else:
            kmer_dict[kmer] = [id_seq]
    return kmer_dict


def search_mates(kmer_dict, sequence, kmer_size):
    """ Give 8 sequences with the most similarities with our sequence
    :param kmer_dict: dictionary of kmers
    :param sequence: sequence
    :param kmer_size: length of kmers
    :return: list of 8 sequences the most similar to sequence
    """

    most_common = Counter(kmer_dict).most_common(kmer_size)
    counter = 0
    result = []
    while True:
        if counter >= kmer_size or counter >= len(most_common) - 1:
            return result
        indexs = most_common[counter][1]
        for index in indexs:
            if not index in result:
                result.append(index)
        counter += 1


def detect_chimera(perc_identity_matrix):
    """ Detect if it is a chimera (1) or not (0)
    :param perc_identity_matrix: matrix
    :return: boolean 1 : it is a chimera, 0 it is not
    """
    perc_identity_matrix = numpy.array(perc_identity_matrix)
    shape = numpy.shape(perc_identity_matrix)
    print(shape)
    moy = 0
    print("perc_identity_matrix: " + str(perc_identity_matrix))
    for i in range(shape[1]):
        ecart_type = numpy.std(perc_identity_matrix[::, i])
        print("ecart_type: " + str(ecart_type))
        moy += ecart_type
    moy = moy / shape[1]
    return moy > 5


def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """
    Return non-chimera sequences with their occurancy
    :param amplicon_file: fasta file
    :param minseqlen: (int) minimal length of a sequence
    :param mincount: (int) minimal number of occurancies
    :param chunk_size: (int) size of each chunk / segment
    :param kmer_size: (int) length of each kmer
    :return: generator of non-chimera sequences with their occurancy [sequence, count]
    """

    sub_sequences = []  # garde en mémoire les sequences et les sub-sequence

    # Step 1
    for sequence, occurancy in dereplication_fulllength(amplicon_file, minseqlen, mincount):
        segments = get_chunks(sequence, chunk_size)
        sub_sequences.append(segments)

        if len(sub_sequences) < 2 :
            continue
        id_seq = 0
        for seq in sub_sequences:
            sequence = seq[0]
            index_sequence = seq[1]
            kmer_dict = {}
            for sub_sequence_index, sub_sequence in enumerate(sequence):
                id_seq += 1
                kmer_dict = get_unique_kmer(kmer_dict, sub_sequence, id_seq, kmer_size)
                mates = []
                for target_seq in range(index_sequence):
                    mates.append(search_mates(kmer_dict,
                                              sequence[target_seq][sub_sequence_index],
                                              kmer_size))

    return print("Non abouti, "
                 "je n'ai pas réussi à comprendre les étapes malgré vos explications")


def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """ Identify OTU
    :param amplicon_file: fasta file
    :param minseqlen: (int) minimal length of a sequence
    :param mincount: (int) minimal number of occurancies
    :param chunk_size: (int) size of each chunk / segment
    :param kmer_size: (int) length of each kmer
    :return: List of OTU with occurancy
    """
    OTU_list = []
    non_chimeral_list = chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size)

    for i, (sequence1, occurancy1) in enumerate(non_chimeral_list):

        for j in range(i + 1, len(non_chimeral_list)):
            sequence2, occurancy2 = non_chimeral_list[j]
            if (sequence1 != sequence2) and \
                    (get_identity([sequence1, sequence2]) > 97) and \
                    (occurancy1 < occurancy2):
                continue
            OTU_list.append(sequence1)

    return OTU_list


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i + width] for i in range(0, len(text), width))


def write_otu(otu_list, output_file):
    """ Write the results on a file
    :param otu_list: list of OTU
    :param output_file: name of output file
    :return: Generate a file named output_file
    """
    with open(output_file, "w") as myfile:
        for counter, otu_element in enumerate(otu_list):
            sequence_otu, occurence_otu = otu_element
            myfile.write(">OTU_{0} occurence:{1}\n".format(counter + 1, occurence_otu))
            myfile.write(fill(sequence_otu) + "\n")
        myfile.close()


def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()


if __name__ == "__main__":
    main()
