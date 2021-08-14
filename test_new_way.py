from typing import Dict

import numpy as np
from Bio import Seq
from numba import jit, njit, vectorize
from main import *

to_int = {
    'A': 0,
    'T': 1,
    'C': 2,
    'G': 3,
    'N': 4
}


def config_sequence(sequences: Dict[str, Seq]) -> Dict:
    result = {}
    for key in sequences:
        sequence = sequences[key]
        sequence_lst = [to_int[x] for x in str(sequence)]
        result[key] = np.array(sequence_lst)
    return result


def hamming_distance_between_array(s1: np.ndarray, s2: np.ndarray, length: int) -> int:
    return length - sum(s1 == s2)


def lcs_hamming(s1: np.ndarray, s2: np.ndarray, k: int, length: int):
    count = 0
    s1_set = np.array([s1[i: i + length] for i in range(len(s1) - length + 1)])
    s2_set = np.array([s2[i: i + length] for i in range(len(s2) - length + 1)])
    for i in range(len(s1_set)):
        for j in range(len(s2_set)):
            ham_distance = hamming_distance_between_array(s1_set[i], s2_set[j], length)
            if ham_distance < k:
                count += 1
    return count


def lcs_between_query_and_targets(query: np.ndarray, targets: Dict[str, np.ndarray], k, length, query_name: str):
    total = 0
    for key in targets.keys():
        total += lcs_hamming(query, targets[key], k, length)
    return total


def reverse_targets_search(query: np.ndarray, targets: Dict[str, Seq], k, length, query_name):
    reversed_targets = {}
    for target_name in targets.keys():
        reversed_targets[target_name] = targets[target_name].reverse_complement()
    reversed_int_targets = config_sequence(reversed_targets)
    return lcs_between_query_and_targets(query, reversed_int_targets, k, length, query_name)


if __name__ == '__main__':
    query = {}  # long non-coding rn
    targets = {}
    reversed_query = {}
    targets_to_query = {}
    query_name = []
    query_file = ''
    query_sequence_file = './ncRNA_genes_fb_2021_02.json'
    chrom_file = './dmel-all-chromosome-r6.40.fasta'
    intron_file = './dmel-all-intron-r6.40.fasta'
    relation_file = ''
    query_name_file = './lncRNAs.txt'
    read_query_name(query_name_file, query_name)
    read_query_from_json(query_sequence_file, query, query_name)
    read_target(chrom_file, targets)
    read_relation_file(relation_file, targets_to_query)

    int_query = config_sequence(query)


    # a = reverse_targets_search(int_query['CR32218'], targets, 1, 15, 'CR32218')