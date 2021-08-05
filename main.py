import json
from typing import Dict, List
from Bio import pairwise2, SeqIO
from Bio.Seq import Seq
import random

# K: number of mismatches allowed in substring
# LENGTH: min length of substring
K_AND_LENGTH = [
    (0, 5),
    (1, 5),
    (0, 6),
    (1, 6),
    (0, 7),
    (1, 7),
    (0, 8),
    (1, 8),
    (2, 8)
]

# bellow are functions, process starts after the line (if __name__ == '__main__')


# under construction (ignore dynamic_lcs for now)
def random_sequence_generator(length: int) -> str:
    s = ''
    for i in range(length):
        r = random.randint(0, 3)
        if r == 0:
            s += 'A'
        elif r == 1:
            s += 'T'
        elif r == 2:
            s += 'C'
        else:
            s += 'G'
    return s


def random_nucleotide_generator() -> str:
    s = ''
    r = random.randint(0, 3)
    if r == 0:
        s = 'A'
    elif r == 1:
        s = 'T'
    elif r == 2:
        s = 'C'
    else:
        s = 'G'
    return s


def randomly_mutate(s: str, k: int):
    new_s = s
    mutated = []
    for i in range(k):
        r = random.randint(0, len(s) - 1)
        if not r in mutated:
            new_s = new_s[: r] + random_nucleotide_generator() + new_s[r + 1:]
            mutated.append(r)
    return new_s


def randomly_mutate_lst(lst: List[str], k):
    new_lst = []
    for item in lst:
        new_lst.append(randomly_mutate(item, k))
    return new_lst


def create_test_targets(lst: List[str], length: int) -> str:
    s = random_sequence_generator(length - len(lst) * len(lst[0]))
    random.shuffle(lst)
    start_idx = len(lst[0])
    change = int(len(s) / len(lst))
    for i in range(len(lst)):
        r = random.randint(start_idx, start_idx + change)
        s = s[0: r] + lst[i] + s[r:]
        start_idx += r
    return s


# helper function for lcs_hamming
def hamming_distance(s1: str, s2: str):
    if len(s1) != len(s2):
        return
    distance = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            distance += 1
    return distance


def lcs_hamming(s1: str, s2: str, k: int, length: int):
    '''
    This lcs use hamming distance as string comparison method.
    '''
    count = 0
    for i in range(0, len(s1) - length + 1):
        for j in range(0, len(s2) - length + 1):
            sub1 = s1[i: i + length]
            sub2 = s2[j: j + length]
            result = hamming_distance(sub1, sub2)
            if result <= k:
                print(str(i) + " " + str(j) + " " + str(length - result) + " " +
                      sub1 + " " + sub2, end="\t")
                count += 1
            else:
                print(str(i) + " " + str(j) + " N N N", end="\t")
        print("")
    print("total matches: " + str(count))
    return


def lcs_hamming_only_matches(s1: str, s2: str, k: int, length: int):
    """
    This lcs use hamming distance as string comparison method.
    """
    count = 0
    for i in range(0, len(s1) - length + 1):
        for j in range(0, len(s2) - length + 1):
            sub1 = s1[i: i + length]
            sub2 = s2[j: j + length]
            result = hamming_distance(sub1, sub2)
            if result <= k:
                # print(str(i) + " " + str(j) + " " + str(result) + " " + sub1 + " " + sub2, end="\n")
                count += 1
    # print("total matches: " + str(count))
    return count


def lcs_hamming_only_matches_0_1(s1: str, s2: str, length: int):
    """
    This lcs use hamming distance as string comparison method.
    """
    count_0 = 0
    count_1 = 0
    for i in range(0, len(s1) - length + 1):
        for j in range(0, len(s2) - length + 1):
            sub1 = s1[i: i + length]
            sub2 = s2[j: j + length]
            result = hamming_distance(sub1, sub2)
            if result <= 0:
                # print(str(i) + " " + str(j) + " " + str(result) + " " + sub1 + " " + sub2, end="\n")
                count_0 += 1
            if result <= 1:
                count_1 += 1
    # print("total matches: " + str(count))
    return [count_0, count_1]

def lcs_pairwise(s1: str, s2: str, k: int, length: int):
    """
    This lcs use pairwise function from biopython as string comparison method.
    """

    for i in range(0, len(s1) - length + 1):
        for j in range(0, len(s2) - length + 1):
            sub1 = s1[i: i + length]
            sub2 = s2[j: j + length]
            result = pairwise2.align.globalxx(sub1, sub2)
            if result[0][2] >= length - k - 1:
                print(str(i) + " " + str(j) + " " + str(length - result[0][2]) + " " +
                      sub1 + " " + sub2, end="\t")
            else:
                print(str(i) + " " + str(j) + " N N N", end="\t")
        print("")


# read from files function bellow. require actual chromosome files to implement them
def read_query_name(file_name: str, query_name: List[str]) -> None:
    file = open(file_name)
    query_name.extend(file.read().split(" "))
    file.close()
    return


def read_query_from_json(file_name: str, query: Dict, query_name: List[str]) -> None:
    file = open(file_name)
    query_data = json.load(file)
    for i in range(len(query_data["data"])):
        ncRNA_name = query_data["data"][i]["gene"]["locusTag"]
        for rna in query_name:
            if rna in ncRNA_name:
                query[rna] = Seq(query_data["data"][i]["sequence"])
    file.close()
    return


def read_target(file_name: str, targets: Dict) -> None:
    for r in SeqIO.parse("dmel-all-intron-r6.40.fasta", "fasta"):
        if r.description.__contains__("loc=Y"):
            targets[r.name] = r.seq


def read_relation_file(relation_file: str, relations: Dict) -> None:
    pass


def targets_query_compare(query: List[str], targets: List[str], k, length):
    for q in query:
        for target in targets:
            lcs_hamming_only_matches(q, target, k, length)


if __name__ == '__main__':


    # 1. initialization
    query = {}  # long non-coding rna
    targets = {}
    reversed_query = {}
    targets_to_query = {}
    query_name = []

    # 2. read data from files
    # this part is empty since files have not been available
    # file names
    query_file = ''
    query_sequence_file = './ncRNA_genes_fb_2021_02.json'
    dna_file = ''
    relation_file = ''
    query_name_file = './lncRNAs.txt'

    # read files
    read_query_name(query_name_file, query_name)
    read_query_from_json(query_sequence_file, query, query_name)
    read_target(dna_file, targets)
    read_relation_file(relation_file, targets_to_query)

    # Here are some data for testing

    # query = {
    #     'rna1': Seq('UUUUUUU'),
    #     'rna2': Seq('AUAUAU'),
    #     'rna3': Seq('AGCGAU')
    # }
    # targets = {
    #     'chromo1': Seq('AACATAATTTTTTTTTT'),
    #     'chromo2': Seq('ATATATATATATATATATATATA'),
    #     'chromo3': Seq('TTTTTTTTAGCGGTTTTTTTTT')
    # }
    # targets_to_query = {
    #     'chromo1': ['rna1'],
    #     'chromo2': ['rna2'],
    #     'chromo3': ['rna3']
    # }



    # 4. use lcs to find common substring
    results = []
    l = [6]
    for length in l:
        for key in query.keys():
            q = query[key]
            result = [0, 0]
            for target in targets.values():
                a = lcs_hamming_only_matches_0_1(str(q), str(target), length)
                result[0] += a[0]
                result[1] += a[1]
            results.append([key, 0, length, result[0]])
            results.append([key, 1, length, result[1]])
            print(results)

    for l in results:
        print(results)
    # # iterate over genes
    # for target in targets.keys():
    #     # iterate over lncRNAs
    #     for i in range(len(query)):
    #         lncRNA = targets_to_query[target][i]
    #         print('query: ' + query[lncRNA])
    #         print('targets: ' + targets[target])
    #         lcs_hamming_only_matches(str(query[lncRNA]), str(targets[target]), K, LENGTH)
    #         print("\n")

    # gene_to_substring and gene_to_lncRNAs are corresponding

    # print('bellow are tests for different distance function\n')
    #
    # print('testing for multiple mismatches using hamming')
    # s1 = "abcdefzzzzzzzzzzhijklmzzzzzzzz"
    # s2 = "xbxdefyyyyyyyyyyhxjxlmyyyyyyyabxdef"
    # print(s1)
    # print(s2)
    # lcs_hamming_only_matches(s1, s2, 2, 6)
    # print("\n")
    #
    # print('testing for no mismatches using pairwise')
    # s3 = 'ATCGATCGAAAAAAAATCGATCG'
    # s4 = 'CCCCCCCCATCGATCGCCCCCCC'
    # print(s3)
    # print(s4)
    # lcs_pairwise(s3, s4, 0, 8)
    #
    # # Here you can define your own sequences for testing
    # s4 = ''
    # s5 = ''
    # lcs_hamming(s4, s5, K, LENGTH)
    # lcs_pairwise(s4, s5, K, LENGTH)

