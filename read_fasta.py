from Bio import SeqIO

if __name__ == '__main__':
    dna = {}

    for r in SeqIO.parse("./dmel-all-chromosome-r6.40.fasta", "fasta"):
        if r.description.__contains__("loc=Y"):
            print(r.description)
            dna[r.name] = r.seq
