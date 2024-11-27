import os
import re 
import gzip
from Bio import SeqIO


def chunkit(genomepath : str, outdir : str) -> list:

    scaffolds = {}
    if genomepath.endswith(".gz"):
        fhd = gzip.open(genomepath, mode='rt')
    else:
        fhd = open(genomepath)

    ctgnames = []
    for rec in SeqIO.parse(fhd, "fasta"):
        scaffolds[rec.id] = rec
        if re.match("chr[1-7]H_", rec.id):
            ctgnames.append(rec.id)

    os.chdir(outdir)
    
    for ctg in ctgnames:
        seq = scaffolds[ctg]
        with open("%s.fa" %(ctg), 'w') as out:
            p = SeqIO.write(seq, out, "fasta")

    unanchored = []
    for ctgid in scaffolds:
        if ctgid not in ctgnames:
            unanchored.append(scaffolds[ctgid])

    if len(unanchored):
        with open("chrUn.fa", 'w') as out:
            p = SeqIO.write(unanchored, out, "fasta")
        ctgnames.append("chrUn")

    return ctgnames
