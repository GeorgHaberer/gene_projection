import re 
import os 
import csv 
from Bio import SeqIO

os.chdir("/mnt/data/bulbosum")

clusters = {}
repcodes = set()
for rec in SeqIO.parse("hbulbosum.evidence.high.cdhit", "fasta"):
    repcodes.add(rec.id)

aktcid = 0
with open("hbulbosum.evidence.high.cdhit.clstr") as fhd:
    for line in fhd:
        if line.startswith(">"):
            aktcid = int(re.split(r"\s", line.strip())[1])
            clusters[aktcid] = []
        else:
            tmp = re.split(r"\s+", line.strip())
            tid = re.sub(r">", "", tmp[2])
            tid = tid[:-3]
            clusters[aktcid].append(tid)

os.chdir("datasets")
with open("hbulb.geneweights.txt", 'w') as out:
    for cid, cluster in clusters.items():
        sz = len(cluster)
        for tid in cluster:
            if tid in repcodes:
                out.write("%s\t%i\n" %(tid, sz))
