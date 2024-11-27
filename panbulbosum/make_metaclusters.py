import re 
import os 
import csv 
from Bio import SeqIO

os.chdir("/mnt/data/bulbosum")

clusters = {}

aktcid = 0
with open("hbulbosum.evidence.high.S10c95.cdhit.clstr") as fhd:
    for line in fhd:
        if line.startswith(">"):
            aktcid = int(re.split(r"\s", line.strip())[1])
            aktcid += 1
            clusters[aktcid] = []
        else:
            tmp = re.split(r"\s+", line.strip())
            tid = re.sub(r">", "", tmp[2])
            tid = tid[:-3]
            clusters[aktcid].append(tid)

os.chdir("datasets")
with open("hbulb.metacluster.txt", 'w') as out:
    for cid, cluster in clusters.items():
        for tid in cluster:
            out.write("%s\t%s\n" %(cid, tid))

os.chdir("/mnt/data/bulbosum/orthotmp")
with open("N0.tsv") as fhd:
    reader = csv.reader(fhd, delimiter='\t')
    header = next(reader)
    with open("hbulb.metacluster.N0s.txt", 'w') as out:
        mid = 1
        for _, _, _, *infos in reader:
            metagroup = []
            for info in infos:
                tids = re.split(", ", info)
                metagroup += tids
            for tid in metagroup:
                if tid.startswith("H.BULB"):
                    out.write("%s\t%s\n" %(mid, tid))
            mid += 1
        
