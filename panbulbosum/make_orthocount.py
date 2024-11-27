import re 
import os 
import csv 


os.chdir("/mnt/data/bulbosum")
os.chdir("orthotmp")
# ---------------------------------------------------------------------------------
genotypes = ["A17", "A40", "A42", "GRA2256_1", "FB19_011_3", ]
haps = [("hap1", "hap2", "hap3", "hap4", "unanchor"), ("hap1", "hap2", "hap3", "hap4", "unanchor"),
        ("hap1", "hap2", "hap3", "hap4", "unanchor"), ("hap1", "hap2", "hap3", "hap4", "unanchor"),
        ("hap1", "hap2", "unanchor")]

colidx = {}
valids = set()
for gt, haplotypes in zip(genotypes, haps):
    for h in haplotypes:
        key = "%s_%s" %(gt, h)
        colidx[key] = 0
        if h != "unanchor":
            valids.add(key)



out = open("/mnt/data/bulbosum/datasets/hbulb.orthocount.txt", 'w')
with open("N0.tsv") as fhd:
    reader = csv.reader(fhd, delimiter='\t')
    header = next(reader)
    for i, tag in enumerate(header):
        if tag in colidx:
            colidx[tag] = i 
    print(header[3:])
    print(valids)

    for _, _, _, *cols in reader:
        tids = []
        cnt = 0
        for c, colname in zip(cols, header[3:]):
            ctids = re.split(r", ", c)
            if len(ctids[0]):
                if colname in valids:
                    cnt += 1
                for t in ctids:
                    tids.append(t)
        cnt = max(1, cnt)
        for tid in tids:
            out.write("%s\t%i\n" %(tid, cnt))
out.close()
