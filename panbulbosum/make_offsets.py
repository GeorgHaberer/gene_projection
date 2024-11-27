import re 
import os 
import csv
import sys
import copy
import gzip
from Bio import SeqIO

sys.path.append("/home/georg/src/pyprojects/bioutils")
from gffparser import parseTaggedGFF
from genedata import GeneLocus, mRNA, Cds
from biosequences import GenomicTemplate, FastaSequence


BULBROOT = "/mnt/data/bulbosum"
RELEASEHOME = os.path.join(BULBROOT, "release_v1")


def getOffset(t):

    tpositions = t.orderedBasePositions()
    cpositions = t.coding.orderedBasePositions()
    tpositions.sort()
    cpositions.sort()
    cmin = cpositions[0]
    cmax = cpositions[-1]
    if t.strand == 0:
        for offset, pos in enumerate(tpositions):
            if pos == cmin:
                break
    else:
        for offset, pos in enumerate(reversed(tpositions)):
            if pos == cmax:
                break
    return offset


# ---------------------------------------------------------------------------------
genotypes = ["A17", "A40", "A42", "GRA2256_1", "FB19_011_3", ]
haps = [("hap1", "hap2", "hap3", "hap4", "unanchor"), ("hap1", "hap2", "hap3", "hap4", "unanchor"),
        ("hap1", "hap2", "hap3", "hap4", "unanchor"), ("hap1", "hap2", "hap3", "hap4", "unanchor"),
        ("hap1", "hap2", "unanchor")]

for gt, haplotypes in zip(genotypes, haps):
    for h in haplotypes:
        print(gt, h)
        OUTDIR = os.path.join(RELEASEHOME, "%s_%s" %(gt, h))
        gff = os.path.join(RELEASEHOME, "%s_%s" %(gt, h), "%s_%s.pgsb.r1.Sep2024.gff3" %(gt, h))
        genes, genome = parseTaggedGFF(gff)
        print(len(genome))
        scaffolds = {}
        genomefile = os.path.join(RELEASEHOME, "%s_%s" %(gt, h), "%s_%s.fa" %(gt, h))
        for rec in SeqIO.parse(genomefile, "fasta"):
            scaffolds[rec.id] = GenomicTemplate(ID=rec.id, seq=str(rec.seq).upper())      

        total = 0
        errors = 0
        out = open(os.path.join(OUTDIR, "%s_%s.offset_transcripts.fa" %(gt, h)), 'w')
        for gene in genome: 
            for tid, t in gene.transcripts.items():
                total += 1
                cod = t.coding
                if cod.minimum() > t.minimum() and cod.maximum() < t.maximum():
                    # we have both 3#- and 5'-UTR -> use transcript directly
                    offset = getOffset(t)
                    transcript = t.getSequence(scaffolds[t.contigID])
                    cds = t.coding.getSequence(scaffolds[t.contigID])
                    a = cds.sequence[:12]
                    b = transcript.sequence[offset:offset+12]
                    if a != b:
                        errors += 1
                        print("original", t.strand)
                        print(a)
                        print(b)
                        print(tid)
                        print()
                    else:
                        trans = FastaSequence(ID=tid, seq=transcript.sequence, descr="CDS:%i:%i" % (offset, len(cds)))
                        out.write("%s\n" % (trans))
                else:
                    # make a new mrna
                    newmrna = mRNA(mid=tid, ctg=t.contigID, strand=t.strand, coords=[], score=t.score, cds=Cds())
                    for ex in t:
                        newmrna.append(copy.deepcopy(ex))
                    newmrna.sort()
                    if cod.minimum() <= t.minimum():
                        if t.strand == 0:
                            x = max(1, newmrna[0].beg-20)
                            newmrna[0].beg = x
                        else:
                            x = max(1, newmrna[0].end - 20)
                            newmrna[0].end = x

                    if cod.maximum() >= t.maximum():
                        scafsize = len(scaffolds[t.contigID])
                        if t.strand == 0:
                            x = min(scafsize, newmrna[-1].end+20)
                            newmrna[-1].end = x
                        else:
                            x = min(scafsize, newmrna[-1].beg+20)
                            newmrna[-1].beg = x

                    for ct in t.coding:
                        newmrna.coding.append(ct)

                    # get offset
                    offset = getOffset(newmrna)
                    transcript = newmrna.getSequence(scaffolds[t.contigID])
                    cds = newmrna.coding.getSequence(scaffolds[t.contigID])
                    a = cds.sequence[:12]
                    b = transcript.sequence[offset:offset+12]
                    if a != b:
                        errors += 1
                        print("adjusted", newmrna.strand)
                        print(a)
                        print(b)
                        print(tid)
                        print()
                    else:
                        trans = FastaSequence(ID=tid, seq=transcript.sequence, descr="CDS:%i:%i" % (offset, len(cds)))
                        out.write("%s\n" % (trans))
        out.close()

        print(total, errors)

        fn = "%s_%s.offsets.txt" %(gt, h)

        with open(os.path.join(OUTDIR, fn), 'w') as out:
            for rec in SeqIO.parse(os.path.join(OUTDIR, "%s_%s.offset_transcripts.fa" %(gt, h)), "fasta"):
                desc = rec.description
                tmp = re.split("\s+", desc)
                tmp2 = re.split(":", tmp[1])
                out.write("%s\t%s\t%s\n" %(rec.id, tmp2[1], tmp2[2]))

