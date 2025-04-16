import re
import os
import csv
import sys
import gzip
from Bio import SeqIO 

sys.path.append($YOUR_PATH_TO_BIOUTILS_MODULE)
from biosequences import FastaSequence, GenomicTemplate
from gffparser import parseTaggedGFF


PROJ_HOME = "DIRECTORY_WITH_PIPELINE_GFFs"
DATA_HOME = "PATH_TO_INPUT_DATA"  # this directory contains pfam & transposon/plastid assignment results
OUTDIR = "PATH_TO_OUTPUT_DIRECTORY"

# first get qualities of input data aka source models
transposons = set()
with open(os.path.join(DATA_HOME, "transposon.codes.txt")) as fhd:
    for line in fhd:
        transposons.add(line.strip())

plastids = set()
with open(os.path.join(DATA_HOME, "plastid.codes.txt")) as fhd:
    for line in fhd:
        plastids.add(line.strip())
      
pfams = {}
with open(os.path.join(DATA_HOME, "domain.codes.txt")) as fhd:
    reader = csv.reader(fhd, delimiter=' ')
    for src, val in reader:
        pfams[src] = float(val)

# these are the genotypes of the project PNWMbarley
genotypes = ["Caribou", "Lightning", "OrNe-1", "Thunder", "Woodie"]   
for genotype in genotypes:

    GFF = os.path.join(PROJ_HOME, "%s.anno.gff.gz" %(genotype))  # GFF of initial preojection
    genes, genome = parseTaggedGFF(GFF)
    finalgenome = []
    for gene in genome:
        src = gene.tags['src']  # the input source is relevant to assign qualities/features to projected model
        rnd = int(gene.tags['round'])
        tid = "%s.m1" %(gene.geneID)
        keep = True
        if rnd == 4 or rnd == 3:
            # for multiple inserted models we keep only those with a PFAM domain
            dom = 100.0
            if src in pfams:
                dom = pfams[src]
            if dom <= 1e-30:
                keep = True
            else:
                keep = False
        if src in transposons:
            gene.tags['tag'] = "transposon_related"
        elif src in plastids:
            gene.tags['tag'] = "plastid_related"
        else:
            gene.tags['tag'] = "projected_gene"
        if keep:
            finalgenome.append(gene)
        
    scaffolds = {}
    fhd = gzip.open("PATH_TO_GENOME_SEQUENCE_OF_%s.fa.gz" %(genotype), mode='rt')
    for rec in SeqIO.parse(fhd, "fasta"):
        scaffolds[rec.id] = GenomicTemplate(ID=rec.id, seq=str(rec.seq).upper())
    fhd.close()
    # now make release files for GFF, cds and proteins
    finalgenome.sort()
    outaa = open(os.path.join(OUTDIR, "%s.prefinal.prot.fa" %(genotype)), 'w')
    outcds = open(os.path.join(OUTDIR, "%s.releaseDec2024.cds.fa" %(genotype)), 'w')
    with open(os.path.join(OUTDIR, "%s.releaseDec2024.gff" %(genotype)), 'w') as out:
        for gene in finalgenome:
            newgid = "HORVU.{0}.PROJ.r1.{1}G{2:08}".format(genotype, gene.contigID, gidcounter)
            newgene = GeneLocus(gid=newgid, ctg=gene.contigID, src="pgsbv1.0")
            newgene.tags['src'] = gene.tags['src']
            newgene.tags['tag'] = gene.tags['tag']
            newgene.tags['method'] = "projection"
            tid = gene.geneID + '.m1'
            t = gene.transcripts[tid]
            newtid = newgid + '.1'
            newmrna = mRNA(mid=newtid, ctg=gene.contigID, strand=t.strand, coords=[], cds=Cds())
            for ex in t:
                newmrna.append(ex)
            for ex in t.coding:
                newmrna.coding.append(ex)
            newgene.addmRNA(newmrna)
            out.write("%s\n" %(newgene.toGFF3(with_exon=True)))
            gidcounter += 1   
            for tid, t in newgene.transcripts.items():
                cds = t.coding.getSequence(scaffolds[t.contigID])
                cds.ID = tid
                outcds.write("%s\n" %(cds))
                aa = cds.translate(readThrough=True)
                prot = FastaSequence(ID=tid, seq=aa)
                outaa.write("%s\n" %(prot))
    outaa.close()
    outcds.close()


