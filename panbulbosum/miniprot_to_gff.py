import re
import csv
import sys
import gzip
import multiprocessing
from Bio import  SeqIO
from Bio.Align import substitution_matrices
from Bio.Align import PairwiseAligner

import biosequences
import genedata
import coords
from biosequences import GenomicTemplate, FastaSequence
from genedata import GeneLocus, mRNA, Cds
from coords import Coords


def transformMiniProtGFF(gff_infile : str, gff_outfile : str = "", genomepath : str = "",
        protfile : str = "", prefix : str = ""):
    proteins = {}
    for rec in SeqIO.parse(protfile, "fasta"):
        proteins[rec.id] = rec

    scaffolds = {}
    if genomepath.endswith(".gz"):
        fhd = gzip.open(genomepath, mode='rt')
    else:
        fhd = open(genomepath)
    for rec in SeqIO.parse(fhd, "fasta"):
        scaffolds[rec.id] = GenomicTemplate(ID=rec.id, seq=str(rec.seq).upper())

    genome = parseMiniProtGFF(gff_infile=gff_infile, scaffolds=scaffolds, proteins=proteins, prefix=prefix)

    with open(gff_outfile, 'w') as out:
        for gene in genome:
            out.write("%s\n" %(gene.toGFF3(with_exon=True)))



def alignBatch(subbatch):
    # define some parameters TODO: provide optional parameters
    matrix = substitution_matrices.load("BLOSUM62")
    gap_open = -10.0
    gap_extend = -0.5
    aligner = PairwiseAligner(mode='global',open_gap_score=gap_open,extend_gap_score=gap_extend,\
                              target_end_gap_score = 0.0,query_end_gap_score = 0.0,\
                              substitution_matrix = matrix)
    scores = {}
    for ngid, seq1, seq2 in subbatch:
        alns = aligner.align(seq1, seq2)
        topali = alns[0]
        scores[ngid] = topali.score
    return scores


def parallel_alignment(batches, nproc=3, chunksize=50000):
    data = [batches[i:i+chunksize] for i in range(0,len(batches),chunksize)]
    results = multiprocessing.Pool(processes=nproc).map(alignBatch, data)
    return results


def parseMiniProtGFF(gff_infile : str, scaffolds : dict, proteins : dict, prefix : str, mincov : float = 0.8) -> list:

    uidcounter = 1

    if gff_infile.endswith(".gz"):
        fhd = gzip.open(gff_infile, mode='rt')
    else:
        fhd = open(gff_infile)
    reader = csv.reader(fhd, delimiter='\t')

    genes = {}
    mid2gid = {}
    oldtid2newtid = {}
    for row in reader:
        if len(row) != 9: continue
        if row[2] == "mRNA":
            ctg = row[0]
            src = row[1]
            tags = re.split(";",row[8])
            tid = ""
            desc = {}
            for tag in tags:
                if tag.startswith("ID="):
                    tid = re.sub("ID=","",tag)
                elif tag.startswith("Target="):
                    tg, complexval = re.split("=",tag)
                    trgval, begval, endval = re.split(" ", complexval)
                    desc['src'] = trgval
                    desc['querybeg'] = begval
                    desc['queryend'] = endval
                else:
                    tg, val = re.split("=",tag)
                    desc[tg] = val
            try:
                sc = float(row[5])
            except:
                sc = 0.0
            gid = "{0}_uid_{1:012}".format(prefix, uidcounter)
            newtid = gid + '.m1'
            oldtid2newtid[tid] = newtid
            locus = GeneLocus(gid=gid,ctg=ctg,src=src,score=sc)
            mstrand = 0 if row[6] == '+' else 1
            for tg,val in desc.items():
                locus.tags[tg] = val
            locus.transcripts[newtid] = mRNA(mid=newtid, ctg=ctg, strand=mstrand, coords=[], cds=Cds())
            genes[gid] = locus
            mid2gid[tid] = gid
            uidcounter += 1

        elif row[2] == "CDS":

            ctg = row[0]
            tags = re.split(";",row[8])
            mid = ""
            for tag in tags:
                if tag.startswith("Parent="):
                    mid = re.sub("Parent=","",tag)

            gid = mid2gid[mid]
            newtid = oldtid2newtid[mid]
            a = int(row[3])
            b = int(row[4])
            exstrand = 0 if row[6] == '+' else 1
            beg = min(a,b) if exstrand == 0 else max(a,b)
            end = max(a,b) if exstrand == 0 else min(a,b)
            ct = Coords(ctg=ctg,beg=beg,end=end,strand=exstrand)
            genes[gid].transcripts[newtid].coding.append(ct)
            genes[gid].transcripts[newtid].append(ct)

    fhd.close()

    # now update scoring and filter weak models
    sequences = []
    genome = []
    for gid, gene in genes.items():
        src = gene.tags['src']
        tid = gid + '.m1'
        tc = gene.transcripts[tid].coding
        cds = tc.getSequence(scaffolds[gene.contigID])
        aa  = cds.translate(readThrough=True)
        srcprot = str(proteins[src].seq)
        if float(len(aa))/float(len(srcprot)) < mincov:
            continue
        sequences.append((gid, aa, srcprot))
        gene.tags['hasstart'] = '1' if aa[0] == 'M' else '0'
        gene.tags['hasstop'] = '1' if aa[-1] == '*' else '0'
        gene.tags['intstop'] = '1' if aa[:-1].find('*') >= 0 else '0'
        genome.append(gene)

    # parallel scoring by global NW alignment to ref proteins
    results = parallel_alignment(sequences)
    scores = {}
    for res in results:
        scores.update(res)

    for m in genome:
        m.score = scores[m.geneID]
    genome.sort()
    return genome
