from __future__ import annotations

import copy
import gzip
import re
import os
import csv
import sys
import multiprocessing
from operator import itemgetter
from Bio import  SeqIO
from Bio.Align import substitution_matrices
from Bio.Align import PairwiseAligner

from genedata import GeneLocus, mRNA, Cds
from coords import Coords
from biosequences import GenomicTemplate, FastaSequence


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
    print("batch aligned")
    return scores


def parallel_alignment(batches, nproc=3, chunksize=50000):
    data = [batches[i:i+chunksize] for i in range(0,len(batches),chunksize)]
    results = multiprocessing.Pool(processes=nproc).map(alignBatch, data)
    return results


def makeCodingByORF(transcriptseq, mrna):

    orfs = transcriptseq.getAllORFs(minLength=30)
    if not len(orfs):
        return None
    tmp = [(len(orf[1]), orf) for orf in orfs]
    tmp.sort(key=itemgetter(0))
    tophit = tmp[-1]
    topstart = tophit[1][0]
    toporf = tmp[-1][1][1]
    for i, letter in enumerate(toporf.sequence):
        if letter == 'M': break
    else:
        i = 0
    atgstart = topstart + 3 * i
    trimup = atgstart - 1
    trimdown = len(transcriptseq) - (topstart - 1 + 3 * len(toporf))
    if trimup > 0:
        pslcoding = mrna.trimBySizer(direction=5, cutsize=trimup)
    else:
        pslcoding = copy.deepcopy(mrna)
    if trimdown > 0:
        pslcoding = pslcoding.trimBySizer(direction=3, cutsize=trimdown)
    return pslcoding


def psl2gff(pslfile : str, scaffolds : dict, offsets : dict, proteins : dict, prefix : str = "",
            mincov : float = 0.7, uidcounter : int = 1) -> tuple[list, int]:

    genome = []
    sequences = []

    with open(pslfile) as fhd:
        reader = csv.reader(fhd, delimiter='\t')
        header = next(reader)
        header = next(reader)
        header = next(reader)
        header = next(reader)
        header = next(reader)
        for counter, row in enumerate(reader):
            matchsize = int(row[0])
            querysize = int(row[10])
            if float(querysize)/float(matchsize) < mincov:
                continue
            ori = row[8]
            src = row[9]
            trg = row[13]
            blksizes = [int(x) for x in re.split(',', row[18])[:-1]]
            questarts = [int(x) for x in re.split(",", row[19])[:-1]]
            trgstarts = [int(x) for x in re.split(",", row[20])[:-1]]
            gid = "{0}_uid_{1:012}".format(prefix, uidcounter)
            gene = GeneLocus(gid=gid, ctg=trg)
            gene.tags['blocksizes'] = re.sub(",$", "", row[18])
            gene.tags['querystart'] = re.sub(",$", "", row[19])
            gene.tags['targetstart'] = re.sub(",$", "", row[20])
            gene.tags['src'] = src

            mid = gid + '.m1'
            strand = 0 if ori == '+' else 1
            mrna = mRNA(mid=mid, ctg=trg, strand=strand, coords=[], score=None, cds=Cds())
            for bs, ts in zip(blksizes, trgstarts):
                a = ts + 1
                b = a + bs - 1
                beg = min(a, b) if strand == 0 else max(a, b)
                end = max(a, b) if strand == 0 else min(a, b)
                ex = Coords(beg=beg, end=end, ctg=trg, strand=strand)
                mrna.append(ex)

            gfftranscript = mrna.getSequence(scaffolds[trg])
            offbeg, offend = offsets[src]
            DERIVED = False
            querymin = min(map(int, questarts))
            pslcoding = None

            if offbeg < querymin:
                # look for longest match
                pslcoding = makeCodingByORF(gfftranscript, mrna)
            else:
                start = offbeg - querymin
                if gfftranscript.sequence[start:start + 3] == "ATG":
                    # take cds coordinates as offset suggests
                    aa = gfftranscript.translate(fromPos=start)
                    trimup = start
                    trimdown = len(gfftranscript) - (start + 3 * len(aa))
                    if trimup > 0:
                        pslcoding = mrna.trimBySizer(direction=5, cutsize=trimup)
                    else:
                        pslcoding = copy.deepcopy(mrna)
                    if trimdown > 0:
                        pslcoding = pslcoding.trimBySizer(direction=3, cutsize=trimdown)
                    DERIVED = True
                else:
                    pslcoding = makeCodingByORF(gfftranscript, mrna)

            if pslcoding is None:
                continue
            for ct in pslcoding:
                mrna.coding.append(ct)
            mrnacds = mrna.coding.getSequence(scaffolds[mrna.contigID])
            aa = mrnacds.translate(readThrough=True)
            srcprot = str(proteins[src].seq)
            if float(len(aa))/float(len(srcprot)) < mincov:
                continue
            sequences.append((gid, aa, srcprot))
            gene.tags['hasstart'] = '1' if aa[0] == 'M' else '0'
            gene.tags['hasstop'] = '1' if aa[-1] == '*' else '0'
            gene.tags['intstop'] = '1' if aa[:-1].find('*') >= 0 else '0'
            gene.tags['origin'] = '1' if DERIVED else '0'
            gene.tags['src'] = src
            gene.addmRNA(mrna)
            genome.append(gene)
            uidcounter += 1

    # parallel scoring by global NW alignment to ref proteins
    results = parallel_alignment(sequences, nproc=10)
    scores = {}
    for res in results:
        scores.update(res)

    for m in genome:
        m.score = scores[m.geneID]

    genome.sort()
    return genome, uidcounter




def main():

    pslpath = sys.argv[1]
    genomepath = sys.argv[2]
    offsetpath = sys.argv[3]
    protfile = sys.argv[4]
    prefix = sys.argv[5]
    outfile = sys.argv[6]

    scaffolds = {}
    if genomepath.endswith(".gz"):
        fhd = gzip.open(genomepath, mode='rt')
    else:
        fhd = open(genomepath)
    for rec in SeqIO.parse(fhd, "fasta"):
        scaffolds[rec.id] = GenomicTemplate(ID=rec.id, seq=str(rec.seq).upper())

    srcoffsets = {}
    with open(offsetpath) as fhd:
        reader = csv.reader(fhd, delimiter='\t')
        for tid, beg, end in reader:
            srcoffsets[tid] = (int(beg), int(end))

    genome = psl2gff(pslpath, scaffolds, srcoffsets, protfile, prefix)
    with open(outfile, 'w') as out:
        for gene in genome:
            if gene.score < 100: continue
            out.write("%s\n" %(gene.toGFF3(with_exon=True)))



if __name__ == '__main__':
    main()
