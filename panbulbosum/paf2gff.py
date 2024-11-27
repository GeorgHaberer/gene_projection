import re
import csv
import sys
import gzip
import copy
import multiprocessing

from Bio import SeqIO
from Bio.Align import substitution_matrices
from Bio.Align import PairwiseAligner

from operator import itemgetter

import biosequences
import genedata
import coords
from biosequences import DnaSequence, FastaSequence, GenomicTemplate
from genedata import GeneLocus, mRNA, Cds
from coords import Coords

from cslongpafentry import extractPafTranscript


#====================================================================
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


def parallel_alignment(batches, nproc=4, chunksize=5000):
    data = [batches[i:i+chunksize] for i in range(0,len(batches),chunksize)]
    results = multiprocessing.Pool(processes=nproc).map(alignBatch, data)
    return results


def makeCodingByORF(paftranscriptseq, entry, mrna):

    orfs = paftranscriptseq.getAllORFs(minLength=30)
    if not len(orfs):
        return None
    tmp = [(len(orf[1]), orf) for orf in orfs]
    tmp.sort(key=itemgetter(0))
    tophit = tmp[-1]
    topstart = tophit[1][0]
    toporf = tmp[-1][1][1]
    toporf.ID = entry.query
    for i, letter in enumerate(toporf.sequence):
        if letter == 'M': break
    else:
        i = 0
    atgstart = topstart + 3 * i
    trimup = atgstart - 1
    trimdown = len(paftranscriptseq) - (topstart - 1 + 3 * len(toporf))
    if trimup > 0:
        pafcoding = mrna.trimBySizer(direction=5, cutsize=trimup)
    else:
        pafcoding = copy.deepcopy(mrna)
    if trimdown > 0:
        pafcoding = pafcoding.trimBySizer(direction=3, cutsize=trimdown)
    return pafcoding


#-----------------------------------------------------------------------
def pafmap2Gff(paffile : str, ctgfile : str, protfile : str, offsetfile : str) -> [list, dict]:

    """transforms PAF formatted transcript mappings <paffilepath>, mapped to scaffold <ctgfilepath>
    and using their corresponding protein sequences <protfilepath> to gff gene models.
    :returns list of gene models"""

    # get protein sequences for scoring
    proteins = {}
    for rec in SeqIO.parse(protfile, "fasta"):
        proteins[rec.id] = rec

    srcoffsets = {}
    with open(offsetfile) as fhd:
        reader = csv.reader(fhd, delimiter='\t')
        for tid, beg, end in reader:
            srcoffsets[tid] = (int(beg), int(end))

    # get sequence of scaffold
    scaffolds = {}
    if ctgfile.endswith(".gz"):
        fhd = gzip.open(ctgfile, mode='rt')
    else:
        fhd = open(ctgfile)
    for rec in SeqIO.parse(fhd, "fasta"):
        scaffolds[rec.id] = GenomicTemplate(ID=rec.id, seq=str(rec.seq).upper())
    fhd.close()

    # do the job
    models = []
    sequences = []
    if paffile.endswith('.gz'):
        fhd = gzip.open(paffile, 'rt')
    else:
        fhd = open(paffile)
    reader = csv.reader(fhd, delimiter='\t')
    uid = 1

    for entry in map(extractPafTranscript, reader):

        pafgid = "uid_{0:010}.{1}.hv{2}".format(uid, entry.query, entry.target)
        paftid = pafgid + '.m1'
        uid += 1
        # tokenize CS data
        data = []
        token = []
        TAGS = ['=', '*', '+', '~', '-']
        for letter in entry.cslong:
            if letter in TAGS:
                if len(token):
                    data.append(token)
                    token = []
            token.append(letter)
        if len(token):
            data.append(token)

        # make genestructure
        structlist = []
        aktexonseq = []
        for token in data:
            tag = token[0]
            datum = token[1:]
            if tag == '=':
                aktexonseq.append(''.join(datum))
            elif tag == '+':
                # sequence in query but not reference! -> skip
                continue
            elif tag == '-':
                # additional sequence in reference
                aktexonseq.append(''.join(datum))
            elif tag == '*':
                # a substitution, first letter is reference base
                aktexonseq.append(datum[0])
            elif tag == '~':
                structlist.append(('e', ''.join(aktexonseq).upper()))
                aktexonseq = []
                intronsize = int(re.sub("[a-zA-Z]", "", ''.join(datum)))
                structlist.append(('i', intronsize))

        structlist.append(('e', ''.join(aktexonseq).upper()))
        tmp = []
        mrna = mRNA(mid=paftid, ctg=entry.target, strand=entry.strand, coords=[], cds=Cds())
        gpos1 = entry.targetmin + 1
        for elmtype, elmdata in structlist:
            if elmtype == 'e':
                tmp.append(elmdata)
                gpos2 = gpos1 + len(elmdata) - 1
                b = min(gpos1, gpos2) if entry.strand == 0 else max(gpos1, gpos2)
                e = max(gpos1, gpos2) if entry.strand == 0 else min(gpos1, gpos2)
                ct = Coords(ctg=entry.target, beg=b, end=e, strand=entry.strand)
                mrna.append(ct)
                gpos1 = gpos2 + 1
            elif elmtype == 'i':
                gpos1 += elmdata
        gfftranscript = mrna.getSequence(scaffolds[entry.target])
        gfftranscript.ID = 'gff:%s' % (mrna.tid)
        paftranscriptseq = DnaSequence(ID=entry.query, seq=''.join(tmp))
        if entry.strand == 1:
            paftranscriptseq = paftranscriptseq.getReverseComplement()
        if paftranscriptseq.sequence != gfftranscript.sequence:
            # very rare case of wild minimap mappings (eg double introns), skip
            continue

        offbeg, offend = srcoffsets[entry.query]
        DERIVED = False

        if offbeg < entry.querymin:
            # look for longest match
            pafcoding = makeCodingByORF(paftranscriptseq, entry, mrna)
        else:
            start = offbeg - entry.querymin
            if paftranscriptseq.sequence[start:start+3] == "ATG":
                # take cds coordinates as offset suggests
                aa = paftranscriptseq.translate(fromPos=start)
                trimup = start
                trimdown = len(paftranscriptseq) - (start + 3 * len(aa))
                if trimup > 0:
                    pafcoding = mrna.trimBySizer(direction=5, cutsize=trimup)
                else:
                    pafcoding = copy.deepcopy(mrna)
                if trimdown > 0:
                    pafcoding = pafcoding.trimBySizer(direction=3, cutsize=trimdown)
                DERIVED = True
            else:
                pafcoding = makeCodingByORF(paftranscriptseq, entry, mrna)

        if pafcoding is None:
            continue

        newgene = GeneLocus(gid=pafgid, ctg=entry.target, src='pgsb.v0', score=0)
        for ct in pafcoding:
            mrna.coding.append(ct)
        mrnacds = mrna.coding.getSequence(scaffolds[mrna.contigID])
        aa = mrnacds.translate(readThrough=True)
        seq = str(proteins[entry.query].seq)
        if len(aa) < 0.8*len(seq):
            continue
        sequences.append((pafgid, aa, seq))
        newgene.tags['hasstart'] = '1' if aa[0] == 'M' else '0'
        newgene.tags['hasstop'] = '1' if aa[-1] == '*' else '0'
        newgene.tags['intstop'] = '1' if aa[:-1].find('*') >= 0 else '0'
        newgene.tags['origin'] = '1' if DERIVED else '0'
        newgene.tags['src'] = entry.query
        newgene.addmRNA(mrna)
        models.append(newgene)

    # parallel scoring by global NW alignment to ref proteins
    results = parallel_alignment(sequences)
    scores = {}
    for res in results:
        scores.update(res)

    for m in models:
        m.score = scores[m.geneID]
    return models
