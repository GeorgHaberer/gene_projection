import os
import re
import sys
import csv
import gzip
import bisect
from Bio import SeqIO
from operator import itemgetter, attrgetter
sys.path.append("/home/pgsb/georg.haberer/src/pyprojects/bioutils")
from gffparser import parseTaggedGFF
from biosequences import GenomicTemplate, FastaSequence

from typing import Callable


def outputAnno(annotation: list, scaffolds: dict, outdir: str, genotype: str, version: int, rndcount: int):
    outaa = open(os.path.join(outdir, "%s.anno.v%i.prot.fasta" % (genotype, version)), 'w')
    outcds = open(os.path.join(outdir, "%s.anno.v%i.cds.fasta" % (genotype, version)), 'w')
    rounds = {}
    for rnd in range(1, rndcount + 1):
        rounds[rnd] = 0
    with open("%s.anno.v%i.gff" % (genotype, version), 'w') as out:
        for anno in annotation:
            rnd = int(anno.tags['round'])
            rounds[rnd] += 1
            out.write("%s\n" % (anno.toGFF3(with_exon=True)))
            tid = anno.geneID + '.m1'
            t = anno.transcripts[tid]
            cds = t.coding.getSequence(scaffolds[anno.contigID])
            cds.ID = tid
            outcds.write("%s\n" % (cds))
            aa = cds.translate(readThrough=True)
            prot = FastaSequence(ID=tid, seq=aa)
            outaa.write("%s\n" % (prot))

    outaa.close()
    outcds.close()
    with open(os.path.join(outdir, "%s.anno.v%i.rounds.txt" % (genotype, version)), 'w') as out:
        for rnd in range(1, rndcount+1):
            out.write("%i\t%i\n" % (rnd, rounds[rnd]))

def complementAnnotation(matches: list, annotation: list, rnd: int, seensrc: set = set(), seenmeta: set = set(),
                         uniqsrc: bool = True, uniqmeta: bool = True, sorter : Callable = attrgetter('score')) \
        -> tuple[list, set, set]:

    matches.sort(key=sorter, reverse=True)
    # if no annotation yet, initialize with top scoring match
    startm = 0
    if not len(annotation):
        m = matches[0]
        m.tags['round'] = "%i" % (rnd)
        annotation.append(m)
        startm = 1
        src = m.tags['src']
        seensrc.add(src)
        cid = m.tags['cluid']
        seenmeta.add(cid)

    # now for each match insert matches by decreasing scores
    # if this match does not overlap any previous match by coding region
    for m in matches[startm:]:
        i = bisect.bisect(annotation, m)
        cmin = m.codingMin()
        cmax = m.codingMax()
        src = m.tags['src']
        cid = m.tags['cluid']

        ########################################################################
        # some filters
        # if we require maximal one match per informant source (unique=True)
        if uniqmeta:
            if cid in seenmeta:
                continue
        if uniqsrc:
            if src in seensrc:
                continue
        # only add non-overlapping matches
        for j in range(max(0, i - 5), min(len(annotation), i + 6)):
            if m.contigID != annotation[j].contigID:
                continue
            cmina = annotation[j].codingMin()
            cmaxa = annotation[j].codingMax()
            if cmin <= cmaxa and cmax >= cmina:
                break
        else:
            m.tags['round'] = "%i" % (rnd)
            annotation.insert(i, m)
            seensrc.add(src)
            seenmeta.add(cid)
    annotation.sort()
    return annotation, seensrc, seenmeta


def getMatches(ANNO_HOME : str, maxscores : dict, weights : dict, domains : dict, ocounts : dict, metamap : dict,
               transposon_codes : set = set(), plastid_codes : set = set(), relcutoff : float = 0.85,
               chrunknown = False, complete : bool = True) -> list:

    os.chdir(ANNO_HOME)
    allmatches = []
    gff_templates = ["psl.gff.gz", "miniprot.gff.gz", "raw.gff.gz"]
    programs = ["psl", "mp", "mm"] 
    if chrunknown == False:
        fafiles = [fafile for fafile in os.listdir('.') if fafile.endswith(".fa")]
        for fn in fafiles:
            ctgid = re.sub(r"\.fa$", "", fn)
            for gffbase, prg in zip(gff_templates, programs):
                os.chdir(ANNO_HOME)
                gff = "%s.%s" %(ctgid, gffbase)
                tmpgenes, tmpgenome = parseTaggedGFF(gff)
                for gene in tmpgenome:
                    # we do some filtering
                    if gene.tags['intstop'] == '1': continue  # we restrict to contiguous ORFs
                    src = gene.tags['src']
                    sc = float(gene.score) / maxscores[src]
                    if sc < relcutoff: continue
                    atg = int(gene.tags['hasstart'])
                    stop = int(gene.tags['hasstop'])
                    if complete:
                        if atg == 0 or stop == 0:
                            continue
                    pfam = domains[src] if src in domains else 100.0  # e-value of 100 for no domain detected
                    repeat = True if src in transposon_codes else False
                    plastid = True if src in plastid_codes else False
                    try:
                        cnt = ocounts[src]
                    except KeyError:
                        cnt = 0
                    w = weights[src]
                    cid = metamap[src]
                    gene.tags['cluid'] = cid
                    tid = gene.geneID + '.m1'
                    t = gene.transcripts[tid]
                    pos = [ex.minimum() for ex in t.coding]
                    pos += [ex.maximum() for ex in t.coding]
                    pos.sort()
                    intronsizes = []
                    if len(pos) == 2:
                        intronsizes.append(0)
                    else:
                        i = 1
                        while i < len(pos)-2:
                            intronsizes.append(pos[i+1]-pos[i])
                            i += 2
                    maxintron = max(intronsizes)
                    sc1 = gene.score*(float(w)/26.0)*(float(cnt)/26.0)
                    sc2 = gene.score*(float(cnt)/26.0)
                    gene.score1 = sc1
                    gene.score2 = sc2
                    gene.prg = 0
                    if prg == "psl":
                        gene.prg = 2
                    elif prg == "mm":
                        gene.prg = 1
                    allmatches.append((gene, atg, stop, pfam, repeat, cnt, w, maxintron, sc, prg, cid, sc1, sc2, plastid))
    
    else:

        for gffbase, prg in zip(gff_templates, programs):
            os.chdir(ANNO_HOME)
            gff = "chrUn.%s" %(gffbase)
            tmpgenes, tmpgenome = parseTaggedGFF(gff)
            for gene in tmpgenome:

                # we do some filtering
                if gene.tags['intstop'] == '1': continue  # we restrict to contiguous ORFs
                atg = int(gene.tags['hasstart'])
                stop = int(gene.tags['hasstop'])
                if atg == 0 or stop == 0:
                    continue  # we do not accept incomplete matches on unanchored scaffolds
                src = gene.tags['src']
                sc = float(gene.score) / maxscores[src]
                if sc < relcutoff: continue
                # fill in remaining stats
                pfam = domains[src] if src in domains else 100.0  # e-value of 100 for no domain detected
                repeat = True if src in transposon_codes else False
                plastid = True if src in plastid_codes else False
                try:
                    cnt = ocounts[src]
                except KeyError:
                    cnt = 0
                w = weights[src]
                cid =  metamap[src]
                gene.tags['cluid'] = cid
                tid = gene.geneID + '.m1'
                t = gene.transcripts[tid]
                pos = [ex.minimum() for ex in t.coding]
                pos += [ex.maximum() for ex in t.coding]
                pos.sort()
                intronsizes = []
                if len(pos) == 2:
                    intronsizes.append(0)
                else:
                    i = 1
                    while i < len(pos)-2:
                        intronsizes.append(pos[i+1]-pos[i])
                        i += 2

                maxintron = max(intronsizes)
                sc1 = gene.score*(float(w)/26.0)*(float(cnt)/26.0)
                sc2 = gene.score*(float(cnt)/26.0)
                gene.score1 = sc1
                gene.score2 = sc2
                gene.prg = 0
                if prg == "psl":
                    gene.prg = 2
                elif prg == "mm":
                    gene.prg = 1
                allmatches.append((gene, atg, stop, pfam, repeat, cnt, w, maxintron, sc, prg, cid, sc1, sc2, plastid))

    return allmatches


# ----------------------------------------------------------------------------

def annoOneLine(genotype : str, genomepath : str):

    BULB_ROOT = "/lustre/groups/pgsb/workspaces/georg.haberer/bulbosum/revision"

    ANNO_HOME = os.path.join(BULB_ROOT, "projections", genotype)
    DATA_HOME = os.path.join(BULB_ROOT, "datasets")

    # read source features
    os.chdir(DATA_HOME)
    maxscores = {}
    transposon_codes = set()
    plastidcodes = set()
    weights = {}
    domains = {}
    ocounts = {}
    metamap = {}

    with open("hbulb.selfhits.hit") as fhd:
        reader = csv.reader(fhd, delimiter='\t')
        for tid, sc in reader:
            maxscores[tid] = float(sc)

    with open("hbulb.transposon.codes") as fhd:
        for line in fhd:
            transposon_codes.add(line.strip())

    with open("hbulb.plastid.codes") as fhd:
        for line in fhd:
            plastidcodes.add(line.strip())
    with open("hbulb.domain.codes.txt") as fhd:
        # generated by : tr -s ' ' < allproteins.cdhit.domtblout | egrep -v ^# | cut -d ' ' -f1,7 > domain.codes.txt
        reader = csv.reader(fhd, delimiter=' ')
        for tid, val in reader:
            try:
                domains[tid] = min(domains[tid], float(val))
            except KeyError:
                domains[tid] = float(val)

    with open("hbulb.geneweights.txt") as fhd:
        reader = csv.reader(fhd, delimiter='\t')
        for tid, cnt in reader:
            weights[tid] = int(cnt)

    with open("hbulb.metacluster.txt") as fhd:
        reader = csv.reader(fhd, delimiter='\t')
        for cluid, tid in reader:
            metamap[tid] = int(cluid)

    with open("hbulb.orthocount.txt") as fhd:
        reader = csv.reader(fhd, delimiter='\t')
        for tid, cnt in reader:
            ocounts[tid] = int(cnt)

    for tid, w in weights.items():
        if tid not in metamap:
            metamap[tid] = 0
        if tid not in ocounts:
            ocounts[tid] = 0

    # read scaffolds for reporting cds/protein sequences of annotations
    scaffolds = {}
    if genomepath.endswith(".gz"):
        fhd = gzip.open(genomepath, mode='rt')
    else:
        fhd = open(genomepath)
    for rec in SeqIO.parse(fhd, "fasta"):
        scaffolds[rec.id] = GenomicTemplate(ID=rec.id, seq=str(rec.seq).upper())

    # read all (contiguous) matches of all programs
    if genotype.endswith("unanchor"):
        cun = True 
    else:
        cun = False    
    allmatches = getMatches(ANNO_HOME=ANNO_HOME, maxscores=maxscores, weights=weights, domains=domains, chrunknown=cun,
                            ocounts=ocounts, metamap=metamap, plastid_codes=plastidcodes, transposon_codes=transposon_codes)
    

    # now run rule based insertion
    #--------------------------------------------------------------
    version = 1
    seen = set()
    meta = set()
    annotation = []
    rndcount = 1
    fsort = attrgetter('score2', 'prg')

    trgmatches = [x[0] for x in allmatches if x[3] <= 0.01 and x[4] == False and x[10] > 0]
    annotation, seen, meta = complementAnnotation(trgmatches, annotation, rndcount, seen, meta, uniqsrc=False,
                                                  uniqmeta=True, sorter=fsort)
    annotation.sort()
    rndcount += 1

    trgmatches = [x[0] for x in allmatches if x[4] == False and x[13] == False and x[10] > 0]
    annotation, seen, meta = complementAnnotation(trgmatches, annotation, rndcount, seen, meta, uniqsrc=False,
                                                  uniqmeta=True, sorter=fsort)
    annotation.sort()
    rndcount += 1

    trgmatches = [x[0] for x in allmatches if x[3] <= 0.01 and x[4] == False and x[13] == False and x[5] >= 6 and x[8] >= 0.9]
    annotation, seen, meta = complementAnnotation(trgmatches, annotation, rndcount, seen, meta, uniqsrc=False,
                                                  uniqmeta=False, sorter=fsort)
    annotation.sort()
    rndcount += 1
    trgmatches = [x[0] for x in allmatches if (x[3] <= 0.01 or x[5] >= 2) and x[4] == False and x[8] >= 0.9]
    annotation, seen, meta = complementAnnotation(trgmatches, annotation, rndcount, seen, meta, uniqsrc=False,
                                                  uniqmeta=False, sorter=fsort)
    annotation.sort()
    rndcount += 1

    trgmatches = [x[0] for x in allmatches if x[4] == True and x[8] >= 0.9]
    annotation, seen, meta = complementAnnotation(trgmatches, annotation, rndcount, seen, meta, uniqsrc=True,
                                                  uniqmeta=True, sorter=fsort)
    annotation.sort()

    outputAnno(annotation, scaffolds, ANNO_HOME, genotype, version, rndcount)


    os.chdir(DATA_HOME)
    ocounts = {}
    metamap = {}

    with open("hbulb.metacluster.N0s.txt") as fhd:
        reader = csv.reader(fhd, delimiter='\t')
        for cluid, tid in reader:
            metamap[tid] = int(cluid)

    with open("hbulb.orthocount.txt") as fhd:
        reader = csv.reader(fhd, delimiter='\t')
        for tid, cnt in reader:
            ocounts[tid] = int(cnt)

    for tid, w in weights.items():
        if tid not in metamap:
            metamap[tid] = 0
        if tid not in ocounts:
            ocounts[tid] = 0

    # read all (contiguous) matches of all programs
    allmatches = getMatches(ANNO_HOME=ANNO_HOME, maxscores=maxscores, weights=weights, domains=domains, chrunknown=cun,
                            ocounts=ocounts, metamap=metamap, plastid_codes=plastidcodes, transposon_codes=transposon_codes)
    

    # now run rule based insertion
    #--------------------------------------------------------------
    version = 2
    seen = set()
    meta = set()
    annotation = []
    rndcount = 1
    fsort = attrgetter('score2', 'prg')

    trgmatches = [x[0] for x in allmatches if x[3] <= 0.01 and x[4] == False and x[10] > 0]
    annotation, seen, meta = complementAnnotation(trgmatches, annotation, rndcount, seen, meta, uniqsrc=False,
                                                  uniqmeta=True, sorter=fsort)
    annotation.sort()
    rndcount += 1

    trgmatches = [x[0] for x in allmatches if x[4] == False and x[13] == False and x[10] > 0]
    annotation, seen, meta = complementAnnotation(trgmatches, annotation, rndcount, seen, meta, uniqsrc=False,
                                                  uniqmeta=True, sorter=fsort)
    annotation.sort()
    rndcount += 1

    trgmatches = [x[0] for x in allmatches if x[3] <= 0.01 and x[4] == False and x[13] == False and x[5] >= 6 and x[8] >= 0.9]
    annotation, seen, meta = complementAnnotation(trgmatches, annotation, rndcount, seen, meta, uniqsrc=False,
                                                  uniqmeta=False, sorter=fsort)
    annotation.sort()
    rndcount += 1

    trgmatches = [x[0] for x in allmatches if (x[3] <= 0.01 or x[5] >= 2) and x[4] == False and x[8] >= 0.9]
    annotation, seen, meta = complementAnnotation(trgmatches, annotation, rndcount, seen, meta, uniqsrc=False,
                                                  uniqmeta=False, sorter=fsort)
    annotation.sort()
    rndcount += 1

    trgmatches = [x[0] for x in allmatches if x[4] == True and x[8] >= 0.9]
    annotation, seen, meta = complementAnnotation(trgmatches, annotation, rndcount, seen, meta, uniqsrc=True,
                                                  uniqmeta=True, sorter=fsort)
    annotation.sort()

    outputAnno(annotation, scaffolds, ANNO_HOME, genotype, version, rndcount)

# ================================================================
if __name__ == '__main__':

    gt = sys.argv[1]
    genomefile = sys.argv[2]
    annoOneLine(genotype=gt, genomepath=genomefile)


