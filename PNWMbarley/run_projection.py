import os
import sys
import bisect
from typing import Callable
from operator import attrgetter

sys.path.append("/home/pgsb/georg.haberer/src/pyprojects/bioutils")
sys.path.append("/home/pgsb/georg.haberer/src/pyprojects/bioutils")
from biosequences import FastaSequence
from mymatches import getMatches, getScaffolds, getSupplementaryInfos



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


# ----------------------------------------------------------------------------

def annoOneLine(genotype: str, genomepath: str):

    WILDHV_ROOT = "/lustre/groups/pgsb/workspaces/georg.haberer/wildbarley"
    JEFF_HOME = "/lustre/groups/pgsb/workspaces/georg.haberer/barley_jeff"
    ANNO_HOME = os.path.join(JEFF_HOME, "projections", genotype)
    DATA_HOME = os.path.join(WILDHV_ROOT, "datasets")

    

    # read source features
    os.chdir(DATA_HOME)
    maxscores, ocounts, metamap, weights, domains, transposon_codes, plastid_codes = getSupplementaryInfos(DATA_HOME)
    scaffolds = getScaffolds(genomepath)
    buscos = set()
    with open("buscoids.txt") as fhd:
        for line in fhd:
            buscos.add(line.strip())

    # get all matches as list for annotation

    # read all (contiguous) matches of all programs
    os.chdir(ANNO_HOME)
    allmatches = getMatches(ANNO_HOME, maxscores=maxscores, weights=weights, domains=domains, ocounts=ocounts,
                            metamap=metamap, transposon_codes=transposon_codes)

    
    
    ######################################################
    # now run rule based insertion

    # --------------------------------------------
    # prefer first trancript-based then protmap for equal scoring matches
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

    trgmatches = [x[0] for x in allmatches if x[4] == False and x[10] > 0]
    annotation, seen, meta = complementAnnotation(trgmatches, annotation, rndcount, seen, meta, uniqsrc=False,
                                                  uniqmeta=True, sorter=fsort)
    annotation.sort()
    rndcount += 1

    trgmatches = [x[0] for x in allmatches if x[3] <= 0.01 and x[4] == False and x[5] >= 6 and x[8] >= 0.9]
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

    return


# ================================================================
if __name__ == '__main__':

    # due to fixed paths (see above) we can simply call this script in the batch queue with the genotype and
    # absolute path of the genome sequence, the rest is done by our hard-coded paths
    gt = sys.argv[1]
    genomefile = sys.argv[2]
    annoOneLine(genotype=gt, genomepath=genomefile)
