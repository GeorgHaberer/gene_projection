from __future__ import annotations

import os
import sys
import csv
import bisect
import argparse
from operator import attrgetter

sys.path.append("/home/georg/PycharmProjects/bioutils")
from gffparser import parseTaggedGFF


# ----------------------------------------------------------------------------
def complementAnnotation(matches : list, annotation : list, rnd : int, seensrc=set(), uniqsrc=True) -> tuple(list, set):
    """inserts from matches non-overlapping matches into current gene annotation list
    :return filled annotation list, updated tracking set of inserted gene IDs"""
    matches.sort(key=attrgetter('score'), reverse=True)
    # if no annotation yet, initialize with top scoring match
    startm = 0
    if not len(annotation):
        m = matches[0]
        m.tags['round'] = "%i" % (rnd)
        annotation.append(m)
        startm = 1
        src = m.tags['src']
        seensrc.add(src)

    # now for each match insert matches by decreasing scores
    # if this match does not overlap any previous match by coding region
    for m in matches[startm:]:

        src = m.tags['src']
        # if we require maximal one match per informant source (unique=True)
        # skip it if we have seen this informant already
        if uniqsrc:
            if src in seensrc:
                continue

        i = bisect.bisect(annotation, m)
        cmin = m.codingMin()
        cmax = m.codingMax()
        # only add non-overlapping matches, overlap of CDS ranges is relevant
        for j in range(max(0, i - 5), min(len(annotation), i + 6)):
            if m.contigID != annotation[j].contigID:
                continue
            cmina = annotation[j].codingMin()
            cmaxa = annotation[j].codingMax()
            if cmin <= cmaxa and cmax >= cmina:
                break
        else:
            m.tags['round'] = "%i" % (rnd)  # each model can be later identified at which filters it has been inserted
            annotation.insert(i, m)
            seensrc.add(src)  # here we keep track which source transcript ids already have been inserted

    return annotation, seensrc


# ----------------------------------------------------------------------------
def annoOneGenotype(genotype : str, projroot: str, pfam_codes : set, plastid_codes : set, transposon_codes : set,
                    maxscores : dict) -> None:

    """gene projection for one genotype: inserts matches into annotation in a stepwise way and based on rules"""

    GT_ROOT = os.path.join(projroot, genotype)
    os.chdir(GT_ROOT)

    # we first process quality bins
    allsourceids = set(list(maxscores.keys()))
    genematch_codes = allsourceids - (transposon_codes | plastidcodes) # only gene, no transposon/plastid-related
    hcmatch_codes = genematch_codes & pfam_codes  # only genes with pfam domain; 'high confidence'

    # second we read in all matches of our minimap2 alignments
    allmatches = []
    raw_gffs = [fn for fn in os.listdir(GT_ROOT) if fn.endswith("raw.gff.gz")]
    for gff in raw_gffs:
        genes, genome = parseTaggedGFF(gff)
        for gene in genome:
            allmatches.append(gene)

    annotation = []
    seensources = set()

    ####################################################################################################
    # in the first two rounds we target complete matches with contiguous ORFs AND pfam hit (hc matches).
    # --------------------------------------------------------------------------------------------------
    # round 1: candidate match with good coverage of source, each source id maximally one insertion
    runde = 1
    completematches = [x for x in allmatches if x.tags['intstop'] == '0' and x.tags['hasstart'] == '1'
                       and x.tags['hasstop'] == '1' and x.score / maxscores[x.tags['src']] >= 0.85]

    trgmatches = [x for x in completematches if x.geneID in hcmatch_codes]
    annotation, seensources = complementAnnotation(trgmatches, annotation, runde, seensrc=seensources, uniqsrc=True)
    annotation.sort()
    # round 2: candidate match with high coverage of source, each source id can be inserted multiple
    #          times -> allow high quality copy number variations
    runde = 2
    completematches = [x for x in allmatches if x.tags['intstop'] == '0' and x.tags['hasstart'] == '1'
                       and x.tags['hasstop'] == '1' and x.score / maxscores[x.tags['src']] >= 0.95]

    trgmatches = [x for x in completematches if x.geneID in hcmatch_codes]
    annotation, seensources = complementAnnotation(trgmatches, annotation, runde, seensrc=seensources, uniqsrc=False)
    annotation.sort()

    ####################################################################################################
    # next we target complete matches with contiguous ORFs and no/opt. pfam hit ('gene' matches).
    # --------------------------------------------------------------------------------------------------
    # round 3: candidate match with good coverage of source, each source id maximally one insertion
    runde = 3
    completematches = [x for x in allmatches if x.tags['intstop'] == '0' and x.tags['hasstart'] == '1'
                       and x.tags['hasstop'] == '1' and x.score / maxscores[x.tags['src']] >= 0.85]
    trgmatches = [x for x in completematches if x.geneID in genematch_codes]
    annotation, seensources = complementAnnotation(trgmatches, annotation, runde, seensrc=seensources, uniqsrc=True)
    annotation.sort()
    # round 4: candidate match with high coverage of source, each source id can be inserted multiple times
    runde = 4
    completematches = [x for x in allmatches if x.tags['intstop'] == '0' and x.tags['hasstart'] == '1'
                       and x.tags['hasstop'] == '1' and x.score / maxscores[x.tags['src']] >= 0.95]
    trgmatches = [x for x in completematches if x.geneID in genematch_codes]
    annotation, seensources = complementAnnotation(trgmatches, annotation, runde, seensrc=seensources, uniqsrc=False)
    annotation.sort()

    # round 5: here we allow for incomplete hc matches for which no source has been yet inserted (unique).
    runde = 5
    trgmatches = [x for x in allmatches if x.tags['intstop'] == '0' and x.geneID in hcmatch_codes and
                  x.score / maxscores[x.tags['src']] >= 0.85]
    annotation, seensources = complementAnnotation(trgmatches, annotation, runde, seensrc=seensources, uniqsrc=True)
    annotation.sort()

    ####################################################################################################
    # lastly we target complete matches with contiguous ORFs for plastid- and transposon-related models
    # --------------------------------------------------------------------------------------------------
    # round 6: plastid genes, each source id complete, contiguous and maximally one insertion
    runde = 6
    # for these two we ask for a score at least 90% of the score of the self-aligned source
    completematches = [x for x in allmatches if x.tags['intstop'] == '0' and x.tags['hasstart'] == '1'
                       and x.tags['hasstop'] == '1' and x.score / maxscores[x.tags['src']] >= 0.9]
    trgmatches = [x for x in completematches if x.geneID in plastid_codes]  # limit to plastids
    annotation, seensources = complementAnnotation(trgmatches, annotation, runde, seensrc=seensources, uniqsrc=True)
    annotation.sort()
    # round 7: transposon-related genes, each source id complete, contiguous and maximally one insertion
    #          latter important to minimize high frequent transfer of high-copy transposons
    runde = 7
    trgmatches = [x for x in completematches if x.geneID in transposon_codes]  # limit to transposons
    annotation, seensources = complementAnnotation(trgmatches, annotation, runde, seensrc=seensources, uniqsrc=True)
    annotation.sort()

    with open("%s.anno.gff" % (genotype), 'w') as out:
        for anno in annotation:
            out.write("%s\n" % (anno.toGFF3(with_exon=True)))
    return



# -----------------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    argParser = argparse.ArgumentParser()
    argParser.add_argument("-r", "--runfile", type=str,
                           help="""absolute path to runlist file defining genotype name to genomepath mapping
                           file format: genotype name <TAB> abs. path to genome fasta sequence""")
    argParser.add_argument("-o", "--projroot", type=str,
                           help="root directory of projection project, previously created from preprocessing run")
    argParser.add_argument("--pfamfile", type=str,
                           help="""file listing transcript ids of source models annotated with a PFAM domain.
                           file format: one transcript id per line""")
    argParser.add_argument("--plastidfile", type=str,
                           help="""file listing transcript ids of source models annotated as plastid-related genes.
                           file format: one transcript id per line""")
    argParser.add_argument("--transposonfile", type=str,
                           help="""file listing gene ids of source models annotated as transposon-related genes.
                           file format: one transcript id per line""")
    argParser.add_argument("--scorefile", type=str,
                           help="""file containing self-scores of transcript source ids.
                           file format: transcript id <TAB> blosum62 score of self-alignment""")


    args = argParser.parse_args()

    runlist = []

    if not os.path.isdir(args.projroot):
        sys.stderr.write("Project root directory %s not found!\n" %(args.projroot))
        sys.exit(2)

    # read accessory information for scoring and filters
    maxscores = {}
    with open(args.scorefile) as fhd:
        reader = csv.reader(fhd, delimiter='\t')
        for tid, sc in reader:
            maxscores[tid] = float(sc)

    plastidcodes = set()
    with open(args.plastidfile) as fhd:
        for line in fhd:
            plastidcodes.add(line.strip())

    transposoncodes = set()
    with open(args.transposonfile) as fhd:
        for line in fhd:
            plastidcodes.add(line.strip())

    pfamcodes = set()
    with open(args.pfamfile) as fhd:
        for line in fhd:
            plastidcodes.add(line.strip())

    # get the whole runlist
    with open(args.runfile) as fhd:
        reader = csv.reader(fhd, delimiter='\t')
        for haplogeno, genomefn in reader:
            runlist.append((haplogeno, genomefn))

    # do the jobs
    for gt, genomefile in runlist:
        annoOneGenotype(genotype=gt, projroot=args.projroot, pfam_codes=pfamcodes, plastid_codes=plastidcodes,
                        transposon_codes=transposoncodes, maxscores=maxscores)



