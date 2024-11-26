import re 
import os 
import csv
import sys 
import bisect

sys.path.append("/home/pgsb/georg.haberer/src/pyprojects/bioutils")
from gffparser import parseTaggedGFF
from genedata import GeneLocus, mRNA, Cds


genotypes = ["A17", "A40", "A42", "GRA2256_1", "FB19_011_3", "FB19_001_1", "FB19_028_3", "FB20_005_1", "FB20_029_7", "PI365428"]
haps = [("hap1", "hap2", "hap3", "hap4", "unanchor"), ("hap1", "hap2", "hap3", "hap4", "unanchor"),
        ("hap1", "hap2", "hap3", "hap4", "unanchor"), ("hap1", "hap2", "hap3", "hap4", "unanchor"),
        ("hap1", "hap2", "unanchor"), ("hap1", "hap2", "hap3", "hap4", "unanchor"),
        ("hap1", "hap2", "hap3", "hap4", "unanchor"), ("hap1", "hap2", "unanchor"),
        ("hap1", "hap2", "unanchor"), ("hap1", "hap2", "unanchor")]


EVIROOT = "/lustre/groups/pgsb/projects/barley/hbulbosum/release_v1"
BULBHOME = "/lustre/groups/pgsb/workspaces/georg.haberer/bulbosum/revision"
PROJROOT = os.path.join(BULBHOME, "projections")
OUTDIR = os.path.join(BULBHOME, "padded")

cutoff = float(sys.argv[1])
version = sys.argv[2]

os.chdir(BULBHOME)
os.chdir("datasets")
domains = {}
with open("hbulb.domain.codes.txt") as fhd:
    reader = csv.reader(fhd, delimiter=' ')
    for tid,val in reader:
        domains[tid] = float(val)

transposons = set()
with open("hbulb.transposon.codes") as fhd:
    reader = csv.reader(fhd, delimiter=' ')
    for line in fhd:
        transposons.add(line.strip())

plastids = set()
with open("hbulb.plastid.codes") as fhd:
    reader = csv.reader(fhd, delimiter=' ')
    for line in fhd:
        plastids.add(line.strip())


print(len(domains))
print(len(transposons))

for gt, haplotypes in zip(genotypes, haps):
   for h in haplotypes:
        trgline = "%s_%s" %(gt, h)
        highgff = os.path.join(EVIROOT, trgline, "%s.pgsb.r1.Sep2024.high.gff3" %(trgline))
        projgff = os.path.join(PROJROOT, trgline, "%s.anno.v1.gff" %(trgline))
        genes, genome = parseTaggedGFF(highgff)
        genome.sort()
        print(trgline, len(genome), sep='\t')
        topid = 0
        for gene in genome:
            mobj = re.search("(?P<lid>\d+$)", gene.geneID)
            if mobj:
                aktid = int(mobj.group('lid'))
                topid = max(topid, aktid)
        
        gidcounter = topid + 10
        addgenes, addgenome = parseTaggedGFF(projgff)
        addgenome.sort()
        for gene in addgenome:
            src = gene.tags['src']
            rnd = int(gene.tags['round'])
            if rnd == 5: continue
            if src in transposons:
                continue
            pfam = 100.0
            try:
                pfam = domains[src]
            except KeyError:
                pass
            if pfam > cutoff:
                continue
            i = bisect.bisect(genome, gene)
            cmin = gene.codingMin()
            cmax = gene.codingMax()
            for j in range(max(0,i-5), min(len(genome), i+6)):
                genej = genome[j]
                if genej.contigID != gene.contigID:
                    continue
                cminj = genej.codingMin()
                cmaxj = genej.codingMax()
                if cmax >= cminj and cmin <= cmaxj:
                    break 
            else:
                # make new geneid
                # H.BULBOSUM.FB19_001_1.r1.1H_1G00000020
                # H.BULBOSUM.FB19_001_1.r1.ctg1022G00001090  : contig_corrected_v1_1022
                # H.BULBOSUM.FB19_011_3.r1.unitig_1038G00000690 : unitig_corrected_v1_1038          
                if gt.endswith("unanchor"):
                    abbctg = re.sub("_corrected_v1", "", gene.contigID)
                    abbctg = re.sub("contig_", "ctg", abbctg)
                else:
                    abbctg = re.sub("chr", "", gene.contigID)
                    newgid = "H.BULBOSUM.{0}.r1.{1}G{2:08}".format(gt, abbctg, gidcounter)
                newgene = GeneLocus(gid=newgid, ctg=gene.contigID, src="PGSB")
                newgene.tags['method'] = "projected"
                thistag = "protein_coding"
                if src in plastids:
                    thistag = "plastid-related"
                elif src in transposons:
                    thistag = "transposon-related"
                newgene.tags['tag'] = thistag
                newgene.tags['confidence'] = "high"
                newgene.tags['Name'] = newgid
                tidcounter = 1
                for oldtid, oldt in gene.transcripts.items():
                    newtid = newgid + '.%i' %(tidcounter)
                    newmrna = mRNA(mid=newtid, ctg=gene.contigID, strand=gene.strand(), coords=[], cds=Cds())
                    newmrna.tags['tag'] = thistag 
                    newmrna.tags['Name'] = newtid
                    newmrna.tags['confidence'] = "high"
                    for ex in oldt:
                        newmrna.append(ex)
                    for ex in oldt.coding:
                        newmrna.coding.append(ex)
                    tidcounter += 1
                    newgene.addmRNA(newmrna)
                genome.insert(i, newgene)
                gidcounter += 10
        
        paddedgff = os.path.join(OUTDIR, "%s.pgsb.v%s.padded.gff" %(trgline, version))
        genome.sort()
        print(trgline, len(genome), sep='\t')
        with open(paddedgff, 'w') as out:
            for gene in genome:
                out.write("%s\n" %(gene.toGFF3(with_exon=True)))
