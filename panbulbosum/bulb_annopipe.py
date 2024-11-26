import re
import os
import sys
import csv
import gzip
import subprocess
import multiprocessing
from Bio import SeqIO

from chunk_genome_bulb import chunkit
from miniprotgff_to_gff import transformMiniProtGFF
from psl2gff import psl2gff
from transformPAF_new import paf_transformer_new

from biosequences import GenomicTemplate, FastaSequence

HBULB_ROOT = "/lustre/groups/pgsb/workspaces/georg.haberer/bulbosum/revision"
OFFSETFILE = os.path.join(HBULB_ROOT, "datasets", "hbulb.offsets.txt")
TFASTA = os.path.join(HBULB_ROOT, "datasets", "hbulb.transcripts.fa")
PFASTA = os.path.join(HBULB_ROOT, "datasets", "hbulbosum.evidence.high.cdhit")
CLEAN_TFASTA = os.path.join(HBULB_ROOT, "datasets", "hbulb.clean.transcripts.fa")
CLEAN_PFASTA = os.path.join(HBULB_ROOT, "datasets", "hbulb.clean.proteins.fa")

MMAP2 = "/home/pgsb/georg.haberer/bin/minimap2/./minimap2"
MINIPROT = "/home/pgsb/georg.haberer/bin/miniprot/./miniprot"
BLAT = "/home/pgsb/georg.haberer/bin/./blat"


def runblat(batch):

    bid, batchproteins, db, gt, pslout = batch
    outdir = os.path.dirname(pslout)
    oocfile = os.path.join(outdir, "%s.11.ooc" %(gt))
    fastapath = os.path.join(outdir, "seq%i.fa" %(bid))
    p = SeqIO.write(batchproteins, fastapath, "fasta")
    CMD = [BLAT, db, fastapath, "-q=dna", "-t=dna", "-ooc=%s" %(oocfile), "-noTrimA", "-fine", "-maxIntron=70000", "-out=psl", pslout]
    subprocess.run(CMD)
    return

def runbatch(genotype, genomefile):

    # first make output directory
    OUTDIR = os.path.join(HBULB_ROOT, "projections", genotype)
    try:
        os.mkdir(OUTDIR)
    except OSError:
        pass

    # chunk genome file, but make only one chrUn if unanchored
    ctgnames = chunkit(genomefile, OUTDIR)

    # preprocess for fast blat queries
    os.chdir(OUTDIR)
    ctgfiles = ["%s.fa" %(ctg) for ctg in ctgnames]
    os.system("cat %s > genome.fasta" %(' '.join(ctgfiles)))
    oocfile = os.path.join(OUTDIR, "%s.11.ooc" %(genotype))
    subprocess.run([BLAT, "genome.fasta", "/dev/null", "/dev/null", "-makeOoc=%s" %(oocfile),  "-repMatch=1024"])

    # now for each chromosome
    tmp_pslfiles = []
    for ctg in ctgnames:

        # first run minimap2 on transcripts
        ctgfile = "%s.fa" %(ctg)
        outfile = "%s.unsorted.paf" %(ctg)
        ctgrun = [MMAP2, '-x', 'splice:hq', '-uf', '-t', '16', '--cs=long', '-o', outfile, ctgfile, TFASTA]
        subprocess.run(ctgrun)
        outsorted = "%s.paf" %(ctg)
        os.system("sort -n -k 12 %s > %s" %(os.path.join(OUTDIR, outfile), os.path.join(OUTDIR, outsorted)))
        subprocess.run(['gzip', outsorted])
        os.system("rm %s" %(os.path.join(OUTDIR, outfile)))
        gffout = os.path.join(OUTDIR, "%s.raw.gff" %(ctg))
        gzipped_outfile = os.path.join(OUTDIR, "%s.gz" %(outsorted))
        paf_transformer_new(gzipped_outfile, ctgfile, PFASTA, gffout, OFFSETFILE)

        subprocess.run(['gzip', gffout])
        # second run miniprot on proteins
        outfile = "%s.mp.original.gff" %(ctg)
        ctgrun = [MINIPROT, '-t', '16', '--gff', '-G', '70k', '--outs', '0.7', '--outc', '0.7', ctgfile, CLEAN_PFASTA]
        # somehow with capturing stdout from multithreaded miniprot slows it down for subprocess.run
        # go with alternative solution os.system (:
        os.system("%s > %s" %(' '.join(ctgrun), outfile))
        subprocess.run(['gzip', outfile])
        outfile += '.gz'
        gffout = os.path.join(OUTDIR, "%s.miniprot.gff" %(ctg))
        uniqprefix = "%s.%s.mp" %(genotype, ctg)
        transformMiniProtGFF(outfile, gffout, ctgfile, CLEAN_PFASTA, uniqprefix)
        subprocess.run(['gzip', gffout])

        # third parallize jobs and run blat
        transcripts = []
        for rec in SeqIO.parse(CLEAN_TFASTA, "fasta"):
            transcripts.append(rec)
        protdict = {}
        for rec in SeqIO.parse(CLEAN_PFASTA, "fasta"):
            protdict[rec.id] = rec

        # # we make 20 batches
        bsize = int(len(transcripts)/20)
        batches = [(i, transcripts[i:i + bsize], os.path.join(OUTDIR, ctgfile), genotype,
                    os.path.join(OUTDIR, "%s.batch%i.psltmp" %(ctg, i))) for i in range(0, len(transcripts), bsize)]
        results = multiprocessing.Pool(processes=20).map(runblat, batches)

        scaffolds = {}
        if ctgfile.endswith(".gz"):
            fhd = gzip.open(ctgfile, mode='rt')
        else:
            fhd = open(ctgfile)
        for rec in SeqIO.parse(fhd, "fasta"):
            scaffolds[rec.id] = GenomicTemplate(ID=rec.id, seq=str(rec.seq).upper())

        srcoffsets = {}
        with open(OFFSETFILE) as fhd:
            reader = csv.reader(fhd, delimiter='\t')
            for tid, beg, end in reader:
                srcoffsets[tid] = (int(beg), int(end))

        pslgenome = []
        uniqprefix = "%s.%s.psl" %(genotype, ctg)
        tmp_pslfiles = []
        uidc = 1
        for _, protbatch, _, _, psltmp in batches:
            tmpgenome, uidc =  psl2gff(psltmp, scaffolds, srcoffsets, protdict, prefix=uniqprefix, uidcounter=uidc)
            tmp_pslfiles.append(psltmp)
            pslgenome += tmpgenome

        pslgenome.sort()
        psloutfile = "%s.psl.gff" %(ctg)
        with open(psloutfile, 'w') as out:
            for gene in pslgenome:
                if gene.score < 100: continue
                out.write("%s\n" %(gene.toGFF3(with_exon=True)))

        subprocess.run(["gzip", psloutfile])
        os.system("cat %s > %s.original.psl" %(' '.join(tmp_pslfiles), ctg))
        subprocess.run(["gzip", "%s.original.psl" %(ctg)])
        for tmppsl in tmp_pslfiles:
            subprocess.run(["rm", tmppsl])


    os.chdir(OUTDIR)
    subprocess.run(["rm", "genome.fasta"])



if __name__ == '__main__':

    gt = sys.argv[1]
    genomefile = sys.argv[2]
    runbatch(genotype=gt, genomefile=genomefile)

