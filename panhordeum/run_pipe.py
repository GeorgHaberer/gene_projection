###########################################################
#                                                         #
#          Main script to run mapping pipeline            #
#                                                         #
###########################################################

import os
import sys
import csv
import argparse
import subprocess

from chunk_genome import chunkit
from transformPAF import paf_transformer


def runbatch(genotype : str, genomefile : str, outdir : str, tfasta: str, pfasta: str, offsetfile: str,
             mmap2 : str):

    # first make output directory
    try:
        os.mkdir(outdir)
    except OSError:
        # TO DO!!!: add option force if overwrite possible, currently overwrites
        # existing results in outdir
        pass

    # chunk genome file, see special modus for unanchored contigs
    ctgnames = chunkit(genomefile, outdir, genotype)

    os.chdir(outdir)  # change to selected output directory
    # now we run our mapping and postprocessing for each chromosome:
    for ctg in ctgnames:
        # run minimap
        ctgfile = "%s.%s.fa" %(genotype, ctg)
        outfile = "%s.%s.unsorted.paf" %(genotype, ctg)
        ctgrun = [mmap2, '-x', 'splice:hq', '-uf', '-t', '8', '--cs=long', '-o', outfile, ctgfile, tfasta]
        subprocess.run(ctgrun)

        # sort by increasing genomic coordinates
        outsorted = "%s.%s.paf" %(genotype, ctg)
        os.system("sort -n -k 8 %s > %s" %(os.path.join(outdir, outfile), os.path.join(outdir, outsorted)))
        subprocess.run(['gzip', outsorted])
        os.system("rm %s" %(os.path.join(outdir, outfile)))

        # transform paf file to gff
        gffout = os.path.join(outdir, "%s.%s.raw.gff" %(genotype, ctg))
        gzipped_outfile = os.path.join(outdir, "%s.gz" %(outsorted))
        paf_transformer(gzipped_outfile, ctgfile, pfasta, gffout, offsetfile)

        # clean up and compress raw gff
        os.system("rm %s" %(ctgfile))
        subprocess.run(['gzip', gffout])


# -----------------------------------------------------------------------------
if __name__ == '__main__':

    argParser = argparse.ArgumentParser()
    argParser.add_argument("-r", "--runfile", type=str,
                           help="""absolute path to runlist file defining genotype name to genomepath mapping
                           file format: genotype name <TAB> abs. path to genome fasta sequence""")
    argParser.add_argument("-o", "--outdir", type=str,
                           help="output directory, created if not present")
    argParser.add_argument("--offsetfile", type=str,
                           help="""file listing coding start in transcripts for each source transcript
                           file format: transcript id <TAB> coding start <TAB> coding length""")
    argParser.add_argument("-t", "--tfasta", type=str,
                           help="abs. path to transcript fasta file to be mapped")
    argParser.add_argument("-p", "--pfasta", type=str,
                           help="abs. path to protein fasta of transcript sources")
    argParser.add_argument("--mmap", type=str,
                           help="path to minimap2 command, ie. your installation path.")
    argParser.add_argument("--batchstart", type=int, default=0,
                           help="start runlist from <batchstart>. Unset (default:0) to start from first entry")
    argParser.add_argument("--batchend", type=int, default=-1,
                           help="run up to <batchend> in runlist. Unset (default:-1) to include all up to last entry. "
                                "Defined as closed interval, ie. inclusive <batchend> entry!")


    args = argParser.parse_args()

    # some sanity checks
    if args.batchstart >= 0 and args.batchend >= 0:
        if args.batchend < args.batchstart:
            sys.stderr.write("Batch end lower than start, please check!\n")
            sys.exit(2)
    runlist = []

    # get the whole runlist
    with open(args.runfile) as fhd:
        reader = csv.reader(fhd, delimiter='\t')
        for haplogeno, genomefn in reader:
            runlist.append((haplogeno, genomefn))

    # do the job within given entry interval
    for gt, genomefile in runlist[args.batchstart:args.batchend+1]:
        runbatch(genotype=gt, genomefile=genomefile, outdir=args.outdir, tfasta=args.tfasta,
                 pfasta=args.pfasta, mmap2=args.mmap, offsetfile=args.offsetfile)
















