from __future__ import annotations

import os
import argparse
from Bio import SeqIO
from Bio.Align import substitution_matrices
from Bio.Align import PairwiseAligner
import multiprocessing


# TODO: currently gap open/ext. penalties and blosum62 matrix fixed, make it variable?
def alignBatch(subbatch : list) -> dict:
    """global needleman wunsch self-alignment of subbatch protein sequences.
    :return dict id -> score"""
    matrix = substitution_matrices.load("BLOSUM62")
    gap_open = -10.0
    gap_extend = -0.5
    aligner = PairwiseAligner(mode='global',open_gap_score=gap_open,extend_gap_score=gap_extend,\
                              target_end_gap_score = 0.0,query_end_gap_score = 0.0,\
                              substitution_matrix = matrix)
    scores = {}
    for seqid, seq in subbatch:
        alns = aligner.align(seq, seq)
        topali = alns[0]
        scores[seqid] = topali.score

    return scores

def parallel_alignment(batches : list, nproc : int = 1, chunksize : int = 5000) -> dict[dict] :
    """<nproc> workers to dispatch alignments in chunks of size <chunksize>
    :return results as dict <id> : <score>"""
    data = [batches[i:i+chunksize] for i in range(0,len(batches), chunksize)]
    results = multiprocessing.Pool(processes=nproc).map(alignBatch, data)
    return results

def runSelfScoring(protfile : str, outfile : str, ncpus : int = 1, chunksize : int = 5000):
    """runs self alignments for entries in protfile and writes out
    <id><TAB><score>"""
    protbatches = []
    for rec in SeqIO.parse(protfile, "fasta"):
        protbatches.append((rec.id, str(rec.seq).upper()))

    results = parallel_alignment(protbatches, nproc=ncpus, chunksize=chunksize)
    with open(os.path.join(outfile), 'w') as out:
        for res in results:
            for gid, sc in res.items():
                out.write("%s\t%.1f\n" %(gid, sc))



# -----------------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    argParser = argparse.ArgumentParser()
    argParser.add_argument("-p", "--protfile", type=str,
                           help="""absolute path to protein fasta file of representative transcripts""")
    argParser.add_argument("-o", "--outfile", type=str,
                           help="abs. path to self-/max-scores of representatives")
    argParser.add_argument("-n", "--ncpus", type=int, default=1,
                           help="""number of parallel cpus for self alignments""")
    argParser.add_argument("--chunksize", type=int, default=5000,
                           help="""number of sequences per parallel process""")


    args = argParser.parse_args()
    runSelfScoring(args.protfile, args.outfile, ncpus=args.ncpus, chunksize=args.chunksize)

