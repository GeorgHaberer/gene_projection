###########################################################
#                                                         #
#           Data container for PAF entry                  #
#                                                         #
###########################################################

from __future__ import annotations
import re

from functools import total_ordering

@total_ordering
class CSlongPafTranscript:
    """container storing minimap mapping of a query id to a target id, also recording strand, cigar
    mapping encoding in long format, query and target minimum and maximum coordinates, and size of
    query and target"""

    def __init__(self, qid : str = "", qsz : int = 0, trg : str = "", tsz : int = 0, qmin : int = 0,
                 qmax : int = 0, tmin : int = 0, tmax : int = 0, ori : int = -1, cigar : str = "") -> None :

        self.query = qid
        self.target = trg
        self.strand = ori
        self.cslong = cigar
        self.querymin = qmin
        self.querymax = qmax
        self.targetmin = tmin
        self.targetmax = tmax
        self.querysize = qsz
        self.targetsize = tsz


    def __eq__(self, other : type[CSlongPafTranscript]) -> bool :

        return (self.target, self.targetmin, self.targetmax) == (other.target, other.targetmin, other.targetmax)


    def __le__(self, other : type[CSlongPafTranscript]) -> bool :

        return (self.target, self.targetmin, self.targetmax) <= (other.target, other.targetmin, other.targetmax)




def extractPafTranscript(row : list) -> CSlongPafTranscript:

    paf = CSlongPafTranscript()
    paf.query = row[0]
    paf.querysize = int(row[1])
    paf.querymin = int(row[2])
    paf.querymax = int(row[3])
    paf.strand = 0 if row[4] == '+' else 1
    paf.target = row[5]
    paf.targetsize = int(row[6])
    paf.targetmin = int(row[7])
    paf.targetmax = int(row[8])
    paf.cslong = re.sub("cs:Z:", "", row[-1])
    return paf

