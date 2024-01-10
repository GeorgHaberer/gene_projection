###########################################################
#                                                         #
#    Data container for one uninterrupted genomic         #
#    sequence element, generally an exon.                 #
#                                                         #
###########################################################

from __future__ import annotations
from typing import Union

from functools import total_ordering

@total_ordering
class Coords:
    """Defines a genetic element in a genome, ie contig id, start and stop. Coordinates
    are inclusive, ie closed interval. Strand is encoded by 0:forward and 1:reverse,
    score is optional."""

    def __init__(self,ctg : str = "", beg : int = -1, end : int = -1,
                 strand : int = -1, score : Union[None, float, int] = None):

        self.contigID = ctg
        self.beg = beg
        self.end = end
        self.strand = strand
        self.score = score


    def size(self) -> int:
        """:return genomic span size [beg,end]"""
        return abs(self.beg-self.end)+1


    def minimum(self) -> int:
        """:return minimum coordinate of beg,end"""
        return min(self.beg,self.end)


    def maximum(self) -> int:
        """:return maximum coordinate of beg,end"""
        return max(self.beg,self.end)


    def overlaps(self, other : Coords) -> bool:
        """:return BOOL for overlap between two coords, one single base sufficient"""
        if self.contigID != other.contigID:
            return False
        if self.minimum() <= other.maximum() and self.maximum() >= other.minimum():
            return True
        return False


    def is_contained_in(self, other : Coords):
        """:return true if self is completely covered by other"""
        if self.contigID != other.contigID:
            return False
        if self.minimum() >= other.minimum() and self.maximum() <= other.maximum():
            return True
        return False


    def __eq__(self, other : Coords) -> bool:
        """two elements are equal if they span the identical genomic region.
        Current implementation doesn't care about score or strand!"""
        return (self.contigID, self.minimum(), self.maximum()) == (other.contigID, other.minimum(),
                                                                   other.maximum())


    def __lt__(self, other : Coords) -> bool:
        """element is lesser if alphanumerical sort of (1) contig ID, or (2) its genomic
        minimum coordinate or (3) its genomic maximum coordinate is less"""
        return (self.contigID, self.minimum(), self.maximum()) < (other.contigID, other.minimum(),
                                                                  other.maximum())


    def __str__(self) -> str:
        """simple print function, mainly diagnostic purposes"""
        return "%s\t%s\t%s\t%s" %(self.contigID, self.strand, self.beg, self.end)


