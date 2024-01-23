#############################################################
#                                                           #
#    Datamodel for piece-wise genetic element consisting    #
#    of several coord tuples, eg an exon-intron structure   #
#                                                           #
#############################################################

from __future__ import annotations
from typing import Union

import sys
import itertools
from operator import attrgetter

from coords import Coords
from biosequences import DnaSequence, GenomicTemplate, FastaSequence

# -------------------------------------------------------
class CoordChain:

    def __init__(self, ctg : str = "", strand : int = -1, coords : list = []):

        self.contigID = ctg
        self.strand = strand
        self._data = []

        if len(coords):
            # right format?
            ins = set([isinstance(ct, Coords) for ct in coords])
            assert len(ins) == 1
            # consistent strands/orientations?
            oris = set([ct.strand for ct in coords])
            if self.strand != -1:
                oris.add(self.strand)
            assert len(oris) == 1
            self.strand = coords[0].strand  # everything consistent, reset/set newly strand
            # are we all on the same pseudomolecule?
            ctgs = set([ct.contigID for ct in coords])
            if self.contigID != "":
                ctgs.add(self.contigID)
            assert len(ctgs) == 1
            # make sure contig is set if coordinates have been provided
            self.contigID = coords[0].contigID
            for ct in coords:
                self._data.append(ct)


    def minimum(self) -> int:
        """:return genomic minimum of coordinates"""
        return min(ct.minimum() for ct in self._data)


    def maximum(self) -> int:
        """:return genomic maximum of coordinates"""
        return max(ct.maximum() for ct in self._data)


    def size(self) -> int:
        """:return spliced size of chain"""
        return sum([ct.size() for ct in self._data])

    def overlaps(self, other : Union[Coords, CoordChain]) -> bool:
        """:return bool if two elements overlap"""
        if self.contigID != other.contigID:
            return False
        if self.minimum() <= other.maximum() and self.maximum() >= other.minimum():
            return True
        return False

    def getSequence(self, scaffold : Union[FastaSequence, GenomicTemplate]) -> DnaSequence:
        """:return (dna)sequence of coordchain for its contig <scaffold>"""
        self.sort()
        r = []
        for el in self:
            x = min(el.beg, el.end)
            y = max(el.beg, el.end)
            r.append(scaffold.sequence[x - 1:y])
        s = DnaSequence(seq=''.join(r))
        if self.strand == 1:
            return s.getReverseComplement()
        else:
            return s


    def trimBySizer(self, direction : int = 5, cutsize : int = 0) -> CoordChain:
        # check vanity of arguments
        if direction not in [3,5]:
            sys.stderr.write("I do not know from what direction you want me to trim!\nArgument <direction> was : %s\n" %(direction))
            raise AttributeError
        if cutsize <= 0:
            sys.stderr.write("Cannot cut %s bp from coordinates!\n" %(cutsize))
            raise AttributeError
        coordstretch = []
        x = cutsize
        for ct in self._data:
            exon = list(range(ct.minimum(),ct.maximum()+1,1))
            coordstretch += exon
        coordstretch.sort()
        if direction == 5:
            # trim from 5'start
            if self.strand == 0:
                while x:
                    coordstretch.pop(0)
                    x -= 1
            elif self.strand == 1:
                while x:
                    coordstretch.pop()
                    x -= 1
        elif direction == 3:
            # trim from 3'-end
            if self.strand == 0:
                while x:
                    coordstretch.pop()
                    x -= 1
            elif self.strand == 1:
                while x:
                    coordstretch.pop(0)
                    x -= 1
        # now construct new coordchain and return to caller
        newchain = CoordChain(ctg=self.contigID,strand=self.strand)
        for k, g in itertools.groupby(enumerate(coordstretch), lambda x: x[0]-x[1]):
            exon = [x[1] for x in g]
            if self.strand == 0:
                ct = Coords(ctg=self.contigID,strand=self.strand,beg=min(exon),end=max(exon))
                newchain.append(ct)
            elif self.strand == 1:
                ct = Coords(ctg=self.contigID,strand=self.strand,beg=max(exon),end=min(exon))
                newchain.append(ct)
        return newchain

    # -----------------------------------------------------------
    # overloading built-ins to simulate features of a python list

    def sort(self) -> None:
        # simulates list sorting
        self._data.sort(key=attrgetter('contigID', 'beg'))


    def append(self, item : Coords, EXPLICIT_ALLOW_OVL : bool = False) -> None:
        """appends new genetic element item <Coords> to chain, simulates list like
        behaviour, ie appends at the end. Make sure that chain is sorted before
        any operation requiring correct order of elements is performed!"""
        assert isinstance(item, Coords)
        # if we add coord item to previously undefined instance coordchain
        # for example instantiated as empty: set strand and contigID
        if self.strand == -1:
            self.strand = item.strand
        if self.contigID == "":
            self.contigID = item.contigID
        for ct in self:
            assert ct.strand == item.strand, "unequal strand at %s,%s, for %s" % (self.contigID, self.strand, item)
            assert ct.contigID == item.contigID, "unequal contig at %s, for %s" % (self.contigID, item)
            if ct.overlaps(item) and not EXPLICIT_ALLOW_OVL:
                raise ValueError
        self._data.append(item)


    def pop(self, index : int = -1) -> Coords:
        """:return Coords element from position index, default to last"""
        ct = self._data[index]
        del self._data[index]
        return ct


    def __len__(self) -> int:
        """:return number of elements, ie exons"""
        return len(self._data)


    def __iter__(self):
        """:return an iterator over the elements/coords"""
        return iter(self._data)


    def __setitem__(self, index : int, item : Coords) -> None:
        """set element at position index"""
        assert isinstance(item, Coords)
        for ct in self:
            if ct.strand != item.strand or ct.contigID != item.contigID:
                raise ValueError
            if ct.overlaps(item):
                raise ValueError

        self._data[index] = item
        if self.strand == -1:
            self.strand = item.strand


    def __getitem__(self, i : int) -> Coords:
        """get element at position i"""
        return self._data[i]


    def __str__(self) -> str:
        """:return simple string representation of chain, mainly for diagnostic purposes"""
        tmp = [self.contigID]
        tmp.append("%i" % (self.strand))
        for ct in self:
            tmp.append("%s" % (ct))
        return '\t'.join(tmp)
