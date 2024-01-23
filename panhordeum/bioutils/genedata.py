##################################################################
#                                                                #
#    Data container to model some attributes and functionality   #
#    of coding, transcript and gene elements                     #
#                                                                #
##################################################################

from __future__ import annotations
from typing import Union

import sys
import itertools
from functools import total_ordering

from bioutils.coordchain import CoordChain
from bioutils.coords import Coords


@total_ordering
class Cds(CoordChain):

    """Data container for a CDS/coding sequence object"""
    def __init__(self, pid : str = "", cid : str = "", ctg : str = "",
                 strand : int = -1, coords : list = []):

        CoordChain.__init__(self, ctg=ctg, strand=strand, coords=coords)
        self.proteinID = pid  # legacy, not really used in pipeline
        self.cdsID = cid

    # overload for sort by genomic position
    def __eq__(self, other : Union[Cds, CoordChain]) -> bool:

        return ((self.contigID, self.minimum(), self.maximum()) ==
                (other.contigID, other.minimum(), other.maximum()))


    def __lt__(self, other : Union[Cds, CoordChain]) -> bool:

        return ((self.contigID, self.minimum(), self.maximum()) <
                (other.contigID, other.minimum(), other.maximum()))


    def getPhaseList(self) -> list:
        """:return a list of frames/phases for each exon. List order corresponds to order
        from min to max genomic positions of exons"""
        frames = []  # stores frames/phases of exons, ordered by min/max order of exons
        self.sort()  # step to ensure min/max order of exons
        frame = 0    # starting frame

        if self.strand == 1:
            for el in reversed(self):
                frames.append(frame)
                remainder = (el.size() - frame) % 3
                if remainder == 0:
                    frame = 0
                elif remainder == 1:
                    frame = 2
                elif remainder == 2:
                    frame = 1
            frames.reverse()
        else:
            for el in self:
                frames.append(frame)
                remainder = (el.size() - frame) % 3
                if remainder == 0:
                    frame = 0
                elif remainder == 1:
                    frame = 2
                elif remainder == 2:
                    frame = 1

        return frames


    def toGFF3(self, parent : str, source : str = '.', CDSscores : bool = True) -> str:
        """transform instance into gff formatted string
        :return formatted text"""
        frames = self.getPhaseList()
        ori = '+' if self.strand == 0 else '-'
        gff = []  # storage for each gff line ~ single element

        for el,frame in zip(self,frames):
            r = []
            r.append("%s" % (self.contigID))
            r.append("%s\tCDS" % (source))
            r.append("%s" % (el.minimum()))
            r.append("%s" % (el.maximum()))

            # we accept numeric or no score type
            if isinstance(el.score, int):
                r.append("%i\t%s\t%s\tParent=%s" % (el.score, ori, frame, parent))
            elif isinstance(el.score, float):
                r.append("%.3f\t%s\t%s\tParent=%s" % (el.score,ori,frame,parent))
            else:
                r.append(".\t%s\t%s\tParent=%s" % (ori,frame,parent))

            gff.append('\t'.join(r))

        return "%s" %('\n'.join(gff))


@total_ordering
class mRNA(CoordChain):

    """Data container for a transcript/mRNA object"""
    def __init__(self, mid : str = "", ctg : str = "", strand : int = -1, coords : list = [],
                 score : Union[None, float, int] = None, cds : Cds = Cds()):

        CoordChain.__init__(self, ctg=ctg, strand=strand, coords=coords)
        self.tid = mid  # element ~ transcript id
        self.score = score
        self.strand = strand
        assert isinstance(cds, Cds)
        self.coding = cds  # coding part of transcript
        if self.strand != -1 and self.coding.strand != -1:
            assert self.strand == self.coding.strand
        self.tags = {}
        self.descriptionLine = ""  # legacy code, better use <tags> attribute to store description

    # overload for sort by genomic position
    def __eq__(self, other : Union[mRNA, CoordChain]) -> bool:

        return ((self.contigID, self.minimum(), self.maximum()) ==
                (other.contigID, other.minimum(), other.maximum()))


    def __lt__(self, other : Union[mRNA, CoordChain]) -> bool:

        return ((self.contigID, self.minimum(), self.maximum()) <
                (other.contigID, other.minimum(), other.maximum()))

    def getUTRs(self) -> tuple[CoordChain, CoordChain]:
        """:return 5'- and 3'-UTRs as coordchains"""
        coordstretch = []
        for ct in self._data:
            exon = list(range(ct.minimum(), ct.maximum()+1,1))
            coordstretch += exon
        coordstretch.sort()
        cmin, cmax = self.coding.minimum(), self.coding.maximum()
        utrmin = []
        utrmax = []
        for pos in coordstretch:
            if pos < cmin:
                utrmin.append(pos)
            else:
                break
        for pos in reversed(coordstretch):
            if pos > cmax:
                utrmax.insert(0, pos)
            else:
                break

        if self.strand == 0:
            utr5 = self._makeChain(utrmin)
            utr3 = self._makeChain(utrmax)
        else:
            utr5 = self._makeChain(utrmax)
            utr3 = self._makeChain(utrmin)
        return utr5, utr3


    def _makeChain(self, stretch : list) -> Union[None, CoordChain]:
        """re-constructs from a list of ordered genomic bases/coordinates the corresponding
        coordchain.
        :return coordchain or None if no coordinates provided"""
        if len(stretch):
            newchain = CoordChain(ctg=self.contigID, strand=self.strand)
            for k, g in itertools.groupby(enumerate(stretch), lambda x: x[0]-x[1]):
                exon = [x[1] for x in g]
                if self.strand == 0:
                    ct = Coords(ctg=self.contigID, strand=self.strand, beg=min(exon), end=max(exon))
                    newchain.append(ct)
                elif self.strand == 1:
                    ct = Coords(ctg=self.contigID, strand=self.strand, beg=max(exon), end=min(exon))
                    newchain.append(ct)
            return newchain
        else:
            return None


    def toGFF3(self, ID : str = "", parent : str = "", source : str = '.',  with_exon : bool = False,
               with_utr : bool = False, CDSscores : bool = True, with_cds : bool = True) -> str:
        """:return gff representation of instance as text string.
        Parent can be manually set or by hierarchical upwards element (usual case),
        element ID should generally be inferred from self.tid, boolean arguments
        include different formatting/output options"""
        gff = []
        ori = '+' if self.strand == 0 else '-'
        r = []
        r.append("%s" % (self.contigID))
        r.append("%s\tmRNA" % (source))
        r.append("%s\t%s" % (self.minimum(), self.maximum()))
        if isinstance(self.score, int):
            r.append("%i\t%s\t." %(self.score, ori))
        elif isinstance(self.score, float):
            r.append("%.3f\t%s\t." %(self.score, ori))
        else:
            r.append(".\t%s\t." %(ori))
        elmID = ID if ID != "" else self.tid
        tags = ["ID=%s;Parent=%s" %(elmID, parent)]
        for tag,value in self.tags.items():
            tags.append("%s=%s" %(tag,value))
        if self.descriptionLine:
            tags.append("Description=%s" %(self.descriptionLine))
        r.append("%s" %(';'.join(tags)))
        gff.append("%s" %('\t'.join(r)))

        if with_utr:
            utr5, utr3 = self.getUTRs()
            if utr5 != None:
                c = 1 if self.strand == 0 else len(self)
                utr5.sort()
                for el in utr5:
                    assert el.strand == self.strand
                    r = []
                    r.append("%s" % (self.contigID))
                    r.append("%s\tfive_prime_UTR" % (source))
                    r.append("%s" % (el.minimum()))
                    r.append("%s" % (el.maximum()))
                    r.append(".\t%s\t.\tID=%s_five_prime_UTR_%i;Parent=%s" % (ori, elmID, c, elmID))
                    c += 1 if self.strand == 0 else -1
                    gff.append('\t'.join(r))

        # print exon elements if asked for
        if with_exon:
            c = 1 if self.strand == 0 else len(self)
            for el in self:
                assert el.strand == self.strand
                r = []
                r.append("%s" % (self.contigID))
                r.append("%s\texon" % (source))
                r.append("%s" % (el.minimum()))
                r.append("%s" % (el.maximum()))
                r.append(".\t%s\t.\tID=%s_exon_%i;Parent=%s" % (ori, elmID, c, elmID))
                c += 1 if self.strand == 0 else -1
                gff.append('\t'.join(r))

        # now get coding region as gff if asked for
        if with_cds:
            if len(self.coding) > 0:
                gff.append("%s" %(self.coding.toGFF3(elmID, source=source, CDSscores=CDSscores)))

        if with_utr:
            if utr3 != None:
                c = 1 if self.strand == 0 else len(self)
                utr3.sort()
                for el in utr3:
                    assert el.strand == self.strand
                    r = []
                    r.append("%s" % (self.contigID))
                    r.append("%s\tthree_prime_UTR" % (source))
                    r.append("%s" % (el.minimum()))
                    r.append("%s" % (el.maximum()))
                    r.append(".\t%s\t.\tID=%s_three_prime_UTR_%i;Parent=%s" % (ori, elmID, c, elmID))
                    c += 1 if self.strand == 0 else -1
                    gff.append('\t'.join(r))

        return "%s" %('\n'.join(gff))


@total_ordering
class GeneLocus:

    """Data container for a gene/locus"""
    def __init__(self, gid : str = "", ctg : str = "", src : str = "",
                 sofa : str ="gene", score : Union[int, float, None] = None):
        self.geneID = gid
        self.contigID = ctg
        self.transcripts = {}  # all alternative transcripts for locus
        self.score = score
        self.source = src
        self.sofa   = sofa
        self.tags = {}
        self.descriptionLine = ""  # as previous classes, legacy code to ensure downwards compatibility


    # overload for sort by genomic position
    def __eq__(self,other):

        return (self.contigID,self.minimum(),self.maximum()) == \
               (other.contigID,other.minimum(),other.maximum())


    def __lt__(self,other):

        return (self.contigID,self.minimum(),self.maximum()) < \
               (other.contigID,other.minimum(),other.maximum())

    def minimum(self) -> int:
        """:return genomic minimum locus position"""
        gmin = 999999999999
        for tid,t in self.transcripts.items():
            gmin = min(gmin,t.minimum())
        return gmin


    def codingMin(self) -> int:
        """:return genomic minimum CDS position"""
        gmin = 999999999999
        for tid,t in self.transcripts.items():
            gmin = min(gmin, t.coding.minimum())
        return gmin

    def maximum(self) -> int:
        """:return genomic maximum locus position"""
        gmax = -1
        for tid,t in self.transcripts.items():
            gmax = max(gmax,t.maximum())
        return gmax


    def codingMax(self) -> int:
        """:return genomic minimum CDS position"""
        gmax = -1
        for tid,t in self.transcripts.items():
            gmax = max(gmax,t.coding.maximum())
        return gmax


    def maxCDSsize(self) -> int:
        """:return maximal spliced size of coding sequences of all associated transcripts"""
        sz = 0
        for tid,t in self.transcripts.items():
            sz = max(sz,t.coding.size())
        return sz


    def overlaps(self, other : Union[GeneLocus, CoordChain]) -> bool:
        """:return bool overlap for testing overlap with <other>, one bp is sufficient"""
        if self.contigID != other.contigID:
            return False
        if self.minimum() <= other.maximum() and self.maximum() >= other.minimum():
            return True
        return False


    def overlap(self, other : Union[GeneLocus, CoordChain]) -> bool:
        """:return bool for overlap with other gene. Legacy code to method overlaps"""
        return self.overlaps(other)


    def codingOverlaps(self, other : GeneLocus) -> bool:
        """:return bool overlap for all coding regions with all coding regions of other gene locus.
        One bp overlap between any cross comparison is sufficient"""
        if self.contigID != other.contigID:
            return False
        for tid1, t1 in self.transcripts.items():
            for tid2, t2 in other.transcripts.items():
                if t1.coding.overlaps(t2.coding):
                    return True
        return False

    def strand(self) -> int:
        """:return strand of gene locus: 0:forward or 1:reverse"""
        strands = set()
        for tid,transcript in self.transcripts.items():
            strands.add(transcript.strand)
        assert len(strands) == 1
        return strands.pop()

    def addmRNA(self, item : mRNA):
        """adds a new mRNA/splice form to locus"""
        assert isinstance(item, mRNA)
        if self.contigID == "":
            self.contigID = item.contigID
        assert self.contigID == item.contigID
        self.transcripts[item.tid] = item

    def basicCoordinates(self) -> tuple[str, int, int, int]:
        """:return most basic description of position of locus in genome.
        Obsolete legacy code"""
        gmax = -1
        gmin = 999999999999
        strand = -1
        ctg    = ""
        for tid,t in self.transcripts.items():
            try:
                gmax = max(gmax,t.maximum())
            except ValueError:
                print(tid)
                print(t)
                sys.exit()
            gmin = min(gmin,t.minimum())
            if strand < 0:
                strand = t.strand
            else:
                assert strand == t.strand
            if not ctg:
                ctg = t.contigID
            else:
                assert ctg == t.contigID

        return ctg, gmin, gmax, strand


    def toGFF3(self, ID : Union[str, None] = None, with_exon : bool = False,
               with_cds : bool = True, with_utr : bool = False, CDSscores : bool = True,
               mrna2report : Union[str, set] = 'all') -> str:
        """:return gff representation of GeneLocus as text string.
        Element <ID> can be manually set or geneID attribute is taken (most common/default use case),
        boolean arguments are transferred downwards to all mRNAs/transcripts.
        <mrna2report> controls which transcripts are reported, set to 'all' for all splice
        variants (most common/default use case)"""
        # first make gene line
        gid = self.geneID if not ID else ID
        gff = []
        ctg, gmin, gmax, strand = self.basicCoordinates()

        ori = '+' if strand == 0 else '-'
        r = []
        r.append("%s" % (self.contigID))
        r.append("%s\t%s" % (self.source, self.sofa))
        r.append("%s\t%s" % (gmin, gmax))
        # score can be numeric or none
        if isinstance(self.score, int):
            r.append("%i\t%s\t." %(self.score, ori))
        elif isinstance(self.score,float):
            r.append("%.3f\t%s\t." %(self.score, ori))
        else:
            r.append(".\t%s\t." %(ori))
        # generate famous column 9 of gff
        tags = ["ID=%s" %(gid)]  # required
        for tag,value in self.tags.items():
            tags.append("%s=%s" %(tag,value))
        if self.descriptionLine:
            tags.append("Description=%s" %(self.descriptionLine))
        r.append("%s" %(';'.join(tags)))
        gff.append('\t'.join(r))

        # now process each mRNA
        mrnas = [t for t in self.transcripts.values()]
        mrnas.sort()
        WITH_CDS = with_cds
        if self.sofa == "pseudogene":
            WITH_CDS = False
        for t in mrnas:
            if mrna2report == 'all':
                gff.append("%s" %(t.toGFF3(ID=t.tid, parent=gid, source=self.source, with_utr=with_utr,
                                           with_exon=with_exon, CDSscores=CDSscores, with_cds=WITH_CDS)))
            else:
                if t.tid in mrna2report:
                   gff.append("%s" %(t.toGFF3(ID=t.tid, parent=gid, source=self.source, with_utr=with_utr,\
                                              with_exon=with_exon, CDSscores=CDSscores, with_cds=WITH_CDS)))

        return "%s" %('\n'.join(gff))
