###########################################################
#                                                         #
#    Datamodels for biological sequences                  #
#                                                         #
###########################################################

# TO DO: should change to biopython's implementation of Seq*objects :)

from __future__ import annotations


# -----------------------------------------------
def complement(letter : str) -> str:
    """ function to complement one nucleotide."""
    if letter == 'A':
        return 'T'
    elif letter == 'C':
        return 'G'
    elif letter == 'G':
        return 'C'
    elif letter == 'T':
        return 'A'
    elif letter == 'N':
        return 'N'
    if letter == 'a':
        return 'T'
    elif letter == 'c':
        return 'G'
    elif letter == 'g':
        return 'C'
    elif letter == 't':
        return 'A'
    elif letter == 'n':
        return 'N'
    elif letter == '-':
        return '-'
    elif letter == '.':
        return '.'
    else:
        return '?'


# --------------------------------------------------------------------
class FastaSequence(object):
    """Datastructure for FASTA formatted biological sequence."""

    def __init__(self, ID : str = "", seq : str = "", descr : str = ""):

        self.ID = ID
        self.sequence = seq
        self.description = descr


    def __getitem__(self, sliced) -> FastaSequence:
        return FastaSequence(ID=self.ID,seq=self.sequence[sliced])


    def __str__(self) -> str:

        if self.description:
            r = [">%s %s" % (self.ID, self.description)]
        else:
            r = [">%s" % (self.ID)]
        for i in range(0, len(self), 80):
            r.append("%s" % (self.sequence[i:i + 80]))
        return '\n'.join(r)


    def __len__(self) -> int:

        return len(self.sequence)




# --------------------------------------------------------------------
class DnaSequence(FastaSequence):
    """Datamodel for DNA sequence"""

    # define class attributes
    Genetic_Code = {'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C', \
                    'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C', \
                    'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*', \
                    'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W', \
                    'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R', \
                    'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R', \
                    'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R', \
                    'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R', \
                    'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S', \
                    'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S', \
                    'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R', \
                    'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R', \
                    'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G', \
                    'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G', \
                    'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G', \
                    'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G', \
                    'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C', \
                    'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C', \
                    'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*', \
                    'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W', \
                    'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R', \
                    'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R', \
                    'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R', \
                    'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R', \
                    'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S', \
                    'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S', \
                    'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R', \
                    'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R', \
                    'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G', \
                    'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G', \
                    'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G', \
                    'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'}

    Complement = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n', \
                  'w': 'w', 's': 's', 'r': 'y', 'y': 'r', 'k': 'm', \
                  'm': 'k', 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', \
                  'N': 'N', 'W': 'W', 'S': 'S', 'R': 'Y', 'Y': 'R', \
                  'K': 'M', 'M': 'K'}


    def __init__(self, ID : str = "", seq : str = "", descr : str = ""):

        FastaSequence.__init__(self, ID, seq, descr)


    def getReverse(self) -> DnaSequence:

        seq = list(self.sequence)
        seq.reverse()
        return DnaSequence(ID="%s" % (self.ID), seq=''.join(seq), descr="%s" % (self.description))


    def getComplement(self) -> DnaSequence:

        compl = []
        for letter in self.sequence:
            try:
                compl.append(DnaSequence.Complement[letter])
            except KeyError:
                compl.append('N')
        return DnaSequence(ID="%s" % (self.ID), seq=''.join(compl), descr="%s" % (self.description))


    def getReverseComplement(self) -> DnaSequence:

        rc = []
        for letter in reversed(self.sequence):
            try:
                rc.append(DnaSequence.Complement[letter])
            except KeyError:
                rc.append('N')
        return DnaSequence(ID="%s" % (self.ID), seq=''.join(rc), descr="%s" % (self.description))


    def translate(self, fromPos : int = 0, readThrough : bool = False) -> str:

        orf = []
        seq = self.sequence.lower()  # make sure mixed upper- and lowercase is handled
        for i in range(fromPos, len(seq), 3):
            codon = seq[i:i + 3]
            if len(codon) < 3: break
            try:
                aa = DnaSequence.Genetic_Code[codon]
            except KeyError:
                aa = 'X'
            orf.append(aa)
            if aa == '*' and not readThrough: break
        return ''.join(orf)


    def getAllORFs(self, minLength : int = 60) -> list:

        assert minLength > 0
        orfs = []
        # for each frame j
        for j in range(3):
            b = j
            while b < len(self.sequence) - minLength + 1:
                orf = self.translate(fromPos=b)
                if len(orf) >= minLength:
                    start = b + 1
                    orfs.append((start, FastaSequence(seq=''.join(orf))))
                b += 3 * len(orf)

        return orfs


# -------------------------------------------------------------------------------
class GenomicTemplate(DnaSequence):
    """Genomic template sequence, derived from DnaSequence Class."""

    def __init__(self, fileLocation : str = "", acc : str = "", ID : str = "",
                 seq : str = "", descr : str = ""):

        DnaSequence.__init__(self, ID, seq, descr)
        self.fileLocation = fileLocation
        self.accession = acc


    def getSequenceByCoords(self, coords : CoordChain) -> DnaSequence:
        # this is kind of obsolete, coordchain has own method to extract its sequence
        # from the genome
        coords.sort()
        r = []
        for el in coords:
            try:
                x = min(el.beg, el.end)
                y = max(el.beg, el.end)
            except AttributeError:
                x = min(el.begin, el.end)
                y = max(el.begin, el.end)
            r.append(self.sequence[x - 1:y])

        s = DnaSequence(seq=''.join(r))
        if coords.getStrand() == 'C':
            return s.getReverseComplement()
        else:
            return s
