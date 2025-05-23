import paf2gff_new


def paf_transformer(paf : str, ctgfn : str, protfn : str, gff : str, offsetfile : str):

    models = paf2gff_new.pafmap2Gff(paffile=paf, ctgfile=ctgfn, protfile=protfn, offsetfile=offsetfile)
    models.sort()
    with open(gff, 'w') as out:
        for gene in models:
            out.write("%s\n" %(gene.toGFF3(with_exon=True)))


if __name__ == '__main__':

    paffile = sys.argv[1]
    ctgfile = sys.argv[2]
    protfile = sys.argv[3]
    gff = sys.argv[4]
    offsetfile = sys.argv[5]
    paf_transformer_new(paffile, ctgfile, protfile, gff, offsetfile)
