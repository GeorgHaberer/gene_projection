This directory contains code applied to generate gene projections of two pan-genome projects. It provides a persistent and detailed description of the pipeline version and parameter settings used. The first project comprised 76 _Hordeum vulgare_ genotypes of which 20 were annotated based on genotype-specific transcriptome data, protein homologies and _de novo_ gene finders. Projections were computed for the residual 56 genotypes. The second project analyzed 11 _Hordeum bulbosum_ genotypes with one annotated by experimental evidences.

> [!NOTE]
> A general description of the workflow of the gene projection pipeline is desribed in detail in the parent directory: https://github.com/GeorgHaberer/gene_projection/edit/main


## _Hordeum vulgare_ projections (pan-barley project)

The projections started with the following input sequences and information available from the evidence-based annotation of the 20 genotypes:
 - transcript sequences of (representatives) annotation of 20 genotypes : _<transcript.all.sources.fa>_
 - start position of start codon in transcripts and CDS size: _<offsets.txt>_
 - protein sequences of annotations of 20 genotypes: _<protein.all.source.fa>_
 - file listing transposon-related genes in above annotations: _<transposon.codes.txt>_
 - file listing plastid-related genes: _<plastid.codes.txt>_


 
 -  
