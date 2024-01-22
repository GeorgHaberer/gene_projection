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

All analysis was performed on a desktop workstation AMD Ryzen 7600 with 48 Gb RAM and 2Tb PCI4.0 SSD harddisk. The following steps were performed to obtain the final projection:

Nearly identical gene models were clustered using cd-hit (v4.8.1)[¹] to obtain a non-redundant set of source transcripts and proteins:
> cd-hit -M 16000 -T 20 -S 4 -c 0.99 -o _nr.protein.cdhit_ -i _protein.all.sources.fa_
>
> non-redundant transcripts were extracted from _nr.protein.cdhit_ to file _nr.transcript.fa_


A local copy of the pfam database was obtained from https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam34.0/ and pfam domains of the non-redundant source models were determined by a conda installation of hmmsearch[²] :
> hmmsearch --domtblout _nr.protein.pfam.domtblout_ --noali -E 0.05 --cpu 8 _Pfam-A.hmm_ _nr.protein.cdhit_
>
> source ids with a pfam domain were extracted from the first column to file _pfam.codes.txt_


Maximal attainable scores for each source model were retrieved from their self-alignments using script scoring_parallel.py in https://github.com/GeorgHaberer/gene_projection/edit/main/geproj_utils:
> python scoring_parallel.py -p _nr.protein.cdhit_ -o _maxscores.txt_ -n 8





 
 -  








# references
[¹] CD-HIT: accelerated for clustering the next generation sequencing data", Limin Fu, Beifang Niu, Zhengwei Zhu, Sitao Wu & Weizhong Li. Bioinformatics, (2012) 28:3150-3152

[²] Pfam: The protein families database in 2021: J. Mistry, S. Chuguransky, L. Williams, M. Qureshi, G.A. Salazar, E.L.L. Sonnhammer, S.C.E. Tosatto, L. Paladin, S. Raj, L.J. Richardson, R.D. Finn, A. Bateman, Nucleic Acids Research (2021) doi: 10.1093/nar/gkaa913 
