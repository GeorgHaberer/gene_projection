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

### a. preprocessing steps

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

Lastly, a file listing all genotype names and the path to their genome sequences was compiled -> _run.list_ 

Preprocessing was then run using the script run_preprocessing.py calling minimap2[³] (https://github.com/lh3/minimap2) as spliced aligner:
> python run_preprocessing.py -r _run.list_ -o _project_outroot_ --offsetfile _offsets.txt_ -t _nr.transcript.fa_ -p _nr.protein.cdhit_ --mmap "path to minimap2 binary"

The final result are _*raw.gff.gz_ files of all matches for each genotype in the respective output directories _project_outroot_/_genotype_name_ . 

### b. projection steps/parameters

Projections were computed for each genotype using the following command:
> python run_geneprojection.py -r _run.list_ -o _proj_outroot_ --pfamfile _pfam.codes.txt_ --transposonfile _transposon.codes.txt_ --plastidfile _plastid.codes.txt_ --scorefile _maxscores.txt_

The rules/filters for each round of insertion are shown in following table, for a description of the workflow and explanation of the rules see parent directory:

| round | rule/features |
| --- | --- |
| 1 | contiguous ORF, start&stop codon, predicted gene, pfam domain, rel.score >= 0.85, unique |
| 2 | contiguous ORF, start&stop codon, predicted gene, pfam domain, rel.score >= 0.95, non-unique |
| 3 | contiguous ORF, start&stop codon, predicted gene, rel.score >= 0.85, unique |
| 4 | contiguous ORF, start&stop codon, predicted gene, rel.score >= 0.95, non-unique |
| 5 | contiguous ORF, predicted gene, pfam domain, rel.score >= 0.95, unique |
| 6 | contiguous ORF, start&stop codon, plastid-related gene, rel.score >= 0.9, unique |
| 7 | contiguous ORF, start&stop codon, transposon-related gene, rel.score >= 0.9, unique |


## _Hordeum bulbosum_ projections (pan-bulbosum project)




 
 -  








# references
[¹] CD-HIT: accelerated for clustering the next generation sequencing data", Limin Fu, Beifang Niu, Zhengwei Zhu, Sitao Wu & Weizhong Li. Bioinformatics, (2012) 28:3150-3152

[²] Pfam: The protein families database in 2021: J. Mistry, S. Chuguransky, L. Williams, M. Qureshi, G.A. Salazar, E.L.L. Sonnhammer, S.C.E. Tosatto, L. Paladin, S. Raj, L.J. Richardson, R.D. Finn, A. Bateman, Nucleic Acids Research (2021) doi: 10.1093/nar/gkaa913 

[³] Li, H. (2021). New strategies to improve minimap2 alignment accuracy. Bioinformatics, 37:4572-4574. doi:10.1093/bioinformatics/btab705


