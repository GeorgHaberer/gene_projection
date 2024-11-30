# Pan-bulbosum project

## 0. Preliminary Remark

  This repository directory describes the gene projection approach applied to the second, revised version of the pan bulbosum genome project. A general overview of the gene projection pipeline is provided in the parent github directory at https://github.com/GeorgHaberer/gene_projection.

## 1. Overview

  In contrast to other genome projects using the gene projection pipeline, in the pan-bulbosum project gene projection were only used to complement for high quality gene models potentially missed in a evidence-based annotation. Figure 1 shows a simplified graphical overview of the workflow, detailed descriptions, scripts and parameters are listed in the following chapters.

  **Fig1. Workflow of the pan-bulbosum projections.**
<p align="center">
  <img src="/images/workflow_scheme.png" width="600">
</p>

<p>
  <br />
</p>

The following tools and modules were used and are required:
- python 3.9 or higher
- biopython v1.81
- bioutils (https://github.com/GeorgHaberer/gene_projection/tree/main/panhordeum/bioutils) either copy scripts in same directory of bulbosum scripts or add to your PYTHONPATH
- CD-hit[[1]](#1)
- orthofinder [[2]](#2)
- minimap2 [[3]](#3)
- miniprot [[4]](#4)
- blat v34 [[5]](#5)
- hmmer/pfam v34 [[6]](#6)
- optional AHRD [[7]](#7)


## 2. Preprocessing data

### 2.1 Non-redundant input (source) models

  Input for the gene projections were evidence-based gene annotations as described in the Materials and Methods section of the manuscript. These proteins were clustered by CD-hit[[1]](#1) to obtain the input source models for the pipeline and to remove (nearly) identical models between evidence-based annotated bulbosum genotypes. This step reduces computational load by maintaining full sensitivity. The script <make_geneweights.py> reports the number of each (nearly) identical gene models per cluster in the file <hbulb.weights.txt>.

> cd-hit -M 16000 -T 20 -S 4 -c 1.0 -o _hbulbosum.evidence.high.cdhit_ -i _protein.all.sources.fa_

> python make_geneweights.py _hbulbosum.evidence.high.cdhit_ _hbulb.weights.txt_

  A second clustering by CD-hit grouped similar source models to so called 'meta-clusters' with the following commands. Meta-cluster information is used in the pipeline to track whether highly similar models have already been integrated in the step-wise projections (see pipeline overview in parent repo directory) to ensure balanced source insertions.

> cd-hit -M 16000 -T 20 -S 10 -c 0.95 -o _hbulbosum.evidence.high.S10c95.cdhit.clstr_ -i _hbulbosum.evidence.high.cdhit_

> python make_metaclusters.py _hbulbosum.evidence.high.S10c95.cdhit.clstr_ _metafile_outputpath_
  
### 2.2 Orthologous source models

  Hierarchical orthogroups were determined for the evidence based genes using orthofinder:

> orthofinder -t 40 -a 8 <bulb.evi.proteins>

  The directory <bulb.evi.proteins> contained for each haplotype of the evidence-based annotated _H.bulbosum_ genotypes the representative high confidence protein models. Hierarchical orthogroups were extracted from the file <N0.tsv> (see orthofinder manual [[2]](#2)). For each source transcript/protein, the number of orthologs was determined by the script make_orthocounts.py.

### 2.3 Additional preprocessing steps

  UTR (5') sizes of source models were computed from the GFF files of the evidence-based annotations using the script make_offsets.py and saved in the <hbulb.offsets.txt> file.
  
  Plastid- and transposon-related genes were identified from the source models running the AHRD pipeline [[7]](#7) and this information was saved in files <plastid.codes> and <transposon.codes>. Based on this gene type information, transposon-derived models were removed from the input protein and transcript IDs and these clean sequence files were separately stored as <hbulb.clean.proteins.fa> and <hbulb.clean.transcripts.fa>. These files were employed as input sources for miniprot and blat alignments to reduce the high complexity of transposon mappings strongly impairing aligner performance.

  Lastly, hmmer [[6]](#6) and PFAM database version 34 searched and identified protein domains in the source model and associated the top evalue to each input protein/transcript.


## 3. De novo projections

  The aim of this project was not to make a _de novo_ annotation of the H. bulbosum genotypes by projections but rather to identify high quality candidates from a compilation of evidence-based gene models to complement their annotations. Because the projection pipeline scores mappings of the input sources relative to their maximal achievable score, ie. alignment of the source protein to itself, self scores of the non-redundant input of section 2 were determined by score_parallel.py:

> python score_parallel.py -i _hbulbosum.evidence.high.cdhit_ -n 5 -o _hbulb.maxscores.txt_
>

  Next, initial mappings of the non-redundant source genes were obtained using bulb_annopipe.py script:

> python bulb_annopipe.py _<genotype>_ _<genome_fasta_file_of_genotype>_
>

  From the initial mappings, a set of projected models as consolidation candidates were generated applying the script anno_projection.py. The pipeline inserts in a step-wise manner non-overlapping initial mappings by quality criteria (see section 2), progressing from higher to lower qualities (see also parent repo directory https://github.com/GeorgHaberer/gene_projection for more details). There were five steps with the following criteria:
  1. No transposon models, pfam domain with an e-value < 0.01, and at least one ortholog
  2. 



## 4. Consolidation




    trgmatches = [x[0] for x in allmatches if x[3] <= 0.01 and x[4] == False and x[10] > 0]
    
    trgmatches = [x[0] for x in allmatches if x[4] == False and x[13] == False and x[10] > 0]
    
    trgmatches = [x[0] for x in allmatches if x[3] <= 0.01 and x[4] == False and x[13] == False and x[5] >= 6 and x[8] >= 0.9]
    
    trgmatches = [x[0] for x in allmatches if (x[3] <= 0.01 or x[5] >= 2) and x[4] == False and x[8] >= 0.9]
    
    trgmatches = [x[0] for x in allmatches if x[4] == True and x[8] >= 0.9]
    



## References

<a id="1">[1]</a> 
Weizhong Li, Adam Godzik, Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences, Bioinformatics, Volume 22, Issue 13, July 2006, Pages 1658–1659, https://doi.org/10.1093/bioinformatics/btl158

<a id="2">[2]</a>
Emms, D.M., Kelly, S. OrthoFinder: solving fundamental biases in whole genome comparisons dramatically improves orthogroup inference accuracy. Genome Biol 16, 157 (2015). https://doi.org/10.1186/s13059-015-0721-2
(see also https://github.com/davidemms/OrthoFinder)

<a id="3">[3]</a>
Li, H. (2021). New strategies to improve minimap2 alignment accuracy. Bioinformatics, 37:4572-4574. doi:10.1093/bioinformatics/btab705
(see also https://github.com/lh3/minimap2)

<a id="4">[4]</a>
Heng Li, Protein-to-genome alignment with miniprot, Bioinformatics, Volume 39, Issue 1, January 2023, btad014, https://doi.org/10.1093/bioinformatics/btad014
(see also https://github.com/lh3/miniprot)

<a id="5">[5]</a>
Kent WJ. BLAT--the BLAST-like alignment tool. Genome Res. 2002 Apr;12(4):656-64. doi: 10.1101/gr.229202. PMID: 11932250; PMCID: PMC187518
(see also http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat)

<a id="6">[6]</a>
Pfam: The protein families database in 2021: J. Mistry, S. Chuguransky, L. Williams, M. Qureshi, G.A. Salazar, E.L.L. Sonnhammer, S.C.E. Tosatto, L. Paladin, S. Raj, L.J. Richardson, R.D. Finn, A. Bateman, Nucleic Acids Research (2021) doi: 10.1093/nar/gkaa913
(see also https://github.com/EddyRivasLab/hmmer)

<a id="7">[7]</a>
AHRD: https://github.com/groupschoof/AHRD




