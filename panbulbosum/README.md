# Pan-bulbosum project

## 0. Preliminary Remark

  This repository directory describes the gene projection approach applied to the second, revised version of the pan bulbosum genome project. A general overview of the gene projection pipeline is provided in the parent github directory at https://github.com/GeorgHaberer/gene_projection.

## 1. Overview

  In contrast to other genome projects using the gene projection pipeline, in the pan-bulbosum project gene projection were only used to complement for high quality gene models potentially missed in a evidence-based annotation. Figure 1 shows a simplified graphical overview of the workflow, detailed descriptions, scripts and parameters are listed in the following chapters.


## 2. Preprocessing data

  Input for the gene projections were evidence-based gene annotations as described in the Materials and Methods section of the manuscript. These proteins were clustered by CD-hit[[1]](#1) to obtain the input source models for the pipeline and to remove (nearly) identical models between evidence-based annotated bulbosum genotypes. 


  A second clustering by CD-hit grouped similar source models to so called 'meta-clusters' with the following commands:

  
  


  In a second 





## References

<a id="1">[1]</a> 
Weizhong Li, Adam Godzik, Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences, Bioinformatics, Volume 22, Issue 13, July 2006, Pages 1658–1659, https://doi.org/10.1093/bioinformatics/btl158

