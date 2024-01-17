# gene_projection

## 0. Preliminary Remark

  This project is under active development and has already been used in several plant (pan-)genome projects[^1] [^2]. To provide a transparent and persistent documentation, description and code for specific versions applied in each project are hosted in sub-directories. A general overview and rationale of the pipeline are provided here at the main site while exact parameter settings and workflows are described in the respective sub-directories.


## 1. Introduction and Motivation

  After the assembly of a genome sequence, determination of its gene content is an essential next task in any genome project and prerequisite for numerous functional and evolutionary downstream analysis. Today's state-of-the-art (structural) gene annotations are generally based on experimental evidences like short and long read RNAseq data. In higher species, however, multiple cell types,  a series of development stages and responses to numerous complex environmental conditions make it nearly impossible to sample the entire transcriptomic space. To address the complete gene space, transcript-based gene models are hence complemented by homology searches with genes of previously annotated species and de novo genefinders to generate a consensus gene annotation. 
  With the recent advances in sequencing technologies, the pan-genome era has started and frequently a few dozens or even hundreds of genotypes or closely related species of one clade are sequenced. Evidence-based gene annotations as described above are very cost- resource- and labor-intensive. Hence, for large-scale pan-genomic studies often only one or a small subset of genotypes is annotated by transcriptomic data. To close this gap and to have an estimate of the gene space of the other genotypes, the gene projection pipeline has been developed. In principle, it is based on homology mappings of the evidence-based gene models to the genome sequences of the unannotated genotypes or species. The aim was to develop a pipeline that:
  + can be run on a single workstation,
  + allows a quick estimation of the gene space in each genotype,
  + provides an accurate estimation of the gene space,
  + should also detect copy number variations (eg. tandemly repeated genes),
  + has a small number of dependencies,
  + has the flexibility to include additional aligners or third-party tools,
  + and is scalable to large pangenome sets.

> [!NOTE] 
> The homology-based gene transfer or annotation also includes several potential limitations a user should be aware of, in particular:
> 
> * genotype-specific exon-intron structures may be insufficiently or incorrectly modelled by the representative genes
> * UTRs rely on the UTRs of the representative gene and may not reflect genotype-specific UTRs
> * alien introgressions eg by horizontal gene transfer is likely not detected.

However, BUSCO analysis and orthologous gene families by orthofinder showed that the gene projection pipeline produces results for projected gene contents that are highly comparable to evidence-based gene annotations.


## 2. Requirements

A fairly powerful workstation (64 bit) with at least 16 Gb RAM, better 32-64 Gb for runs using multiple cores is recommended. You may have to adjust the number of parallel cpus in the scripts to your resources. For very large pan-genome projects, or to gain additional speed-up, the analysis could be run in a batch queue. 

The following software is required:
- python >=3.8 (https://www.python.org/)
- biopython >=1.8 (https://biopython.org)

For the preprocessing steps, several third-party tools should be available: a sequence clustering software, an annotation tool and at least splice-aware aligner. Here are the tools that have been used in several genome projects, of course any can be supplemented by your favorite clustering/annotation/mapping tool. 

- minimap2 (https://github.com/lh3/minimap2)
- cd-hit (https://github.com/weizhongli/cdhit/tree/master)
- hmmer/hmmsearch (https://github.com/EddyRivasLab/hmmer)
- AHRD (https://github.com/groupschoof/AHRD)

Currently only minimap2 is supported in the repository but I am planning to provide additional integration for the following aligners:
- miniprot (https://github.com/lh3/miniprot)
- gmap (https://github.com/juliangehring/GMAP-GSNAP/tree/master)
- blat (https://github.com/djhshih/blat; conda or UCSC genome hub https://hgdownload.cse.ucsc.edu/admin/exe/)


## 3. Principle Workflow

The rationale of the pipeline is pretty simple: it maps the non-redundant transcripts of a previous annotation (the **source models**; preferably from a set of evidence-based annotations of closely related genotypes or species) to the genome sequences of the **target** genotypes/species and estimates from these matches the gene content by applying user-defined rules. There are two major steps, first preprocessing and mapping the input sequences, and second a step-wise and rule-based selection of non-overlapping matches. The following section provides a graphical scheme and explanation of the two steps.

### 3.1. Preprocessing and Mapping



**Fig1. Pipeline preprocessing steps.**
<p align="center">
  <img src="/images/preprocess_overview.png" width="600">
</p>


### 3.2. Annotation by Projections

blalballa


**Fig 2. Annotation by stepwise match insertion.**
<p align="center">
  <img src="/images/projection_insert_scheme.png" width="600">
</p>

blalbla


## 4. Citation

A publication of the pipeline is in progress, in the meantime, if you use code from this repository, please refer to it as github URL, and include one of the following references:

[^1]: Walkowiak, S, Gao, L, Monat, C, Haberer G, et al. Multiple wheat genomes reveal global variation in modern breeding. Nature 588, 277–283 (2020). https://doi.org/10.1038/s41586-020-2961-x

[^2]: Jayakodi M, Padmarasu S, Haberer G, et al. The barley pan-genome reveals the hidden legacy of mutation breeding. Nature. 2020 Dec;588(7837):284-289. https://doi: 10.1038/s41586-020-2947-8.












