# gene_projection

## 0. Preliminary Remark

  This project is under active development and has already been used in several plant (pan-)genome projects. To provide a transparent and persistent documentation, description and code for specific versions applied in each project are hosted in sub-directories. A general overview and rationale of the pipeline are provided here at the main site while exact parameter settings and workflows are described in the respective sub-directories.


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

The homology-based gene transfer or annotation also includes several potential limitations a user should be aware of, in particular:

* genotype-specific exon-intron structures may be insufficiently or incorrectly modelled by the representative genes
* UTRs rely on the UTRs of the representative gene and may not reflect genotype-specific UTRs
* alien introgressions eg by horizontal gene transfer is likely not detected.

However, BUSCO analysis and orthologous gene families by orthofinder showed that the gene projection pipeline produces results for projected gene contents that are highly comparable to evidence-based gene annotations.

## 2. Requirements

A fairly powerful workstation (64 bit) with at least 16 Gb RAM, better 32-64 Gb for runs using multiple cores is recommended. You may have to adjust the number of parallel cpus in the scripts to your resources. For very large pan-genome projects, or to gain additional speed-up, the analysis could be run in a batch queue. 

The following software is required:
- python >=3.8 (https://www.python.org/)
- biopython >=1.8 (https://biopython.org)


In addition, one of the following alignment tools should be installed (only minimap2 currently supported in the repository):

- minimap2 (https://github.com/lh3/minimap2)

- miniprot (https://github.com/lh3/miniprot)
- gmap (https://github.com/juliangehring/GMAP-GSNAP/tree/master)
- blat (https://github.com/djhshih/blat; conda or UCSC genome hub https://hgdownload.cse.ucsc.edu/admin/exe/)

## 3. Principle Workflow


blabla

**Fig1. Pipeline preprocessing steps.**
<p align="center">
  <img src="/images/preprocess_overview.png" width="600">
</p>


blalballa


**Fig 2. Annotation by stepwise match insertion.**
<p align="center">
  <img src="/images/projection_insert_scheme.png" width="600">
</p>

blalbla


## 4. Citation











