---
title: 'methometR: Extracting methylation information from Optical Genome Maps'
author: "Surajit Bhattacharya,Hayk Barsheghyan,Seth Berger and Eric Vilain"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{nanotatoR}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

# Introduction
Short-read bisulfite sequencing (**SRBS**) and Illumina Epic methylation arrays (**IEMA**) have helped in identifying epigenetic signatures in research and clinic. Although, due to technical limitations, they do not provide long-range haplotype specific methylation states, but rather detect signals that are averaged for specific genomic positions. These limitations can be alleviated with a novel dual-label optical genome mapping (DL-OGM) technology for detection epigenetic changes. The method relies on differential labeling of high molecular weight DNA. Long DNA molecules are nicked with BspQI endonuclease and labeled with red fluorescent nucleotides, followed by treatment with M.TaqI methyltransferase that attaches green fluorescent cofactor onto non-methylated (hypomethylated) CpGs in TCGA sequences throughout the genome. Though the technology, identifies and visually characterizes regions of hypomethylations, no tool is available that performs quantification of the signals, to differentiate between an affected and a control allele/condition.

We have built a R based tool, ***methometR***, that quantifies/normalizes the methylation signals to evaluate differential methylation patterns between cases and controls.

***methometR*** is built using R, and ggplots is used for visualization. ***methometR*** functions in 3 steps. First, extract and map green labels to reference genome, to identify hypomethylated labels, across the whole genome. Next, quantify /normalize the methylation levels of a region based on coverage of molecules that span the region, Finally, visualization of methylation levels, in the form of barplots, across user specified region range, to differentiate patterns between cases and controls.

#Package Installation

***methometR*** is currently available from the GitHub repository. Installation method is as follows:
```{r eval=FALSE}
library("devtools")
devtools::install_github("VilainLab/methometR")
```

```{r eval=TRUE}
library("methometR")
```

The ***methometR*** package is compatible with R versions ≥ 4.1.


#Package Functionalities

Given a list of molecule and contig infromation map it helps to identify the methylation patterns in the OGM maps.

##Modifying Molecule map

```{r eval=FALSE}
xmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1/merged_smaps/SampContigMolecule.xmap", package="methometR")
 cmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1/merged_smaps/SampMolecule_q.cmap", package="methometR")
 modcmap <- readingCmap(cmap)
 modxmap <- readingXmap(xmap)
 modMolcmap <- modmolcmap(molcmap = modcmap,
   xmapdata = modxmap,
   cntNick = 1,
   cntMeth=2)
```
##Mapping Nick position on the contigs to the refference


```{r eval=FALSE}
refcontigXmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1_sv/merged_smaps/ContigRef.xmap", package="methometR")
 refCmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1_sv/merged_smaps/hg19ref_r.cmap", package="methometR")
 Contigqcmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1_sv/merged_smaps/SampContig_q.cmap", package="methometR")
 refxmapdat <- readingXmap(refcontigXmap)
 refcmapdat <- readingCmap(refcmap)
 contigcmapdat <- readingCmap(Contigqcmap)
 nickRef<-nickReference(refxmap = refxmapdat, refcmap = refcmapdat, 
     contigcmap = contigcmapdat, contigID = 6701, 
    returnMethod =c("dataFrame"),  chrom = 4)
```

##Mapping Nick position on the contigs to the refference


#References


