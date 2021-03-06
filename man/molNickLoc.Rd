% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methylationAnalysis_Step1.r
\name{molNickLoc}
\alias{molNickLoc}
\title{Mapping nick label location from contig to molecule}
\usage{
molNickLoc(nickrefLoc, molcmap, xmapdata, contigID, cntNick, cntMeth)
}
\arguments{
\item{nickrefLoc}{character. Path to Bionano compressed file.}

\item{molcmap}{character. molecule cmap data.}

\item{xmapdata}{character. contig-molecule xmap data.}

\item{contigID}{character. query contig ID.}

\item{cntNick}{integer. Nick label value.}

\item{cntMeth}{integer. Meth label value.}
}
\value{
dataframe containing the molecule file, with the nick position
}
\description{
Mapping nick label location from contig to molecule
}
\examples{
refcontigXmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1_sv/merged_smaps/ContigRef.xmap", package="methometR")
refCmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1_sv/merged_smaps/hg19ref_r.cmap", package="methometR")
Contigqcmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1_sv/merged_smaps/SampContig_q.cmap", package="methometR")
refxmapdat <- readingXmap(refcontigXmap)
refcmapdat <- readingCmap(refCmap)
contigcmapdat <- readingCmap(Contigqcmap)
nickRef<-nickReference(refxmap = refxmapdat, refcmap = refcmapdat, 
    contigcmap = contigcmapdat, contigID = 6701, 
   returnMethod =c("dataFrame"),  chrom = 4)
xmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1/merged_smaps/SampContigMolecule.xmap", package="methometR")
molcmap <- system.file("extdata",  "Cmapdir/output/contigs/exp_refineFinal1/merged_smaps/SampMolecule_q.cmap", package="methometR")
contigcmap <- system.file("extdata",  "Cmapdir/output/contigs/exp_refineFinal1/merged_smaps/SampContig_r.cmap", package="methometR")
contigID = 6701; cntNick = 1; cntMeth=2
modxmap <- readingXmap(xmap)
cmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1/merged_smaps/SampMolecule_q.cmap", package="methometR")
modcmap <- readingCmap(cmap) 
modMolcmap <- modmolcmap(molcmap = modcmap,
  xmapdata = modxmap,
  cntNick = 1,
  cntMeth=2)
contigmodcmap <- readingXmap(contigcmap)
modmolcmapNick <-molNickLoc(nickrefLoc = nickRef,
   molcmap = modMolcmap,
   xmapdata = modxmap, 
   contigID = contigID, cntNick = cntNick, 
   cntMeth = cntMeth)
}
