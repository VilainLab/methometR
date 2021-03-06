% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methylationAnalysis_Step1.r
\name{molMethylRefLoc}
\alias{molMethylRefLoc}
\title{Liftover of molecule methylation and nick  contig position to reference coordinates}
\usage{
molMethylRefLoc(
  molcmap,
  outputType = c("dataframe", "text"),
  cntMeth = 1,
  cntNick = 2
)
}
\arguments{
\item{molcmap}{character. molecule cmap data.}

\item{outputType}{character. Output method. dataframe or text.}

\item{cntMeth}{integer. Meth label value.}

\item{cntNick}{integer. Nick label value.}
}
\value{
dataframe containing the methyl label positions mapped to the reference.
}
\description{
Liftover of molecule methylation and nick  contig position to reference coordinates
}
\examples{
xmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1/merged_smaps/SampContigMolecule.xmap", package="methometR")
cmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1/merged_smaps/SampMolecule_q.cmap", package="methometR")
modcmap <- readingCmap(cmap)
modxmap <- readingXmap(xmap)
modMolcmap <- modmolcmap(molcmap = modcmap,
  xmapdata = modxmap,
  cntNick = 1,
  cntMeth=2)
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
contigmodcmap <- readingXmap(contigcmap)
modmolcmapNick <-molNickLoc(nickrefLoc = nickRef,
   molcmap = modMolcmap,
   xmapdata = modxmap, 
   contigID = contigID, cntNick = cntNick, 
   cntMeth = cntMeth)
MethylRefLoc<-molMethylRefLoc(molcmap = modmolcmapNick,
   outputType = c("dataframe"),
   cntNick = cntNick, 
   cntMeth = cntMeth)
}
