% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scrHLAmatrix.R
\name{HLA_Matrix}
\alias{HLA_Matrix}
\title{Deduping and Correcting scrHLAtag counts and Creating Seurat-compatible Matrices}
\usage{
HLA_Matrix(
  reads,
  seu,
  hla_recip = character(),
  hla_donor = character(),
  QC_mm2 = TRUE,
  res_conflict_per_gene = TRUE,
  LD_correct = TRUE,
  remove_alleles = character(),
  s1_belowmax = 0.75,
  AS_belowmax = 0.85,
  NM_thresh = 15,
  de_thresh = 0.015,
  parallelize = FALSE,
  CB_rev_com = FALSE,
  return_stats = FALSE,
  Ct = 0
)
}
\arguments{
\item{reads}{is the scrHLAtag count file including columns for CB, UMI, and HLA alleles (https://github.com/furlan-lab/scrHLAtag).}

\item{seu}{is the Seurat object associated with the scrHLAtag count file (https://satijalab.org/seurat/index.html).}

\item{hla_recip}{is a character list of recipient-specific HLA alleles if known; default is an empty character vector.}

\item{hla_donor}{is a character list of donor-specific HLA alleles if known; default is an empty character vector.}

\item{QC_mm2}{is a logical, called TRUE if removing low quality reads based on minimap2 tags is desired.}

\item{res_conflict_per_gene}{is a logical, called TRUE if resolving per-HLA genotype conflicts is desired, with the assumption that each Cell can have no more than 2 alleles of the same HLA gene.}

\item{LD_correct}{is a logical, called TRUE if Linkage Diseqilibrium (LD) correction in the HLA-DR locus is desired, with the assumption of a very strong LD between certain DRB1 allele families and the DRB2, DRB3, DRB4, DRB5, DRB6, DRB7, DRB8, and DRB9 loci.}

\item{remove_alleles}{is a character list of HLA alleles to be removed from the count file if desired; default is an empty character vector.}

\item{s1_belowmax}{is a proportion (0 to 1) of the maximum value (best quality) of the minimap2 's1' tag above which the quality of the read is acceptable; default at 0.75 of the max s1 score.}

\item{AS_belowmax}{is a proportion (0 to 1) of the maximum value (best quality) of the minimap2 'AS' tag above which the quality of the read is acceptable; default at 0.85 of the max AS score.}

\item{NM_thresh}{is the number of mismatches and gaps in the minimap2 alignment at or below which the quality of the read is acceptable; default is 15.}

\item{de_thresh}{is the gap-compressed per-base sequence divergence in the minimap2 alignment at or below which the quality of the read is acceptable; the number is between 0 and 1, and default is 0.015.}

\item{parallelize}{is a logical, called TRUE if using parallel processing (multi-threading) is desired; default is FALSE.}

\item{CB_rev_com}{is a logical, called TRUE if the need to obtained the reverse complement of Cell Barcodes (CBs) is desired; default is FALSE.}

\item{return_stats}{is a logical, when TRUE returns step-by-step read statistics and UMI duplication rate in a list of dataframes and plot, in addition to the Seurat-compatible count matrix; will require additional computations which may noticeably slow down the function; defualt is FALSE.}

\item{Ct}{is the count threshold for the PCR copies of UMIs to retain; default is 0.}
}
\value{
an Assay type matrix
}
\description{
Deduping and Correcting scrHLAtag counts and Creating Seurat-compatible Matrices
}
\examples{
samples <- c("AML_101_BM", "AML_101_34")
mol_info <- c("molecule_info_gene.txt.gz", "molecule_info_mRNA.txt.gz")
for (i in 1:length(mol_info)){
  dl<-lapply(samples, function(sample){
    d<-read.table(file.path("path/to/scrHLAtag/out/files", sample,
                            mol_info[i]), header = F, sep=" ", fill = T) 
    d$V1<-paste0(sample, "_", d$V1, "-1")
    colnames(d)<-c("name","CB", "nb", "UMI", "gene", "query_len","start", "mapq", "cigar", "NM", "AS", "s1", "de", "seq")
    d$samp <- sample
    d
  })
  ctsu<-do.call(rbind,dl)
  rm(dl)
  cts[[str_sub(strsplit(mol_info[i], "\\\\.")[[1]][1], start= -4)]] <- ctsu
  rm(ctsu)
}
HLA_Matrix(reads = cts[["mRNA"]], seu = your_Seurat_obj, hla_recip = c("A*24:02:01", "DRB1*04:01:01", "DRB4*01:03:02"), hla_donor = c("A*33:03:01", "B*42:01:01"))
}
