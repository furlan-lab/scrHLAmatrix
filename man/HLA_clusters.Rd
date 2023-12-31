% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scrHLAmatrix.R
\name{HLA_clusters}
\alias{HLA_clusters}
\title{Getting raw scrHLAtag counts and analyzing distribution of HLA alleles per Cell Barcodes in UMAP space}
\usage{
HLA_clusters(
  reads,
  k = 2,
  seu = NULL,
  CB_rev_com = FALSE,
  geno_metadata_id = NULL,
  hla_with_counts_above = 0,
  CBs_with_counts_above = 0,
  match_CB_with_seu = TRUE,
  umap_first_n_PCs = 25,
  QC_mm2 = TRUE,
  s1_belowmax = 0.75,
  AS_belowmax = 0.85,
  NM_thresh = 15,
  de_thresh = 0.015,
  parallelize = FALSE,
  pt_size = 0.5,
  ...
)
}
\arguments{
\item{reads}{is the scrHLAtag count file including columns for CB, UMI, and HLA alleles (https://github.com/furlan-lab/scrHLAtag).}

\item{k}{is the number of clusters to partition the UMAP space into, e.g. the number of entities or genotypes you 'think' there might be in your captured sample.}

\item{seu}{is the Seurat object associated with the scrHLAtag count file (https://satijalab.org/seurat/index.html).}

\item{CB_rev_com}{is a logical, called TRUE if the need to obtained the reverse complement of Cell Barcodes (CBs) is desired; default is FALSE.}

\item{geno_metadata_id}{is a character, the column ID of the Seurat metadata designated to distinguish genotypes, if this information is available. 'NULL' by default or when genotyping information is not available.}

\item{hla_with_counts_above}{is the number of total reads accross CBs at or above which an HLA allele is retained in the matrix.}

\item{CBs_with_counts_above}{is the number of total reads accross HLA alleles at or above which a CB is retained in the matrix. Note: 'princomp()' can only be used with at least as many units (CBs) as variables (HLAs), thus the function will make sure that number of CBs is equal or more than available HLA alleles in the matrix.}

\item{match_CB_with_seu}{is a logical, called TRUE if filtering CBs in the scrHLAtag count file with matching ones in the Seurat object is desired.}

\item{umap_first_n_PCs}{is a number representing the first 'n' principal components to input into the umap function. Note: current version only supports default parameters of the 'umap()' function (more info on: https://cran.r-project.org/web/packages/uwot/index.html)}

\item{QC_mm2}{is a logical, called TRUE if removing low quality reads based on minimap2 tags is desired.}

\item{s1_belowmax}{is a proportion (0 to 1) of the maximum value (best quality) of the minimap2 's1' tag above which the quality of the read is acceptable; default at 0.75 of the max s1 score.}

\item{AS_belowmax}{is a proportion (0 to 1) of the maximum value (best quality) of the minimap2 'AS' tag above which the quality of the read is acceptable; default at 0.85 of the max AS score.}

\item{NM_thresh}{is the number of mismatches and gaps in the minimap2 alignment at or below which the quality of the read is acceptable; default is 15.}

\item{de_thresh}{is the gap-compressed per-base sequence divergence in the minimap2 alignment at or below which the quality of the read is acceptable; the number is between 0 and 1, and default is 0.015.}

\item{parallelize}{is a logical, called TRUE if using parallel processing (multi-threading) is desired; default is FALSE.}

\item{pt_size}{is a number, the size of the geometric point displayed by ggplot2.}

\item{...}{arguments passed onto uwot::umap()}
}
\value{
a large list containing DataFrame with UMAP coordinates and ggplot of HLA clusters
}
\description{
Getting raw scrHLAtag counts and analyzing distribution of HLA alleles per Cell Barcodes in UMAP space
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
HLA_clusters(reads = cts[["mRNA"]], k = 2, seu = your_Seurat_Obj, geno_metadata_id = "geno", hla_with_counts_above = 5, CBs_with_counts_above = 35)
}
