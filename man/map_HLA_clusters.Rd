% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scrHLAmatrix.R
\name{map_HLA_clusters}
\alias{map_HLA_clusters}
\title{Mapping the HLA clusters found by the 'HLA_clusters()' function back into the scrHLAtag count files}
\usage{
map_HLA_clusters(
  reads.list,
  cluster_coordinates,
  CB_rev_com = FALSE,
  parallelize = FALSE
)
}
\arguments{
\item{reads.list}{is a scrHLAtag count file (or a list of scrHLAtag count files) including columns for CB, UMI, and HLA alleles (https://github.com/furlan-lab/scrHLAtag).}

\item{cluster_coordinates}{is the UMAP coordinates dataframe with HLA clustering information (found by the 'HLA_clusters()' function) from which clusters are extracted and mapped into the scrHLAtag count files by matching Cell Barcodes. Currently the barcode format supported is: SAMPLE_AATGCTTGGTCCATTA-1}

\item{CB_rev_com}{is a logical, called TRUE if the need to obtained the reverse complement of Cell Barcodes (CBs) is desired; default is FALSE.}

\item{parallelize}{is a logical, called TRUE if using parallel processing (multi-threading) is desired; default is FALSE.}
}
\value{
a large list containing scrHLAtag count file(s) including columns for CB, UMI, and HLA alleles, with the addition of HLA_clusters
}
\description{
Mapping the HLA clusters found by the 'HLA_clusters()' function back into the scrHLAtag count files
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
map_HLA_clusters(reads.list = cts, k = 2, cluster_coordinates = UMAP_dataframe_from_HLA_clusters_function)
}
