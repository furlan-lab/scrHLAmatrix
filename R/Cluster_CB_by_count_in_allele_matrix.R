#' Getting raw scrHLAtag counts and analyzing distribution of HLA alleles per Cell Barcodes in UMAP space
#' 
#' @param reads  is the scrHLAtag count file including columns for CB, UMI, and HLA alleles (https://github.com/furlan-lab/scrHLAtag).
#' @param k  is the number of clusters to partition the UMAP space into, e.g. the number of entities or genotypes you 'think' there might be in your captured sample.
#' @param seu  is the Seurat object associated with the scrHLAtag count file (https://satijalab.org/seurat/index.html).
#' @param CB_rev_com  is a logical, called \code{TRUE} if the need to obtained the reverse complement of Cell Barcodes (CBs) is desired; default is \code{FALSE}. 
#' @param geno_metadata_id  is a character, the column ID of the Seurat metadata designated to distinguish genotypes, if this information is available. \code{NULL} by default or when genotyping information is not available. 
#' @param hla_with_counts_above  is the number of total reads accross CBs at or above which an HLA allele is retained in the matrix.
#' @param CBs_with_counts_above  is the number of total reads accross HLA alleles at or above which a CB is retained in the matrix. Note: \code{stats::princomp()} can only be used with at least as many units (CBs) as variables (HLAs), thus the function will make sure that number of CBs is equal or more than available HLA alleles in the matrix.
#' @param match_CB_with_seu  is a logical, called \code{TRUE} if filtering CBs in the scrHLAtag count file with matching ones in the Seurat object is desired. 
#' @param method  is the graph-based clustering method to be used for partitioning cells based on their HLA count patterns. The choice is between a Connectivity-based method: \code{"hclust"} or \code{"umap_hclust"}, a Centroid-based method: \code{"kmeans"}, or a Distribution based method: \code{"mclust"} (for Gaussian Mixture Model). We found \code{"hclust"} had the best ability to separate allogeneic entities so we set it as Default. In some cases, applying hierarchical clustering directly on UMAP coordinates gives good allogeneic entity separation so we provide the \code{"umap_hclust"} option.
#' @param QC_mm2  is a logical, called \code{TRUE} if removing low quality reads based on minimap2 tags is desired.
#' @param s1_percent_pass_score  is a percentage (\code{0} to \code{100}) cuttoff from the maximum score (best quality) of the minimap2 's1' tag, which a read needs to acheive to pass as acceptable; default at \code{80} and becomes less inclusive if value increases.
#' @param AS_percent_pass_score  is a percentage (\code{0} to \code{100}) cuttoff from the maximum score (best quality) of the minimap2 'AS' tag, which a read needs to acheive to pass as acceptable; default at \code{80} and becomes less inclusive if value increases.
#' @param NM_thresh  is the number of mismatches and gaps in the minimap2 alignment at or below which the quality of the read is acceptable; default is \code{15}.
#' @param de_thresh  is the gap-compressed per-base sequence divergence in the minimap2 alignment at or below which the quality of the read is acceptable; the number is between \code{0} and \code{1}, and default is \code{0.01}.
#' @param hclust_algorithm  applies to \code{stats::hclust()}; is the agglomeration algorithm to be used. Values include \code{"ward.D"}, \code{"ward.D2"}, \code{"single"}, \code{"complete"}, \code{"average"}, \code{"mcquitty"}, \code{"median"}, or \code{"centroid"}; for more information: \code{?stats::hclust}.
#' @param kmeans_algorithm  applies to \code{stats::kmeans()}; is the k-means algorithm to be used. Values include \code{"Hartigan-Wong"}, \code{"Lloyd"}, \code{"Forgy"}, or \code{"MacQueen"}; for more information: \code{?stats::kmeans}.
#' @param parallelize  is a logical, called \code{TRUE} if using parallel processing (multi-threading) is desired; default is \code{FALSE}.
#' @param pt_size  is a number, the size of the geometric point displayed by ggplot2. 
#' @param return_heavy  is a logical, if \code{TRUE} it also returns the now processed scrHLAtag count file (minimap2 QCed, CB reverse comp'ed, etc..) but the returned object is significantly heavier; default is \code{FALSE}. 
#' @param ...  arguments passed onto \code{uwot::umap()}.
#' @import stringr
#' @import pbmcapply
#' @import parallel
#' @import Matrix
#' @import magrittr
#' @import htmltools
#' @import ggplot2
#' @import Seurat
#' @import dplyr
#' @return a large list containing DataFrame with UMAP coordinates and ggplot of HLA clusters
#' @examples
#' samples <- c("AML_101_BM", "AML_101_34")
#' mol_info <- c("molecule_info_gene.txt.gz", "molecule_info_mRNA.txt.gz")
#' cts <- list()
#' for (i in 1:length(mol_info)){
#'   dl<-lapply(samples, function(sample){
#'     d<-read.table(file.path("path/to/scrHLAtag/out/files", sample,
#'                             mol_info[i]), header = F, sep=" ", fill = T) 
#'     d$V1<-paste0(sample, "_", d$V1, "-1")
#'     colnames(d)<-c("name","CB", "nb", "UMI", "gene", "query_len","start", "mapq", "cigar", "NM", "AS", "s1", "de", "seq")
#'     d$samp <- sample
#'     d
#'   })
#'   ctsu<-do.call(rbind,dl)
#'   rm(dl)
#'   cts[[str_sub(strsplit(mol_info[i], "\\.")[[1]][1], start= -4)]] <- ctsu
#'   rm(ctsu)
#' }
#' HLA_umap <- HLA_clusters(reads = cts[["mRNA"]], k = 2, seu = your_Seurat_Obj, geno_metadata_id = "geno", hla_with_counts_above = 5, CBs_with_counts_above = 35)
#' @export

HLA_clusters <- function(reads, k = 2, seu = NULL, CB_rev_com = FALSE, geno_metadata_id = NULL, 
                         hla_with_counts_above = 0, CBs_with_counts_above = 50, match_CB_with_seu = TRUE, method = "hclust", 
                         QC_mm2 = TRUE, s1_percent_pass_score = 80, AS_percent_pass_score = 80, NM_thresh = 15, de_thresh = 0.01, 
                         hclust_algorithm = "complete", kmeans_algorithm = "Hartigan-Wong",
                         parallelize = FALSE, pt_size = 0.5, return_heavy = FALSE, ...) {
  if (!requireNamespace("mclust", quietly = TRUE)) { stop("Package 'mclust' needed for this function to work. Please install it.", call. = FALSE) }
  if (!"package:mclust" %in% search()) {library(mclust)}
  ## parallelize
  if (parallelize) {
    multi_thread <- parallel::detectCores()
    message(cat("\nMulti-threading! Available cores: ", parallel::detectCores(), sep = ""))
  } else {
    multi_thread <- 1
  }
  ## Remove low quality reads based on minimap2 tags
  if (QC_mm2) {
    message(cat("\nRemoving low quality reads based on minimap2 tags"))
    reads <- split(data.table::setDT(reads), by = "gene")
    reads <- parallel::mclapply(reads, data.table::setDF, mc.cores = multi_thread)
    reads <- pbmcapply::pbmclapply(reads, function(df){df[df$s1 > (s1_percent_pass_score/100)*max(df$s1) & df$AS > (AS_percent_pass_score/100)*max(df$AS) & df$NM <= NM_thresh & df$de <= de_thresh,]}, mc.cores = multi_thread)
    reads <- data.table::rbindlist(reads)
    row.names(reads)<-NULL
  } 
  ## Reverse Complement the CB
  if (CB_rev_com) {
    message(cat("\nConverting Cell Barcodes to their reverse complements"))
    # reads$CB <- pbmcapply::pbmclapply(reads$CB, function(x) as.character(Biostrings::reverseComplement(DNAString(x))), mc.cores = multi_thread) %>% unlist() # slow
    reads$CB <- pbmcapply::pbmclapply(reads$CB, function(x) intToUtf8(rev(utf8ToInt(chartr('ATGC', 'TACG', x)))), mc.cores = multi_thread) %>% unlist()        # fast
  }  
  ## the 3-field level resolution of HLA is actually 2-field for MICA and MICB. must fix this glitch:
  special <- "[-_*|?.+$^]"
  reads[c("hla", "leftover")] <- stringr::str_split_fixed(reads$gene, special, 2)
  reads$leftover <- NULL
  reads$a <- sapply(reads$gene, function(x) strsplit(x, ":")[[1]][1])
  reads$b <- sapply(reads$gene, function(x) strsplit(x, ":")[[1]][2])
  reads$c <- sapply(reads$gene, function(x) strsplit(x, ":")[[1]][3])
  reads$a <- ifelse(is.na(reads$b), reads$a, paste0(reads$a, ":", reads$b))
  mic <- table(reads$a) %>% names() %>% grep("MIC", ., value = T)
  for (x in 1:length(mic)){
    reads$a <- ifelse(reads$a %in% mic[x], 
                      paste0(reads$a, ":", table(reads[reads$a == mic[x],]$c) %>% which.max() %>% names()), 
                      reads$a)
  }
  reads$a <- ifelse(endsWith(reads$a, ":"), substr(reads$a, 1, nchar(reads$a)-1), reads$a)
  reads$gene <- ifelse(reads$hla %in% c("MICA", "MICB"), reads$a, reads$gene)
  reads$a <- NULL
  reads$b <- NULL
  reads$c <- NULL
  rm(mic)
  ## Matrix formation
  alleles <- unique(reads$gene) %>% sort()
  reads$seu_barcode <- paste0(reads$samp,"_",reads$CB,"-1")
  message(cat("\nCreating an HLA Count Matrix"))
  reads <- split(data.table::setDT(reads), by = "seu_barcode")
  reads <- parallel::mclapply(reads, data.table::setDF, mc.cores = multi_thread)
  HLA.matrix <- matrix(0, nrow = length(alleles), ncol = length(reads), dimnames = list(alleles, names(reads)))
  pb <- pbmcapply::progressBar(min = 0, max = length(reads), style = "ETA", char = "=")
  for (i in 1:length(reads)) {
    counts <- table(reads[[i]]$gene)
    HLA.matrix[, i] <- counts[alleles]
    setTxtProgressBar(pb, i)
  }
  close(pb)
  ## do we keep or jettison the 'reads'? 
  if (return_heavy) {
    reads <- data.table::rbindlist(reads)
    row.names(reads)<-NULL
  } else {
    reads <- NULL
  }
  HLA.matrix[is.na(HLA.matrix)]<-0
  #HLA.matrix<-Matrix(HLA.matrix,sparse = T)
  ## Matching with Seurat colnames
  if (is.null(seu)) {
    part_HLA<- HLA.matrix
  } else {
    if (class(seu) == "Seurat") {
      message(cat("\nObject of class 'Seurat' detected"))
      message(cat(crayon::green("Note: "), "Currently the Seurat Barcode (i.e. colnames or Cells) supported format is: SAMPLE_AATGCTTGGTCCATTA-1", sep = ""))
      if(match_CB_with_seu) {
        part_HLA<- HLA.matrix[,colnames(HLA.matrix) %in% Cells(seu)]
      } else {
        part_HLA<- HLA.matrix
      }
    } else {
      stop("Single-cell dataset container must be of class 'Seurat'")
    }
  }
  ## removing HLA alleles with low counts overall
  r <- hla_with_counts_above
  part_HLA <- part_HLA[which(rowSums(part_HLA)>=r),]
  ## removing cell barcodes (CBs) with low counts overall
  n <- CBs_with_counts_above
  if (dim(part_HLA[,which(colSums(part_HLA)>=n)])[2] < dim(part_HLA)[1]) {
    part_HLA <- part_HLA[,order(-colSums(part_HLA))[1:dim(part_HLA)[1]]]
  } else {
    part_HLA <- part_HLA[,which(colSums(part_HLA)>=n)]
  }
  message(cat("\nRunning PCA (from 'stats') and UMAP (from 'uwot')"))
  ## normalize by size factor
  part_HLA <- prop.table(part_HLA, margin = 2)
  ## log-transform
  part_HLA <- log10(part_HLA+0.01)
  ## princomp from stats
  pc <- stats::princomp(t(part_HLA))
  pcv<-as.data.frame(pc$scores)
  elbow <- barplot((pc$sdev^2/sum(pc$sdev^2))[1:100])
  umat<-pcv[,1:floor(0.8*ncol(pcv))] %>% as.matrix()
  umapout<-uwot::umap(umat, verbose = TRUE, ...) # fix `batch = TRUE, seed = 1985` (for example) for more consistent umaps.
  colnames(umapout)<-c("umap1", "umap2")
  umapout <- as.data.frame(umapout)
  if (!is.null(seu) & !is.null(geno_metadata_id)){
    umapout<-cbind(umapout, seu@meta.data[match(rownames(umapout), colnames(seu)),])
    g0 <- ggplot(umapout, aes(x=umap1, y=umap2, color=!!sym(geno_metadata_id)))+geom_point(size=pt_size)#+scale_color_manual(values=pals::glasbey())+theme_bw()
  }
  if (!(length(method) == 1L && method %in% c("hclust", "kmeans", "mclust", "umap_hclust"))) {
    message(cat("\nArgument `method` should be one of 'hclust', 'kmeans', 'mclust', or 'umap_hclust'. Using default."))
    method <- "hclust"
  } %>% suppressWarnings()
  if (method == "hclust") {
    message(cat("\nConnectivity-based Clustering: ", crayon::red("Hierarchical Clustering"), " (agglomerative) on PCA space...", sep = ""))
    if (hclust_algorithm == "centroid") {exp <- 2} else {exp <- 1}
    humapout <- stats::hclust(dist(umat[, 1:min(2000, ncol(umat))])^exp, method = hclust_algorithm)
    umapout$hla_clusters <- stats::cutree(humapout, k = k)
  } else if (method == "kmeans") {
    message(cat("\nCentroid-based Clustering: ", crayon::red("k-means Clustering"), " on PCA space...", sep = ""))
    humapout <- stats::kmeans(umat, centers = k, algorithm = kmeans_algorithm)
    umapout$hla_clusters <- humapout$cluster
  } else if (method == "mclust") {
    message(cat("\nDistribution-based Clustering: ", crayon::red("HGaussian Mixture Model Clustering"), " on PCA space...", sep = ""))
    humapout <- mclust::Mclust(umat[, 1:min(50, ncol(umat))], G = k) # beyond 50 cols, the algorithm takes forever
    umapout$hla_clusters <- humapout$classification
  } else if (method == "umap_hclust") {
    message(cat("\nConnectivity-based Clustering: ", crayon::red("Hierarchical Clustering"), " (agglomerative) on UMAP coordinates...", sep = ""))
    if (hclust_method == "centroid") {exp <- 2} else {exp <- 1}
    humapout <- stats::hclust(dist(as.matrix(umapout[,1:2]))^exp, method = hclust_method) 
    umapout$hla_clusters <- stats::cutree(humapout, k = k)
  }
  umapout$hla_clusters <- as.factor(umapout$hla_clusters)
  g <- ggplot(umapout, aes(x=umap1, y=umap2, color=hla_clusters))+geom_point(size=pt_size)#+scale_color_manual(values=pals::glasbey())+theme_bw()
  #message(cat("\nDone!!"))
  if (!is.null(seu) & !is.null(geno_metadata_id)){
    return(list(UMAP_coordinates = umapout, HLA_clusters_on_umap = g, genotype_on_umap = g0, reads = reads, top80_PC = umat))
  } else {
    return(list(UMAP_coordinates = umapout, HLA_clusters_on_umap = g, reads = reads, top80_PC = umat))
  }
}


#' Mapping the HLA clusters found by the 'HLA_clusters()' function back into the scrHLAtag count files
#'
#' @param reads.list  is a scrHLAtag count file (or a list of scrHLAtag count files) including columns for CB, UMI, and HLA alleles (https://github.com/furlan-lab/scrHLAtag).
#' @param cluster_coordinates  is the UMAP coordinates dataframe with HLA clustering information (found by the \code{scrHLAmatrix::HLA_clusters()} function) from which clusters are extracted and mapped into the scrHLAtag count files by matching Cell Barcodes. Currently the barcode format supported is: SAMPLE_AATGCTTGGTCCATTA-1
#' @param CB_rev_com  is a logical, called \code{TRUE} if the need to obtained the reverse complement of Cell Barcodes (CBs) is desired; default is \code{FALSE}.
#' @param additional_col_id  is a character; if mapping of any additional column ID of the Seurat metadata collected in the \code{cluster_coordinates} dataframe, to the scrHLAtag count files is desired. \code{NULL} by default.
#' @param parallelize  is a logical, called \code{TRUE} if using parallel processing (multi-threading) is desired; default is \code{FALSE}.
#' @import stringr
#' @return a large list containing scrHLAtag count file(s) including columns for CB, UMI, and HLA alleles, with the addition of HLA_clusters
#' @examples
#' samples <- c("AML_101_BM", "AML_101_34")
#' mol_info <- c("molecule_info_gene.txt.gz", "molecule_info_mRNA.txt.gz")
#' cts <- list()
#' for (i in 1:length(mol_info)){
#'   dl<-lapply(samples, function(sample){
#'     d<-read.table(file.path("path/to/scrHLAtag/out/files", sample,
#'                             mol_info[i]), header = F, sep=" ", fill = T) 
#'     d$V1<-paste0(sample, "_", d$V1, "-1")
#'     colnames(d)<-c("name","CB", "nb", "UMI", "gene", "query_len","start", "mapq", "cigar", "NM", "AS", "s1", "de", "seq")
#'     d$samp <- sample
#'     d
#'   })
#'   ctsu<-do.call(rbind,dl)
#'   rm(dl)
#'   cts[[str_sub(strsplit(mol_info[i], "\\.")[[1]][1], start= -4)]] <- ctsu
#'   rm(ctsu)
#' }
#' cts <- map_HLA_clusters(reads.list = cts, cluster_coordinates = UMAP_dataframe_from_HLA_clusters_function, additional_col_id = "geno")
#' @export

map_HLA_clusters <- function(reads.list, cluster_coordinates, CB_rev_com = FALSE, additional_col_id = NULL, parallelize = FALSE) {
  ## parallelize
  if (parallelize) {
    multi_thread <- parallel::detectCores()
    # message(cat("\nMulti-threading! Available cores: ", parallel::detectCores(), sep = ""))
  } else {
    multi_thread <- 1
  }  
  if (any(class(reads.list) %in% c("data.frame", "data.table"))) {
    reads.list <- list(reads.list=reads.list)
  } else {
    if (!(sapply(reads.list, function(x) class(x)) %in% list("data.frame", "data.table", c("data.table", "data.frame"), c("data.frame", "data.table")) %>% all()) | class(reads.list) != "list") {
      stop("The class of 'reads.list' should be either a 'data.frame', a 'data.table', or a 'list' of the formers")
    }
  }
  if (CB_rev_com) {
    for (i in 1:length(reads.list)) {
      if (!is.null(additional_col_id)) {
        reads.list[[i]][, additional_col_id] <-  NA
        reads.list[[i]][, additional_col_id] <- cluster_coordinates[, additional_col_id][
          match(paste0(reads.list[[i]]$samp,"_", parallel::mclapply(reads.list[[i]]$CB, function(x) intToUtf8(rev(utf8ToInt(chartr('ATGC', 'TACG', x)))), mc.cores = multi_thread) %>% unlist(),"-1"), 
                rownames(cluster_coordinates))
        ]
      }
      reads.list[[i]]$hla_clusters <- NA
      reads.list[[i]]$hla_clusters <- cluster_coordinates$hla_clusters[
        match(paste0(reads.list[[i]]$samp,"_", parallel::mclapply(reads.list[[i]]$CB, function(x) intToUtf8(rev(utf8ToInt(chartr('ATGC', 'TACG', x)))), mc.cores = multi_thread) %>% unlist(),"-1"), 
              rownames(cluster_coordinates))
      ]
    }
  } else {
    for (i in 1:length(reads.list)) {
      if (!is.null(additional_col_id)) {
        reads.list[[i]][, additional_col_id] <-  NA
        reads.list[[i]][, additional_col_id] <- cluster_coordinates[, additional_col_id][
          match(paste0(reads.list[[i]]$samp,"_",reads.list[[i]]$CB,"-1"), 
                rownames(cluster_coordinates))
        ]
      }
      reads.list[[i]]$hla_clusters <- NA
      reads.list[[i]]$hla_clusters <- cluster_coordinates$hla_clusters[
        match(paste0(reads.list[[i]]$samp,"_",reads.list[[i]]$CB,"-1"), 
              rownames(cluster_coordinates))
      ]
    }
  }
  if (length(reads.list) == 1) {
    reads.list <- reads.list[[1]]
    return(reads.list)
  } else {
    return(reads.list)
  }
}