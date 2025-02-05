#' Getting raw scrHLAtag counts and analyzing distribution of HLA alleles per Cell Barcodes in UMAP space
#' 
#' @param reads  is the scrHLAtag count file including columns for CB, UMI, and HLA alleles (\url{https://github.com/furlan-lab/scrHLAtag}).
#' @param k  can be \code{NULL} or a fixed number of cluster counts to partition the datapoints into, e.g. the number of entities or genotypes you \emph{think} there might be in your captured sample. If \code{NULL}, each clustering method will automatically determine \code{k} clusters on its own (will not work for \code{"hclust"} and \code{"kmeans"}, which need predefined \code{k}s).
#' @param seu  is the Seurat object associated with the scrHLAtag count file (\url{https://satijalab.org/seurat/index.html}).
#' @param CB_rev_com  logical, called \code{TRUE} if the need to obtain the reverse complement of Cell Barcodes (CBs) is desired; default is \code{FALSE}. 
#' @param geno_metadata_id  a character, the column ID of the Seurat metadata designated to distinguish genotypes, if this information is available. \code{NULL} by default or when genotyping information is not available. 
#' @param hla_with_counts_above  number of total reads accross CBs at or above which an HLA allele is retained in the matrix.
#' @param CBs_with_counts_above  number of total reads accross HLA alleles at or above which a CB is retained in the matrix. Note: \code{stats::princomp()} can only be used with at least as many units (CBs) as variables (HLAs), thus the function will make sure that number of CBs is equal or more than available HLA alleles in the matrix.
#' @param match_CB_with_seu  logical, called \code{TRUE} if filtering CBs in the scrHLAtag count file with matching ones in the Seurat object is desired. 
#' @param method  the name(s) of the graph-based clustering method(s) to be used for partitioning cells based on their HLA count patterns. The choice is between one or a combination of a Community structure detection method: \code{"leiden"}, a Density-based method: \code{"dbscan"}, a Connectivity-based method: \code{"hclust"}, a Centroid-based method: \code{"kmeans"}, and a Distribution-based method: \code{"gmm"} (for Gaussian Mixture Model). The methods are run with their respective Default parameters. Some of those methods may predict \emph{true} allogeneic entities with better accuracy than others; as we cannot know a priori which is the best method, we propose the method: \code{"consensus"}, which groups cells in the same cluster if they agree on membership in a majority of methods, otherwise they are unclassified (\code{NA}s). Note: \code{"consensus"} is equivalent to \code{c("leiden", "hclust", "kmeans", "gmm", "dbscan", "consensus")}.
#' @param n_PCs  the number of top principal components to retain in downstream clustering and umap analyses; default is \code{50} or the top 80 percent of PCs, whichever is smaller.
#' @param dbscan_minPts  only works for the  \code{"dbscan"} method: number of minimum points required in the epsilon neighborhood radius (\code{eps}) of core points. While the other methods require 1 parameter (e.g., \code{k}), \code{"dbscan"} requires 2: \code{eps} and \code{minPts}. To acheive desired \code{k} clusters, a range of \code{eps} parameter is tested against a fixed \code{minPts}, which is provided here. Default at \code{30}, but can be adjusted higher or lower depending on how small and 'clumped' an allogeneic entity is suspected to be. 
#' @param QC_mm2  logical, called \code{TRUE} if removing low quality reads based on minimap2 tags is desired.
#' @param s1_percent_pass_score  percentage, \code{0} to \code{100}, cuttoff from the maximum score (best quality) of the minimap2 's1' tag, which a read needs to acheive to pass as acceptable; default at \code{80} and becomes less inclusive if value increases.
#' @param AS_percent_pass_score  percentage, \code{0} to \code{100}, cuttoff from the maximum score (best quality) of the minimap2 'AS' tag, which a read needs to acheive to pass as acceptable; default at \code{80} and becomes less inclusive if value increases.
#' @param NM_thresh  the number of mismatches and gaps in the minimap2 alignment at or below which the quality of the read is acceptable; default is \code{15}.
#' @param de_thresh  the gap-compressed per-base sequence divergence in the minimap2 alignment at or below which the quality of the read is acceptable; the number is between \code{0} and \code{1}, and default is \code{0.01}.
#' @param parallelize  logical, called \code{TRUE} if using parallel processing (multi-threading) is desired; default is \code{FALSE}.
#' @param pt_size   a number, the size of the geometric point displayed by ggplot2. 
#' @param return_heavy   logical, if \code{TRUE} it also returns the now processed scrHLAtag count file (minimap2 QCed, CB reverse comp'ed, etc..) which is usually a heavy object memory-wise; default is \code{FALSE}. 
#' @param seed   numeric (or \code{NULL}), to set seed (or not) in the environment for reproducibility
#' @param suppress_plots  called \code{TRUE} to suppress plots from appearing while running the function.
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
#' dirs_path <- "path/to/scrHLAtag/out/files"
#' dirs<-list.dirs(path=dirs_path, full.names = T, recursive = F)
#' dirs<- lapply(dirs, list.dirs, recursive = F) %>% unlist
#' dirs<- lapply(dirs, dir, pattern = "unguided_hla_align_corrected", recursive = F, full.names = T) %>% unlist
#' dirnames <- c("AML_101_BM", "AML_101_34", "TN_BM", "TN_34") # this is how the samples were organized in the directories
#' ## Load the counts files
#' cts <- HLA_load(directories = dirs, dir_names = dirnames, seu = your_Seurat_obj)
#' ## Process those count files
#' HLA_umap <- HLA_clusters(reads = cts[["mRNA"]], k = 2, seu = your_Seurat_Obj, geno_metadata_id = "geno", hla_with_counts_above = 5, CBs_with_counts_above = 35)
#' @export

HLA_clusters <- function(reads, k = 2, seu = NULL, CB_rev_com = FALSE, geno_metadata_id = NULL, 
                         hla_with_counts_above = 0, CBs_with_counts_above = 25, match_CB_with_seu = TRUE, 
                         method = "consensus", n_PCs = 50, dbscan_minPts = 30,
                         QC_mm2 = TRUE, s1_percent_pass_score = 80, AS_percent_pass_score = 80, NM_thresh = 15, de_thresh = 0.01, 
                         parallelize = FALSE, pt_size = 0.5, return_heavy = FALSE, seed = NULL, suppress_plots = FALSE, ...) {
  if (!requireNamespace("mclust", quietly = TRUE)) { stop("Package 'mclust' needed for this function to work. Please install it.", call. = FALSE) }
  if (!requireNamespace("dbscan", quietly = TRUE)) { stop("Package 'dbscan' needed for this function to work. Please install it.", call. = FALSE) }
  if (!requireNamespace("FNN", quietly = TRUE)) { stop("Package 'FNN' needed for this function to work. Please install it.", call. = FALSE) }
  if (!all(sapply(c("leiden", "igraph", "reticulate"), requireNamespace, quietly = TRUE))) { stop("Install 'leiden' and associated 'igraph' and 'reticulate' pakages for this function to work.\nMoreover, install the python dependencies in R (if you haven't already):\n  reticulate::install_python(version = '<version>') #example '3.8.2'\n  reticulate::py_install('python-igraph')\n  reticulate::py_install('leidenalg', forge = TRUE)\n  reticulate::py_config()", call. = FALSE) }
  if (!is.null(k)) if (abs(k)!=as.integer(k) | k==0) { stop("'k' must be a whole positive number or 'NULL'.", call. = FALSE)}
  if (is.null(k) & !any(method %in% c("leiden", "dbscan", "gmm"))) { stop("'k' cannot be 'NULL' while using methods 'hclust', 'kmeans', or 'consensus'.", call. = FALSE)}
  if (!"package:mclust" %in% search()) {suppressPackageStartupMessages({library(mclust)})}
  if (!is.null(seed)) set.seed(seed)
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
  reads$seu_barcode <- stringr::str_c(reads$samp, reads$id_cb_separator, reads$CB, reads$id_cb_suffix)
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
    if ("Seurat" %in% class(seu)) {
      if (match_CB_with_seu) {
        part_HLA<- HLA.matrix[,colnames(HLA.matrix) %in% Cells(seu)]
        if (ncol(part_HLA) == 0) { stop("Seurat Barcodes did not match any of the CBs in the scrHLAtag counts object. Make sure `seu` and `reads` objects are related.", call. = FALSE) }
      } else {
        part_HLA<- HLA.matrix
      }
    } else {
      stop("Single-cell dataset container (in argument 'seu') must be of class 'Seurat'", call. = FALSE)
    }
  }
  ## removing HLA alleles with low counts overall
  r <- hla_with_counts_above
  if (dim(part_HLA[which(rowSums(part_HLA)>=r), , drop = F])[1] < 2) { # must keep at least 2 rows in the matrix
    part_HLA <- part_HLA[order(-rowSums(part_HLA))[1:2], , drop = F]
  } else {
    part_HLA <- part_HLA[which(rowSums(part_HLA)>=r), , drop = F]
  }
  ## removing cell barcodes (CBs) with low counts overall
  n <- CBs_with_counts_above
  if (dim(part_HLA[,which(colSums(part_HLA)>=n), drop = F])[2] < dim(part_HLA)[1]) {
    part_HLA <- part_HLA[,order(-colSums(part_HLA))[1:dim(part_HLA)[1]], drop = F]
  } else {
    part_HLA <- part_HLA[,which(colSums(part_HLA)>=n), drop = F]
  }
  message(cat("\nRunning PCA (from 'stats') and UMAP (from 'uwot')"))
  ## normalize by size factor
  part_HLA <- prop.table(part_HLA, margin = 2)
  ## log-transform
  part_HLA <- log10(part_HLA+0.01)
  ## princomp from stats
  pc <- stats::princomp(t(part_HLA))
  pcv<-as.data.frame(pc$scores)
  if (!suppress_plots) elbow <- barplot((pc$sdev^2/sum(pc$sdev^2))[1:100])
  umat<-pcv[,1:floor(min(n_PCs, 0.8*ncol(pcv)))] %>% as.matrix()
  umapout<-uwot::umap(umat, verbose = TRUE, batch = TRUE, seed = 1985, ...) #batch and seed are fixed to promote consistency and repeatability
  colnames(umapout)<-c("umap1", "umap2")
  umapout <- as.data.frame(umapout)
  if (!is.null(seu) & !is.null(geno_metadata_id)){
    umapout<-cbind(umapout, seu@meta.data[match(rownames(umapout), colnames(seu)),])
    g0 <- ggplot(umapout, aes(x=umap1, y=umap2, color=!!sym(geno_metadata_id)))+geom_point(size=pt_size)#+scale_color_manual(values=pals::glasbey())+theme_bw()
  } else {
    g0 <- NULL
  }
  if (!any(method %in% c("leiden", "hclust", "kmeans", "gmm", "dbscan", "consensus"))) {
    method <- "consensus"
    message(cat("\nArgument `method` should be one or a combination of 'leiden' 'dbscan', 'hclust', 'kmeans', 'gmm', or 'consensus'. Defaulting to: '", method, "'", sep = ""))
  } %>% suppressWarnings()
  message(cat("\nGraph-based Clustering:"))
  if (any(c("leiden", "consensus") %in% method)) {
    message(cat(crayon::red(format(Sys.time(), "%H:%M:%S"), "- Community detection (leiden) on PCA space"), sep = ""))
    # Create a k-nearest neighbors graph, then convert to igraph
    make_knn_graph <- function(data, k) {
      dist_matrix <- as.matrix(stats::dist(data))
      knn_graph <- matrix(0, nrow = nrow(data), ncol = nrow(data))
      for (i in 1:nrow(data)) {
        neighbors <- order(dist_matrix[i, ])[2:(k + 1)]
        knn_graph[i, neighbors] <- 1
        knn_graph[neighbors, i] <- 1
      }
      return(knn_graph)
    }
    knn_graph <- make_knn_graph(as.matrix(umat[, 1:min(50, ncol(umat))]), 15)
    ig <- igraph::graph_from_adjacency_matrix(knn_graph, mode = "undirected")
    if (is.null(k)) {
      leiden_res <- suppressWarnings(leiden::leiden(ig, resolution_parameter = 1, seed = 1985)) 
      umapout$hla_clusters  <- as.factor(leiden_res)
      umapout$clust_leiden  <- as.factor(leiden_res)
    } else {
      find_optimal_res <- function(data, clust) {
        x <- seq(from = 0.01, to = 3, length.out = max(8, parallel::detectCores()))
        Res <- pbmcapply::pbmclapply(x, function(res) {
          leiden_res <- suppressWarnings(leiden::leiden(data, resolution_parameter = res, seed = 1985))
          add <- data.frame(res = res, clust_num = length(unique(as.factor(leiden_res))))
          return(add)
        }, mc.cores = parallel::detectCores()) # always attempt to parallelize
        Res <- do.call(rbind, Res)
        R <- Res[which(Res$clust_num == clust), ]
        if (nrow(R) == 0) { 
          # second attempt
          up <- min(Res$clust_num[Res$clust_num > clust]) %>% suppressWarnings()
          dn <- max(Res$clust_num[Res$clust_num < clust]) %>% suppressWarnings()
          if (!(up %in% c(-Inf, Inf)) & !(dn %in% c(-Inf, Inf))) {
            R <- Res[which(Res$clust_num %in% dn:up), ]
            x <- seq(from = min(R$res), to = max(R$res), length.out = max(8, parallel::detectCores()))
            Res <- pbmcapply::pbmclapply(x, function(res) {
              leiden_res <- suppressWarnings(leiden::leiden(data, resolution_parameter = res, seed = 1985))
              add <- data.frame(res = res, clust_num = length(unique(as.factor(leiden_res))))
              return(add)
            }, mc.cores = parallel::detectCores()) # always attempt to parallelize
            Res <- do.call(rbind, Res)
            R <- Res[which(Res$clust_num == clust), ]
          }
        }
        if (nrow(R) == 0) {return(NA)} else {return(R[which(R$clust_num == clust), ]$res %>% max())}
      }
      optimal_res <- find_optimal_res(ig, clust = k)
      if (!is.na(optimal_res)) {
        leiden_res <- suppressWarnings(leiden::leiden(ig, resolution_parameter = optimal_res, seed = 1985))
        umapout$hla_clusters  <- as.factor(leiden_res)
        umapout$clust_leiden  <- as.factor(leiden_res)
      } else {
        umapout$hla_clusters  <- NA
        umapout$clust_leiden  <- NA
        method <- "consensus"
        message(cat(crayon::red("could not compute 'leiden' with desired number of clusters, attempting 'dbscan', 'hclust', 'kmeans', 'gmm', and their 'consensus'"), sep = ""))
      }
    }
  }
  if (any(c("dbscan", "consensus") %in% method)) {
    message(cat(crayon::red(format(Sys.time(), "%H:%M:%S"), "- Density-based (dbscan) on UMAP coordinates"), sep = ""))
    if (is.null(k)) {
      dbscan_knn <- dbscan::kNNdist(as.matrix(umapout[,1:2]), k =  (dbscan_minPts-1))
      dbscan_knn_sort <- sort(dbscan_knn)
      infl <- c(FALSE, diff(diff(dbscan_knn_sort)>0)!=0) #get inflection points
      infl_select <- dbscan_knn_sort[intersect(names(dbscan_knn_sort), names(infl[which(infl == T)]))] %>% unname()
      optimal_eps <- infl_select[which.max(infl_select)]
      # # the optimal eps is at the 'elbow' (inflection point)
      # plot_df <- data.frame(x = 1:length(dbscan_knn_sort), y = dbscan_knn_sort)
      # ggplot(plot_df, aes(x = x, y = y)) +
      #   geom_line() +
      #   labs(title = "kNN Distance Plot", x = "Points sorted by distance", y = "kNN Distance") +
      #   geom_hline(yintercept = optimal_eps, color = "red", linetype = "dashed")
      dbscan_result <- dbscan::dbscan(as.matrix(umapout[,1:2]), eps = optimal_eps, minPts = mpts)
      umapout$hla_clusters  <- dbscan_result$cluster
      umapout$clust_dbscan <- dbscan_result$cluster
      # De-noise dbscan: reassign points with the noise cluster '0' to the classification of their nearest non-0 neighbor
      if (any(umapout$clust_dbscan == 0)) {
        nearest_neighbors <- FNN::get.knnx(umapout[umapout$clust_dbscan != 0, c("umap1", "umap2")], 
                                           umapout[which(umapout$clust_dbscan == 0), c("umap1", "umap2")], 
                                           k = 1)
        umapout$clust_dbscan[which(umapout$clust_dbscan == 0)] <- umapout$clust_dbscan[umapout$clust_dbscan != 0][nearest_neighbors$nn.index]
      }
    } else {
      find_optimal_eps <- function(data, desired_clusters, minPts = 30) {
        x <- c(0.0001 %o% 1.0903^(1:125))
        Eps <- pbmcapply::pbmclapply(x, function(eps) {
          dbscan_result <- dbscan::dbscan(data, eps = eps, minPts = minPts)
          num_clusters <- length(unique(dbscan_result$cluster)) - (sum(dbscan_result$cluster == 0) > 0) # excluding noise cluster (0)
          if (num_clusters == desired_clusters) {
            add <- data.frame(eps = eps, n_clust0= length(dbscan_result$cluster[dbscan_result$cluster == 0]), p_clust0 = sum(dbscan_result$cluster == 0)/length(dbscan_result$cluster))
          } else {
            add <- data.frame(eps = numeric(0), n_clust0 = numeric(0), p_clust0 = numeric(0))
          }
          return(add)
        }, mc.cores = parallel::detectCores()) # always attempt to parallelize
        Eps <- do.call(rbind, Eps)
        # fail if, despite having desired number of clusters, the 'noise' cluster contained more than 20% of data points
        if (nrow(Eps) == 0) {return(NA)} else if (Eps[which.min(Eps$n_clust0),]$p_clust0 > 0.2) {return(NA)} else {return(Eps[which.min(Eps$n_clust0),]$eps)}
      }
      mpts <- dbscan_minPts
      optimal_eps <- find_optimal_eps(as.matrix(umapout[,1:2]), desired_clusters = k, minPts = mpts)
      if (!is.na(optimal_eps)) {
        dbscan_result <- dbscan::dbscan(as.matrix(umapout[,1:2]), eps = optimal_eps, minPts = mpts)
        umapout$hla_clusters  <- dbscan_result$cluster #keep this 'noisy' here for now.
        umapout$clust_dbscan  <- dbscan_result$cluster
        # De-noise dbscan: reassign points with the noise cluster '0' to the classification of their nearest non-0 neighbor
        if (any(umapout$clust_dbscan == 0)) {
          nearest_neighbors <- FNN::get.knnx(umapout[umapout$clust_dbscan != 0, c("umap1", "umap2")], 
                                             umapout[which(umapout$clust_dbscan == 0), c("umap1", "umap2")], 
                                             k = 1)
          umapout$clust_dbscan[which(umapout$clust_dbscan == 0)] <- umapout$clust_dbscan[umapout$clust_dbscan != 0][nearest_neighbors$nn.index]
        }
      } else {
        umapout$hla_clusters  <- NA
        umapout$clust_dbscan  <- NA
        method <- "consensus"
        message(cat(crayon::red("could not compute 'dbscan' with desired number of clusters, attempting 'hclust', 'kmeans', 'gmm', and their 'consensus'"), sep = ""))
      }
    }
  }  
  if (any(c("hclust", "consensus") %in% method) & !is.null(k)) {
    message(cat(crayon::red(format(Sys.time(), "%H:%M:%S"), "- Hierarchical (agglomerative) on PCA space"), sep = ""))
    humapout <- stats::hclust(stats::dist(umat[, 1:min(250, ncol(umat))])) #limiting allowable number of PCs to 250 to prevent the algorithm from being needlessly slow
    umapout$hla_clusters <- stats::cutree(humapout, k = k)
    umapout$clust_hclust <- stats::cutree(humapout, k = k)
  }  
  if (any(c("kmeans", "consensus") %in% method) & !is.null(k)) {
    message(cat(crayon::red(format(Sys.time(), "%H:%M:%S"), "- k-means on PCA space"), sep = ""))
    humapout <- stats::kmeans(umat, centers = k)
    umapout$hla_clusters <- humapout$cluster
    umapout$clust_kmeans <- humapout$cluster
  }  
  if (any(c("gmm", "consensus") %in% method)) {
    message(cat(crayon::red(format(Sys.time(), "%H:%M:%S"), "- Gaussian Mixture Model on PCA space"), sep = ""))
    humapout <- mclust::Mclust(umat[, 1:min(25, ncol(umat))], G = k) # beyond 25 cols, the algorithm takes forever
    umapout$hla_clusters <- humapout$classification
    umapout$clust_gmm   <-  humapout$classification
  }  
  if (sum(method %in% c("leiden", "hclust", "kmeans", "gmm", "dbscan", "consensus"), na.rm = TRUE) >= 3L | "consensus" %in% method) {
    message(cat(crayon::red(format(Sys.time(), "%H:%M:%S"), "- Finding Consensus: redraw clusters with cells agreeing on membership in a majority of methods"), sep = ""))
    metamat <- umapout[, grepl("^clust_", names(umapout))]
    metamat <- metamat[, !apply(metamat, 2, function(x) all(is.na(x)))] # make sure there are no cols entirely NAs
    ## match the cluster naming (1, 2, 3, etc.) as consistent as possible between the different algorithms outputs
    ## i.e. cells in cluster 1 in one algorithm are also most likely in cluster 1 from another algorithm
    sub_metamat <- list()
    for (i in 2:(ncol(metamat))) {
      sub_df <- metamat %>% select(names(metamat)[1], names(metamat)[i])
      sub_metamat[[length(sub_metamat) + 1]] <- sub_df
    }
    sub_metamat <- lapply(sub_metamat, function(sub_df){
      # rename cols for function consistency, but keep old names for later
      old.names <- colnames(sub_df)
      colnames(sub_df) <- c("ref_idx", "q_idx")
      # define permutations. be carefull with big numbers, it follows a factorial growth
      perms <- combinat::permn(1:length(unique(sub_df$ref_idx)))
      # apply permutations
      apply_perm <- function(query, permutation) {
        sapply(query, function(x) permutation[x])
      }
      try_perms <- lapply(perms, function(p) {
        sub_df %>% mutate(q_idx = apply_perm(q_idx, p))
      })
      results <- lapply(try_perms, function(d) {
        d %>% group_by(ref_idx) %>% summarise(concurrence = mean(q_idx == ref_idx)) %>% as.data.frame()
      })
      # the winner is one which maximises the concurrence
      winner <- lapply(results, function(s) {
        score <- s$concurrence %>% sum()
      }) %>% unlist() %>% which.max()
      sub_df$q_idx <- try_perms[[winner]]$q_idx
      colnames(sub_df) <- old.names
      return(sub_df)
    })
    metamat <- cbind(select(sub_metamat[[1]], 1), lapply(sub_metamat, function(s){s <- select(s, 2)}) %>% do.call(cbind, .))
    # next, find consensus among clusters
    get_consensus <- function(row) {
      freq <- table(row)
      max_freq <- max(freq)
      consensus_value <- paste(names(freq[freq == max_freq]), collapse = "/")
      probability <- max_freq / length(row)
      return(list(clust_cons = consensus_value, prob = probability))
    }
    consensus_res <- apply(metamat, 1, get_consensus)
    cons <- cbind(metamat, consensus = sapply(consensus_res, function(x) x$clust_cons), 
                  probability = sapply(consensus_res, function(x) x$prob))
    umapout$hla_clust_cons <- paste0(cons$consensus, " (prob = ", cons$probability, ")")
    unique_clust <- unique(sapply(consensus_res, function(x) x$clust_cons))
    if (length(unique_clust[!stringr::str_detect(unique_clust, "/")]) > 0) {
      unique_clust <- unique_clust[!stringr::str_detect(unique_clust, "/")]
    }
    umapout$hla_clusters <- NA
    # register only if it most probably belongs to a unique cluster and that probability is > 0.65
    umapout$hla_clusters <- ifelse(cons$consensus %in% unique_clust & cons$probability > 0.65, cons$consensus, NA)
    # if the number of consensus clusters is smaller than is supposed to (smaller than k), re-attempt to register with less stringency, i.e. with probability > 0.5
    if (length(na.omit(unique(umapout$hla_clusters))) < k) {
      umapout$hla_clusters <- ifelse(cons$consensus %in% unique_clust & cons$probability > 0.5, cons$consensus, NA)
    }
    g1 <- ggplot(umapout, aes(x=umap1, y=umap2, color=hla_clust_cons))+geom_point(size=pt_size)
  } else {g1 <- NULL}
  umapout$hla_clusters <- as.factor(umapout$hla_clusters)
  message(cat(crayon::red(format(Sys.time(), "%H:%M:%S"), "- Predicting clusters finished"), sep = ""))
  g <- ggplot(umapout, aes(x=umap1, y=umap2, color=hla_clusters))+geom_point(size=pt_size)#+scale_color_manual(values=pals::glasbey())+theme_bw()
  detach("package:mclust", unload = TRUE) #to avoid inter-package object masking and conflicts
  return(list(UMAP_coordinates = umapout, HLA_clusters_on_umap = g, genotype_on_umap = g0, consensus_efficiency_on_umap = g1, reads = reads, top80_PC = umat))
}


#' Mapping the HLA clusters found by the 'HLA_clusters()' function back into the scrHLAtag count files
#'
#' @param reads.list  is a scrHLAtag count file (or a list of scrHLAtag count files) including columns for CB, UMI, and HLA alleles (\url{https://github.com/furlan-lab/scrHLAtag}).
#' @param cluster_coordinates  is the UMAP coordinates dataframe with HLA clustering information (found by the \code{scrHLAmatrix::HLA_clusters()} function) from which clusters are extracted and mapped into the scrHLAtag count files by matching Cell Barcodes.
#' @param CB_rev_com  logical, called \code{TRUE} if the need to obtain the reverse complement of Cell Barcodes (CBs) is desired; default is \code{FALSE}.
#' @param additional_col_id  is a character; if mapping of any additional column ID of the Seurat metadata collected in the \code{cluster_coordinates} dataframe, to the scrHLAtag count files is desired. \code{NULL} by default.
#' @param parallelize  logical, called \code{TRUE} if using parallel processing (multi-threading) is desired; default is \code{FALSE}.
#' @import stringr
#' @return a large list containing scrHLAtag count file(s) including columns for CB, UMI, and HLA alleles, with the addition of HLA_clusters
#' @examples
#' dirs_path <- "path/to/scrHLAtag/out/files"
#' dirs<-list.dirs(path=dirs_path, full.names = T, recursive = F)
#' dirs<- lapply(dirs, list.dirs, recursive = F) %>% unlist
#' dirs<- lapply(dirs, dir, pattern = "unguided_hla_align_corrected", recursive = F, full.names = T) %>% unlist
#' dirnames <- c("AML_101_BM", "AML_101_34", "TN_BM", "TN_34") # this is how the samples were organized in the directories
#' ## Load the counts files
#' cts <- HLA_load(directories = dirs, dir_names = dirnames, seu = your_Seurat_obj)
#' ## Process those count files
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
      stop("The class of 'reads.list' should be either a 'data.frame', a 'data.table', or a 'list' of the formers", call. = FALSE)
    }
  } 
  if (CB_rev_com) {
    for (i in 1:length(reads.list)) {
      if (!is.null(additional_col_id)) {
        reads.list[[i]][, additional_col_id] <-  NA
        reads.list[[i]][, additional_col_id] <- cluster_coordinates[, additional_col_id][
          match(paste0(reads.list[[i]]$samp,reads.list[[i]]$id_cb_separator, parallel::mclapply(reads.list[[i]]$CB, function(x) intToUtf8(rev(utf8ToInt(chartr('ATGC', 'TACG', x)))), mc.cores = multi_thread) %>% unlist(),reads.list[[i]]$id_cb_suffix), 
                rownames(cluster_coordinates))
        ]
      }
      reads.list[[i]]$hla_clusters <- NA
      reads.list[[i]]$hla_clusters <- cluster_coordinates$hla_clusters[
        match(paste0(reads.list[[i]]$samp,reads.list[[i]]$id_cb_separator, parallel::mclapply(reads.list[[i]]$CB, function(x) intToUtf8(rev(utf8ToInt(chartr('ATGC', 'TACG', x)))), mc.cores = multi_thread) %>% unlist(),reads.list[[i]]$id_cb_suffix), 
              rownames(cluster_coordinates))
      ]
    }
  } else {
    for (i in 1:length(reads.list)) {
      if (!is.null(additional_col_id)) {
        reads.list[[i]][, additional_col_id] <-  NA
        reads.list[[i]][, additional_col_id] <- cluster_coordinates[, additional_col_id][
          match(paste0(reads.list[[i]]$samp, reads.list[[i]]$id_cb_separator, reads.list[[i]]$CB, reads.list[[i]]$id_cb_suffix), 
                rownames(cluster_coordinates))
        ]
      }
      reads.list[[i]]$hla_clusters <- NA
      reads.list[[i]]$hla_clusters <- cluster_coordinates$hla_clusters[
        match(paste0(reads.list[[i]]$samp, reads.list[[i]]$id_cb_separator, reads.list[[i]]$CB, reads.list[[i]]$id_cb_suffix), 
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