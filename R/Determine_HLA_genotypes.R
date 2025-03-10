#' A comprehensive function to extract the top HLA alleles from the scrHLAtag count files (based on the most counts per Cell Barcode)
#' 
#' @param reads_1  is the primary scrHLAtag count file (1 of 2 files containing either the mRNA molecular info or the genomic (gene) molecular info). It includes columns for CB, UMI, and HLA alleles (\url{https://github.com/furlan-lab/scrHLAtag}).
#' @param reads_2  is the secondary scrHLAtag count file (the alternative file vs. the one designated in '\code{reads_1}' argument). It includes columns for CB, UMI, and HLA alleles (\url{https://github.com/furlan-lab/scrHLAtag}). Default is \code{NULL}, in which case it will not be able to count alternative aligment.
#' @param k  can be \code{NULL} or a positive integer setting the number of cluster counts to partition the datapoints into, e.g. the number of entities or genotypes you \emph{think} there might be in your captured sample. If \code{NULL}, each clustering method will automatically determine \code{k} clusters on its own (will not work for \code{"hclust"} and \code{"kmeans"}, which need predefined \code{k}s).
#' @param seu  is the Seurat object associated with the scrHLAtag count file (\url{https://satijalab.org/seurat/index.html}), and entered here if matching CBs in count file with Seurat colnames is desired.
#' @param CB_rev_com  a logical, called \code{TRUE} if the need to obtain the reverse complement of Cell Barcodes (CBs) is desired; default is \code{FALSE}. 
#' @param hla_with_counts_above  number of total reads accross CBs at or above which an HLA allele is retained in the matrix.
#' @param CBs_with_counts_above  number of total reads accross HLA alleles at or above which a CB is retained in the matrix. Note: at present, the function will make sure that number of CBs is equal or more than available HLA alleles in the matrix.
#' @param match_CB_with_seu  a logical, called \code{TRUE} if filtering CBs in the scrHLAtag count file with matching ones in the Seurat object is desired. 
#' @param clust_method  the name of the graph-based clustering method to be used for partitioning cells based on their HLA count patterns. The choice is between a Community structure detection method: \code{"leiden"}, a Density-based method: \code{"dbscan"}, a Connectivity-based method: \code{"hclust"}, a Centroid-based method: \code{"kmeans"}, or a Distribution-based method: \code{"gmm"} (for Gaussian Mixture Model). The methods are run with their respective Default parameters. Some of those methods may predict \emph{true} allogeneic entities with better accuracy than others; as we cannot know a priori which is the best method, we propose the method: \code{"consensus"}, which groups cells in the same cluster if they agree on membership in > 65 percent of methods, otherwise they are unclassified (\code{NA}s).
#' @param n_PCs  the number of top principal components to retain in downstream clustering and umap analyses; default is \code{50} or the top 80 percent of PCs, whichever is smaller.
#' @param dbscan_minPts  only works for the  \code{"dbscan"} clust_method: number of minimum points required in the epsilon neighborhood radius (\code{eps}) of core points. While the other methods require 1 parameter (e.g., \code{k}), \code{"dbscan"} requires 2: \code{eps} and \code{minPts}. To acheive desired \code{k} clusters, a range of \code{eps} parameter is tested against a fixed \code{minPts}, which is provided here. Default at \code{30}, but can be adjusted higher or lower depending on how small and \emph{clumped} an allogeneic entity is suspected to be. 
#' @param top_cumulative_frac  fraction, \code{0} to \code{1}, of total counts for a particular HLA locus, where the highest-ranking alleles (or genotypes) cumulatively accounting for that fraction of the total, shall be retained in the candidate list of alleles; default at \code{0.85}.
#' @param bulk_to_perCB_threshold  threshold of number of uniquely mapped HLA alleles in the primary count file \code{reads_1} above which listing the top alleles uses the Pseudo-Bulk Algorithm and below which it uses the Per Single-Cell Algorithm. Default is \code{2500}.
#' @param allowed_alleles_per_cell  a numeric (single integer or range) determining the minimum and maximum number of highest ranking allele genotypes per cell to keep if such number is beyond those limits after filtering by fraction; default is \code{c(1, 200)}, usefull in the early scrHLAtag iterations to give minimap2 lots of room to align; once you are ready to finalize the top HLA allele list, you can try \code{c(1, 2)} if you assume a cell can have a min of 1 allele (homozygous) and a max of 2 (heterozygous).
#' @param stringent_mode  a logical, when called \code{TRUE}, the algorithm detects when the final iteration is near (unique alleles in read file is equal or less than 200 per allogeneic entity); thus, getting top alleles becomes more stringent, with \code{allowed_alleles_per_cell} switching to \code{c(1, 2)}. This argument, however, will be ignored if the \emph{Pseudo-Bulk} algorithm is running instead of the \emph{Per Single-Cell} algorithm.
#' @param correct_alleles  a logical. Minimap2 of scrHLAtag will preferentially map reads that are in reality \emph{DPA1*02:02:02}, \emph{A*03:01:01}, \emph{B*13:02:01}, \emph{C*02:02:02}, or \emph{C*04:01:01}, to the presumed erroneous \emph{DPA1*02:38Q}, \emph{A*03:437Q}, \emph{B*13:123Q}, \emph{C*02:205Q}, or \emph{C*04:09N/C*04:61N}, respectively. When called \code{TRUE}, the algorithm will replace the unlikely allele(s) with their putative \emph{correct} version(s). Will work if \code{stringent_mode} is \code{TRUE}.
#' @param field_resolution  integer fron \code{1} to \code{3}, to select the HLA nomenclature level of Field resolution, where \code{1}, \code{2}, or \code{3} will take into consideration the first, the first two, or the first three field(s) of HLA designation; default is \code{3}.
#' @param QC_mm2  a logical, called \code{TRUE} if removing low quality reads based on minimap2 tags is desired.
#' @param s1_percent_pass_score  percentage, \code{0} to \code{100}, cuttoff from the maximum score (best quality) of the minimap2 's1' tag, which a read needs to acheive to pass as acceptable; default at \code{80} and becomes less inclusive if value increases.
#' @param AS_percent_pass_score  percentage, \code{0} to \code{100}, cuttoff from the maximum score (best quality) of the minimap2 'AS' tag, which a read needs to acheive to pass as acceptable; default at \code{80} and becomes less inclusive if value increases.
#' @param NM_thresh  the number of mismatches and gaps in the minimap2 alignment at or below which the quality of the read is acceptable; default is \code{15}.
#' @param de_thresh  the gap-compressed per-base sequence divergence in the minimap2 alignment at or below which the quality of the read is acceptable; the number is between \code{0} and \code{1}, and default is \code{0.01}.
#' @param parallelize  a logical, called \code{TRUE} if using parallel processing (multi-threading) is desired; default is \code{TRUE}.
#' @param umap_spread for \code{uwot::umap()}; effective scale of embedded points determining how clustered/clumped the embedded points are.
#' @param umap_min_dist for \code{uwot::umap()}; effective minimum distance between embedded points. Smaller values will result in a more clustered/clumped embedding.
#' @param umap_repulsion for \code{uwot::umap()}; weighting applied to negative samples in low dimensional embedding optimization.
#' @param suppress_plots  called \code{TRUE} to suppress plots from appearing while running the function.
#' @param seed  numeric (or \code{NULL}), to set seed (or not) in the environment for reproducibility
#' @param ... other arguments passed onto \code{uwot::umap()}.
#' @import stringr
#' @import pbmcapply
#' @import parallel
#' @import Matrix
#' @import magrittr
#' @import htmltools
#' @import Seurat
#' @import dplyr
#' @return a Vector of the top HLA alleles in the count files (in terms of reads per Cell Barcode).
#' @examples
#' dirs_path <- "path/to/scrHLAtag/out/files"
#' dirs<-list.dirs(path=dirs_path, full.names = T, recursive = F)
#' dirs<- lapply(dirs, list.dirs, recursive = F) %>% unlist
#' dirs<- lapply(dirs, dir, pattern = "unguided_hla_align_corrected", recursive = F, full.names = T) %>% unlist
#' dirnames <- c("AML_101_BM", "AML_101_34", "TN_BM", "TN_34") # this is how the samples were organized in the directories
#' ## Load the counts files
#' cts <- HLA_load(directories = dirs, dir_names = dirnames, seu = your_Seurat_obj)
#' ## Process those count files
#' top_alleles <- Top_HLA_list(reads_1 = cts[["mRNA"]], reads_2 = cts[["gene"]], seu = your_Seurat_Obj)
#' # 
#' # Note: the function is optimized to choose whether to find top HLA alleles
#' # using a Pseudo-Bulk algorithm (large number of uniquely mapped alleles) or 
#' # a Single-Cell algorithm (number of alleles is fewer than 2500).
#' # 
#' # Note: if for a particular HLA, the alleles with the most counts are in a tie 
#' # between 3 or more alleles in a particular Cell Barcode, we cannot know which 
#' # ones are the top two alleles, so that CB is not counted. This is similar if  
#' # there were no counts for that allele (all zeros).
#' @export

Top_HLA_list <- function(reads_1, reads_2 = NULL, k = 1, seu = NULL, CB_rev_com = FALSE,
                         hla_with_counts_above = 0, CBs_with_counts_above = 25, match_CB_with_seu = TRUE, 
                         clust_method = "consensus", n_PCs = 50, dbscan_minPts = 30,
                         QC_mm2 = TRUE, s1_percent_pass_score = 80, AS_percent_pass_score = 80, NM_thresh = 15, de_thresh = 0.01,
                         top_cumulative_frac = 0.85, bulk_to_perCB_threshold = 2500,
                         allowed_alleles_per_cell = c(1, 200), stringent_mode = TRUE, correct_alleles = TRUE,                         
                         field_resolution = 3, parallelize = TRUE, 
                         umap_spread = 5, umap_min_dist = 0.001, umap_repulsion = 0.001, seed = NULL, suppress_plots = TRUE, ...) {
  s <- Sys.time()
  # creating HLA umap clusters
  HLA_umap_clusters <- scrHLAmatrix::HLA_clusters(reads = reads_1, 
                                                  k = k, 
                                                  seu = seu, 
                                                  CB_rev_com = CB_rev_com,  # TRUE for 3prime 10x on pacbio
                                                  match_CB_with_seu = match_CB_with_seu,
                                                  hla_with_counts_above = hla_with_counts_above, 
                                                  CBs_with_counts_above = CBs_with_counts_above,
                                                  QC_mm2 = QC_mm2, s1_percent_pass_score = s1_percent_pass_score, 
                                                  AS_percent_pass_score = AS_percent_pass_score, NM_thresh = NM_thresh, de_thresh = de_thresh,
                                                  return_heavy = TRUE, n_PCs = n_PCs,
                                                  method = clust_method, dbscan_minPts = dbscan_minPts,
                                                  spread = umap_spread, min_dist = umap_min_dist, repulsion_strength = umap_repulsion, seed = seed, suppress_plots = suppress_plots, ...)
  if (!suppress_plots) print(HLA_umap_clusters[[2]]+scale_color_manual(values=pals::glasbey())+theme_classic())
  if ((reads_1$gene %>% unique() %>% length()) > bulk_to_perCB_threshold) {
    message(cat("\nReads count file shows greater than ", bulk_to_perCB_threshold, " uniquely mapped HLA alleles; extracting top alleles using the Pseudo-Bulk Algorithm", sep = ""))
    reads_1 <- map_HLA_clusters(reads.list = reads_1, HLA_umap_clusters[[1]], CB_rev_com = CB_rev_com)
    if (!is.null(reads_2)) {reads_2 <- map_HLA_clusters(reads.list = reads_2, HLA_umap_clusters[[1]], CB_rev_com = CB_rev_com)}
    cluster <- levels(reads_1$hla_clusters)
    # alt_align <- rep(c(TRUE, FALSE), length.out = length(cluster))
    reads.list <- list(reads_1 = reads_1, reads_2 = reads_2)
    reads.list <- reads.list[!sapply(reads.list, is.null)]
    top_alleles_HLA <- c()
    for (i in seq_along(cluster)) {
      top_alle <- c()
      for (j in seq_along(reads.list)) {
        top_al <- scrHLAmatrix:::Top_HLA_list_bulk(reads_1 = reads.list[[j]][reads.list[[j]]$hla_clusters %in% cluster[i], ], 
                                                   frac = top_cumulative_frac,
                                                   insert_pop_most_freq = T,   # recommended to be TRUE for first iteration
                                                   use_alt_align_ABC = F)
        top_alle <- c(top_alle, top_al)
      }
      top_alleles_HLA <- c(top_alleles_HLA, top_alle)
    }
    top_alleles_HLA <- top_alleles_HLA %>% unique() %>% sort()
  } else {
    reads_1 <- scrHLAmatrix::map_HLA_clusters(reads.list = HLA_umap_clusters[["reads"]], HLA_umap_clusters[[1]], CB_rev_com = F)
    cl <- unique(reads_1$hla_clusters)
    cl <- cl[!is.na(cl)]
    # is this the last round?
    reads <- lapply(cl, function(x) {
      tmp <- reads_1[reads_1$hla_clusters %in% x,]
      return(tmp)
    })
    if (parallelize) {multi_thread <- parallel::detectCores()} else {multi_thread <- 1}
    last_round <- mclapply(1:length(reads), function(m) {
      sub_reads <- with(reads[[m]], split(reads[[m]], list(hla = hla)))
      last_one <- lapply(1:length(sub_reads), function(o) {
        tmp <- length(unique(sub_reads[[o]]$gene))>3
        return(tmp)
      })
      last_one <- !any(last_one)
      return(last_one)
    }, mc.cores = multi_thread)  
    rm(reads)
    # get number of unique alleles
    unique_alleles <- lapply(cl, function(x) {
      tmp <- reads_1[reads_1$hla_clusters %in% x,]
      return(tmp)
    })
    unique_alleles <- lapply(unique_alleles, function(x) {
      tmp <- x$gene %>% unique() %>% length()
      return(tmp)
    }) %>% unlist()
    if (all(unique_alleles <= length(cl)*100) & stringent_mode) {
      if (all(last_round)) {
        message(cat("\nReads count file shows ", length(cl)*100, " or fewer uniquely mapped HLA alleles per allogeneic entity;\nStringent mode active: ", crayon::green("extracting the FINAL list of HLA alleles using the Per Single-Cell algorithm"), sep = ""))
        allowed_alleles_per_cell <- c(1, 2)
        } else {
        message(cat("\nReads count file shows ", length(cl)*100, " or fewer uniquely mapped HLA alleles per allogeneic entity;\nStringent mode active: ", crayon::green("extracting a refined list of top HLA alleles using the Per Single-Cell algorithm (not final yet!)"), sep = ""))
        allowed_alleles_per_cell <- c(1, 2)
        }
    } else {
      message(cat("\nReads count file shows ", bulk_to_perCB_threshold, " or fewer uniquely mapped HLA alleles; extracting top alleles using the Per Single-Cell Algorithm", sep = ""))
    }
    top_alleles_HLA <- scrHLAmatrix:::Top_HLA_list_byCB_preprocessed(reads = reads_1,
                                                                     seu = seu,
                                                                     match_CB_with_seu = match_CB_with_seu,
                                                                     hla_with_counts_above = hla_with_counts_above,
                                                                     CBs_with_counts_above = CBs_with_counts_above,
                                                                     frac = top_cumulative_frac,
                                                                     allowed_alleles_per_cell = allowed_alleles_per_cell,
                                                                     field_resolution = field_resolution,
                                                                     parallelize = parallelize)
    if (all(unique_alleles <= length(cl)*100) & stringent_mode) {
      problematic_alleles <- c("DPA1*02:38Q", "A*03:437Q", "B*13:123Q", "C*02:205Q", "C*04:09N", "C*04:61N")
      names(problematic_alleles) <- c("DPA1*02:02:02", "A*03:01:01", "B*13:02:01", "C*02:02:02", "C*04:01:01", "C*04:01:01")
      is_the_allele_correct <- problematic_alleles[which(problematic_alleles %in% top_alleles_HLA)]
      names(is_the_allele_correct) <- names(problematic_alleles)[which(problematic_alleles %in% is_the_allele_correct)]
      for (x in seq_along(top_alleles_HLA[which(top_alleles_HLA %in% is_the_allele_correct)])) {
        message(cat(crayon::red("\nWarning: "), 
                    crayon::bgWhite(" ", is_the_allele_correct[x], " "),
                    crayon::red(" detected in final list. \nMake sure the correct allele is not "),
                    crayon::bgWhite(" ", names(is_the_allele_correct)[x], " "),
                    "\n  The mRNA (i.e. cDNA) reference IMGT/HLA sequence of the rare allele ", is_the_allele_correct[x], 
                    "\n  is more extended/complete at the 3' end than the similar but more common allele ", names(is_the_allele_correct)[x], ".",
                    "\n  Minimap2 of scrHLAtag will preferentially map ", names(is_the_allele_correct)[x], 
                    " reads to the ", is_the_allele_correct[x], " ref.", 
                    sep = ""))
      }
      if (correct_alleles & length(is_the_allele_correct) > 0) {
        message(cat("\nReplacing problematic allele(s) with their alternative version(s). \n  Note: 'correct_alleles = FALSE' to turn off this feature."))
        top_alleles_HLA <- sapply(top_alleles_HLA, function(x) if (x %in% is_the_allele_correct) names(is_the_allele_correct)[is_the_allele_correct == x] else x)
        names(top_alleles_HLA) <- NULL
        top_alleles_HLA <- top_alleles_HLA %>% sort()
      }
    } 
  }
  e <- difftime(Sys.time(), s, units = "sec") %>% as.numeric() %>% abs()
  message(cat("\nDone!! (runtime: ", format(as.POSIXlt(e, origin = "1970-01-01", tz = "UTC"), "%H:%M:%S", tz = "UTC"), ")", sep = ""))
  return(top_alleles_HLA)
}