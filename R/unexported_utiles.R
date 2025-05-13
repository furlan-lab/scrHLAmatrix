#' Extracting the top HLA alleles from the scrHLAtag count files in a Single-Cell approach (based on the most counts per Cell Barcode)
#' 
#' @param reads  is the scrHLAtag count file including columns for CB, UMI, and HLA alleles (https://github.com/furlan-lab/scrHLAtag).
#' @param cell_data_obj  is a cell dataset object associated with the scrHLAtag count file. Currently the function is compatible with Seurat (\url{https://satijalab.org/seurat/index.html}).
#' @param CB_rev_com  is a logical, called TRUE if the need to obtained the reverse complement of Cell Barcodes (CBs) is desired; default is FALSE. 
#' @param hla_with_counts_above  is the number of total reads accross CBs at or above which an HLA allele is retained in the matrix.
#' @param CBs_with_counts_above  is the number of total reads accross HLA alleles at or above which a CB is retained in the matrix. Note: at present, the function will make sure that number of CBs is equal or more than available HLA alleles in the matrix.
#' @param match_CB_with_obj  is a logical, called TRUE if filtering CBs in the scrHLAtag count file with matching ones in the cell dataset object is desired. 
#' @param frac  is the fraction (0 to 1) of total reads for a particular HLA gene, which incorporates the highest ranking alleles of that gene in terms of number of reads; default at 0.85 .
#' @param allowed_alleles_per_cell  is a numeric (single or range) determining the minimum and maximum number of highest ranking allele genotypes per cell to keep if such number is beyond those limits when filtering by fraction 'frac'; default is c(1, 200), usefull in the early scrHLAtag iterations to give minimap2 lots of room to align; once you are ready to finalize the top HLA allele list, you can try c(1, 2) if you assume a cell can have a min of 1 allele (homozygous) and a max of 2 (heterozygous).
#' @param flag_bad_alleles  is a logical, called TRUE if instead of removing alleles with imperfect alignment scores (when finalizing the top HLA allele list), those alleles are flagged and returned separately.
#' @param field_resolution  is a numeric, to select the HLA nomenclature level of Field resolution, where 1, 2, or 3 will take into consideration the first, the first two, or the first three field(s) of HLA designation; default is 3.
#' @param QC_mm2  is a logical, called TRUE if removing low quality reads based on minimap2 tags is desired.
#' @param s1_percent_pass_score  is a percentage (\code{0} to \code{100}) cuttoff from the maximum score (best quality) of the minimap2 's1' tag, which a read needs to acheive to pass as acceptable; default at \code{80} and becomes less inclusive if value increases.
#' @param AS_percent_pass_score  is a percentage (\code{0} to \code{100}) cuttoff from the maximum score (best quality) of the minimap2 'AS' tag, which a read needs to acheive to pass as acceptable; default at \code{80} and becomes less inclusive if value increases.
#' @param NM_thresh  is the number of mismatches and gaps in the minimap2 alignment at or below which the quality of the read is acceptable; default is 15.
#' @param de_thresh  is the gap-compressed per-base sequence divergence in the minimap2 alignment at or below which the quality of the read is acceptable; the number is between 0 and 1, and default is 0.015.
#' @param parallelize  is a logical, called TRUE if using parallel processing (multi-threading) is desired; default is TRUE.
#' @import stringr
#' @import pbmcapply
#' @import parallel
#' @import Matrix
#' @import magrittr
#' @import htmltools
#' @import Seurat
#' @import dplyr
#' @import kSamples
#' @return a Vector of the top HLA alleles in the count files (in terms of reads per Cell Barcode).
#' @noRd

Top_HLA_list_byCB <- function(reads, cell_data_obj = NULL, CB_rev_com = FALSE, hla_with_counts_above = 5, CBs_with_counts_above = 15, match_CB_with_obj = TRUE, frac = 0.85, allowed_alleles_per_cell = c(1, 200), flag_bad_alleles = TRUE, field_resolution = 3, QC_mm2 = TRUE, s1_percent_pass_score = 80, AS_percent_pass_score = 80, NM_thresh = 15, de_thresh = 0.01, parallelize = TRUE) {
  seu <- cell_data_obj
  cell_data_obj <- NULL
  ## parallelize
  if (parallelize) {
    multi_thread <- parallel::detectCores()
    message(cat("\nMulti-threading! Available cores: ", parallel::detectCores()))
  } else {
    multi_thread <- 1
  }
  ## Remove low quality reads based on minimap2 tags
  if (QC_mm2) {
    message(cat("\nRemoving low quality reads based on minimap2 tags"))
    reads <- with(reads, split(reads, list(gene=gene)))
    reads <- pbmclapply(reads, function(df){df[df$s1 > (s1_percent_pass_score/100)*max(df$s1) & df$AS > (AS_percent_pass_score/100)*max(df$AS) & df$NM <= NM_thresh & df$de <= de_thresh,]}, mc.cores = multi_thread)
    reads <- data.table::rbindlist(reads)
    row.names(reads)<-NULL
    reads <- setDF(reads)
  } 
  ## Reverse Complement the CB
  if (CB_rev_com) {
    message(cat("\nConverting Cell Barcodes to their reverse complements"))
    # reads$CB <- pbmclapply(reads$CB, function(x) as.character(Biostrings::reverseComplement(DNAString(x))), mc.cores = multi_thread) %>% unlist() # slow
    reads$CB <- pbmclapply(reads$CB, function(x) intToUtf8(rev(utf8ToInt(chartr('ATGC', 'TACG', x)))), mc.cores = multi_thread) %>% unlist()        # fast
  }  
  ## the 3-field level resolution of HLA is actually 2-field for MICA and MICB. must fix this glitch:
  special <- "[-_*|?.+$^]"
  reads$hla <- stringr::str_split_fixed(reads$gene, special, 2)[ ,1]
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
  ## Selecting HLA designation level of field resolution.
  if (field_resolution == 2) {
    reads$a <- sapply(reads$gene, function(x) strsplit(x, ":")[[1]][1])
    reads$b <- sapply(reads$gene, function(x) strsplit(x, ":")[[1]][2])
    reads$a <- ifelse(is.na(reads$b), reads$a, paste0(reads$a, ":", reads$b))
    reads$gene <- reads$a
    reads$a <- NULL
    reads$b <- NULL
  }
  if (field_resolution == 1) {
    reads$a <- sapply(reads$gene, function(x) strsplit(x, ":")[[1]][1])
    reads$gene <- reads$a
    reads$a <- NULL
  }
  if (!(field_resolution %in% c(1:3))) {
    message(cat(crayon::red("Note: HLA field resolution must be 1, 2, or 3 to take into consideration the first, the first two, or the first three field(s) of HLA designation.\nKeeping field resolution at 3 by default.")))
  }
  ## Matrix formation
  message(cat("\nCreating HLA Count Matrix(ces)"))
  if ("hla_clusters" %in% colnames(reads)) {
    cl <- unique(reads$hla_clusters)
    cl <- cl[!is.na(cl)]
    message(cat("  Clusters generated from HLA distribution per Cell Barcode in UMAP space detected! Number of Clusters: ", length(cl)))
    reads <- lapply(cl, function(x) {
      tmp <- reads[reads$hla_clusters %in% x,]
      return(tmp)
    })
  } else {
    reads <- list(reads = reads)
    message(cat(crayon::green("Recommended:"), "To refine your results, first analyze distribution of alleles per Cell Barcodes in UMAP space using 'HLA_clusters()', then map the generated HLA Clusters back to your count data using 'map_HLA_clusters()'."))
  }  
  # is this the last round?
  last_round <- mclapply(1:length(reads), function(m) {
    sub_reads <- with(reads[[m]], split(reads[[m]], list(hla = hla)))
    last_one <- lapply(1:length(sub_reads), function(o) {
      tmp <- length(unique(sub_reads[[o]]$gene))>3*length(cl)
      return(tmp)
    })
    last_one <- !any(last_one)
    return(last_one)
  }, mc.cores = multi_thread)  
  # identifying potentially bad alleles
  reads <- mclapply(1:length(reads), function(m) {
    reads[[m]]$hla_field1 <- str_split_fixed(reads[[m]]$gene, ":", 4)[,1]
    return(reads[[m]])
  }, mc.cores = multi_thread)
  # if this is the final top HLA list, try capture potentially bad alleles (based on mm2 metrics) in a list
  if (!(min(allowed_alleles_per_cell) == 1 && max(allowed_alleles_per_cell) <= 3)) {
    bad_alleles <- mclapply(1:length(reads), function(m) { # `m` is an unused arg
      tmp <- character()
      return(tmp)
    }, mc.cores = multi_thread)
    names(bad_alleles) <- paste0("allo_cluster_", seq_along(bad_alleles))
  } else {
    # in preparation for the upcoming Anderson Darling tests, skip the test if the tested vector is constant across the categories (otherwise it will return an error):
    safe_ad_pval <- function(var, data) {
      y <- data[[var]]
      if (length(unique(y)) < 2) {
        return(1)    # give it the max p-value = 1 which will give -log10(1) = 0
      }
      ad_res <- kSamples::ad.test(stats::as.formula(sprintf("%s ~ gene", var)), data = data)
      return(ad_res$ad[1, 3])
    }
    bad_alleles <- mclapply(1:length(reads), function(m) {
      kw <- data.frame(matrix(ncol = 6, nrow = 0, dimnames = list(NULL, c("ad_s1", "ad_AS", "ad_NM", "ad_de", "sum_ad", "best")))) 
      bad <- character()
      for (j in 1:length(unique(reads[[m]]$hla_field1))) {
        sub_ctsu <- reads[[m]][reads[[m]]$hla_field1 == unique(reads[[m]]$hla_field1)[j],]
        max_s1 <- max(sub_ctsu$s1)
        max_AS <- max(sub_ctsu$AS)
        sc <- data.frame(matrix(ncol = length(colnames(sub_ctsu)), nrow = 0, dimnames = list(NULL, colnames(sub_ctsu)))) 
        for (i in 1:length(unique(sub_ctsu$gene))){
          if (nrow(sub_ctsu[sub_ctsu$gene == unique(sub_ctsu$gene)[i],]) > 1000) {
            sct <- sub_ctsu[sub_ctsu$gene == unique(sub_ctsu$gene)[i],][sample(nrow(sub_ctsu[sub_ctsu$gene == unique(sub_ctsu$gene)[i],]), 1000),]
            sct$s1 <- sct$s1*(max_s1/max(sct$s1)) #normalizing the s1 values among the gene alleles
            sct$AS <- sct$AS*(max_AS/max(sct$AS)) #normalizing the AS values among the gene alleles
          } else {
            sct <- sub_ctsu[sub_ctsu$gene == unique(sub_ctsu$gene)[i],]
            sct$s1 <- sct$s1*(max_s1/max(sct$s1)) #normalizing the s1 values among the gene alleles
            sct$AS <- sct$AS*(max_AS/max(sct$AS)) #normalizing the AS values among the gene alleles
          }
          sc <- rbind(sc, sct)
          rm(sct)
        }
        sub_ctsu <- sc
        rm(sc)
        if (length(unique(sub_ctsu$gene))==1){
          k <- data.frame(ad_s1 = 0,
                          ad_AS = 0,
                          ad_NM = 0,
                          ad_de = 0,
                          sum_ad= 0,
                          best=NA)
        } else {
          k <- data.frame(ad_s1 = trunc(-log10(safe_ad_pval("s1",  sub_ctsu))*10^2)/10^2,
                          ad_AS = trunc(-log10(safe_ad_pval("AS",  sub_ctsu))*10^2)/10^2,
                          ad_NM = trunc(-log10(safe_ad_pval("NM",  sub_ctsu))*10^2)/10^2,
                          ad_de = trunc(-log10(safe_ad_pval("de",  sub_ctsu))*10^2)/10^2,
                          sum_ad= 0,
                          best=NA)
        }
        row.names(k) <- unique(reads[[m]]$hla_field1)[j]
        k[1:4][ is.na(k[1:4]) ] <- 0
        k$sum_ad <- rowSums(k[,1:4])
        best <- c()
        if (k$ad_s1 > 10 & k$sum_ad > 30) {
          st <- group_by(sub_ctsu, gene) %>%
            summarise(
              count = n(),
              mean = mean(s1, na.rm = TRUE),
              sd = sd(s1, na.rm = TRUE),
              q25 = quantile(s1, probs = .25, na.rm =TRUE),
              median = median(s1, na.rm = TRUE),
              q75 = quantile(s1, probs = .75, na.rm =TRUE),
              IQR = IQR(s1, na.rm = TRUE)
            )
          star <- st$gene[which.max(st$mean)]
          other_alleles <- setdiff(unique(sub_ctsu$gene), star)
          # Run pairwise AD tests between best and others to identify the ones that rank up next to the best, then capture them in the list
          pairwise_results <- lapply(other_alleles, function(lvl) {
            subset_data <- subset(sub_ctsu, gene %in% c(star, lvl))
            test_result <- trunc(-log10(safe_ad_pval("s1",  sub_ctsu))*10^2)/10^2
            data.frame(
              comparison_versus = lvl,
              mlog10Pval = test_result
            )
          })
          pairwise_results <- do.call(rbind, pairwise_results)
          best <- c(best, star, pairwise_results$comparison_versus[which(pairwise_results$mlog10Pval < 10)])
        }
        if (k$ad_AS > 10 & k$sum_ad > 30) {
          st <- group_by(sub_ctsu, gene) %>%
            summarise(
              count = n(),
              mean = mean(AS, na.rm = TRUE),
              sd = sd(AS, na.rm = TRUE),
              q25 = quantile(AS, probs = .25, na.rm =TRUE),
              median = median(AS, na.rm = TRUE),
              q75 = quantile(AS, probs = .75, na.rm =TRUE),
              IQR = IQR(AS, na.rm = TRUE)
            )
          star <- st$gene[which.max(st$mean)]
          other_alleles <- setdiff(unique(sub_ctsu$gene), star)
          # Run pairwise AD tests between best and others to identify the ones that rank up next to the best, then capture them in the list
          pairwise_results <- lapply(other_alleles, function(lvl) {
            subset_data <- subset(sub_ctsu, gene %in% c(star, lvl))
            test_result <- trunc(-log10(safe_ad_pval("AS",  sub_ctsu))*10^2)/10^2
            data.frame(
              comparison_versus = lvl,
              mlog10Pval = test_result
            )
          })
          pairwise_results <- do.call(rbind, pairwise_results)
          best <- c(best, star, pairwise_results$comparison_versus[which(pairwise_results$mlog10Pval < 10)])
        }
        if (k$ad_NM > 10 & k$sum_ad > 30) {
          st <- group_by(sub_ctsu, gene) %>%
            summarise(
              count = n(),
              mean = mean(NM, na.rm = TRUE),
              sd = sd(NM, na.rm = TRUE),
              q25 = quantile(NM, probs = .25, na.rm =TRUE),
              median = median(NM, na.rm = TRUE),
              q75 = quantile(NM, probs = .75, na.rm =TRUE),
              IQR = IQR(NM, na.rm = TRUE)
            )
          star <- st$gene[which.min(st$mean)]
          other_alleles <- setdiff(unique(sub_ctsu$gene), star)
          # Run pairwise AD tests between best and others to identify the ones that rank up next to the best, then capture them in the list
          pairwise_results <- lapply(other_alleles, function(lvl) {
            subset_data <- subset(sub_ctsu, gene %in% c(star, lvl))
            test_result <- trunc(-log10(safe_ad_pval("NM",  sub_ctsu))*10^2)/10^2
            data.frame(
              comparison_versus = lvl,
              mlog10Pval = test_result
            )
          })
          pairwise_results <- do.call(rbind, pairwise_results)
          best <- c(best, star, pairwise_results$comparison_versus[which(pairwise_results$mlog10Pval < 10)])
        }
        if (k$ad_de > 10 & k$sum_ad > 30) {
          st <- group_by(sub_ctsu, gene) %>%
            summarise(
              count = n(),
              mean = mean(de, na.rm = TRUE),
              sd = sd(de, na.rm = TRUE),
              q25 = quantile(de, probs = .25, na.rm =TRUE),
              median = median(de, na.rm = TRUE),
              q75 = quantile(de, probs = .75, na.rm =TRUE),
              IQR = IQR(de, na.rm = TRUE)
            )
          star <- st$gene[which.min(st$mean)]
          other_alleles <- setdiff(unique(sub_ctsu$gene), star)
          # Run pairwise AD tests between best and others to identify the ones that rank up next to the best, then capture them in the list
          pairwise_results <- lapply(other_alleles, function(lvl) {
            subset_data <- subset(sub_ctsu, gene %in% c(star, lvl))
            test_result <- trunc(-log10(safe_ad_pval("de",  sub_ctsu))*10^2)/10^2
            data.frame(
              comparison_versus = lvl,
              mlog10Pval = test_result
            )
          })
          pairwise_results <- do.call(rbind, pairwise_results)
          best <- c(best, star, pairwise_results$comparison_versus[which(pairwise_results$mlog10Pval < 10)])
        }
        best <- names(table(best))[which(table(best) >=2 )] # each allele has to be best in at least 2 of the 4 metrics
        if (sum(k[,c(1:4)] > 10) < 2) {best <- NA} # the rule is that 2 of the log Pvals should be greater than 10
        if (!last_round[[m]]) {best <- NA} # if this is not the last round, no need to get best quality alleles just yet
        k$best <- paste(best, collapse = "_")
        if (!is.null(best) && !anyNA(best)) {bad <- c(bad, setdiff(unique(st$gene), best))} #if there is a best allele, capture the bad ones in a list
        bad <- sort(bad)
        kw <- rbind(kw, k)
        kw <- kw[order(row.names(kw)), ] #`kw` is a dataframe that is not returned, only useful for debug
        rm(k)
      }
      return(bad)
    }, mc.cores = multi_thread)
    names(bad_alleles) <- paste0("allo_cluster_", seq_along(bad_alleles))
  } 
  # creating as many count matrices as there are elements (an element for each HLA umap "cluster") in the "reads" list
  matrices <- mclapply(1:length(reads), function(m) {
    alleles <- unique(reads[[m]]$gene) %>% sort()
    reads[[m]]$seu_barcode <- paste0(reads[[m]]$samp, reads[[m]]$id_cb_separator, reads[[m]]$CB, reads[[m]]$id_cb_suffix)
    reads[[m]] <- with(reads[[m]], split(reads[[m]], list(seu_barcode=seu_barcode)))
    HLA.matrix <- matrix(0, nrow = length(alleles), ncol = length(reads[[m]]), dimnames = list(alleles, names(reads[[m]])))
    pb <- pbmcapply::progressBar(min = 0, max = length(reads[[m]]), style = "ETA", char = "=")
    for (i in 1:length(reads[[m]])) {
      counts <- table(reads[[m]][[i]]$gene)
      HLA.matrix[, i] <- counts[alleles]
      setTxtProgressBar(pb, i)
    }
    close(pb)
    HLA.matrix[is.na(HLA.matrix)]<-0
    #HLA.matrix<-Matrix(HLA.matrix,sparse = T)
    ## Matching with Seurat colnames
    if (is.null(seu)) {
      part_HLA<- HLA.matrix
    } else {
      if ("Seurat" %in% class(seu)) {
        if(match_CB_with_obj) {
          part_HLA<- HLA.matrix[,colnames(HLA.matrix) %in% Cells(seu)]
        } else {
          part_HLA<- HLA.matrix
        }
      } else {
        stop("Single-cell dataset container (in argument 'cell_data_obj') of class 'Seurat' are currently compatible", call. = FALSE)
      }
    }
    ## removing HLA alleles with low counts overall
    r <- hla_with_counts_above
    part_HLA <- part_HLA[which(rowSums(part_HLA)>=r), , drop = F]
    ## removing cell barcodes (CBs) with low counts overall
    find_n <- function(part_HLA, n) {
      while (n > 0 && dim(part_HLA[, which(colSums(part_HLA) >= n), drop = FALSE])[2] < dim(part_HLA)[1]) {
        n <- n - 1
      }
      return(n)  # Return the last valid n before the inequality becomes false
    }
    n <- find_n(part_HLA = part_HLA, n = CBs_with_counts_above)
    part_HLA <- part_HLA[,which(colSums(part_HLA)>=n), drop = F]
    return(part_HLA)
  }, mc.cores = 1) # more efficient if running on 1 CPU
  rm(reads)
  #-------------------------------------------------------------------------------
  message(cat("\nCounting Top Two alleles per HLA gene per Cell Barcode"))
  # creating as many top2list dataframes as there are elements in the "matrices" list (example, 2 "top2list" dataframes if there was 2 HLA "clusters")
  top2list <- mclapply(1:length(matrices), function(m) {
    matrices[[m]] <- matrices[[m]] %>% as.data.frame()
    special <- "[-_*|?.+$^]"
    matrices[[m]]$hla <- stringr::str_split_fixed(row.names(matrices[[m]]), special, 2)[ ,1]
    matrices[[m]] <- with(matrices[[m]], split(matrices[[m]], list(hla=hla)))
    top2cb <- pbmclapply(1:length(matrices[[m]]), function(j) {
      top2tab <- mclapply(1:(ncol(matrices[[m]][[j]])-1), function(i) {
        colval <- matrices[[m]][[j]][,i][matrices[[m]][[j]][,i]!=0] # column values that are greater than zero
        rnames <- row.names(matrices[[m]][[j]])[matrices[[m]][[j]][,i]!=0] # only row.names with column values greater than zero
        names(colval) <- rnames
        cname <- colnames(matrices[[m]][[j]])[i] # the name of the CB
        max1 <- names(which(colval==max(colval))) %>% suppressWarnings() 
        max2 <- names(which(colval[colval!=max(colval)]==max(colval[colval!=max(colval)]))) %>% suppressWarnings()
        if (length(max1) > 2 | length(max1) == 0) {
          # if more than 2 alleles show up as having the top value of reads, we don't know which 2 of them are the correct ones, so we leave them empty
          # if there's no counts for that gene (length(max1) is 0, and so is length(max2) automatically), then also we leave empty
          max1 <- NULL
          max2 <- NULL
        } else {
          if (length(max1) == 2) {
            # if 2 alleles are in a tie with having top value of reads, we have our 2 alleles per cell!
            max1 <- max1 %>% sort() %>% paste(., collapse = "_")
            max2 <- NULL
          } else {
            if (length(max1) == 1 & length(max2) == 1) {
              max1 <- c(max1, max2) %>% sort() %>% paste(., collapse = "_")
              max2 <- NULL
            } else {
              # if the second most top values of reads include a number of alleles different than 1, we don't know which one of them is the second allele, so we leave empty
              max2 <- NA
              max1 <- c(max1, max2) %>% paste(., collapse = "_")
              max2 <- NULL
            }
          }
          maxt <- data.frame("seu_barcode"=cname, "toptwo"=max1, "hla"=matrices[[m]][[j]]$hla[1])
          return(maxt)
        }
      }, mc.cores = 1) # more efficient if inner loop always running on a single CPU (mc.cores = 1)
      return(do.call("rbind", top2tab))
    }, mc.cores = multi_thread)
    top2cb <- do.call("rbind", top2cb)
    return(top2cb)
  }, mc.cores = 1) # more efficient if running on 1 CPU
  ## count the alleles in the CBs and prepare them in a list
  message(cat("\nWriting List of Top HLA alleles"))
  # creating as many top_a_list character lists as there are elements in the "top2list" list (example, 2 "top_a_list" lists if there was 2 HLA "clusters")
  top_a_list <- pbmclapply(1:length(top2list), function(m) {
    top_a <- c()
    a<-unique(top2list[[m]]$hla) %>% sort()
    for (j in 1:length(a)){
      t<-as.data.table(table(top2list[[m]][top2list[[m]]$hla == a[j],]$toptwo))
      t$toptwo <- t$V1
      t <- aggregate(N ~ toptwo, data = t, FUN = sum)
      t<- t[order(-t$N, t$toptwo),]
      row.names(t)<- NULL
      t$hlagene <- a[j]
      # order by N
      t<- t[order(t$N),]
      t$csum <- cumsum(t$N)
      # ----------------------------------------------------------------
      # include alleles covering a specific fraction 'frac' of all reads
      if (!is.numeric(allowed_alleles_per_cell)) allowed_alleles_per_cell <- c(1, 200)
      min_alleles_keep <- min(allowed_alleles_per_cell)
      max_alleles_keep <- max(allowed_alleles_per_cell)
      if (length(t$toptwo) <= min_alleles_keep) {
        top <- t[order(-t$N),]
      } else {
        if (length(t$toptwo[t$csum/sum(t$N) > (1-frac)]) <= min_alleles_keep) {
          top <- t[order(-t$N),][1:min_alleles_keep,]
        } else {
          top <- t[t$csum/sum(t$N) > (1-frac),]
          top <- top[order(-top$N),]
        }
      }
      if (nrow(top) > max_alleles_keep) top <- top[1:max_alleles_keep,]
      top <- top$toptwo
      top_a <- c(top_a, top)
    }
    top_a <- sapply(top_a, function(x) strsplit(x, "_")[[1]]) %>% as.character()
    top_a <- unique(top_a)
    top_a <- top_a[top_a != "NA"]
    top_a <- top_a %>% sort()
    if (flag_bad_alleles) {
      remove_a <- character() #if flag is true, dont remove
    } else {
      remove_a <- bad_alleles[[m]]
    }
    top_a <- setdiff(top_a, remove_a)
    return(top_a)
  }, mc.cores = multi_thread)
  top_a_list <- do.call("c", top_a_list)
  top_a_list <- top_a_list %>% unique() %>% sort()
  flagged_alleles <- "alleles imperfectly associated with reads that are worth deeper examination will be listed here in later iterations"
  if (min(allowed_alleles_per_cell) == 1 && max(allowed_alleles_per_cell) <= 3 && any(lapply(bad_alleles, function(l) length(l)>0)) ) {
    flagged_alleles <- bad_alleles
    message(cat("Reads associated with these allele(s) may require a closer look at their sequence(s):"))
    for(cluster in names(bad_alleles)) {
      allele_boxes <- paste0(crayon::bgWhite(" ", bad_alleles[[cluster]], " "), collapse = " ") # wrap each allele in a white box
      message(cat(paste0("in ", crayon::bgWhite(" ", cluster, " "), " >> ", allele_boxes)))
    }
  } 
  message(cat("\nDone!!"))
  return(list(top_alleles = top_a_list, flagged_alleles = flagged_alleles))
}


#' Extracting the top HLA alleles from a PREPROCESSED scrHLAtag count files in a Single-Cell approach (based on the most counts per Cell Barcode)
#' 
#' @param reads  is the scrHLAtag count file including columns for CB, UMI, and HLA alleles (https://github.com/furlan-lab/scrHLAtag).
#' @param cell_data_obj  is a cell dataset object associated with the scrHLAtag count file. Currently the function is compatible with Seurat (\url{https://satijalab.org/seurat/index.html}).
#' @param hla_with_counts_above  is the number of total reads accross CBs at or above which an HLA allele is retained in the matrix.
#' @param CBs_with_counts_above  is the number of total reads accross HLA alleles at or above which a CB is retained in the matrix. Note: at present, the function will make sure that number of CBs is equal or more than available HLA alleles in the matrix.
#' @param match_CB_with_obj  is a logical, called TRUE if filtering CBs in the scrHLAtag count file with matching ones in the cell dataset object is desired. 
#' @param frac  is the fraction (0 to 1) of total reads for a particular HLA gene, which incorporates the highest ranking alleles of that gene in terms of number of reads; default at 0.85 .
#' @param allowed_alleles_per_cell  is a numeric (single or range) determining the minimum and maximum number of highest ranking allele genotypes per cell to keep if such number is beyond those limits when filtering by fraction 'frac'; default is c(1, 200), usefull in the early scrHLAtag iterations to give minimap2 lots of room to align; once you are ready to finalize the top HLA allele list, you can try c(1, 2) if you assume a cell can have a min of 1 allele (homozygous) and a max of 2 (heterozygous).
#' @param flag_bad_alleles  is a logical, called TRUE if instead of removing alleles with imperfect alignment scores (when finalizing the top HLA allele list), those alleles are flagged and returned separately.
#' @param field_resolution  is a numeric, to select the HLA nomenclature level of Field resolution, where 1, 2, or 3 will take into consideration the first, the first two, or the first three field(s) of HLA designation; default is 3.
#' @param parallelize  is a logical, called TRUE if using parallel processing (multi-threading) is desired; default is TRUE.
#' @import stringr
#' @import pbmcapply
#' @import parallel
#' @import Matrix
#' @import magrittr
#' @import htmltools
#' @import Seurat
#' @import dplyr
#' @import kSamples
#' @return a Vector of the top HLA alleles in the count files (in terms of reads per Cell Barcode).
#' @noRd

Top_HLA_list_byCB_preprocessed <- function(reads, cell_data_obj = NULL, hla_with_counts_above = 5, CBs_with_counts_above = 15, match_CB_with_obj = TRUE, frac = 0.85, allowed_alleles_per_cell = c(1, 200), flag_bad_alleles = TRUE, field_resolution = 3, parallelize = TRUE) {
  seu <- cell_data_obj
  cell_data_obj <- NULL
  ## parallelize
  if (parallelize) {
    multi_thread <- parallel::detectCores()
    # message(cat("\nMulti-threading! Available cores: ", parallel::detectCores()))
  } else {
    multi_thread <- 1
  }
  ## assuming, minimap2 quality control, CB reverse complement (if necessary), and 3 field MICA/B correction has been performed elswhere, then: 
  ## Selecting HLA designation level of field resolution.
  if (field_resolution == 2) {
    reads$a <- sapply(reads$gene, function(x) strsplit(x, ":")[[1]][1])
    reads$b <- sapply(reads$gene, function(x) strsplit(x, ":")[[1]][2])
    reads$a <- ifelse(is.na(reads$b), reads$a, paste0(reads$a, ":", reads$b))
    reads$gene <- reads$a
    reads$a <- NULL
    reads$b <- NULL
  }
  if (field_resolution == 1) {
    reads$a <- sapply(reads$gene, function(x) strsplit(x, ":")[[1]][1])
    reads$gene <- reads$a
    reads$a <- NULL
  }
  if (!(field_resolution %in% c(1:3))) {
    message(cat(crayon::red("Note: HLA field resolution must be 1, 2, or 3 to take into consideration the first, the first two, or the first three field(s) of HLA designation.\nKeeping field resolution at 3 by default.")))
  }
  ## Matrix formation
  message(cat("\nCreating HLA Count Matrix(ces)"))
  if ("hla_clusters" %in% colnames(reads)) {
    cl <- unique(reads$hla_clusters)
    cl <- cl[!is.na(cl)]
    # message(cat("  Clusters generated from HLA distribution per Cell Barcode in UMAP space detected! Number of Clusters: ", length(cl)))
    reads <- lapply(cl, function(x) {
      tmp <- reads[reads$hla_clusters %in% x,]
      return(tmp)
    })
  } else {
    reads <- list(reads = reads)
    message(cat(crayon::green("Recommended:"), "To refine your results, first analyze distribution of alleles per Cell Barcodes in UMAP space using 'HLA_clusters()', then map the generated HLA Clusters back to your count data using 'map_HLA_clusters()'."))
  }
  # is this the last round?
  last_round <- mclapply(1:length(reads), function(m) {
    sub_reads <- with(reads[[m]], split(reads[[m]], list(hla = hla)))
    last_one <- lapply(1:length(sub_reads), function(o) {
      tmp <- length(unique(sub_reads[[o]]$gene))>3*length(cl)
      return(tmp)
    })
    last_one <- !any(last_one)
    return(last_one)
  }, mc.cores = multi_thread)  
  # identifying potentially bad alleles
  reads <- mclapply(1:length(reads), function(m) {
    reads[[m]]$hla_field1 <- str_split_fixed(reads[[m]]$gene, ":", 4)[,1]
    return(reads[[m]])
  }, mc.cores = multi_thread)
  # if this is the final top HLA list, try capture potentially bad alleles (based on mm2 metrics) in a list
  if (!(min(allowed_alleles_per_cell) == 1 && max(allowed_alleles_per_cell) <= 3)) {
    bad_alleles <- mclapply(1:length(reads), function(m) { # `m` is an unused arg
      tmp <- character()
      return(tmp)
    }, mc.cores = multi_thread)
    names(bad_alleles) <- paste0("allo_cluster_", seq_along(bad_alleles))
  } else {
    # in preparation for the upcoming Anderson Darling tests, skip the test if the tested vector is constant across the categories (otherwise it will return an error):
    safe_ad_pval <- function(var, data) {
      y <- data[[var]]
      if (length(unique(y)) < 2) {
        return(1)    # give it the max p-value = 1 which will give -log10(1) = 0
      }
      ad_res <- kSamples::ad.test(stats::as.formula(sprintf("%s ~ gene", var)), data = data)
      return(ad_res$ad[1, 3])
    }
    bad_alleles <- mclapply(1:length(reads), function(m) {
      kw <- data.frame(matrix(ncol = 6, nrow = 0, dimnames = list(NULL, c("ad_s1", "ad_AS", "ad_NM", "ad_de", "sum_ad", "best")))) 
      bad <- character()
      for (j in 1:length(unique(reads[[m]]$hla_field1))) {
        sub_ctsu <- reads[[m]][reads[[m]]$hla_field1 == unique(reads[[m]]$hla_field1)[j],]
        max_s1 <- max(sub_ctsu$s1)
        max_AS <- max(sub_ctsu$AS)
        sc <- data.frame(matrix(ncol = length(colnames(sub_ctsu)), nrow = 0, dimnames = list(NULL, colnames(sub_ctsu)))) 
        for (i in 1:length(unique(sub_ctsu$gene))){
          if (nrow(sub_ctsu[sub_ctsu$gene == unique(sub_ctsu$gene)[i],]) > 1000) {
            sct <- sub_ctsu[sub_ctsu$gene == unique(sub_ctsu$gene)[i],][sample(nrow(sub_ctsu[sub_ctsu$gene == unique(sub_ctsu$gene)[i],]), 1000),]
            sct$s1 <- sct$s1*(max_s1/max(sct$s1)) #normalizing the s1 values among the gene alleles
            sct$AS <- sct$AS*(max_AS/max(sct$AS)) #normalizing the AS values among the gene alleles
          } else {
            sct <- sub_ctsu[sub_ctsu$gene == unique(sub_ctsu$gene)[i],]
            sct$s1 <- sct$s1*(max_s1/max(sct$s1)) #normalizing the s1 values among the gene alleles
            sct$AS <- sct$AS*(max_AS/max(sct$AS)) #normalizing the AS values among the gene alleles
          }
          sc <- rbind(sc, sct)
          rm(sct)
        }
        sub_ctsu <- sc
        rm(sc)
        if (length(unique(sub_ctsu$gene))==1){
          k <- data.frame(ad_s1 = 0,
                          ad_AS = 0,
                          ad_NM = 0,
                          ad_de = 0,
                          sum_ad= 0,
                          best=NA)
        } else {
          k <- data.frame(ad_s1 = trunc(-log10(safe_ad_pval("s1",  sub_ctsu))*10^2)/10^2,
                          ad_AS = trunc(-log10(safe_ad_pval("AS",  sub_ctsu))*10^2)/10^2,
                          ad_NM = trunc(-log10(safe_ad_pval("NM",  sub_ctsu))*10^2)/10^2,
                          ad_de = trunc(-log10(safe_ad_pval("de",  sub_ctsu))*10^2)/10^2,
                          sum_ad= 0,
                          best=NA)
        }
        row.names(k) <- unique(reads[[m]]$hla_field1)[j]
        k[1:4][ is.na(k[1:4]) ] <- 0
        k$sum_ad <- rowSums(k[,1:4])
        best <- c()
        if (k$ad_s1 > 10 & k$sum_ad > 30) {
          st <- group_by(sub_ctsu, gene) %>%
            summarise(
              count = n(),
              mean = mean(s1, na.rm = TRUE),
              sd = sd(s1, na.rm = TRUE),
              q25 = quantile(s1, probs = .25, na.rm =TRUE),
              median = median(s1, na.rm = TRUE),
              q75 = quantile(s1, probs = .75, na.rm =TRUE),
              IQR = IQR(s1, na.rm = TRUE)
            )
          star <- st$gene[which.max(st$mean)]
          other_alleles <- setdiff(unique(sub_ctsu$gene), star)
          # Run pairwise AD tests between best and others to identify the ones that rank up next to the best, then capture them in the list
          pairwise_results <- lapply(other_alleles, function(lvl) {
            subset_data <- subset(sub_ctsu, gene %in% c(star, lvl))
            test_result <- trunc(-log10(safe_ad_pval("s1",  sub_ctsu))*10^2)/10^2
            data.frame(
              comparison_versus = lvl,
              mlog10Pval = test_result
            )
          })
          pairwise_results <- do.call(rbind, pairwise_results)
          best <- c(best, star, pairwise_results$comparison_versus[which(pairwise_results$mlog10Pval < 10)])
        }
        if (k$ad_AS > 10 & k$sum_ad > 30) {
          st <- group_by(sub_ctsu, gene) %>%
            summarise(
              count = n(),
              mean = mean(AS, na.rm = TRUE),
              sd = sd(AS, na.rm = TRUE),
              q25 = quantile(AS, probs = .25, na.rm =TRUE),
              median = median(AS, na.rm = TRUE),
              q75 = quantile(AS, probs = .75, na.rm =TRUE),
              IQR = IQR(AS, na.rm = TRUE)
            )
          star <- st$gene[which.max(st$mean)]
          other_alleles <- setdiff(unique(sub_ctsu$gene), star)
          # Run pairwise AD tests between best and others to identify the ones that rank up next to the best, then capture them in the list
          pairwise_results <- lapply(other_alleles, function(lvl) {
            subset_data <- subset(sub_ctsu, gene %in% c(star, lvl))
            test_result <- trunc(-log10(safe_ad_pval("AS",  sub_ctsu))*10^2)/10^2
            data.frame(
              comparison_versus = lvl,
              mlog10Pval = test_result
            )
          })
          pairwise_results <- do.call(rbind, pairwise_results)
          best <- c(best, star, pairwise_results$comparison_versus[which(pairwise_results$mlog10Pval < 10)])
        }
        if (k$ad_NM > 10 & k$sum_ad > 30) {
          st <- group_by(sub_ctsu, gene) %>%
            summarise(
              count = n(),
              mean = mean(NM, na.rm = TRUE),
              sd = sd(NM, na.rm = TRUE),
              q25 = quantile(NM, probs = .25, na.rm =TRUE),
              median = median(NM, na.rm = TRUE),
              q75 = quantile(NM, probs = .75, na.rm =TRUE),
              IQR = IQR(NM, na.rm = TRUE)
            )
          star <- st$gene[which.min(st$mean)]
          other_alleles <- setdiff(unique(sub_ctsu$gene), star)
          # Run pairwise AD tests between best and others to identify the ones that rank up next to the best, then capture them in the list
          pairwise_results <- lapply(other_alleles, function(lvl) {
            subset_data <- subset(sub_ctsu, gene %in% c(star, lvl))
            test_result <- trunc(-log10(safe_ad_pval("NM",  sub_ctsu))*10^2)/10^2
            data.frame(
              comparison_versus = lvl,
              mlog10Pval = test_result
            )
          })
          pairwise_results <- do.call(rbind, pairwise_results)
          best <- c(best, star, pairwise_results$comparison_versus[which(pairwise_results$mlog10Pval < 10)])
        }
        if (k$ad_de > 10 & k$sum_ad > 30) {
          st <- group_by(sub_ctsu, gene) %>%
            summarise(
              count = n(),
              mean = mean(de, na.rm = TRUE),
              sd = sd(de, na.rm = TRUE),
              q25 = quantile(de, probs = .25, na.rm =TRUE),
              median = median(de, na.rm = TRUE),
              q75 = quantile(de, probs = .75, na.rm =TRUE),
              IQR = IQR(de, na.rm = TRUE)
            )
          star <- st$gene[which.min(st$mean)]
          other_alleles <- setdiff(unique(sub_ctsu$gene), star)
          # Run pairwise AD tests between best and others to identify the ones that rank up next to the best, then capture them in the list
          pairwise_results <- lapply(other_alleles, function(lvl) {
            subset_data <- subset(sub_ctsu, gene %in% c(star, lvl))
            test_result <- trunc(-log10(safe_ad_pval("de",  sub_ctsu))*10^2)/10^2
            data.frame(
              comparison_versus = lvl,
              mlog10Pval = test_result
            )
          })
          pairwise_results <- do.call(rbind, pairwise_results)
          best <- c(best, star, pairwise_results$comparison_versus[which(pairwise_results$mlog10Pval < 10)])
        }
        best <- names(table(best))[which(table(best) >=2 )] # each allele has to be best in at least 2 of the 4 metrics
        if (sum(k[,c(1:4)] > 10) < 2) {best <- NA} # the rule is that 2 of the log Pvals should be greater than 10
        if (!last_round[[m]]) {best <- NA} # if this is not the last round, no need to get best quality alleles just yet
        k$best <- paste(best, collapse = "_")
        if (!is.null(best) && !anyNA(best)) {bad <- c(bad, setdiff(unique(st$gene), best))} #if there is a best allele, capture the bad ones in a list
        bad <- sort(bad)
        kw <- rbind(kw, k)
        kw <- kw[order(row.names(kw)), ] #`kw` is a dataframe that is not returned, only useful for debug
        rm(k)
      }
      return(bad)
    }, mc.cores = multi_thread)
    names(bad_alleles) <- paste0("allo_cluster_", seq_along(bad_alleles))
  } 
  # creating as many count matrices as there are elements (an element for each HLA umap "cluster") in the "reads" list
  matrices <- mclapply(1:length(reads), function(m) {
    alleles <- unique(reads[[m]]$gene) %>% sort()
    reads[[m]]$seu_barcode <- paste0(reads[[m]]$samp, reads[[m]]$id_cb_separator, reads[[m]]$CB, reads[[m]]$id_cb_suffix)
    reads[[m]] <- with(reads[[m]], split(reads[[m]], list(seu_barcode=seu_barcode)))
    HLA.matrix <- matrix(0, nrow = length(alleles), ncol = length(reads[[m]]), dimnames = list(alleles, names(reads[[m]])))
    pb <- pbmcapply::progressBar(min = 0, max = length(reads[[m]]), style = "ETA", char = "=")
    for (i in 1:length(reads[[m]])) {
      counts <- table(reads[[m]][[i]]$gene)
      HLA.matrix[, i] <- counts[alleles]
      setTxtProgressBar(pb, i)
    }
    close(pb)
    HLA.matrix[is.na(HLA.matrix)]<-0
    #HLA.matrix<-Matrix(HLA.matrix,sparse = T)
    ## Matching with Seurat colnames
    if (is.null(seu)) {
      part_HLA<- HLA.matrix
    } else {
      if ("Seurat" %in% class(seu)) {
        if(match_CB_with_obj) {
          part_HLA<- HLA.matrix[,colnames(HLA.matrix) %in% Cells(seu)]
        } else {
          part_HLA<- HLA.matrix
        }
      } else {
        stop("Single-cell dataset container (in argument 'cell_data_obj') of class 'Seurat' are currently compatible", call. = FALSE)
      }
    }
    ## removing HLA alleles with low counts overall
    r <- hla_with_counts_above
    part_HLA <- part_HLA[which(rowSums(part_HLA)>=r), , drop = F]
    ## removing cell barcodes (CBs) with low counts overall
    n <- CBs_with_counts_above
    find_n <- function(part_HLA, n) {
      while (n > 0 && dim(part_HLA[, which(colSums(part_HLA) >= n), drop = FALSE])[2] < dim(part_HLA)[1]) {
        n <- n - 1
      }
      return(n)  # Return the last valid n before the inequality becomes false
    }
    n <- find_n(part_HLA = part_HLA, n = CBs_with_counts_above)
    part_HLA <- part_HLA[,which(colSums(part_HLA)>=n), drop = F]
    return(part_HLA)
  }, mc.cores = 1) # more efficient if running on 1 CPU
  rm(reads)
  #-------------------------------------------------------------------------------
  message(cat("\nCounting Top Two alleles per HLA gene per Cell Barcode"))
  # creating as many top2list dataframes as there are elements in the "matrices" list (example, 2 "top2list" dataframes if there was 2 HLA "clusters")
  top2list <- mclapply(1:length(matrices), function(m) {
    matrices[[m]] <- matrices[[m]] %>% as.data.frame()
    special <- "[-_*|?.+$^]"
    matrices[[m]]$hla <- stringr::str_split_fixed(row.names(matrices[[m]]), special, 2)[ ,1]
    matrices[[m]] <- with(matrices[[m]], split(matrices[[m]], list(hla=hla)))
    top2cb <- pbmclapply(1:length(matrices[[m]]), function(j) {
      top2tab <- mclapply(1:(ncol(matrices[[m]][[j]])-1), function(i) {
        colval <- matrices[[m]][[j]][,i][matrices[[m]][[j]][,i]!=0] # column values that are greater than zero
        rnames <- row.names(matrices[[m]][[j]])[matrices[[m]][[j]][,i]!=0] # only row.names with column values greater than zero
        names(colval) <- rnames
        cname <- colnames(matrices[[m]][[j]])[i] # the name of the CB
        max1 <- names(which(colval==max(colval))) %>% suppressWarnings() 
        max2 <- names(which(colval[colval!=max(colval)]==max(colval[colval!=max(colval)]))) %>% suppressWarnings()
        if (length(max1) > 2 | length(max1) == 0) {
          # if more than 2 alleles show up as having the top value of reads, we don't know which 2 of them are the correct ones, so we leave them empty
          # if there's no counts for that gene (length(max1) is 0, and so is length(max2) automatically), then also we leave empty
          max1 <- NULL
          max2 <- NULL
        } else {
          if (length(max1) == 2) {
            # if 2 alleles are in a tie with having top value of reads, we have our 2 alleles per cell!
            max1 <- max1 %>% sort() %>% paste(., collapse = "_")
            max2 <- NULL
          } else {
            if (length(max1) == 1 & length(max2) == 1) {
              max1 <- c(max1, max2) %>% sort() %>% paste(., collapse = "_")
              max2 <- NULL
            } else {
              # if the second most top values of reads include a number of alleles different than 1, we don't know which one of them is the second allele, so we leave empty
              max2 <- NA
              max1 <- c(max1, max2) %>% paste(., collapse = "_")
              max2 <- NULL
            }
          }
          maxt <- data.frame("seu_barcode"=cname, "toptwo"=max1, "hla"=matrices[[m]][[j]]$hla[1])
          return(maxt)
        }
      }, mc.cores = 1) # more efficient if inner loop always running on a single CPU (mc.cores = 1)
      return(do.call("rbind", top2tab))
    }, mc.cores = multi_thread)
    top2cb <- do.call("rbind", top2cb)
    return(top2cb)
  }, mc.cores = 1) # more efficient if running on 1 CPU
  ## count the alleles in the CBs and prepare them in a list
  message(cat("\nWriting List of Top HLA alleles"))
  # creating as many top_a_list character lists as there are elements in the "top2list" list (example, 2 "top_a_list" lists if there was 2 HLA "clusters")
  top_a_list <- pbmclapply(1:length(top2list), function(m) {
    top_a <- c()
    a<-unique(top2list[[m]]$hla) %>% sort()
    for (j in 1:length(a)){
      t<-as.data.table(table(top2list[[m]][top2list[[m]]$hla == a[j],]$toptwo))
      t$toptwo <- t$V1
      t <- aggregate(N ~ toptwo, data = t, FUN = sum)
      t<- t[order(-t$N, t$toptwo),]
      row.names(t)<- NULL
      t$hlagene <- a[j]
      # order by N
      t<- t[order(t$N),]
      t$csum <- cumsum(t$N)
      # ----------------------------------------------------------------
      # include alleles covering a specific fraction 'frac' of all reads
      if (!is.numeric(allowed_alleles_per_cell)) allowed_alleles_per_cell <- c(1, 200)
      min_alleles_keep <- min(allowed_alleles_per_cell)
      max_alleles_keep <- max(allowed_alleles_per_cell)
      if (length(t$toptwo) <= min_alleles_keep) {
        top <- t[order(-t$N),]
      } else {
        if (length(t$toptwo[t$csum/sum(t$N) > (1-frac)]) <= min_alleles_keep) {
          top <- t[order(-t$N),][1:min_alleles_keep,]
        } else {
          top <- t[t$csum/sum(t$N) > (1-frac),]
          top <- top[order(-top$N),]
        }
      }
      if (nrow(top) > max_alleles_keep) top <- top[1:max_alleles_keep,]
      top <- top$toptwo
      top_a <- c(top_a, top)
    }
    top_a <- sapply(top_a, function(x) strsplit(x, "_")[[1]]) %>% as.character()
    top_a <- unique(top_a)
    top_a <- top_a[top_a != "NA"]
    top_a <- top_a %>% sort()
    if (flag_bad_alleles) {
      remove_a <- character() #if flag is true, dont remove
    } else {
      remove_a <- bad_alleles[[m]]
    }
    top_a <- setdiff(top_a, remove_a)
    return(top_a)
  }, mc.cores = multi_thread)
  top_a_list <- do.call("c", top_a_list)
  top_a_list <- top_a_list %>% unique() %>% sort()
  flagged_alleles <- "alleles imperfectly associated with reads that are wort deeper examination will be listed here in later iterations"
  if (min(allowed_alleles_per_cell) == 1 && max(allowed_alleles_per_cell) <= 3 && any(lapply(bad_alleles, function(l) length(l)>0)) ) {
    flagged_alleles <- bad_alleles
    message(cat("Reads associated with these allele(s) may require a closer look at their sequence(s):"))
    for(cluster in names(bad_alleles)) {
      allele_boxes <- paste0(crayon::bgWhite(" ", bad_alleles[[cluster]], " "), collapse = " ") # wrap each allele in a white box
      message(cat(paste0("in ", crayon::bgWhite(" ", cluster, " "), " >> ", allele_boxes)))
    }
  } 
  # message(cat("\nDone!!"))
  return(list(top_alleles = top_a_list, flagged_alleles = flagged_alleles))
}




#' Extracting the top HLA alleles from the scrHLAtag count files in a Pseudo-Bulk approach
#' 
#' @param reads_1  is the primary scrHLAtag count file (1 of 2 files containing either the mRNA molecular info or the genomic (gene) molecular info). It includes columns for CB, UMI, and HLA alleles (https://github.com/furlan-lab/scrHLAtag).
#' @param reads_2  is the secondary scrHLAtag count file (the alternative file vs. the one designated in 'reads_1' argument). It includes columns for CB, UMI, and HLA alleles (https://github.com/furlan-lab/scrHLAtag). Default is NULL, in which case it will not be able to count alternative aligment and argument 'use_alt_align_ABC' will become irrelevant.
#' @param frac  is the fraction (0 to 1) of total reads for a particular HLA gene, which incorporates the highest ranking alleles of that gene in terms of number of reads; default at 0.85 .
#' @param min_alleles_keep  is a numeric representing the minimum number of highest ranking alleles to keep despite filtering by fraction 'frac'; default is 2.
#' @param min_reads_per_gene  is a numeric representing minimum number of total reads per HLA gene (including all its alleles) below which the gene is filtered out; default is 20. 
#' @param insert_pop_most_freq  is a logical, whether to to include the HLA allele with the highest frequency in the human population despite low reads in the count files; default is TRUE.
#' @param use_alt_align_ABC  is a logical, whether to use the count file from the alternative alignment (rather than the primary alignment) to count reads for the HLA-A, -B, and -C genes. It was observed in some cases that using genomic alignments has better accuracy in predicting genotype versus mRNA alignments (not the case for Class-II and other HLA genes); default is FALSE.
#' @import stringr
#' @import cowplot
#' @import magrittr
#' @import htmltools
#' @import ggplot2
#' @import dplyr
#' @return a Vector of the top HLA alleles in the count files (in terms of reads).
#' @noRd

Top_HLA_list_bulk <- function(reads_1, reads_2 = NULL, frac = 0.85, min_alleles_keep = 2, min_reads_per_gene = 20, insert_pop_most_freq = TRUE, use_alt_align_ABC = FALSE){
  if (is.null(reads_2)) {
    reads_2 <- reads_1
  }
  # extract the HLA genes that appear in the reads
  special <- "[_*|?.+$^]"
  reads_2$gene0 <- gsub(special, "-", reads_2$gene)
  reads_2$hla <- stringr::str_split_fixed(reads_2$gene, special, 2)[ ,1]
  reads_1$gene0 <- gsub(special, "-", reads_1$gene)
  reads_1$hla <- stringr::str_split_fixed(reads_1$gene, special, 2)[ ,1]
  # get the HLA genes in a list
  cts_abc <- unique(reads_2$hla)[order(unique(reads_2$hla))]
  if (use_alt_align_ABC) {
    if (length(cts_abc[which(cts_abc %in% c("A", "B", "C"))]) != 3) {
      stop("the secondary molecule info count file does not contain alleles belonging to all of HLA-A, -B, and -C", call. = FALSE)
    }
  }
  cts_notabc <- unique(reads_1$hla)[order(unique(reads_1$hla))]
  top_alleles <- c()
  for (j in 1:length(cts_notabc)){
    if (!use_alt_align_ABC) {
      t<-as.data.table(table(reads_1[reads_1$hla == cts_notabc[j],]$gene))
    } else {
      if (cts_notabc[j] %in% c("A", "B", "C")) {
        t<-as.data.table(table(reads_2[reads_2$hla == cts_notabc[j],]$gene))
      } else {
        t<-as.data.table(table(reads_1[reads_1$hla == cts_notabc[j],]$gene))
      }
    }
    # t$twofield <- sapply(t$V1, function(x) strsplit(x, ":")[[1]][1:2] %>% paste0(collapse=":"))
    t$twofield <- t$V1
    t <- aggregate(N ~ twofield, data = t, FUN = sum)
    t<- t[order(-t$N, t$twofield),]
    row.names(t)<- NULL
    # the 3-field level resolution of HLA mRNA is actually 2-field for MICA and MICB. must fix this glitch:
    if (strsplit(t$twofield, "\\*")[[1]][1] %in% c("MICA", "MICB") & !all(is.na(sapply(t$twofield, function(x) strsplit(x, ":")[[1]][3])))) {
      t$a <- sapply(t$twofield, function(x) strsplit(x, ":")[[1]][1])
      t$c <- sapply(t$twofield, function(x) strsplit(x, ":")[[1]][2])
      t$a <- ifelse(is.na(t$c), t$a, paste0(t$a, ":", t$c))
      t$b <- sapply(t$twofield, function(x) strsplit(x, ":")[[1]][3])
      y <- aggregate(b ~ a, data = t, function(x) x[order(x)][1])
      t <- aggregate(N ~ a, data = t, FUN = sum)
      t$a <- ifelse(is.na(y$b[match(t$a, y$a)]), t$a, paste0(t$a, ":", y$b[match(t$a, y$a)]) )
      colnames(t)[1] <- "twofield"
      rm(y)
      t<- t[order(-t$N, t$twofield),]
      row.names(t)<- NULL
    }
    freq<-scrHLAmatrix::freq 
    freq<- freq[complete.cases(freq$allele_freq),]
    sapply(freq, class)
    freq$sample_size <- gsub(",", "", freq$sample_size)
    freq$sample_size <- as.numeric(freq$sample_size)
    freq$allele_pop <- freq$allele_freq*freq$sample_size*2
    freq$allele_pop <- round(freq$allele_pop, digits = 0)
    t$fscore <- NA
    t$hlagene <- cts_notabc[j]
    for (i in 1:nrow(t)) {
      if (nrow(freq[freq$allele == t[i,]$twofield,]) != 0) {
        agg<-cbind(
          aggregate(sample_size ~ allele, data = freq[freq$allele == t[i,]$twofield,], FUN = mean),
          aggregate(allele_pop ~ allele, data = freq[freq$allele == t[i,]$twofield,], FUN = mean)[,2]
        )
        agg$mean_freq<-agg[,3]/(agg[,2]*2)
        t[i,3] <- agg$mean_freq
      }
    }
    t$onefield <- sapply(t$twofield, function(x) strsplit(x, ":")[[1]][1])
    # drop out the genes with fewer than x amount of reads
    if (sum(t$N) < min_reads_per_gene) {
      t <- t[0,]
    }
    # order by N
    t<- t[order(t$N),]
    t$csum <- cumsum(t$N)
    # ----------------------------------------------------------------
    # include alleles covering a specific fraction 'frac' of all reads
    if (length(t$twofield) <= min_alleles_keep) {
      top <- t$twofield
    } else {
      if (length(t$twofield[t$csum/sum(t$N) > (1-frac)]) <= min_alleles_keep) {
        top <- t[order(-t$N),][1:min_alleles_keep,]$twofield
      } else {
        top <- t$twofield[t$csum/sum(t$N) > (1-frac)]
      }
    }
    # include the most common allele in the population even if not in the top alleles by fraction of all reads (from previous step)
    if (insert_pop_most_freq) {
      u <- unique(t$onefield[t$csum/sum(t$N) > (1-frac)])
      spl <- with(t, split(t, list(onefield=onefield)))
      spl <- spl[names(spl) %in% u]
      if (length(spl) > 0) {
        for (i in 1:length(spl)) {
          if (length(spl[[i]]$fscore[!is.na(spl[[i]]$fscore)]) > 0) {
            top <- c(top, spl[[i]]$twofield[spl[[i]]$fscore == max(spl[[i]]$fscore[!is.na(spl[[i]]$fscore)])])
          }
        }
      }
    }
    top <- top[!is.na(top)] #remove NAs that were introduced in previous step
    top <- unique(top) #remove any duplicates
    top_alleles <- c(top_alleles, top)
  }
  top_alleles<-top_alleles[order(top_alleles)]
  return(top_alleles)
}