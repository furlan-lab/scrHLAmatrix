#' Getting the top two alleles (for each HLA gene) with the most counts per Cell Barcode
#' 
#' @param reads  is the scrHLAtag count file including columns for CB, UMI, and HLA alleles (\url{https://github.com/furlan-lab/scrHLAtag}).
#' @param seu  is the Seurat object associated with the scrHLAtag count file (\url{https://satijalab.org/seurat/index.html}), and entered here if matching CBs in count file with Seurat colnames is desired.
#' @param CB_rev_com  logical, called \code{TRUE} if the need to obtain the reverse complement of Cell Barcodes (CBs) is desired; default is \code{FALSE}. 
#' @param hla_with_counts_above  the number of total reads accross CBs at or above which an HLA allele is retained in the matrix.
#' @param CBs_with_counts_above  the number of total reads accross HLA alleles at or above which a CB is retained in the matrix. Note: at present, the function will make sure that number of CBs is equal or more than available HLA alleles in the matrix.
#' @param match_CB_with_seu  logical, called \code{TRUE} if filtering CBs in the scrHLAtag count file with matching ones in the Seurat object is desired. 
#' @param cluster_index  a numeric or integer, representing the index of the HLA cluster present in the counts data (as previously analyzed by \code{HLA_clusters()} and mapped back to the data by \code{map_HLA_clusters()}). It allows to subset the visualization based on the selected HLA cluster. \code{NULL} is the default in which case the counts data will not be subsetted.
#' @param top_hla  the number of top HLA alleles (i.e. with the highest number of reads) per HLA gene to display in the plot; default is \code{10}.
#' @param field_resolution  integer fron \code{1} to \code{3}, to select the HLA nomenclature level of Field resolution, where \code{1}, \code{2}, or \code{3} will take into consideration the first, the first two, or the first three field(s) of HLA designation; default is \code{3}.
#' @param QC_mm2  logical, called \code{TRUE} if removing low quality reads based on minimap2 tags is desired.
#' @param s1_percent_pass_score  percentage, \code{0} to \code{100}, cuttoff from the maximum score (best quality) of the minimap2 's1' tag, which a read needs to acheive to pass as acceptable; default at \code{80} and becomes less inclusive if value increases.
#' @param AS_percent_pass_score  percentage, \code{0} to \code{100}, cuttoff from the maximum score (best quality) of the minimap2 'AS' tag, which a read needs to acheive to pass as acceptable; default at \code{80} and becomes less inclusive if value increases.
#' @param NM_thresh  number of mismatches and gaps in the minimap2 alignment at or below which the quality of the read is acceptable; default is \code{15}.
#' @param de_thresh  gap-compressed per-base sequence divergence in the minimap2 alignment at or below which the quality of the read is acceptable; the number is between \code{0} and \code{1}, and default is \code{0.01}.
#' @param default_theme  logical, called \code{TRUE} by default to output a 'plot_grid' of all HLA alleles in the count file in a single figure; if \code{FALSE} the output is a list of plots without a specified theme.
#' @param return_genotype_data  logical, called \code{TRUE} if the HLA genotypes of each CB, ordered in a dataframe, is desired in addition to the plots; \code{FALSE} by default. 
#' @param parallelize  logical, called \code{TRUE} if using parallel processing (multi-threading) is desired; default is \code{TRUE}.
#' @import stringr
#' @import cowplot
#' @import pbmcapply
#' @import parallel
#' @import Matrix
#' @import magrittr
#' @import htmltools
#' @import ggplot2
#' @import Seurat
#' @import dplyr
#' @return ggplot of HLA genotypes per Cell Barcodes and the number of times they occur
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
#' top_per_cb_plot <- Top_HLA_plot_byCB(reads = cts[["mRNA"]], seu = your_Seurat_Obj, hla_with_counts_above = 5, CBs_with_counts_above = 35)
#' top_per_cb_plot %>% gridExtra::grid.arrange()
#' # 
#' # Note: if for a particular HLA, the alleles with the most counts are in a tie 
#' # between 3 or more alleles in a particular Cell Barcode, we cannot know which 
#' # ones are the top two alleles, so that CB is not counted. This is similar if  
#' # there were no counts for that allele (all zeros).
#' @export

Top_HLA_plot_byCB <- function(reads, seu = NULL, CB_rev_com = FALSE, hla_with_counts_above = 5, CBs_with_counts_above = 15, match_CB_with_seu = TRUE, cluster_index = NULL, top_hla = 10, field_resolution = 3, QC_mm2 = TRUE, s1_percent_pass_score = 80, AS_percent_pass_score = 80, NM_thresh = 15, de_thresh = 0.01, default_theme = TRUE, return_genotype_data = FALSE, parallelize = TRUE) {
  ## parallelize
  if (parallelize) {
    multi_thread <- parallel::detectCores()
    message(cat("\nMulti-threading! Available cores: ", parallel::detectCores()))
  } else {
    multi_thread <- 1
  }
  ## if HLA clusters based on distribution on UMAP was analyzed, and visualizing top alleles per cluster is desired
  if (!is.null(cluster_index)) {
    if (!(is.numeric(cluster_index) | is.integer(cluster_index))) {stop("'cluster_index' must be numeric or integer", call. = FALSE)}
    if ("hla_clusters" %in% colnames(reads)) {
      cl <- unique(reads$hla_clusters) %>% as.factor()
      #cl <- cl[!is.na(cl)]
      idx <- cluster_index
      message(cat("\nClusters based on HLA genotype paterns in reads count object! Number of Clusters: ", length(levels(cl))))
      reads <- reads[reads$hla_clusters %in% levels(cl)[idx],]
    } else {
      message(cat(crayon::red("Warning: Colname 'hla_clusters' not detected in reads count object."), "\n  Did you analyze distribution of alleles per Cell Barcodes using ", crayon::bgWhite(" HLA_clusters() "), ", \n  then map the generated HLA Clusters back to your counts data object using ", crayon::bgWhite(" map_HLA_clusters() "), "?", sep = ""))
    }
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
    warning("HLA field resolution must be 1, 2, or 3 to take into consideration the first, the first two, or the first three field(s) of HLA designation.\nKeeping field resolution at 3 by default.")
  }
  ## Matrix formation
  alleles <- unique(reads$gene) %>% sort()
  reads$seu_barcode <- paste0(reads$samp,"_",reads$CB,"-1")
  reads <- split(data.table::setDT(reads), by = "seu_barcode")
  reads <- parallel::mclapply(reads, data.table::setDF, mc.cores = multi_thread)
  # matrix
  HLA.matrix <- matrix(0, nrow = length(alleles), ncol = length(reads), dimnames = list(alleles, names(reads)))
  message(cat("\nCreating an HLA Count Matrix"))
  pb <- pbmcapply::progressBar(min = 0, max = length(reads), style = "ETA", char = "=")
  for (i in 1:length(reads)) {
    counts <- table(reads[[i]]$gene)
    HLA.matrix[, i] <- counts[alleles]
    setTxtProgressBar(pb, i)
  }
  close(pb)
  rm(reads)
  HLA.matrix[is.na(HLA.matrix)]<-0
  #HLA.matrix<-Matrix(HLA.matrix,sparse = T)
  ## Matching with Seurat colnames
  if (is.null(seu)) {
    part_HLA<- HLA.matrix
  } else {
    if ("Seurat" %in% class(seu)) {
      message(cat("\nObject of class 'Seurat' detected"))
      message(cat(crayon::green("Note: "), "Currently the Seurat Barcode (i.e. colnames or Cells) supported format is: SAMPLE_AATGCTTGGTCCATTA-1", sep = ""))
      if(match_CB_with_seu) {
        part_HLA<- HLA.matrix[,colnames(HLA.matrix) %in% Cells(seu)]
      } else {
        part_HLA<- HLA.matrix
      }
    } else {
      stop("Single-cell dataset container (in argument 'seu') must be of class 'Seurat'", call. = FALSE)
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
  #-------------------------------------------------------------------------------
  message(cat("\nCounting Top Two alleles per HLA gene per Cell Barcode"))
  part_HLA <- part_HLA %>% as.data.frame()
  special <- "[-_*|?.+$^]"
  part_HLA[c("hla", "leftover")] <- stringr::str_split_fixed(row.names(part_HLA), special, 2)
  part_HLA$leftover <- NULL
  part_HLA <- with(part_HLA, split(part_HLA, list(hla=hla)))
  top2cb <- pbmclapply(1:length(part_HLA), function(j) {
    top2tab <- mclapply(1:(ncol(part_HLA[[j]])-1), function(i) {
      colval <- part_HLA[[j]][,i][part_HLA[[j]][,i]!=0] # column values that are greater than zero
      rnames <- row.names(part_HLA[[j]])[part_HLA[[j]][,i]!=0] # only row.names with column values greater than zero
      names(colval) <- rnames
      cname <- colnames(part_HLA[[j]])[i] # the name of the CB
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
        maxt <- data.frame("seu_barcode"=cname, "toptwo"=max1, "hla"=part_HLA[[j]]$hla[1])
        return(maxt)
      }
    }, mc.cores = 1) # more efficient if inner loop always running on a single CPU (mc.cores = 1)
    return(do.call("rbind", top2tab))
  }, mc.cores = multi_thread)
  top2cb <- do.call("rbind", top2cb)
  
  ## count the alleles in the CBs and prepare them for ggplot2
  message(cat("\nDrawing Plots"))
  tab <- data.frame("toptwo"=character(0), "N"=integer(0), "hlagene"=character(0), "csum"=integer(0))
  a<-unique(top2cb$hla) %>% sort()
  for (j in 1:length(a)){
    t<-as.data.table(table(top2cb[top2cb$hla == a[j],]$toptwo))
    t$toptwo <- t$V1
    t <- aggregate(N ~ toptwo, data = t, FUN = sum)
    t<- t[order(-t$N, t$toptwo),]
    row.names(t)<- NULL
    t$hlagene <- a[j]
    # order by N
    t<- t[order(t$N),]
    t$csum <- cumsum(t$N)
    #----------------------------------------------------------------------
    # prepare table for plots, just get the top 'n' alleles with most reads
    t<- t[order(-t$N, t$toptwo),]
    if (nrow(t)>top_hla){
      t<- t[1:top_hla,]
    }
    tab <- rbind(tab, t)
    tab <- tab[order(tab$hlagene),]
    row.names(tab) <- NULL
  }
  if (default_theme) {
    plots <- c()
    pb <- pbmcapply::progressBar(min = 0, max = length(unique(tab$hlagene)), style = "ETA", char = "=")
    for (j in 1:length(unique(tab$hlagene))) {
      t <- tab[tab$hlagene == unique(tab$hlagene)[j],]
      t <- t[order(-t$N),]
      row.names(t)<- NULL
      t$toptwo <- paste0(str_pad(row.names(t), 2, pad = "0"), "_",t$toptwo)
      g <- ggplot(t, aes(x= toptwo, y= N))+
        geom_bar(stat = 'identity')+
        xlab("Cell Barcodes (n)")+ 
        ylab("Top Genotype Combinations")+
        scale_x_discrete(label=function(x) sub('...', '', x))+
        theme(text = element_text(size = 9),legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
      plots<- c(plots,list(g))
      setTxtProgressBar(pb, j)
    }
    close(pb)
    pl <- do.call("plot_grid", c(plots, align = "hv", ncol=floor(sqrt(length(plots)))))
    y.grob <- grid::textGrob("Cell Barcodes (n)", gp=grid::gpar(fontsize=9), rot=90)
    x.grob <- grid::textGrob("Top Genotype Combinations", gp=grid::gpar(fontsize=9))
    pl0 <- gridExtra::grid.arrange(arrangeGrob(pl, left = y.grob, bottom = x.grob))
    message(cat("\nDone!! you can display your plots by calling gridExtra::grid.arrange()"))
    if (return_genotype_data) {
      return(list(Plots = pl0, Genotype_data_per_CB = top2cb))
    } else {
      return(pl0)
    }
  } else {
    plots <- c()
    pb <- pbmcapply::progressBar(min = 0, max = length(unique(tab$hlagene)), style = "ETA", char = "=")
    for (j in 1:length(unique(tab$hlagene))) {
      t <- tab[tab$hlagene == unique(tab$hlagene)[j],]
      t <- t[order(-t$N),]
      row.names(t)<- NULL
      t$toptwo <- paste0(str_pad(row.names(t), 2, pad = "0"), "_",t$toptwo)
      g <- ggplot(t, aes(x= toptwo, y= N))+
        geom_bar(stat = 'identity')+
        xlab("Cell Barcodes (n)")+ 
        ylab("Top Genotype Combinations")+
        scale_x_discrete(label=function(x) sub('...', '', x))+
        theme(legend.position = "none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
      plots<- c(plots,list(g))
      setTxtProgressBar(pb, j)
    }
    close(pb)
    message(cat("\nDone!! you can display your plots by calling gridExtra::grid.arrange()"))
    if (return_genotype_data) {
      return(list(Plots = plots, Genotype_data_per_CB = top2cb))
    } else {
      return(plots)
    }
  }
}





#' Plotting the top HLA alleles extracted in a Pseudo-Bulk approach from the scrHLAtag count files 
#' 
#' @param reads_1  is the primary scrHLAtag count file (1 of 2 files containing either the mRNA molecular info or the genomic (gene) molecular info). It includes columns for CB, UMI, and HLA alleles (\url{https://github.com/furlan-lab/scrHLAtag}).
#' @param reads_2  is the secondary scrHLAtag count file (the alternative file vs. the one designated in \code{reads_1} argument). It includes columns for CB, UMI, and HLA alleles (\url{https://github.com/furlan-lab/scrHLAtag}). Default is \code{NULL}, in which case it will not be able to count alternative aligment and argument \code{use_alt_align_ABC} will become irrelevant.
#' @param cluster_index  a numeric or integer, representing the index of the HLA cluster present in the counts data (as previously analyzed by \code{HLA_clusters()} and mapped back to the data by \code{map_HLA_clusters()}). It allows to subset the visualization based on the selected HLA cluster. \code{NULL} is the default in which case the counts data will not be subsetted.
#' @param top_hla  numeric, representing the number of top HLA alleles (i.e. with the highest number of reads) per HLA gene to display in the plot; default is \code{10}.
#' @param min_reads_per_gene  numeric representing minimum number of total reads per HLA gene (including all its alleles) below which the gene is filtered out; default is \code{20}. 
#' @param use_alt_align_ABC  logical, whether to use the count file from the alternative alignment (rather than the primary alignment) to count reads for the HLA-A, -B, and -C genes. It was observed in some cases that using genomic alignments has better accuracy in predicting genotype versus mRNA alignments (not the case for Class-II and other HLA genes); default is \code{FALSE}.
#' @param color_pal  is a character list of colors to visualize HLA polulation frequencies when available. When \code{color_pal} is not provided (\code{NULL}), it defaults to \code{viridis::viridis(n = 10, option = "C")}.
#' @import stringr
#' @import cowplot
#' @import magrittr
#' @import htmltools
#' @import ggplot2
#' @import dplyr
#' @return a Vector of the top HLA alleles in the count files (in terms of reads).
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
#' Top_HLA_plot_bulk(reads_1 = cts[["mRNA"]], reads_2 = cts[["gene"]], use_alt_align_ABC = TRUE)
#' @export

Top_HLA_plot_bulk <- function(reads_1, reads_2 = NULL, cluster_index = NULL, top_hla = 10, min_reads_per_gene = 20, use_alt_align_ABC = FALSE, color_pal = NULL){
  if (is.null(reads_2)) {
    if (!is.null(cluster_index)) {
      if (!(is.numeric(cluster_index) | is.integer(cluster_index))) {stop("'cluster_index' must be numeric or integer", call. = FALSE)}
      if ("hla_clusters" %in% colnames(reads_1)) {
        cl <- unique(reads_1$hla_clusters) %>% as.factor()
        idx <- cluster_index
        message(cat("\nClusters based on HLA genotype paterns in reads count object! Number of Clusters: ", length(levels(cl))))
        reads_1 <- reads_1[reads_1$hla_clusters %in% levels(cl)[idx],]
      } else {
        message(cat(crayon::red("Warning: Colname 'hla_clusters' not detected in reads count object."), "\n  Did you analyze distribution of alleles per Cell Barcodes in UMAP space using ", crayon::bgWhite(" HLA_clusters() "), ", \n  then map the generated HLA Clusters back to your counts data object using ", crayon::bgWhite(" map_HLA_clusters() "), "?", sep = ""))
      }
    } 
    reads_2 <- reads_1
    message(cat(crayon::red("Warning: "), "A secondary reads count file does not seem to be included (reads_2 = NULL).\n  The function will run but the argument 'use_alt_align_ABC' will be irrelevant.", sep = ""))
    #warning("The molecule_info_gene.txt.gz count file does not seem to be included. The function will run but the argument 'use_alt_align_ABC' will be irrelevant.")
  } else {
    if (!is.null(cluster_index)) {
      if (!(is.numeric(cluster_index) | is.integer(cluster_index))) {stop("'cluster_index' must be numeric or integer", call. = FALSE)}
      if ("hla_clusters" %in% colnames(reads_1)) {
        cl <- unique(reads_1$hla_clusters) %>% as.factor()
        idx <- cluster_index
        message(cat("Clusters based on HLA genotype paterns in",   crayon::bgWhite(" first "), "reads count object! Number of Clusters: ", length(levels(cl))))
        reads_1 <- reads_1[reads_1$hla_clusters %in% levels(cl)[idx],]
      } else {
        message(cat(crayon::red("Warning: Colname 'hla_clusters' not detected in first reads count object."), "\n  Did you analyze distribution of alleles per Cell Barcodes in UMAP space using ", crayon::bgWhite(" HLA_clusters() "), ", \n  then map the generated HLA Clusters back to your counts data object using ", crayon::bgWhite(" map_HLA_clusters() "), "?", sep = ""))
      }
      if ("hla_clusters" %in% colnames(reads_2)) {
        cl <- unique(reads_2$hla_clusters) %>% as.factor()
        idx <- cluster_index
        message(cat("Clusters based on HLA genotype paterns in",   crayon::bgWhite(" second "), "reads count object! Number of Clusters: ", length(levels(cl))))
        reads_2 <- reads_2[reads_2$hla_clusters %in% levels(cl)[idx],]
      } else {
        message(cat(crayon::red("Warning: Colname 'hla_clusters' not detected in second reads count object."), "\n  Did you analyze distribution of alleles per Cell Barcodes in UMAP space using ", crayon::bgWhite(" HLA_clusters() "), ", \n  then map the generated HLA Clusters back to your counts data object using ", crayon::bgWhite(" map_HLA_clusters() "), "?", sep = ""))
      }
    }
  } 
  # extract the HLA genes that appear in the reads
  special <- "[_*|?.+$^]"
  reads_2$gene0 <- gsub(special, "-", reads_2$gene)
  reads_2[c("hla", "leftover")] <- stringr::str_split_fixed(reads_2$gene, special, 2)
  reads_2$leftover <- NULL
  reads_1$gene0 <- gsub(special, "-", reads_1$gene)
  reads_1[c("hla", "leftover")] <- stringr::str_split_fixed(reads_1$gene, special, 2)
  reads_1$leftover <- NULL
  # get the HLA genes in a list
  cts_abc <- unique(reads_2$hla)[order(unique(reads_2$hla))]
  if (use_alt_align_ABC) {
    if (length(cts_abc[which(cts_abc %in% c("A", "B", "C"))]) != 3) {
      stop("the secondary molecule info count file does not contain alleles belonging to all of HLA-A, -B, and -C", call. = FALSE)
    }
  }
  cts_notabc <- unique(reads_1$hla)[order(unique(reads_1$hla))]
  tab <- data.frame("twofield"=character(0), "N"=integer(0), "fscore"=numeric(0), "hlagene"=character(0), "onefield"=character(0), "csum"=integer(0))
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
    if (strsplit(t$twofield, "\\*")[[1]][1] %in% c("MICA", "MICB")) {
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
    #---------------------------------------------------------------------
    # prepare table for plots, just get the top 10 alleles with most reads
    t<- t[order(-t$N, -t$fscore, t$twofield),]
    if (nrow(t)>top_hla){
      t<- t[1:top_hla,]
    }
    tab <- rbind(tab, t)
    tab <- tab[order(tab$hlagene),]
    row.names(tab) <- NULL
  }
  if (is.null(color_pal)){
    color_pal <- viridis::viridis(n = 10, option = "C")
  }
  message(cat("\nDrawing Plots..."))
  plots <- c()
  for (j in 1:length(unique(tab$hlagene))) {
    t <- tab[tab$hlagene == unique(tab$hlagene)[j],]
    t <- t[order(-t$N),]
    row.names(t)<- NULL
    t$twofield <- paste0(str_pad(row.names(t), 2, pad = "0"), "_",t$twofield)
    g <- ggplot(t, aes(x= twofield, y= N, fill= fscore))+
      geom_bar(stat = 'identity')+
      scale_fill_gradientn(limits = c(0,1), colours = rev(color_pal), na.value = "grey35")+
      xlab("Reads (n)")+ 
      ylab("Top 10 alleles")+
      scale_x_discrete(label=function(x) sub('...', '', x))+
      theme(text = element_text(size = 9),legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    plots<- c(plots,list(g))
  }
  pl <- do.call("plot_grid", c(plots, align = "hv", ncol=floor(sqrt(length(plots)))))
  leg <- get_legend(ggplot(tab, aes(x= twofield, y= N, fill= fscore, color=""))+
                      geom_bar(stat = 'identity')+
                      scale_fill_gradientn(limits = c(0,1), colours = rev(color_pal), breaks=c(0,0.25,0.5,0.75,1),labels=c(0.00,0.25,0.50,0.75,1.00), na.value = "grey35")+
                      scale_color_manual(values=NA, na.value = NA)+
                      labs(fill="Population\nFreq", color = "No data")+
                      #guides(colour=guide_legend(override.aes=list(colour="grey70")))+
                      theme(text = element_text(size = 8), legend.key.size = unit(4, 'mm'), legend.title = element_text(size=6))
  )
  pl2 <- plot_grid(pl, gridExtra::grid.arrange(leg), rel_widths = c(10,1)) 
  print(pl2)
  message(cat("\nDone!! "))
  return(pl2)
}

