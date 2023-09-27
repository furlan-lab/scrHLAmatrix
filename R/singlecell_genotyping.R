#' Getting the top two alleles (for each HLA gene) with the most counts per Cell Barcode
#' 
#' @param reads  is the scrHLAtag count file including columns for CB, UMI, and HLA alleles (https://github.com/furlan-lab/scrHLAtag).
#' @param seu  is the Seurat object associated with the scrHLAtag count file (https://satijalab.org/seurat/index.html).
#' @param CB_rev_com  is a logical, called TRUE if the need to obtained the reverse complement of Cell Barcodes (CBs) is desired; default is FALSE. 
#' @param hla_with_counts_above  is the number of total reads accross CBs at or above which an HLA allele is retained in the matrix.
#' @param CBs_with_counts_above  is the number of total reads accross HLA alleles at or above which a CB is retained in the matrix. Note: at present, the function will make sure that number of CBs is equal or more than available HLA alleles in the matrix.
#' @param match_CB_with_seu  is a logical, called TRUE if filtering CBs in the scrHLAtag count file with matching ones in the Seurat object is desired. 
#' @param cluster_index  is a numeric, representing the index of the HLA cluster present in the counts data (as previously analyzed by 'HLA_clusters()' and mapped back to the data by 'map_HLA_clusters()'). It allows to subset the visualization based on the selected HLA cluster. 'NULL' is the default in which case the counts data will not be subsetted.
#' @param top_hla  is a numeric, representing the number of top HLA alleles (i.e. with the highest number of reads) per HLA gene to display in the plot; default is 10.
#' @param field_resolution  is a numeric, to select the HLA nomenclature level of Field resolution, where 1, 2, or 3 will take into consideration the first, the first two, or the first three field(s) of HLA designation; default is 3.
#' @param QC_mm2  is a logical, called TRUE if removing low quality reads based on minimap2 tags is desired.
#' @param s1_belowmax  is a proportion (0 to 1) of the maximum value (best quality) of the minimap2 's1' tag above which the quality of the read is acceptable; default at 0.75 of the max s1 score.
#' @param AS_belowmax  is a proportion (0 to 1) of the maximum value (best quality) of the minimap2 'AS' tag above which the quality of the read is acceptable; default at 0.85 of the max AS score.
#' @param NM_thresh  is the number of mismatches and gaps in the minimap2 alignment at or below which the quality of the read is acceptable; default is 15.
#' @param de_thresh  is the gap-compressed per-base sequence divergence in the minimap2 alignment at or below which the quality of the read is acceptable; the number is between 0 and 1, and default is 0.015.
#' @param default_theme  is a logical, called TRUE by default to output a 'plot_grid' of all HLA alleles in the count file in a single figure; if FALSE the output is a list of plots without a specified theme.
#' @param return_genotype_data  is a logical, called TRUE if the HLA genotypes of each CB, ordered in a dataframe, is desired in addition to the plots; FALSE by default. 
#' @param parallelize  is a logical, called TRUE if using parallel processing (multi-threading) is desired; default is TRUE.
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
#' HLA_alleles_per_CB(reads = cts[["mRNA"]], seu = your_Seurat_Obj, hla_with_counts_above = 5, CBs_with_counts_above = 35)
#' # 
#' # Note: if for a particular HLA, the alleles with the most counts are in a tie 
#' # between 3 or more alleles in a particular Cell Barcode, we cannot know which 
#' # ones are the top two alleles, so that CB is not counted. This is similar if  
#' # there were no counts for that allele (all zeros).
#' @export

HLA_alleles_per_CB <- function(reads, seu = NULL, CB_rev_com = FALSE, hla_with_counts_above = 0, CBs_with_counts_above = 0, match_CB_with_seu = TRUE, cluster_index = NULL, top_hla = 10, field_resolution = 3, QC_mm2 = TRUE, s1_belowmax = 0.75, AS_belowmax = 0.85, NM_thresh = 15, de_thresh = 0.015, default_theme = TRUE, return_genotype_data = FALSE, parallelize = TRUE) {
  ## parallelize
  if (parallelize) {
    multi_thread <- parallel::detectCores()
    message(cat("\nMulti-threading! Available cores: ", parallel::detectCores()))
  } else {
    multi_thread <- 1
  }
  ## if HLA clusters based on distribution on UMAP was analyzed, and visualizing top alleles per cluster is desired
  if (!is.null(cluster_index)) {
    if ("hla_clusters" %in% colnames(reads)) {
      cl <- unique(reads$hla_clusters)
      cl <- cl[!is.na(cl)]
      idx <- cluster_index
      message(cat("\nClusters generated from HLA distribution per Cell Barcode in UMAP space detected! Number of Clusters: ", length(cl)))
      reads <- reads[reads$hla_clusters %in% cl[idx],]
    } else {
      warning("  Colname 'hla_clusters' not detected in counts data.\n  Did you analyze distribution of alleles per Cell Barcodes in UMAP space using 'HLA_clusters()', \n  then map the generated HLA Clusters back to your counts data using 'map_HLA_clusters()'?")
    }
  }  
  ## Remove low quality reads based on minimap2 tags
  if (QC_mm2) {
    message(cat("\nRemoving low quality reads based on minimap2 tags"))
    reads <- with(reads, split(reads, list(gene=gene)))
    reads <- pbmclapply(reads, function(df){df[df$s1 > s1_belowmax*max(df$s1) & df$AS > AS_belowmax*max(df$AS) & df$NM <= NM_thresh & df$de <= de_thresh,]}, mc.cores = multi_thread)
    reads <-  do.call("rbind", reads)
    row.names(reads)<-NULL
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
  reads <- with(reads, split(reads, list(seu_barcode=seu_barcode)))
  # matrix
  HLA.matrix <- matrix(0, nrow = length(alleles), ncol = length(reads), dimnames = list(alleles, names(reads)))
  message(cat("\nCreating an HLA Count Matrix"))
  pb <- txtProgressBar(min = 0, max = length(reads), style = 3, char = "=")
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
    if (class(seu) == "Seurat") {
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
    }
    
    pl <- do.call("plot_grid", c(plots, align = "hv", ncol=floor(sqrt(length(plots)))))
    y.grob <- grid::textGrob("Cell Barcodes (n)", gp=grid::gpar(fontsize=9), rot=90)
    x.grob <- grid::textGrob("Top Genotype Combinations", gp=grid::gpar(fontsize=9))
    pl0 <- grid.arrange(arrangeGrob(pl, left = y.grob, bottom = x.grob))
    message(cat("\nDone!!"))
    if (return_genotype_data) {
      return(list(Plots = pl0, Genotype_data_per_CB = top2cb))
    } else {
      return(pl0)
    }
  } else {
    plots <- c()
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
    }
    message(cat("\nDone!!"))
    if (return_genotype_data) {
      return(list(Plots = plots, Genotype_data_per_CB = top2cb))
    } else {
      return(plots)
    }
  }
}


#' Extracting the top HLA alleles from the scrHLAtag count files in a Single-Cell approach (based on the most counts per Cell Barcode)
#' 
#' @param reads  is the scrHLAtag count file including columns for CB, UMI, and HLA alleles (https://github.com/furlan-lab/scrHLAtag).
#' @param seu  is the Seurat object associated with the scrHLAtag count file (https://satijalab.org/seurat/index.html).
#' @param CB_rev_com  is a logical, called TRUE if the need to obtained the reverse complement of Cell Barcodes (CBs) is desired; default is FALSE. 
#' @param hla_with_counts_above  is the number of total reads accross CBs at or above which an HLA allele is retained in the matrix.
#' @param CBs_with_counts_above  is the number of total reads accross HLA alleles at or above which a CB is retained in the matrix. Note: at present, the function will make sure that number of CBs is equal or more than available HLA alleles in the matrix.
#' @param match_CB_with_seu  is a logical, called TRUE if filtering CBs in the scrHLAtag count file with matching ones in the Seurat object is desired. 
#' @param frac  is the fraction (0 to 1) of total reads for a particular HLA gene, which incorporates the highest ranking alleles of that gene in terms of number of reads; default at 0.85 .
#' @param min_alleles_keep  is a numeric representing the minimum number of highest ranking alleles to keep despite filtering by fraction 'frac'; default is 1.
#' @param field_resolution  is a numeric, to select the HLA nomenclature level of Field resolution, where 1, 2, or 3 will take into consideration the first, the first two, or the first three field(s) of HLA designation; default is 3.
#' @param QC_mm2  is a logical, called TRUE if removing low quality reads based on minimap2 tags is desired.
#' @param s1_belowmax  is a proportion (0 to 1) of the maximum value (best quality) of the minimap2 's1' tag above which the quality of the read is acceptable; default at 0.75 of the max s1 score.
#' @param AS_belowmax  is a proportion (0 to 1) of the maximum value (best quality) of the minimap2 'AS' tag above which the quality of the read is acceptable; default at 0.85 of the max AS score.
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
#' @return a Vector of the top HLA alleles in the count files (in terms of reads per Cell Barcode).
#' @examples
#' samples <- c("AML_101_BM", "AML_101_34")
#' mol_info <- c("molecule_info_gene.txt.gz", "molecule_info_mRNA.txt.gz")
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
#' Top_HLA_list(reads = cts[["mRNA"]], seu = your_Seurat_Obj, hla_with_counts_above = 5, CBs_with_counts_above = 35, frac = 0.9, min_alleles_keep = 2)
#' # 
#' # Note: if for a particular HLA, the alleles with the most counts are in a tie 
#' # between 3 or more alleles in a particular Cell Barcode, we cannot know which 
#' # ones are the top two alleles, so that CB is not counted. This is similar if  
#' # there were no counts for that allele (all zeros).
#' @export

Top_HLA_list <- function(reads, seu = NULL, CB_rev_com = FALSE, hla_with_counts_above = 0, CBs_with_counts_above = 0, match_CB_with_seu = TRUE, frac = 0.85, min_alleles_keep = 1, field_resolution = 3, QC_mm2 = TRUE, s1_belowmax = 0.75, AS_belowmax = 0.85, NM_thresh = 15, de_thresh = 0.015, parallelize = TRUE) {
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
    reads <- pbmclapply(reads, function(df){df[df$s1 > s1_belowmax*max(df$s1) & df$AS > AS_belowmax*max(df$AS) & df$NM <= NM_thresh & df$de <= de_thresh,]}, mc.cores = multi_thread)
    reads <-  do.call("rbind", reads)
    row.names(reads)<-NULL
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
    message(cat("  Suggestion: To refine your results, first analyze distribution of alleles per Cell Barcodes in UMAP space using 'HLA_clusters()', \n  then map the generated HLA Clusters back to your count data using 'map_HLA_clusters()'."))
  }  
  matrices <- mclapply(1:length(reads), function(m) {
    alleles <- unique(reads[[m]]$gene) %>% sort()
    reads[[m]]$seu_barcode <- paste0(reads[[m]]$samp,"_",reads[[m]]$CB,"-1")
    reads[[m]] <- with(reads[[m]], split(reads[[m]], list(seu_barcode=seu_barcode)))
    HLA.matrix <- matrix(0, nrow = length(alleles), ncol = length(reads[[m]]), dimnames = list(alleles, names(reads[[m]])))
    pb <- txtProgressBar(min = 0, max = length(reads[[m]]), style = 3, char = "=")
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
      if (class(seu) == "Seurat") {
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
    return(part_HLA)
  }, mc.cores = 1) # more efficient if running on 1 CPU
  rm(reads)
  #-------------------------------------------------------------------------------
  message(cat("\nCounting Top Two alleles per HLA gene per Cell Barcode"))
  top2list <- mclapply(1:length(matrices), function(m) {
    matrices[[m]] <- matrices[[m]] %>% as.data.frame()
    special <- "[-_*|?.+$^]"
    matrices[[m]][c("hla", "leftover")] <- stringr::str_split_fixed(row.names(matrices[[m]]), special, 2)
    matrices[[m]]$leftover <- NULL
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
      if (length(t$toptwo) <= min_alleles_keep) {
        top <- t$toptwo
      } else {
        if (length(t$toptwo[t$csum/sum(t$N) > (1-frac)]) <= min_alleles_keep) {
          top <- t[order(-t$N),][1:min_alleles_keep,]$toptwo
        } else {
          top <- t$toptwo[t$csum/sum(t$N) > (1-frac)]
        }
      }
      top_a <- c(top_a, top)
    }
    top_a <- sapply(top_a, function(x) strsplit(x, "_")[[1]]) %>% as.character()
    top_a <- unique(top_a)
    top_a <- top_a[top_a != "NA"]
    top_a <- top_a %>% sort()
    return(top_a)
  }, mc.cores = multi_thread)
  top_a_list <- do.call("c", top_a_list)
  top_a_list <- top_a_list %>% unique() %>% sort()
  message(cat("\nDone!!"))
  return(top_a_list)
}
