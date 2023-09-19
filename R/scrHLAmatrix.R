#' Deduping and Correcting scrHLAtag counts and Creating Seurat-compatible Matrices
#' 
#' @param reads  is the scrHLAtag count file including columns for CB, UMI, and HLA alleles (https://github.com/furlan-lab/scrHLAtag).
#' @param seu  is the Seurat object associated with the scrHLAtag count file (https://satijalab.org/seurat/index.html).
#' @param hla_recip  is a character list of recipient-specific HLA alleles if known; default is an empty character vector.
#' @param hla_donor  is a character list of donor-specific HLA alleles if known; default is an empty character vector.
#' @param QC_mm2  is a logical, called TRUE if removing low quality reads based on minimap2 tags is desired.
#' @param s1_belowmax  is a proportion (0 to 1) of the maximum value (best quality) of the minimap2 's1' tag above which the quality of the read is acceptable; default at 0.75 of the max s1 score.
#' @param AS_belowmax  is a proportion (0 to 1) of the maximum value (best quality) of the minimap2 'AS' tag above which the quality of the read is acceptable; default at 0.85 of the max AS score.
#' @param NM_thresh  is the number of mismatches and gaps in the minimap2 alignment at or below which the quality of the read is acceptable; default is 15.
#' @param de_thresh  is the gap-compressed per-base sequence divergence in the minimap2 alignment at or below which the quality of the read is acceptable; the number is between 0 and 1, and default is 0.015.
#' @param parallelize  is a logical, called TRUE if using parallel processing (multi-threading) is desired; default is TRUE.
#' @param CB_rev_com  is a logical, called TRUE if the need to obtained the reverse complement of Cell Barcodes (CBs) is desired; default is FALSE. 
#' @param Ct  is the count threshold for the PCR copies of UMIs to retain; default is 0.
#' @import stringr
#' @import pbmcapply
#' @import parallel
#' @import Matrix
#' @import magrittr
#' @import htmltools
#' @import ggplot2
#' @import Seurat
#' @import dplyr
#' @return an Assay type matrix
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
#' HLA_Matrix(reads = cts[["mRNA"]], seu = your_Seurat_obj, hla_recip = c("A*24:02:01", "DRB1*04:01:01", "DRB4*01:03:02"), hla_donor = c("A*33:03:01", "B*42:01:01"))
#' @export

HLA_Matrix <- function(reads, seu, hla_recip = character(), hla_donor = character(), QC_mm2 = TRUE, s1_belowmax = 0.75, AS_belowmax = 0.85, NM_thresh = 15, de_thresh = 0.015, parallelize = TRUE, CB_rev_com = FALSE, Ct = 0) {
  ## check Seurat object
  if (class(seu) != "Seurat") {
    stop("Single-cell dataset container must be of class 'Seurat'")
  }
  ## parallelize
  if (parallelize) {
    multi_thread <- parallel::detectCores()
    message(cat("\nMulti-threading! Available cores: ", parallel::detectCores()))
  } else {
    multi_thread <- 1
  }
  ## add the dash into the HLA allele as it is the one accepted in the names of the Seurat assay features
  special <- "[_*|?.+$^]"
  if (any(grepl(special, reads$gene))) {
    reads$gene0 <- gsub(special, "-", reads$gene)
    reads[c("hla", "leftover")] <- str_split_fixed(reads$gene, special, 2)
    reads$leftover <- NULL
    message(cat("\nAvailable reads per gene"))
    print(table(reads$hla, useNA = "ifany"))
  } else if (all(grepl("-", reads$gene))){
    reads$gene0 <- reads$gene
    reads[c("hla", "leftover")] <- str_split_fixed(reads$gene, "-", 2)
    reads$leftover <- NULL
    message(cat("\nAvailable reads per gene"))
    print(table(reads$hla, useNA = "ifany"))
  } else {
    stop("The HLA allele column is unrecognizable or has incorrect format. \nMake sure gene and allele are separated by standard nomenclature asterisk (or other special character).")
  }
  reads$cbumi <- paste0(reads$CB, ":", reads$UMI)
  ## check format of 'hla_recip' and 'hla_donor'
  if (any(grepl(special, hla_recip))) {
    hla_recip <- gsub(special, "-", hla_recip)
  } else if (all(grepl("-", hla_recip))) {
    hla_recip <- hla_recip
  } else {
    stop("Incorrect format for recipient-defined and/or donor-defined HLA alleles. \nMake sure gene and allele are separated by standard nomenclature asterisk (or other special character).")
  }
  if (any(grepl(special, hla_donor))) {
    hla_donor <- gsub(special, "-", hla_donor)
  } else if (all(grepl("-", hla_donor))) {
    hla_donor <- hla_donor
  } else {
    stop("Incorrect format for recipient-defined and/or donor-defined HLA alleles. \nMake sure gene and allele are separated by standard nomenclature asterisk (or other special character).")
  }
  ## check class of 'Ct'
  if (class(Ct) != "integer"){
    if (class(Ct) != "numeric"){
      stop("Count threshold 'Ct' must be numeric or integer")
    }
  }
  ## check class of 'NM_thresh'
  if (class(NM_thresh) != "integer"){
    if (class(NM_thresh) != "numeric"){
      stop("NM tag threshold 'NM_thresh' must be numeric or integer")
    }
  }
  ## check 's1_belowmax', 'AS_belowmax', and 'de_thresh' thresholds are between 0 and 1
  if (!all(between(c(s1_belowmax, AS_belowmax, de_thresh), 0, 1))){
    stop("'s1_belowmax', 'AS_belowmax', and 'de_thresh' should be numerics between 0 and 1")
  }
  ## check relevant col names exist
  if (!all(c("CB", "UMI", "gene", "samp") %in% colnames(reads))){
    stop("scrHLAtag output 'reads' dataframe must at least contain the columns 'CB', 'UMI', 'gene' (with HLA alleles), and 'samp' (matching the sample names in the corresponding Seurat object)")
  }
  ## Reverse Complement the CB
  if (CB_rev_com) {
    message(cat("\nConverting Cell Barcodes to their reverse complements"))
    # reads$CB <- pbmclapply(reads$CB, function(x) as.character(Biostrings::reverseComplement(DNAString(x))), mc.cores = multi_thread) %>% unlist() # slow
    reads$CB <- pbmclapply(reads$CB, function(x) intToUtf8(rev(utf8ToInt(chartr('ATGC', 'TACG', x)))), mc.cores = multi_thread) %>% unlist()        # fast
  } 
  reads$seu_barcode <- paste0(reads$samp,"_",reads$CB,"-1")
  message(cat("\nProportions of Cell Barcodes found (TRUE) or not found (FALSE) in the Seurat object colnames: "))
  print(reads$seu_barcode %in% colnames(seu) %>% table() / dim(reads)[1])
  ## Remove low quality reads based on minimap2 tags
  if (QC_mm2) {
    message(cat("\nRemoving low quality reads based on minimap2 tags"))
    reads <- with(reads, split(reads, list(gene0=gene0)))
    reads <- pbmclapply(reads, function(df){df[df$s1 > s1_belowmax*max(df$s1) & df$AS > AS_belowmax*max(df$AS) & df$NM <= NM_thresh & df$de <= de_thresh,]}, mc.cores = multi_thread)
    reads <-  do.call("rbind", reads)
    row.names(reads)<-NULL
  } else {
    reads <- reads[order(reads$gene0), ]
  }
  ## see if more than 1 allele are present per umi at a time
  # count all the problematic CB:UMIs for which a molecular swap is suspected
  reads$mol_swap <- NA
  reads$mol_swap <- as.factor(reads$mol_swap)
  reads$class_swap <- NA
  reads$class_swap <- as.factor(reads$class_swap)
  alleles <- unique(reads$gene0)  
  # split,   this is computationally heavy (about 10min for 10M rows)
  reads <- with(reads, split(reads, list(cbumi=cbumi))) 
  # message(cat("\nEstimating UMI duplication rate"))
  # n_umi<-pbmclapply(reads, nrow, mc.cores = multi_thread) %>% unlist()
  # umi_counts<- data.frame(n_umi)
  # umi_counts$dummy <- 1
  # umi_counts <- umi_counts[order(umi_counts$n_umi),]
  # row.names(umi_counts)<- NULL
  # g <- ggplot(umi_counts, aes(x= as.numeric(row.names(umi_counts)), y=n_umi, group= dummy))+
  #   #geom_smooth(size=2, method = "gam")+
  #   geom_line()+
  #   scale_y_log10(name = "PCR copies per UMI")+
  #   scale_x_continuous(name = "Rank", n.breaks = 8) +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  # print(g)
  message(cat("\nEstimating molecular swap"))
  pb <- txtProgressBar(min = 0, max = length(reads), style = 3, char = "=")
  for(j in 1:length(reads)){
    reads[[j]]$mol_swap <- ifelse(length(unique(reads[[j]]$gene)) > 1, 
                                      reads[[j]]$mol_swap <- "yes",
                                      reads[[j]]$mol_swap <- "no")
    setTxtProgressBar(pb, j)
  }
  close(pb)
  cl1<-c("A" , "B" , "C", "E" , "F" , "G")
  cl2<-c("DRA", "DRB1", "DRB3", "DRB4", "DRB5", "DQA1", "DQA2", "DQB1", "DQB2", "DPA1", "DPA2", "DPB1", "DPB2", "DMA", "DMB", "DOA", "DOB")
  for(j in 1:length(reads)){
    reads[[j]]$class_swap <- ifelse(any(reads[[j]]$hla %in% cl1) & any(reads[[j]]$hla %in% cl2), 
                                        reads[[j]]$class_swap <- "yes",
                                        reads[[j]]$class_swap <- "no")
    setTxtProgressBar(pb, j)
  }
  close(pb)
  # count the molecular swap rate 
  mol_swap_rate_per_read <- length(which(sapply(reads, function(df) "yes" %in% df$mol_swap))) / length(which(sapply(reads, function(df) nrow(df)>1)))
  mol_swap_rate_per_read
  class_swap_rate_per_read<-length(which(sapply(reads, function(df) "yes" %in% df$class_swap))) / length(which(sapply(reads, function(df) nrow(df)>1)))
  class_swap_rate_per_read
  intraclass_swap_rate<- mol_swap_rate_per_read - class_swap_rate_per_read
  intraclass_swap_rate
  message(cat("\n  intra-class molecular swap rate per UMI ", 
            format(round(100*intraclass_swap_rate, 3), nsmall = 1),
            "% of Cells\n  inter-class molecular swap rate per UMI ",
            format(round(100*class_swap_rate_per_read, 3), nsmall = 1),
            "% of Cells"))  
  ## Count Threshold 'Ct' above which the count of UMI copies becomes acceptable  
  if (Ct > 0) {
    reads <- reads[sapply(reads, nrow) > Ct]
  }
  ## Function to keep the most occurring HLA when more than 1 HLA is present per cb:umi
  keep_one <- function(df) {
    n_hla <- table(df$gene0)
    if (length(n_hla) > 1) {
      most_hla <- names(which.max(n_hla))
      most_hla_count <- n_hla[most_hla]
      if (sum(n_hla == most_hla_count) > 1) {
        # If there is a tie, remove the entire df
        df <- NULL
      } else {
        # Filter the df to keep only rows with the most occurring hla
        df <- df[df$gene0 == most_hla, ]
      }
    }
    return(df)
  }
  # Applying the function
  message(cat("\nCorrecting Molecular Swap: keeping the reads per UMI with the HLA allele occuring the most"))
  reads <- pbmclapply(reads, keep_one, mc.cores = multi_thread)
  ## Performing Dedup
  message(cat("\nPerforming Dedup on UMIs: removing PCR duplicates"))
  reads <- pbmclapply(reads, function(df){df[1,]}, mc.cores = multi_thread)
  reads <-  do.call("rbind", reads)
  row.names(reads)<-NULL
  ## Clean-up HLA conflicts per CB
  message(cat("\nEstimating Genotype Conflicts: assuming a cell cannot have both recipient and donor-origin HLA allele"))
  # remove obsolete cols and add Seurat barcode
  reads$mol_swap <- NULL
  reads$class_swap <- NULL
  reads$seu_barcode <- paste0(reads$samp,"_",reads$CB,"-1")
  reads$hla_conflict <- NA
  reads$hla_conflict <- as.factor(reads$hla_conflict)
  # split by Seurat barcode
  reads <- with(reads, split(reads, list(seu_barcode=seu_barcode)))
  # detect HLA conflicts (i.e. donor-spec and recipient-spec HLA in the same barcode)
  pb <- txtProgressBar(min = 0, max = length(reads), style = 3, char = "=")
  for(j in 1:length(reads)){
    reads[[j]]$hla_conflict <- ifelse(any(reads[[j]]$gene0 %in% hla_recip) & any(reads[[j]]$gene0 %in% hla_donor),
                                             reads[[j]]$hla_conflict <- "yes",
                                             reads[[j]]$hla_conflict <- "no")
    setTxtProgressBar(pb, j)
  }
  close(pb)
  hla_conflict_rate <- length(which(sapply(reads, function(df) "yes" %in% df$hla_conflict))) / length(reads)
  message(cat("  Conflicting HLA (both donor and recipient HLA within the same Cell Barcode) affects ", 
            format(round(100*hla_conflict_rate, 3), nsmall = 1),
            "% of Cells"))
  # function to clean-up HLA conflicts by Seurat barcode
  keep_two <- function(df, recip, donor) {
    if (any(df$gene0 %in% recip) & any(df$gene0 %in% donor)){
      # determine the count difference between recipient-specific HLA-associated UMIs & donor-specific HLA-associated UMIs
      diff <- sum(table(df[which(df$gene0 %in% donor),]$gene0))-sum(table(df[which(df$gene0 %in% recip),]$gene0))
      if (diff < 0){
        # if more recipient-specific UMIs, remove donor UMIs
        df <- df[-which(df$gene0 %in% donor),]
      } else if (diff > 0) {
        # if more donor-specific UMIs, remove recipient UMIs
        df <- df[-which(df$gene0 %in% recip),]
      } else {
        # in case of tie, remove UMIs with both recipient and donor alleles, leaving only the shared alleles
        df <- df[-which(df$gene0 %in% donor),]
        df <- df[-which(df$gene0 %in% recip),]
      }
    }
    return(df)
  }
  # keep two alleles per HLA gene!
  if (hla_conflict_rate == 0) {
    message(cat("\nResolving Genotype Conflicts: no conflicts to resolve"))
  } else {
    message(cat("\nResolving Genotype Conflicts: keeping recipient-origin or donor-origin HLA alleles per Cell, whichever are the most occuring"))
    reads <- pbmclapply(reads, keep_two, recip = hla_recip, donor = hla_donor, mc.cores = multi_thread)
  }
  ## Matrix formation
  # matrix
  HLA.matrix <- matrix(0, nrow = length(alleles), ncol = length(reads), dimnames = list(alleles, names(reads)))
  message(cat("\nCreating the HLA Count Matrix compatible with the Seurat object"))
  pb <- txtProgressBar(min = 0, max = length(reads), style = 3, char = "=")
  for (i in 1:length(reads)) {
    counts <- table(reads[[i]]$gene0)
    HLA.matrix[, i] <- counts[alleles]
    setTxtProgressBar(pb, i)
  }
  close(pb)
  # rm(reads)
  # merging with Seurat cell names
  HLA.matrix<-as.data.frame(t(HLA.matrix))
  HLA.matrix <- HLA.matrix[rownames(HLA.matrix) %in% colnames(seu),]
  cells <- data.frame(row.names = colnames(seu))
  HLA.matrix<- merge(HLA.matrix, cells, by= "row.names", all=T)
  rownames(HLA.matrix) <- HLA.matrix$Row.names
  HLA.matrix <- subset(HLA.matrix, select = -1)
  HLA.matrix <- HLA.matrix[order(rownames(HLA.matrix)),]
  HLA.matrix[is.na(HLA.matrix)]<-0
  HLA.matrix<-as.matrix(t(HLA.matrix))
  HLA.matrix<-Matrix(HLA.matrix,sparse = T)
  HLA <- CreateAssayObject(counts = HLA.matrix)
  message(cat("\nDone!!"))
  return(HLA)
}



#' Extracting the top HLA alleles from the scrHLAtag count files
#' 
#' @param reads_1  is the primary scrHLAtag count file (1 of 2 files containing either the mRNA molecular info or the genomic (gene) molecular info). It includes columns for CB, UMI, and HLA alleles (https://github.com/furlan-lab/scrHLAtag).
#' @param reads_2  is the secondary scrHLAtag count file (the alternative file vs. the one designated in 'reads_1' argument). It includes columns for CB, UMI, and HLA alleles (https://github.com/furlan-lab/scrHLAtag). Default is NULL, in which case it will not be able to count alternative aligment and argument 'use_alt_align_ABC' will become irrelevant.
#' @param frac  is the fraction (0 to 1) of total reads for a particular HLA gene, which incorporates the highest ranking alleles of that gene in terms of number of reads; default at 0.75 .
#' @param min_alleles_keep  is a numeric representing the minimum number of highest ranking alleles to keep despite filtering by fraction 'frac'; default is 5.
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
#' Top_HLA_list(reads_1 = cts[["mRNA"]], reads_2 = cts[["gene"]], frac = 0.8, min_alleles_keep = 2, use_alt_align_ABC = TRUE)
#' @export

Top_HLA_list <- function(reads_1, reads_2 = NULL, frac = 0.75, min_alleles_keep = 5, min_reads_per_gene = 20, insert_pop_most_freq = TRUE, use_alt_align_ABC = FALSE){
  if (is.null(reads_2)) {
    reads_2 <- reads_1
  }
  # extract the HLA genes that appear in the reads
  special <- "[_*|?.+$^]"
  reads_2$gene0 <- gsub(special, "-", reads_2$gene)
  reads_2[c("hla", "leftover")] <- str_split_fixed(reads_2$gene, special, 2)
  reads_2$leftover <- NULL
  reads_1$gene0 <- gsub(special, "-", reads_1$gene)
  reads_1[c("hla", "leftover")] <- str_split_fixed(reads_1$gene, special, 2)
  reads_1$leftover <- NULL
  # get the HLA genes in a list
  cts_abc <- unique(reads_2$hla)[order(unique(reads_2$hla))]
  if (use_alt_align_ABC) {
    if (length(cts_abc[which(cts_abc %in% c("A", "B", "C"))]) != 3) {
      stop("the secondary molecule info count file does not contain alleles belonging to all of HLA-A, -B, and -C")
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



#' Plotting the top HLA alleles from the scrHLAtag count files
#' 
#' @param reads_1  is the primary scrHLAtag count file (1 of 2 files containing either the mRNA molecular info or the genomic (gene) molecular info). It includes columns for CB, UMI, and HLA alleles (https://github.com/furlan-lab/scrHLAtag).
#' @param reads_2  is the secondary scrHLAtag count file (the alternative file vs. the one designated in 'reads_1' argument). It includes columns for CB, UMI, and HLA alleles (https://github.com/furlan-lab/scrHLAtag). Default is NULL, in which case it will not be able to count alternative aligment and argument 'use_alt_align_ABC' will become irrelevant.
#' @param top_hla  is a numeric, representing the number of top HLA alleles (i.e. with the highest number of reads) per HLA gene to display in the plot; default is 10.
#' @param min_reads_per_gene  is a numeric representing minimum number of total reads per HLA gene (including all its alleles) below which the gene is filtered out; default is 20. 
#' @param use_alt_align_ABC  is a logical, whether to use the count file from the alternative alignment (rather than the primary alignment) to count reads for the HLA-A, -B, and -C genes. It was observed in some cases that using genomic alignments has better accuracy in predicting genotype versus mRNA alignments (not the case for Class-II and other HLA genes); default is FALSE.
#' @param color_pal  is a character list of colors to visualize HLA polulation frequencies when available. When 'color_pal' is not provided (NULL), it defaults to viridis::viridis(n = 10, option = "C").
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
#' Top_HLA_plot(reads_1 = cts[["mRNA"]], reads_2 = cts[["gene"]], use_alt_align_ABC = TRUE)
#' @export

Top_HLA_plot <- function(reads_1, reads_2 = NULL, top_hla = 10, min_reads_per_gene = 20, use_alt_align_ABC = FALSE, color_pal = NULL){
  if (is.null(reads_2)) {
    reads_2 <- reads_1
    warning("The molecule_info_gene.txt.gz count file does not seem to be included. The function will run but the argument 'use_alt_align_ABC' will be irrelevant.")
  }  
  # extract the HLA genes that appear in the reads
  special <- "[_*|?.+$^]"
  reads_2$gene0 <- gsub(special, "-", reads_2$gene)
  reads_2[c("hla", "leftover")] <- str_split_fixed(reads_2$gene, special, 2)
  reads_2$leftover <- NULL
  reads_1$gene0 <- gsub(special, "-", reads_1$gene)
  reads_1[c("hla", "leftover")] <- str_split_fixed(reads_1$gene, special, 2)
  reads_1$leftover <- NULL
  # get the HLA genes in a list
  cts_abc <- unique(reads_2$hla)[order(unique(reads_2$hla))]
  if (use_alt_align_ABC) {
    if (length(cts_abc[which(cts_abc %in% c("A", "B", "C"))]) != 3) {
      stop("the secondary molecule info count file does not contain alleles belonging to all of HLA-A, -B, and -C")
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
  return(pl2)
}



#' Getting raw scrHLAtag counts and analyzing distribution of HLA alleles per Cell Barcodes in UMAP space
#' 
#' @param reads  is the scrHLAtag count file including columns for CB, UMI, and HLA alleles (https://github.com/furlan-lab/scrHLAtag).
#' @param k  is the number of clusters to partition the UMAP space into, e.g. the number of entities or genotypes you 'think' there might be in your captured sample.
#' @param seu  is the Seurat object associated with the scrHLAtag count file (https://satijalab.org/seurat/index.html).
#' @param CB_rev_com  is a logical, called TRUE if the need to obtained the reverse complement of Cell Barcodes (CBs) is desired; default is FALSE. 
#' @param geno_metadata_id  is a character, the column ID of the Seurat metadata designated to distinguish genotypes, if this information is available. 'NULL' by default or when genotyping information is not available. 
#' @param hla_with_counts_above  is the number of total reads accross CBs at or above which an HLA allele is retained in the matrix.
#' @param CBs_with_counts_above  is the number of total reads accross HLA alleles at or above which a CB is retained in the matrix. Note: 'princomp()' can only be used with at least as many units (CBs) as variables (HLAs), thus the function will make sure that number of CBs is equal or more than available HLA alleles in the matrix.
#' @param match_CB_with_seu  is a logical, called TRUE if filtering CBs in the scrHLAtag count file with matching ones in the Seurat object is desired. 
#' @param umap_first_n_PCs  is a number representing the first 'n' principal components to input into the umap function. Note: current version only supports default parameters of the 'umap()' function (more info on: https://cran.r-project.org/web/packages/uwot/index.html)
#' @param QC_mm2  is a logical, called TRUE if removing low quality reads based on minimap2 tags is desired.
#' @param s1_belowmax  is a proportion (0 to 1) of the maximum value (best quality) of the minimap2 's1' tag above which the quality of the read is acceptable; default at 0.75 of the max s1 score.
#' @param AS_belowmax  is a proportion (0 to 1) of the maximum value (best quality) of the minimap2 'AS' tag above which the quality of the read is acceptable; default at 0.85 of the max AS score.
#' @param NM_thresh  is the number of mismatches and gaps in the minimap2 alignment at or below which the quality of the read is acceptable; default is 15.
#' @param de_thresh  is the gap-compressed per-base sequence divergence in the minimap2 alignment at or below which the quality of the read is acceptable; the number is between 0 and 1, and default is 0.015.
#' @param parallelize  is a logical, called TRUE if using parallel processing (multi-threading) is desired; default is TRUE.
#' @param pt_size  is a number, the size of the geometric point displayed by ggplot2. 
#' @param ...  arguments passed onto uwot::umap()
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
#' HLA_clusters(reads = cts[["mRNA"]], k = 2, seu = your_Seurat_Obj, geno_metadata_id = "geno", hla_with_counts_above = 1, CBs_with_counts_above = 40)
#' @export

HLA_clusters <- function(reads, k = 2, seu = NULL, CB_rev_com = FALSE, geno_metadata_id = NULL, hla_with_counts_above = 0, CBs_with_counts_above = 0, match_CB_with_seu = TRUE, umap_first_n_PCs = 25, QC_mm2 = TRUE, s1_belowmax = 0.75, AS_belowmax = 0.85, NM_thresh = 15, de_thresh = 0.015, parallelize = TRUE, pt_size = 0.5, ...) {
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
  alleles <- unique(reads$gene)
  reads$seu_barcode <- paste0(reads$samp,"_",reads$CB,"-1")
  reads <- with(reads, split(reads, list(seu_barcode=seu_barcode)))
  ## Matrix formation
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
  message(cat("\nRunning PCA (from 'stats') and UMAP (from 'uwot')"))
  ## normalize by size factor
  part_HLA <- prop.table(part_HLA, margin = 2)
  ## log-transform
  part_HLA <- log10(part_HLA+0.01)
  ## princomp from stats
  pc <- stats::princomp(t(part_HLA))
  pcv<-as.data.frame(pc$scores)
  elbow <- barplot((pc$sdev^2/sum(pc$sdev^2))[1:100])
  # pcv<-cbind(pcv, seu@meta.data[match(row.names(pcv), colnames(seu)),])
  # ggplot(pcv, aes(x=Comp.1, y=Comp.2, color=logUMI))+geom_point(size=0.25)+scale_color_viridis_b()+theme_bw()
  # ggplot(pcv, aes(x=Comp.1, y=Comp.2, color=geno))+geom_point(size=0.5)+scale_color_manual(values=pals::glasbey())+theme_bw()
  umat<-pcv[,1:umap_first_n_PCs] %>% as.matrix()
  umapout<-uwot::umap(umat, verbose=T, ...)
  colnames(umapout)<-c("umap1", "umap2")
  umapout <- as.data.frame(umapout)
  if (!is.null(seu)){
    if (!is.null(geno_metadata_id)) {
      umapout<-cbind(umapout, seu@meta.data[match(rownames(umapout), colnames(seu)),])
      g0 <- ggplot(umapout, aes(x=umap1, y=umap2, color=!!sym(geno_metadata_id)))+geom_point(size=pt_size)#+scale_color_manual(values=pals::glasbey())+theme_bw()
    }
  }
  # ggplot(umapout, aes(x=umap1, y=umap2, color=celltype))+geom_point(size=0.25)+scale_color_manual(values=pals::glasbey())+theme_bw()
  # ggplot(umapout, aes(x=umap1, y=umap2, color=geno))+geom_point(size=0.5)+scale_color_manual(values=pals::glasbey())+theme_bw()
  # ggplot(umapout, aes(x=umap1, y=umap2, color=logUMI))+geom_point(size=0.25)+scale_color_viridis_b()+theme_bw()
  message(cat("\nClustering on the UMAP space using Hierarchical Clustering (from 'stats')"))
  humapout <- stats::hclust(dist(as.matrix(umapout[,1:2]))) 
  umapout$hla_clusters <- stats::cutree(humapout, k = k)
  umapout$hla_clusters <- as.factor(umapout$hla_clusters)
  g <- ggplot(umapout, aes(x=umap1, y=umap2, color=hla_clusters))+geom_point(size=pt_size)#+scale_color_manual(values=pals::glasbey())+theme_bw()
  message(cat("\nDone!!"))
  if (!is.null(seu)){
    return(list(UMAP_coordinates = umapout, HLA_clusters_on_umap = g, genotype_on_umap = g0))
  } else {
    return(list(UMAP_coordinates = umapout, HLA_clusters_on_umap = g))
  }
}




#' Mapping the HLA clusters found by the 'HLA_clusters()' function back into the scrHLAtag count files
#'
#' @param reads.list  is a list of scrHLAtag count file(s) including columns for CB, UMI, and HLA alleles (https://github.com/furlan-lab/scrHLAtag).
#' @param cluster_coordinates  is the UMAP coordinates dataframe with HLA clustering information (found by the 'HLA_clusters()' function) from which clusters are extracted and mapped into the scrHLAtag count files by matching Cell Barcodes. Currently the barcode format supported is: SAMPLE_AATGCTTGGTCCATTA-1
#' @param CB_rev_com  is a logical, called TRUE if the need to obtained the reverse complement of Cell Barcodes (CBs) is desired; default is FALSE.
#' @param parallelize  is a logical, called TRUE if using parallel processing (multi-threading) is desired; default is TRUE.
#' @import stringr
#' @return a large list containing scrHLAtag count file(s) including columns for CB, UMI, and HLA alleles, with the addition of HLA_clusters
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
#' map_HLA_clusters(reads.list = cts, k = 2, cluster_coordinates = UMAP_dataframe_from_HLA_clusters_function)
#' @export

map_HLA_clusters <- function(reads.list, cluster_coordinates, CB_rev_com = FALSE, parallelize = TRUE) {
  ## parallelize
  if (parallelize) {
    multi_thread <- parallel::detectCores()
    # message(cat("\nMulti-threading! Available cores: ", parallel::detectCores()))
  } else {
    multi_thread <- 1
  }  
  if (CB_rev_com) {
    message(cat("\nConverting Cell Barcodes to their reverse complements"))
    for (i in 1:length(reads.list)) {
      reads.list[[i]]$hla_clusters <- NA
      reads.list[[i]]$hla_clusters <- cluster_coordinates$hla_clusters[
        match(paste0(reads.list[[i]]$samp,"_",pbmclapply(reads.list[[i]]$CB, function(x) intToUtf8(rev(utf8ToInt(chartr('ATGC', 'TACG', x)))), mc.cores = multi_thread) %>% unlist(),"-1"), 
              rownames(cluster_coordinates))
      ]
    }
  } else {
    for (i in 1:length(reads.list)) {
      reads.list[[i]]$hla_clusters <- NA
      reads.list[[i]]$hla_clusters <- cluster_coordinates$hla_clusters[
        match(paste0(reads.list[[i]]$samp,"_",reads.list[[i]]$CB,"-1"), 
              rownames(cluster_coordinates))
      ]
    }
  }
  return(reads.list)
}

