#' Deduping and Correcting scrHLAtag counts and Creating Seurat-compatible Matrices
#' 
#' @param cts  is the scrHLAtag count file including columns for CB, UMI, and HLA alleles (https://github.com/furlan-lab/scrHLAtag).
#' @param seu  is the Seurat object associated with the scrHLAtag count file (https://satijalab.org/seurat/index.html).
#' @param hla_recip  is a character list of recipient-specific HLA alleles if known; default is an empty character vector.
#' @param hla_donor  is a character list of donor-specific HLA alleles if known; default is an empty character vector.
#' @param s1_belowmax  is a proportion (0 to 1) of the maximum value (best quality) of the minimap2 's1' tag above which the quality of the read is acceptable; default at 0.75 of the max s1 score.
#' @param AS_belowmax  is a proportion (0 to 1) of the maximum value (best quality) of the minimap2 'AS' tag above which the quality of the read is acceptable; default at 0.85 of the max AS score.
#' @param NM_thresh  is the number of mismatches and gaps in the minimap2 alignment at or below which the quality of the read is acceptable; default is 15.
#' @param de_thresh  is the gap-compressed per-base sequence divergence in the minimap2 alignment at or below which the quality of the read is acceptable; the number is between 0 and 1, and default is 0.015.
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
#' HLA_Matrix(cts = cts[["mRNA"]], seu = your_Seurat_obj, hla_recip = c("A*24:02:01", "DRB1*04:01:01", "DRB4*01:03:02"), hla_donor = c("A*33:03:01", "B*42:01:01"))
#' @export


HLA_Matrix <- function(cts, seu, hla_recip = character(), hla_donor = character(), s1_belowmax = 0.75, AS_belowmax = 0.85, NM_thresh = 15, de_thresh = 0.015, Ct = 0) {
  ## check Seurat object
  if (class(seu) != "Seurat") {
    stop("Single-cell dataset container must be of class 'Seurat'")
  }
  ## add the dash into the HLA allele as it is the one accepted in the names of the Seurat assay features
  special <- "[_*|?.+$^]"
  if (any(grepl(special, cts$gene))) {
    cts$gene0 <- gsub(special, "-", cts$gene)
    cts[c("hla", "leftover")] <- str_split_fixed(cts$gene, special, 2)
    cts$leftover <- NULL
    print("Available reads per gene", quote=F)
    print(table(cts$hla, useNA = "ifany"))
  } else if (all(grepl("-", cts$gene))){
    cts$gene0 <- cts$gene
    cts[c("hla", "leftover")] <- str_split_fixed(cts$gene, "-", 2)
    cts$leftover <- NULL
    print("Available reads per gene", quote=F)
    print(table(cts$hla, useNA = "ifany"))
  } else {
    stop("The HLA allele column is unrecognizable or has incorrect format. \nMake sure gene and allele are separated by standard nomenclature asterisk (or other special character).")
  }
  cts$cbumi <- paste0(cts$CB, ":", cts$UMI)
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
  if (!all(c("CB", "UMI", "gene", "samp") %in% colnames(cts))){
    stop("scrHLAtag output 'cts' dataframe must at least contain the columns 'CB', 'UMI', 'gene' (with HLA alleles), and 'samp' (matching the sample names in the corresponding Seurat object)")
  }
  ## Remove low quality reads based on minimap2 tags
  print("1/6 - Removing low quality reads based on minimap2 tags", quote=F)
  cts.split <- with(cts, split(cts, list(gene0=gene0)))
  cts.split <- pbmclapply(cts.split, function(df){df[df$s1 > s1_belowmax*max(df$s1) & df$AS > AS_belowmax*max(df$AS) & df$NM <= NM_thresh & df$de <= de_thresh,]}, mc.cores = parallel::detectCores())
  cts <-  do.call("rbind", cts.split)
  row.names(cts)<-NULL
  rm(cts.split)

  ## see if more than 1 allele are present per umi at a time
  # count all the problematic CB:UMIs for which a molecular swap is suspected
  cts$mol_swap <- NA
  cts$mol_swap <- as.factor(cts$mol_swap)
  cts$class_swap <- NA
  cts$class_swap <- as.factor(cts$class_swap)
  # split,   this is computationally heavy (about 10min for 10M rows)
  cts.split <- with(cts, split(cts, list(cbumi=cbumi))) 
  # print("1/6 - Estimating UMI duplication rate", quote=F)
  # n_umi<-pbmclapply(cts.split, nrow, mc.cores = parallel::detectCores()) %>% unlist()
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
  print("2/6 - Estimating molecular swap", quote=F)
  pb <- txtProgressBar(min = 0, max = length(cts.split), style = 3, char = "=")
  for(j in 1:length(cts.split)){
    cts.split[[j]]$mol_swap <- ifelse(length(unique(cts.split[[j]]$gene)) > 1, 
                                      cts.split[[j]]$mol_swap <- "yes",
                                      cts.split[[j]]$mol_swap <- "no")
    setTxtProgressBar(pb, j)
  }
  close(pb)
  unique(cts$hla)
  cl1<-c("B" , "A" , "C")
  cl2<-c("DPA1", "DPB1", "DRB4", "DRB1", "DQA1", "DRB3", "DQB1")
  any(cts$hla %in% cl1)
  for(j in 1:length(cts.split)){
    cts.split[[j]]$class_swap <- ifelse(any(cts.split[[j]]$hla %in% cl1) & any(cts.split[[j]]$hla %in% cl2), 
                                        cts.split[[j]]$class_swap <- "yes",
                                        cts.split[[j]]$class_swap <- "no")
    setTxtProgressBar(pb, j)
  }
  close(pb)
  # count the molecular swap rate 
  mol_swap_rate_per_read <- length(which(sapply(cts.split, function(df) "yes" %in% df$mol_swap))) / length(which(sapply(cts.split, function(df) nrow(df)>1)))
  mol_swap_rate_per_read
  class_swap_rate_per_read<-length(which(sapply(cts.split, function(df) "yes" %in% df$class_swap))) / length(which(sapply(cts.split, function(df) nrow(df)>1)))
  class_swap_rate_per_read
  intraclass_swap_rate<- mol_swap_rate_per_read - class_swap_rate_per_read
  intraclass_swap_rate
  print(cat("\nintra-class molecular swap rate per UMI ", 
            format(round(100*intraclass_swap_rate, 3), nsmall = 1),
            "% of Cells\ninter-class molecular swap rate per UMI ",
            format(round(100*class_swap_rate_per_read, 3), nsmall = 1),
            "% of Cells\n"), quote = F)
  alleles <- unique(cts$gene0)
  
  ## Count Threshold 'Ct' above which the count of UMI copies becomes acceptable  
  if (Ct == 0) {
    cts.split.ct <- cts.split
    rm(cts.split)
  } else {
    cts.split.ct <- cts.split[sapply(cts.split, nrow) > Ct]
    rm(cts.split)
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
  print("3/6 - Correcting Molecular Swap: keeping the reads per UMI with the HLA allele occuring the most", quote = F)
  cts.split.ct <- pbmclapply(cts.split.ct, keep_one, mc.cores = parallel::detectCores())
  
  ## Performing Dedup
  print("4/6 - Performing Dedup on UMIs: removing PCR duplicates", quote = F)
  cts.fltr.dedup <- pbmclapply(cts.split.ct, function(df){df[1,]}, mc.cores = parallel::detectCores())
  cts.fltr.dedup <-  do.call("rbind", cts.fltr.dedup)
  row.names(cts.fltr.dedup)<-NULL
  rm(cts.split.ct)
  
  ## Clean-up HLA conflicts per CB
  # remove obsolete cols and add Seurat barcode
  cts.fltr.dedup$mol_swap <- NULL
  cts.fltr.dedup$class_swap <- NULL
  cts.fltr.dedup$seu_barcode <- paste0(cts.fltr.dedup$samp,"_",cts.fltr.dedup$CB,"-1")
  cts.fltr.dedup$hla_conflict <- NA
  cts.fltr.dedup$hla_conflict <- as.factor(cts.fltr.dedup$hla_conflict)
  # split by Seurat barcode
  cts.dedup.cb <- with(cts.fltr.dedup, split(cts.fltr.dedup, list(seu_barcode=seu_barcode)))
  rm(cts.fltr.dedup)
  # detect HLA conflicts (i.e. donor-spec and recipient-spec HLA in the same barcode)
  pb <- txtProgressBar(min = 0, max = length(cts.dedup.cb), style = 3, char = "=")
  for(j in 1:length(cts.dedup.cb)){
    cts.dedup.cb[[j]]$hla_conflict <- ifelse(any(cts.dedup.cb[[j]]$gene0 %in% hla_recip) & any(cts.dedup.cb[[j]]$gene0 %in% hla_donor),
                                             cts.dedup.cb[[j]]$hla_conflict <- "yes",
                                             cts.dedup.cb[[j]]$hla_conflict <- "no")
    setTxtProgressBar(pb, j)
  }
  close(pb)
  hla_conflict_rate <- length(which(sapply(cts.dedup.cb, function(df) "yes" %in% df$hla_conflict))) / length(cts.dedup.cb)
  print(cat("Conflicting HLA (both donor and recipient HLA within the same Cell Barcode) affects ", 
            format(round(100*hla_conflict_rate, 1), nsmall = 1),
            "% of Cells\n"), quote = F)
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
  print("5/6 - Conflict Correction: assuming a cell cannot have both recipient and donor-origin HLA allele, keeping only the most occuring", quote = F)
  cts.dedup.cb <- pbmclapply(cts.dedup.cb, keep_two, recip = hla_recip, donor = hla_donor, mc.cores = parallel::detectCores())
  
  ## Matrix formation
  # matrix
  HLA.matrix <- matrix(0, nrow = length(alleles), ncol = length(cts.dedup.cb), dimnames = list(alleles, names(cts.dedup.cb)))
  print("6/6 - Creating the HLA Count Matrix compatible with the Seurat object", quote = F)
  pb <- txtProgressBar(min = 0, max = length(cts.dedup.cb), style = 3, char = "=")
  for (i in 1:length(cts.dedup.cb)) {
    counts <- table(cts.dedup.cb[[i]]$gene0)
    HLA.matrix[, i] <- counts[alleles]
    setTxtProgressBar(pb, i)
  }
  close(pb)
  rm(cts.dedup.cb)
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
  return(HLA)
}


#' Extracting the top HLA alleles from the scrHLAtag count files
#' 
#' @param cts_1  is the primary scrHLAtag count file (1 of 2 files containing either the mRNA molecular info or the genomic (gene) molecular info). It includes columns for CB, UMI, and HLA alleles (https://github.com/furlan-lab/scrHLAtag).
#' @param cts_2  is the secondary scrHLAtag count file (the alternative file vs. the one designated in 'cts_1' argument). It includes columns for CB, UMI, and HLA alleles (https://github.com/furlan-lab/scrHLAtag). Default is NULL, in which case it will not be able to count alternative aligment and argument 'use_alt_align_ABC' will become irrelevant.
#' @param frac  is the fraction (0 to 1) of total reads for a particular HLA gene, which incorporates the highest ranking alleles of that gene in terms of number of reads; default at 0.65 .
#' @param min_alleles_keep  is a numeric representing the minimum number of highest ranking alleles to keep despite filtering by fraction 'frac'; default is 5.
#' @param min_reads_per_gene  is a numeric representing minimum number of total reads per HLA gene (including all its alleles) below which the gene is filtered out; default is 200. 
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
#' Top_HLA_list(cts_1 = cts[["mRNA"]], cts_2 = cts[["gene"]], frac = 0.8, min_alleles_keep = 2, use_alt_align_ABC = TRUE)
#' @export

Top_HLA_list <- function(cts_1, cts_2 = NULL, frac = 0.65, min_alleles_keep = 5, min_reads_per_gene = 200, insert_pop_most_freq = TRUE, use_alt_align_ABC = FALSE){
  if (is.null(cts_2)) {
    cts_2 <- cts_1
    warning("The molecule_info_gene.txt.gz count file does not seem to be included. The function will run but the argument 'use_alt_align_ABC' will be irrelevant.")
  }
  # extract the HLA genes that appear in the reads
  special <- "[_*|?.+$^]"
  cts_2$gene0 <- gsub(special, "-", cts_2$gene)
  cts_2[c("hla", "leftover")] <- str_split_fixed(cts_2$gene, special, 2)
  cts_2$leftover <- NULL
  cts_1$gene0 <- gsub(special, "-", cts_1$gene)
  cts_1[c("hla", "leftover")] <- str_split_fixed(cts_1$gene, special, 2)
  cts_1$leftover <- NULL
  # get the HLA genes in a list
  cts_abc <- unique(cts_2$hla)[order(unique(cts_2$hla))]
  if (use_alt_align_ABC) {
    if (length(cts_abc[which(cts_abc %in% c("A", "B", "C"))]) != 3) {
      stop("the molecule_info_gene.txt.gz file does not contain alleles belonging to all of HLA-A, -B, and -C")
    }
  }
  cts_notabc <- unique(cts_1$hla)[order(unique(cts_1$hla))]
  top_alleles <- c()
  for (j in 1:length(cts_notabc)){
    if (!use_alt_align_ABC) {
      t<-as.data.table(table(cts_1[cts_1$hla == cts_notabc[j],]$gene))
    } else {
      if (cts_notabc[j] %in% c("A", "B", "C")) {
        t<-as.data.table(table(cts_2[cts_2$hla == cts_notabc[j],]$gene))
      } else {
        t<-as.data.table(table(cts_1[cts_1$hla == cts_notabc[j],]$gene))
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
    if (length(unique(t$twofield)) <= min_alleles_keep) {
      top <- t$twofield
    } else {
      top <- t$twofield[t$csum/sum(t$N) > (1-frac)]
    }
    # include the most common allele in the population even if not in the top alleles by fraction of all reads (from previous step)
    if (insert_pop_most_freq) {
      u <- unique(t$onefield[t$csum/sum(t$N) > (1-frac)])
      spl <- with(t, split(t, list(onefield=onefield)))
      spl <- spl[names(spl) %in% u]
      if (length(spl) > 0) {
        for (i in 1:length(spl)) {
          top <- c(top, spl[[i]]$twofield[spl[[i]]$fscore == max(spl[[i]]$fscore[!is.na(spl[[i]]$fscore)])])
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
#' @param cts_1  is the primary scrHLAtag count file (1 of 2 files containing either the mRNA molecular info or the genomic (gene) molecular info). It includes columns for CB, UMI, and HLA alleles (https://github.com/furlan-lab/scrHLAtag).
#' @param cts_2  is the secondary scrHLAtag count file (the alternative file vs. the one designated in 'cts_1' argument). It includes columns for CB, UMI, and HLA alleles (https://github.com/furlan-lab/scrHLAtag). Default is NULL, in which case it will not be able to count alternative aligment and argument 'use_alt_align_ABC' will become irrelevant.
#' @param top_hla  is a numeric, representing the number of top HLA alleles (i.e. with the highest number of reads) per HLA gene to display in the plot; default is 10.
#' @param min_reads_per_gene  is a numeric representing minimum number of total reads per HLA gene (including all its alleles) below which the gene is filtered out; default is 200. 
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
#' Top_HLA_plot(cts_1 = cts[["mRNA"]], cts_2 = cts[["gene"]], use_alt_align_ABC = TRUE)
#' @export

Top_HLA_plot <- function(cts_1, cts_2 = NULL, top_hla = 10, min_reads_per_gene = 200, use_alt_align_ABC = FALSE, color_pal = NULL){
  if (is.null(cts_2)) {
    cts_2 <- cts_1
    warning("The molecule_info_gene.txt.gz count file does not seem to be included. The function will run but the argument 'use_alt_align_ABC' will be irrelevant.")
  }  
  # extract the HLA genes that appear in the reads
  special <- "[_*|?.+$^]"
  cts_2$gene0 <- gsub(special, "-", cts_2$gene)
  cts_2[c("hla", "leftover")] <- str_split_fixed(cts_2$gene, special, 2)
  cts_2$leftover <- NULL
  cts_1$gene0 <- gsub(special, "-", cts_1$gene)
  cts_1[c("hla", "leftover")] <- str_split_fixed(cts_1$gene, special, 2)
  cts_1$leftover <- NULL
  # get the HLA genes in a list
  cts_abc <- unique(cts_2$hla)[order(unique(cts_2$hla))]
  if (use_alt_align_ABC) {
    if (length(cts_abc[which(cts_abc %in% c("A", "B", "C"))]) != 3) {
      stop("the molecule_info_gene.txt.gz file does not contain alleles belonging to all of HLA-A, -B, and -C")
    }
  }
  cts_notabc <- unique(cts_1$hla)[order(unique(cts_1$hla))]
  tab <- data.frame("twofield"=character(0), "N"=integer(0), "fscore"=numeric(0), "hlagene"=character(0), "onefield"=character(0), "csum"=integer(0))
  for (j in 1:length(cts_notabc)){
    if (!use_alt_align_ABC) {
      t<-as.data.table(table(cts_1[cts_1$hla == cts_notabc[j],]$gene))
    } else {
      if (cts_notabc[j] %in% c("A", "B", "C")) {
        t<-as.data.table(table(cts_2[cts_2$hla == cts_notabc[j],]$gene))
      } else {
        t<-as.data.table(table(cts_1[cts_1$hla == cts_notabc[j],]$gene))
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
      scale_fill_gradientn(colours = rev(color_pal), na.value = "grey70")+
      xlab("Reads (n)")+ 
      ylab("Top 10 alleles")+
      theme(text = element_text(size = 9),legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    plots<- c(plots,list(g))
  }
  pl <- do.call("plot_grid", c(plots, align = "hv", ncol=floor(sqrt(length(plots)))))
  return(pl)
}
