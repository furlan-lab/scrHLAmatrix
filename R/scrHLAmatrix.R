#' Deduping and Correcting scrHLAtag counts and Creating Seurat-compatible Matrices
#' 
#' @param reads  is the scrHLAtag count file including columns for CB, UMI, and HLA alleles (https://github.com/furlan-lab/scrHLAtag).
#' @param seu  is the Seurat object associated with the scrHLAtag count file (https://satijalab.org/seurat/index.html).
#' @param hla_recip  is a character list of recipient-specific HLA alleles if known; default is an empty character vector.
#' @param hla_donor  is a character list of donor-specific HLA alleles if known; default is an empty character vector.
#' @param QC_mm2  is a logical, called TRUE if removing low quality reads based on minimap2 tags is desired.
#' @param res_conflict_per_gene  is a logical, called TRUE if resolving per-HLA genotype conflicts is desired, with the assumption that each Cell can have no more than 2 alleles of the same HLA gene.
#' @param LD_correct  is a logical, called TRUE if Linkage Diseqilibrium (LD) correction in the HLA-DR locus is desired, with the assumption of a very strong LD between certain DRB1 allele families and the DRB2, DRB3, DRB4, DRB5, DRB6, DRB7, DRB8, and DRB9 loci.
#' @param remove_alleles  is a character list of HLA alleles to be removed from the count file if desired; default is an empty character vector.
#' @param s1_belowmax  is a proportion (0 to 1) of the maximum value (best quality) of the minimap2 's1' tag above which the quality of the read is acceptable; default at 0.75 of the max s1 score.
#' @param AS_belowmax  is a proportion (0 to 1) of the maximum value (best quality) of the minimap2 'AS' tag above which the quality of the read is acceptable; default at 0.85 of the max AS score.
#' @param NM_thresh  is the number of mismatches and gaps in the minimap2 alignment at or below which the quality of the read is acceptable; default is 15.
#' @param de_thresh  is the gap-compressed per-base sequence divergence in the minimap2 alignment at or below which the quality of the read is acceptable; the number is between 0 and 1, and default is 0.015.
#' @param parallelize  is a logical, called TRUE if using parallel processing (multi-threading) is desired; default is FALSE.
#' @param CB_rev_com  is a logical, called TRUE if the need to obtained the reverse complement of Cell Barcodes (CBs) is desired; default is FALSE. 
#' @param return_stats  is a logical, when TRUE returns step-by-step read statistics and UMI duplication rate in a list of dataframes and plot, in addition to the Seurat-compatible count matrix; will require additional computations which may noticeably slow down the function; defualt is FALSE.
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

HLA_Matrix <- function(reads, seu, hla_recip = character(), hla_donor = character(), QC_mm2 = TRUE, res_conflict_per_gene = TRUE, LD_correct = TRUE, remove_alleles = character(), s1_belowmax = 0.80, AS_belowmax = 0.80, NM_thresh = 15, de_thresh = 0.01, parallelize = FALSE, CB_rev_com = FALSE, return_stats = FALSE) {
  s <- Sys.time()
  #message(cat(format(s, "%F %H:%M:%S")))
  n_reads <- nrow(reads)
  if (n_reads > 1e+06) {
    message(cat("\nLarge scrHLAtag count file detected (", n_reads, " rows counted); expect resource-demanding processing and long run times", sep = ""))
    message(cat(crayon::green("Note: "), "Multi-threading errors had been more frequently experienced while running large count files; it is recommended to maintain 'parallelize = FALSE'", sep = ""))
  }
  ## check Seurat object
  if (class(seu) == "Seurat") {
    message(cat("\nObject of class 'Seurat' detected"))
    message(cat(crayon::green("Note: "), "Currently the Seurat Barcode (i.e. colnames or Cells) supported format is: SAMPLE_AATGCTTGGTCCATTA-1", sep = ""))
  } else {
    stop("Single-cell dataset container must be of class 'Seurat'")
  }
  ## parallelize
  if (parallelize) {
    multi_thread <- parallel::detectCores()
    message(cat("\nMulti-threading! Available cores: ", parallel::detectCores(), sep = ""))
  } else {
    multi_thread <- 1
  }
  ## check relevant col names exist
  if (!all(c("CB", "UMI", "gene", "samp") %in% colnames(reads))){
    stop("scrHLAtag output 'reads' dataframe must at least contain the columns 'CB', 'UMI', 'gene' (with HLA alleles), and 'samp' (matching the sample names in the corresponding Seurat object)")
  }
  ## jettison unused columns
  reads <- reads[, colnames(reads) %in% c("CB", "UMI", "gene", "NM", "AS", "s1", "de", "samp"), drop = F]
  ## add the dash into the HLA allele as it is the one accepted in the names of the Seurat assay features
  special <- "[_*|?.+$^]"
  if (any(grepl(special, reads$gene))) {
    reads$gene0 <- gsub(special, "-", reads$gene)
    reads <- reads %>% mutate(gene1 = strsplit(gene0, "-") %>% sapply(., function(x) c(x[1], substr(x[2], 1, 2))) %>% apply(., 2, function(x) paste(x, collapse = "-")))
    reads[c("hla", "leftover")] <- stringr::str_split_fixed(reads$gene, special, 2)
    reads$leftover <- NULL
    reads$gene <- NULL
  } else if (all(grepl("-", reads$gene))){
    reads$gene0 <- reads$gene
    reads <- reads %>% mutate(gene1 = strsplit(gene0, "-") %>% sapply(., function(x) c(x[1], substr(x[2], 1, 2))) %>% apply(., 2, function(x) paste(x, collapse = "-")))
    reads[c("hla", "leftover")] <- stringr::str_split_fixed(reads$gene, "-", 2)
    reads$leftover <- NULL
    reads$gene <- NULL
  } else {
    stop("The HLA allele column is unrecognizable or has incorrect format. \nMake sure gene and allele are separated by standard nomenclature asterisk (or other special character)")
  }
  reads$cbumi <- stringr::str_c(reads$samp, ":", reads$CB, ":", reads$UMI)
  reads$UMI <- NULL # no longer needed 
  ## check format of 'hla_recip' and 'hla_donor'
  if (any(grepl(special, hla_recip))) {
    hla_recip <- gsub(special, "-", hla_recip)
  } else if (all(grepl("-", hla_recip))) {
    hla_recip <- hla_recip
  } else {
    stop("Incorrect format for recipient-defined and/or donor-defined HLA alleles. \nMake sure gene and allele are separated by standard nomenclature asterisk (or other special character)")
  }
  if (any(grepl(special, hla_donor))) {
    hla_donor <- gsub(special, "-", hla_donor)
  } else if (all(grepl("-", hla_donor))) {
    hla_donor <- hla_donor
  } else {
    stop("Incorrect format for recipient-defined and/or donor-defined HLA alleles. \nMake sure gene and allele are separated by standard nomenclature asterisk (or other special character)")
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
  ## Reverse Complement the CB
  if (CB_rev_com) {
    message(cat("\nConverting Cell Barcodes to their reverse complements"))
    # reads$CB <- pbmcapply::pbmclapply(reads$CB, function(x) as.character(Biostrings::reverseComplement(DNAString(x))), mc.cores = multi_thread) %>% unlist() # slow
    reads$CB <- pbmcapply::pbmclapply(reads$CB, function(x) intToUtf8(rev(utf8ToInt(chartr('ATGC', 'TACG', x)))), mc.cores = multi_thread) %>% unlist()        # fast
  } 
  ## Estimating number of reads and number of CBs
  reads$seu_barcode <- stringr::str_c(reads$samp,"_",reads$CB,"-1")
  reads$samp <- NULL # no longer needed
  stats_df <- data.frame(reads = numeric(), reads_found = numeric(), cbs = numeric(), cbs_found = numeric(), cbs_seu = numeric(), cb_seu_match_rate = numeric(), step = character())
  reads_in_seu <- as.numeric(reads$seu_barcode %in% colnames(seu) %>% table())[2]
  ureads_in_seu <- as.numeric(unique(reads$seu_barcode) %in% colnames(seu) %>% table())[2]
  l_ureads_cb <- length(unique(reads$seu_barcode))
  n_reads <- nrow(reads)
  n_cells <- length(colnames(seu))
  message(cat("\nReads in count file: ", n_reads, 
              ", including ", reads_in_seu,
              " (", format(round(100*(reads_in_seu/n_reads), 2), nsmall = 1),
              "%) belonging to Cells found in Seurat object", sep = ""))
  message(cat("Unique Cell Barcodes (CB): ", l_ureads_cb, 
              ", including ", 
              ureads_in_seu, 
              " (", format(round(100*(ureads_in_seu/l_ureads_cb), 2), nsmall = 1), 
              "%) matching the ", n_cells, " Cells in Seurat object (match rate ",
              format(round(100*(ureads_in_seu/n_cells), 2), nsmall = 1), "%)", sep = ""))
  if (return_stats) {
  stats_df <- rbind(stats_df,
                    data.frame(reads = n_reads, 
                               reads_found = reads_in_seu, 
                               cbs = l_ureads_cb, 
                               cbs_found = ureads_in_seu, 
                               cbs_seu = n_cells, 
                               cb_seu_match_rate = ureads_in_seu/n_cells, 
                               step = "0_raw_reads"))
  }
  message(cat("Available reads per gene:"))
  print(table(reads$hla, useNA = "ifany"))
  if (return_stats) {
    message(cat("\nEstimating UMI duplication rate"))
    reads <- split(data.table::setDT(reads), by = "cbumi")
    reads <- parallel::mclapply(reads, data.table::setDF, mc.cores = multi_thread)
    n_umi <- pbmcapply::pbmclapply(reads, nrow, mc.cores = multi_thread) %>% unlist()
    umi_counts<- data.frame(n_umi)
    umi_counts$dummy <- 1 #had to add this dummy var for the code to work, removed later; can also use drop = F in next line
    umi_counts <- umi_counts[order(umi_counts$n_umi),]
    row.names(umi_counts)<- NULL
    umi_counts$rank <- row.names(umi_counts) %>% as.numeric()
    umi_counts$dummy <- NULL
    g <- ggplot(umi_counts, aes(x= rank, y=n_umi))+
      #geom_smooth(size=2, method = "gam")+
      geom_line()+
      scale_y_log10(name = "PCR copies per UMI")+
      scale_x_continuous(name = "Rank (nth UMI)", n.breaks = 8) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    print(g)  
    reads <- data.table::rbindlist(reads)
    row.names(reads)<-NULL
  }
  ## Remove low quality reads based on minimap2 tags
  if (QC_mm2 & (!all(c("NM", "AS", "s1", "de") %in% colnames(reads))) ) {
    message(cat("\nCannot remove low quality reads based on minimap2 tags: scrHLAtag output 'reads' dataframe must contain the mm2 columns 'NM', 'AS', 's1', and 'de'; ", crayon::red("continuing without mm2 QC"), sep = ""))
    QC_mm2 <- FALSE
  }
  if (QC_mm2) {
    message(cat("\nRemoving low quality reads based on minimap2 tags"))   
    reads <- split(data.table::setDT(reads), by = "gene0")   
    reads <- pbmcapply::pbmclapply(reads, function(df){df[df$s1 > s1_belowmax*max(df$s1) & df$AS > AS_belowmax*max(df$AS) & df$NM <= NM_thresh & df$de <= de_thresh,]}, mc.cores = multi_thread)    
    reads <- data.table::rbindlist(reads) 
    row.names(reads)<-NULL
    if (return_stats) {
      reads_in_seu <- as.numeric(reads$seu_barcode %in% colnames(seu) %>% table())[2]
      ureads_in_seu <- as.numeric(unique(reads$seu_barcode) %in% colnames(seu) %>% table())[2]
      l_ureads_cb <- length(unique(reads$seu_barcode))
      n_reads <- nrow(reads)
      n_cells <- length(colnames(seu))
      message(cat("  Reads remaining: ", n_reads, 
                  ", including ", reads_in_seu,
                  " (", format(round(100*(reads_in_seu/n_reads), 2), nsmall = 1),
                  "%) belonging to Cells found in Seurat object", sep = ""))
      message(cat("  CBs remaining: ", l_ureads_cb, 
                  ", including ", 
                  ureads_in_seu, 
                  " (", format(round(100*(ureads_in_seu/l_ureads_cb), 2), nsmall = 1), 
                  "%) matching the ", n_cells, " Cells in Seurat object (match rate ",
                  format(round(100*(ureads_in_seu/n_cells), 2), nsmall = 1), "%)", sep = ""))
      stats_df <- rbind(stats_df,
                        data.frame(reads = n_reads, 
                                   reads_found = reads_in_seu, 
                                   cbs = l_ureads_cb, 
                                   cbs_found = ureads_in_seu, 
                                   cbs_seu = n_cells, 
                                   cb_seu_match_rate = ureads_in_seu/n_cells, 
                                   step = "1_qc_mm2"))
    }
  } else {
    reads <- reads[order(reads$gene0), ]
  }
  ## jettison unused mm2 columns
  reads <- reads[, !colnames(reads) %in% c("NM", "AS", "s1", "de")]
  ## see if more than 1 allele are present per umi at a time
  # count all the problematic CB:UMIs for which a molecular swap is suspected
  reads$mol_swap <- NA
  reads$mol_swap <- as.factor(reads$mol_swap)
  reads$class_swap <- NA
  reads$class_swap <- as.factor(reads$class_swap)  
  # split
  message(cat("\nEstimating Molecular Swap (excluding unduplicated UMIs where molecular swap cannot be estimated)"))   
  reads <- split(data.table::setDT(reads), by = "cbumi") 
  pb <- pbmcapply::progressBar(min = 0, max = length(reads), style = "ETA", char = "=")
  for(j in 1:length(reads)){
    reads[[j]] <- data.table::setDF(reads[[j]])
    reads[[j]]$mol_swap <- ifelse(length(unique(reads[[j]]$gene0)) > 1, 
                                      reads[[j]]$mol_swap <- "yes",
                                      reads[[j]]$mol_swap <- "no")
    setTxtProgressBar(pb, j)
  }
  close(pb)
  cl1<-c("A" , "B" , "C", "E" , "F" , "G", "H", "K", "Y", "J", "L", "N", "P", "S", "T", "U", "V", "W")
  cl2<-c("DRA", "DRB1", "DRB2", "DRB3", "DRB4", "DRB5", "DRB6", "DRB7", "DRB8", "DRB9", "DQA1", "DQA2", "DQB1", "DQB2", "DPA1", "DPA2", "DPB1", "DPB2", "DMA", "DMB", "DOA", "DOB")
  cl0<-c("HFE", "MICA", "MICB", "TAP1", "TAP2")
  for(j in 1:length(reads)){
    reads[[j]]$class_swap <- ifelse((any(reads[[j]]$hla %in% cl1) & any(reads[[j]]$hla %in% cl2)) | (any(reads[[j]]$hla %in% cl1) & any(reads[[j]]$hla %in% cl0)) | (any(reads[[j]]$hla %in% cl2) & any(reads[[j]]$hla %in% cl0)), 
                                        reads[[j]]$class_swap <- "yes",
                                        reads[[j]]$class_swap <- "no")
    setTxtProgressBar(pb, j)
  }
  close(pb)
  # count the molecular swap rate 
  mol_swap_rate_per_read <- length(which(sapply(reads, function(df) "yes" %in% df$mol_swap))) / length(which(sapply(reads, function(df) nrow(df)>1)))
  class_swap_rate_per_read<-length(which(sapply(reads, function(df) "yes" %in% df$class_swap))) / length(which(sapply(reads, function(df) nrow(df)>1)))
  intraclass_swap_rate<- mol_swap_rate_per_read - class_swap_rate_per_read
  message(cat("\n  intra-class molecular swap rate per UMI ", 
    format(round(100*intraclass_swap_rate, 3), nsmall = 1),
    "% of Cells\n  inter-class molecular swap rate per UMI ",
    format(round(100*class_swap_rate_per_read, 3), nsmall = 1),
    "% of Cells", sep = ""))  
  ## Function to keep the most occurring HLA when more than 1 HLA is present per cb:umi
  keep_one <- function(df) {
    n_hla <- table(df$gene0)
    if (length(n_hla) > 1) {
      most_hla <- names(which(n_hla==max(n_hla)))
      if (length(most_hla) > 1) {
        # If there is a tie, remove the entire df
        df <- NULL
      } else {
        # Filter the df to keep only rows with the most occurring hla
        df <- df[df$gene0 %in% most_hla, ]
      }
    }
    return(df)
  }
  # Applying the function
  message(cat("\nCorrecting Molecular Swap: keeping the reads per UMI with the HLA allele having the highest statistical probability"))
  reads <- pbmcapply::pbmclapply(reads, keep_one, mc.cores = multi_thread) 
  reads <- Filter(Negate(is.null), reads)
  if (return_stats) {
    reads <- data.table::rbindlist(reads)      
    row.names(reads)<-NULL
    reads_in_seu  <- as.numeric(reads$seu_barcode %in% colnames(seu) %>% table())[2]
    ureads_in_seu <- as.numeric(unique(reads$seu_barcode) %in% colnames(seu) %>% table())[2]
    l_ureads_cb  <- length(unique(reads$seu_barcode))
    n_reads <- nrow(reads)
    n_cells <- length(colnames(seu))
    message(cat("  Reads remaining: ", n_reads, 
                ", including ", reads_in_seu,
                " (", format(round(100*(reads_in_seu/n_reads), 2), nsmall = 1),
                "%) belonging to Cells found in Seurat object", sep = ""))
    message(cat("  CBs remaining: ", l_ureads_cb, 
                ", including ", 
                ureads_in_seu, 
                " (", format(round(100*(ureads_in_seu/l_ureads_cb), 2), nsmall = 1), 
                "%) matching the ", n_cells, " Cells in Seurat object (match rate ",
                format(round(100*(ureads_in_seu/n_cells), 2), nsmall = 1), "%)", sep = ""))
    stats_df <- rbind(stats_df,
                      data.frame(reads = n_reads, 
                                 reads_found = reads_in_seu, 
                                 cbs = l_ureads_cb, 
                                 cbs_found = ureads_in_seu, 
                                 cbs_seu = n_cells, 
                                 cb_seu_match_rate = ureads_in_seu/n_cells, 
                                 step = "2_mol_swap"))     
    reads <- split(data.table::setDT(reads), by = "cbumi") 
  }
  ## Performing Dedup
  message(cat("\nPerforming Dedup on UMIs: removing PCR duplicates"))
  reads <- pbmcapply::pbmclapply(reads, function(df){
    df <- data.table::setDF(df)
    df <- df[1,]
    return(df)
  }, mc.cores = multi_thread)
  reads <- data.table::rbindlist(reads)  
  row.names(reads)<-NULL
  if (return_stats) {
    reads_in_seu <- as.numeric(reads$seu_barcode %in% colnames(seu) %>% table())[2]
    ureads_in_seu <- as.numeric(unique(reads$seu_barcode) %in% colnames(seu) %>% table())[2]
    l_ureads_cb <- length(unique(reads$seu_barcode))
    n_reads <- nrow(reads)
    n_cells <- length(colnames(seu))
    message(cat("  Reads remaining after Dedup: ", n_reads, 
                ", including ", reads_in_seu,
                " (", format(round(100*(reads_in_seu/n_reads), 2), nsmall = 1),
                "%) belonging to Cells found in Seurat object", sep = ""))
    message(cat("  CBs remaining: ", l_ureads_cb, 
                ", including ", 
                ureads_in_seu, 
                " (", format(round(100*(ureads_in_seu/l_ureads_cb), 2), nsmall = 1), 
                "%) matching the ", n_cells, " Cells in Seurat object (match rate ",
                format(round(100*(ureads_in_seu/n_cells), 2), nsmall = 1), "%)", sep = ""))
    stats_df <- rbind(stats_df,
                      data.frame(reads = n_reads, 
                                 reads_found = reads_in_seu, 
                                 cbs = l_ureads_cb, 
                                 cbs_found = ureads_in_seu, 
                                 cbs_seu = n_cells, 
                                 cb_seu_match_rate = ureads_in_seu/n_cells, 
                                 step = "3_dedup"))
  }
  # remove obsolete cols and add Seurat barcode
  reads$mol_swap <- NULL
  reads$class_swap <- NULL
  reads$hla_conflict <- NA
  reads$hla_conflict <- as.factor(reads$hla_conflict)
  ## Remove undesired HLa alleles
  if (length(remove_alleles) > 0) {
    message(cat("\nRemoving the following alleles from the counts file:", remove_alleles, sep = " "))
    remove_alleles <- remove_alleles %>% gsub(special, "-", .)
    if (length(which(reads$gene0 %in% remove_alleles)) > 0) {
      reads <- reads[-which(reads$gene0 %in% remove_alleles),]
    }
    if (return_stats) {
      reads_in_seu <- as.numeric(reads$seu_barcode %in% colnames(seu) %>% table())[2]
      ureads_in_seu <- as.numeric(unique(reads$seu_barcode) %in% colnames(seu) %>% table())[2]
      l_ureads_cb <- length(unique(reads$seu_barcode))
      n_reads <- nrow(reads)
      n_cells <- length(colnames(seu))
      message(cat("  Reads remaining: ", n_reads, 
                  ", including ", reads_in_seu,
                  " (", format(round(100*(reads_in_seu/n_reads), 2), nsmall = 1),
                  "%) belonging to Cells found in Seurat object", sep = ""))
      message(cat("  CBs remaining: ", l_ureads_cb, 
                  ", including ", 
                  ureads_in_seu, 
                  " (", format(round(100*(ureads_in_seu/l_ureads_cb), 2), nsmall = 1), 
                  "%) matching the ", n_cells, " Cells in Seurat object (match rate ",
                  format(round(100*(ureads_in_seu/n_cells), 2), nsmall = 1), "%)", sep = ""))
      stats_df <- rbind(stats_df,
                        data.frame(reads = n_reads, 
                                   reads_found = reads_in_seu, 
                                   cbs = l_ureads_cb, 
                                   cbs_found = ureads_in_seu, 
                                   cbs_seu = n_cells, 
                                   cb_seu_match_rate = ureads_in_seu/n_cells, 
                                   step = "4_rmv_alleles"))
    }
  }
  ## Clean-up HLA conflicts per CB
  message(cat("\nDonor-v-Recipient Genotype Conflict: assuming a cell cannot have both recipient and donor-origin HLA allele"))
  # split by Seurat barcode
  reads <- split(data.table::setDT(reads), by = "seu_barcode")
  # detect HLA conflicts (i.e. donor-spec and recipient-spec HLA in the same barcode)
  pb <- pbmcapply::progressBar(min = 0, max = length(reads), style = "ETA", char = "=")
  for(j in 1:length(reads)){
    reads[[j]] <- data.table::setDF(reads[[j]])
    reads[[j]]$hla_conflict <- ifelse(any(reads[[j]]$gene0 %in% hla_recip) & any(reads[[j]]$gene0 %in% hla_donor),
                                             reads[[j]]$hla_conflict <- "yes",
                                             reads[[j]]$hla_conflict <- "no")
    setTxtProgressBar(pb, j)
  }
  close(pb)
  hla_conflict_rate <- length(which(sapply(reads, function(df) "yes" %in% df$hla_conflict))) / length(reads)
  message(cat("  Both donor and recipient HLA within the same Cell Barcode in ", 
            format(round(100*hla_conflict_rate, 3), nsmall = 1),
            "% of Cells", sep = ""))
  # function to clean-up HLA conflicts by Seurat barcode
  remove_conflict <- function(df, recip, donor) {
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
  if (hla_conflict_rate == 0) {
    message(cat("\nResolving Donor-v-Recipient Genotype Conflicts: no conflicts to resolve"))
    if (return_stats) {
      reads <- data.table::rbindlist(reads)
      row.names(reads)<-NULL
    }
  } else {
    message(cat("\nResolving Donor-v-Recipient Genotype Conflicts:\n  keeping either donor-specific or recipient specific HLA-associated UMIs, based on their count difference per Cell"))
    reads <- pbmcapply::pbmclapply(reads, remove_conflict, recip = hla_recip, donor = hla_donor, mc.cores = multi_thread)
    if (return_stats) {
      reads <- data.table::rbindlist(reads)
      row.names(reads)<-NULL
      reads_in_seu  <- as.numeric(reads$seu_barcode %in% colnames(seu) %>% table())[2]
      ureads_in_seu <- as.numeric(unique(reads$seu_barcode) %in% colnames(seu) %>% table())[2]
      l_ureads_cb  <- length(unique(reads$seu_barcode))
      n_reads <- nrow(reads)
      n_cells <- length(colnames(seu))
      message(cat("  Reads remaining: ", n_reads, 
                  ", including ", reads_in_seu,
                  " (", format(round(100*(reads_in_seu/n_reads), 2), nsmall = 1),
                  "%) belonging to Cells found in Seurat object", sep = ""))
      message(cat("  CBs remaining: ", l_ureads_cb, 
                  ", including ", 
                  ureads_in_seu, 
                  " (", format(round(100*(ureads_in_seu/l_ureads_cb), 2), nsmall = 1), 
                  "%) matching the ", n_cells, " Cells in Seurat object (match rate ",
                  format(round(100*(ureads_in_seu/n_cells), 2), nsmall = 1), "%)", sep = ""))
      stats_df <- rbind(stats_df,
                        data.frame(reads = n_reads, 
                                   reads_found = reads_in_seu, 
                                   cbs = l_ureads_cb, 
                                   cbs_found = ureads_in_seu, 
                                   cbs_seu = n_cells, 
                                   cb_seu_match_rate = ureads_in_seu/n_cells, 
                                   step = "5_dn_rp_conflict"))
    }
  }
  ## Resolving per gene conflicts
  if (res_conflict_per_gene) {
    if (!return_stats) {
      reads <- data.table::rbindlist(reads)
      row.names(reads)<-NULL
    }
    message(cat("\nResolving per-HLA Genotype Conflicts: assuming each cell has a max of 2 genotypes per HLA gene and keeping those with the most counts"))
    reads$cb_hla <- paste0(reads$CB,"_",reads$hla)
    reads$CB <- NULL # no longer needed
    reads <- split(data.table::setDT(reads), by = "cb_hla")
    keep_two <- function(df) {
      df <- data.table::setDF(df)
      n_hla <- table(df$gene0)
      if (length(n_hla) > 2) {
        max1 <- names(which(n_hla==max(n_hla))) %>% suppressWarnings()
        max2 <- names(which(n_hla[n_hla!=max(n_hla)]==max(n_hla[n_hla!=max(n_hla)]))) %>% suppressWarnings()
        if (length(max1) > 2) {
          # If there is a tie with 3 or more alleles, remove the entire df
          df <- NULL
        } else {
          if (length(max1) == 2) {
            # if 2 alleles are in a tie, we have our 2 alleles per cell!
            df <- df[df$gene0 %in% max1, ]
          } else {
            if (length(max1) == 1 & length(max2) == 1) {
              df <- df[df$gene0 %in% c(max1, max2), ]
            } else {
              # if the second most common allele is a list of 2 or more alleles, or zero alleles, we won't know which to pick, so leave out
              df <- df[df$gene0 %in% max1, ]
            }
          }
        }
      }
      return(df)
    }
    # Applying the function
    reads <- pbmcapply::pbmclapply(reads, keep_two, mc.cores = multi_thread)
    reads <- data.table::rbindlist(reads)
    row.names(reads)<-NULL
    if (return_stats) {
      reads_in_seu  <- as.numeric(reads$seu_barcode %in% colnames(seu) %>% table())[2]
      ureads_in_seu <- as.numeric(unique(reads$seu_barcode) %in% colnames(seu) %>% table())[2]
      l_ureads_cb  <- length(unique(reads$seu_barcode))
      n_reads <- nrow(reads)
      n_cells <- length(colnames(seu))
      message(cat("  Reads remaining: ", n_reads, 
                  ", including ", reads_in_seu,
                  " (", format(round(100*(reads_in_seu/n_reads), 2), nsmall = 1),
                  "%) belonging to Cells found in Seurat object", sep = ""))
      message(cat("  CBs remaining: ", l_ureads_cb, 
                  ", including ", 
                  ureads_in_seu, 
                  " (", format(round(100*(ureads_in_seu/l_ureads_cb), 2), nsmall = 1), 
                  "%) matching the ", n_cells, " Cells in Seurat object (match rate ",
                  format(round(100*(ureads_in_seu/n_cells), 2), nsmall = 1), "%)", sep = ""))
      stats_df <- rbind(stats_df,
                        data.frame(reads = n_reads, 
                                   reads_found = reads_in_seu, 
                                   cbs = l_ureads_cb, 
                                   cbs_found = ureads_in_seu, 
                                   cbs_seu = n_cells, 
                                   cb_seu_match_rate = ureads_in_seu/n_cells,
                                   step = "6_per_gene_conflict"))
    }
    reads <- split(data.table::setDT(reads), by = "seu_barcode")
  }
  ## Linkage Diseqilibrium correction in the DR region
  if (LD_correct) {
    message(cat("\nLinkage Disequilibrium Correction in the HLA-DR locus: assuming strong LD in
  the DR1  subregion haplotype: DRB1*01 and *10 in LD with DRB6 and DRB9
  the DR51 subregion haplotype: DRB1*15 and *16 in LD with DRB6, DRB5, and DRB9
  the DR52 subregion haplotype: DRB1*03, *11, *12, *13, and *14 in LD with DRB2, DRB3, and DRB9
  the DR8  subregion haplotype: DRB1*08 in LD with DRB9
  the DR53 subregion haplotype: DRB1*04, *07, and *09 in LD with DRB7, DRB8, DRB4, and DRB9"))
    LD <- function(df) {
      df <- data.table::setDF(df)
      ld <- list(
        DR1 = c("DRB1-01", "DRB1-10"),
        DR51 = c("DRB1-15", "DRB1-16"),
        DR52 = c("DRB1-03", "DRB1-11", "DRB1-12", "DRB1-13", "DRB1-14"),
        DR8 = c("DRB1-08"),
        DR53 = c("DRB1-04", "DRB1-07", "DRB1-09")
      )
      rmv <- list(
        DR1 = c("DRB2", "DRB3", "DRB4", "DRB5", "DRB7", "DRB8"), # haplotypes DRB1*01/*10 in LD with DRB6 and DRB9
        DR51 = c("DRB2", "DRB3", "DRB4", "DRB7", "DRB8"), # haplotypes DRB1*15/*16 in LD with DRB6, DRB5, and DRB9
        DR52 = c("DRB4", "DRB5", "DRB6", "DRB7", "DRB8"), # haplotypes DRB1*03/*11/*12/*13/*14 in LD with DRB2, DRB3, and DRB9
        DR8 = c("DRB2", "DRB3", "DRB4", "DRB5", "DRB6", "DRB7", "DRB8"), # haplotypes DRB1*08 in LD with DRB9
        DR53 =  c("DRB2", "DRB3", "DRB5", "DRB6") # haplotypes DRB1*04/*07/*09 in LD with DRB7, DRB8, DRB4, and DRB9
      )
      if (any(df$gene1 %in% ld[["DR1"]]))  {haplo <- c(DR1 = T)} else {haplo <- c(DR1 = F)}
      if (any(df$gene1 %in% ld[["DR51"]])) {haplo <- c(haplo, DR51 = T)} else {haplo <- c(haplo, DR51 = F)}
      if (any(df$gene1 %in% ld[["DR52"]])) {haplo <- c(haplo, DR52 = T)} else {haplo <- c(haplo, DR52 = F)}
      if (any(df$gene1 %in% ld[["DR8"]]))  {haplo <- c(haplo, DR8 = T)} else {haplo <- c(haplo, DR8 = F)}
      if (any(df$gene1 %in% ld[["DR53"]])) {haplo <- c(haplo, DR53 = T)} else {haplo <- c(haplo, DR53 = F)}
      for (i in 1:length(rmv)) {
        if (!haplo[i]) {rmv[[i]] <- NA}
      }
      rmv <- rmv[!is.na(rmv)]
      rmv <- Reduce(intersect, rmv)
      if (length(which(df$hla %in% rmv)) > 0) {df <- df[-which(df$hla %in% rmv),]}
      return(df)
    }    
    reads <- pbmcapply::pbmclapply(reads, LD, mc.cores = multi_thread)
    if (return_stats) {
      reads <- data.table::rbindlist(reads)
      row.names(reads)<-NULL
      reads_in_seu  <- as.numeric(reads$seu_barcode %in% colnames(seu) %>% table())[2]
      ureads_in_seu <- as.numeric(unique(reads$seu_barcode) %in% colnames(seu) %>% table())[2]
      l_ureads_cb  <- length(unique(reads$seu_barcode))
      n_reads <- nrow(reads)
      n_cells <- length(colnames(seu))
      message(cat("  Reads remaining: ", n_reads, 
                  ", including ", reads_in_seu,
                  " (", format(round(100*(reads_in_seu/n_reads), 2), nsmall = 1),
                  "%) belonging to Cells found in Seurat object", sep = ""))
      message(cat("  CBs remaining: ", l_ureads_cb, 
                  ", including ", 
                  ureads_in_seu, 
                  " (", format(round(100*(ureads_in_seu/l_ureads_cb), 2), nsmall = 1), 
                  "%) matching the ", n_cells, " Cells in Seurat object (match rate ",
                  format(round(100*(ureads_in_seu/n_cells), 2), nsmall = 1), "%)", sep = ""))
      stats_df <- rbind(stats_df,
                        data.frame(reads = n_reads, 
                                   reads_found = reads_in_seu, 
                                   cbs = l_ureads_cb, 
                                   cbs_found = ureads_in_seu, 
                                   cbs_seu = n_cells, 
                                   cb_seu_match_rate = ureads_in_seu/n_cells,
                                   step = "7_ld_correct"))
      reads <- split(data.table::setDT(reads), by = "seu_barcode")
      reads <- parallel::mclapply(reads, data.table::setDF, mc.cores = multi_thread)
    }
  }
  ## Matrix formation
  # matrix
  message(cat("\nCreating the HLA Count Matrix compatible with the Seurat object"))
  alleles <- mclapply(1:length(reads), function(x) reads[[x]]$gene0, mc.cores = multi_thread) %>% unlist() %>% na.omit() %>% unique() %>% sort()
  HLA.matrix <- matrix(0, nrow = length(alleles), ncol = length(reads), dimnames = list(alleles, names(reads)))
  pb <- pbmcapply::progressBar(min = 0, max = length(reads), style = "ETA", char = "=")
  for (i in 1:length(reads)) {
    counts <- table(reads[[i]]$gene0)
    HLA.matrix[, i] <- counts[alleles]
    setTxtProgressBar(pb, i)
  }
  close(pb)
  rm(reads)
  # merging with Seurat cell names
  HLA.matrix <- as.data.frame(t(HLA.matrix))
  HLA.matrix <- HLA.matrix[rownames(HLA.matrix) %in% colnames(seu),]
  cells <- data.frame(row.names = colnames(seu))
  HLA.matrix<- merge(HLA.matrix, cells, by= "row.names", all=T)
  rownames(HLA.matrix) <- HLA.matrix$Row.names
  HLA.matrix <- subset(HLA.matrix, select = -1)
  HLA.matrix <- HLA.matrix[order(rownames(HLA.matrix)),]
  HLA.matrix[is.na(HLA.matrix)] <- 0
  HLA.matrix <- HLA.matrix[,which(colSums(HLA.matrix)>0)]
  HLA.matrix <- as.matrix(t(HLA.matrix))
  HLA.matrix <- Matrix(HLA.matrix,sparse = T)
  HLA <- CreateAssayObject(counts = HLA.matrix)
  #message(cat("\nDone!!  ", format(Sys.time(), "%F %H:%M:%S"), " (runtime: ", difftime(Sys.time(), s, units = "min") %>% as.numeric() %>% abs(), " min)", sep = ""))
  e <- difftime(Sys.time(), s, units = "sec") %>% as.numeric() %>% abs()
  message(cat("\nDone!! (runtime: ", format(as.POSIXlt(e, origin = "1970-01-01", tz = "UTC"), "%H:%M:%S", tz = "UTC"), ")", sep = ""))
  if (return_stats) {
    return(list(matrix = HLA, per_step_stats = stats_df, umi_dupl_counts = umi_counts, umi_dupl_plot = g))
  } else {
    return(HLA)
  }
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
#' @param parallelize  is a logical, called TRUE if using parallel processing (multi-threading) is desired; default is FALSE.
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
#' HLA_clusters(reads = cts[["mRNA"]], k = 2, seu = your_Seurat_Obj, geno_metadata_id = "geno", hla_with_counts_above = 5, CBs_with_counts_above = 35)
#' @export

HLA_clusters <- function(reads, k = 2, seu = NULL, CB_rev_com = FALSE, geno_metadata_id = NULL, hla_with_counts_above = 0, CBs_with_counts_above = 0, match_CB_with_seu = TRUE, umap_first_n_PCs = 25, QC_mm2 = TRUE, s1_belowmax = 0.8, AS_belowmax = 0.8, NM_thresh = 15, de_thresh = 0.01, parallelize = FALSE, pt_size = 0.5, ...) {
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
    reads <- pbmcapply::pbmclapply(reads, function(df){df[df$s1 > s1_belowmax*max(df$s1) & df$AS > AS_belowmax*max(df$AS) & df$NM <= NM_thresh & df$de <= de_thresh,]}, mc.cores = multi_thread)
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
  rm(reads)
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
  # pcv<-cbind(pcv, seu@meta.data[match(row.names(pcv), colnames(seu)),])
  # ggplot(pcv, aes(x=Comp.1, y=Comp.2, color=logUMI))+geom_point(size=0.25)+scale_color_viridis_b()+theme_bw()
  # ggplot(pcv, aes(x=Comp.1, y=Comp.2, color=geno))+geom_point(size=0.5)+scale_color_manual(values=pals::glasbey())+theme_bw()
  umat<-pcv[,1:umap_first_n_PCs] %>% as.matrix()
  umapout<-uwot::umap(umat, verbose=T, ...)
  colnames(umapout)<-c("umap1", "umap2")
  umapout <- as.data.frame(umapout)
  if (!is.null(seu) & !is.null(geno_metadata_id)){
    umapout<-cbind(umapout, seu@meta.data[match(rownames(umapout), colnames(seu)),])
    g0 <- ggplot(umapout, aes(x=umap1, y=umap2, color=!!sym(geno_metadata_id)))+geom_point(size=pt_size)#+scale_color_manual(values=pals::glasbey())+theme_bw()
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
  if (!is.null(seu) & !is.null(geno_metadata_id)){
    return(list(UMAP_coordinates = umapout, HLA_clusters_on_umap = g, genotype_on_umap = g0))
  } else {
    return(list(UMAP_coordinates = umapout, HLA_clusters_on_umap = g))
  }
}


#' Mapping the HLA clusters found by the 'HLA_clusters()' function back into the scrHLAtag count files
#'
#' @param reads.list  is a scrHLAtag count file (or a list of scrHLAtag count files) including columns for CB, UMI, and HLA alleles (https://github.com/furlan-lab/scrHLAtag).
#' @param cluster_coordinates  is the UMAP coordinates dataframe with HLA clustering information (found by the 'HLA_clusters()' function) from which clusters are extracted and mapped into the scrHLAtag count files by matching Cell Barcodes. Currently the barcode format supported is: SAMPLE_AATGCTTGGTCCATTA-1
#' @param CB_rev_com  is a logical, called TRUE if the need to obtained the reverse complement of Cell Barcodes (CBs) is desired; default is FALSE.
#' @param additional_col_id  is a character; if mapping of any additional column ID of the Seurat metadata collected in the 'cluster_coordinates' dataframe, to the scrHLAtag count files is desired. 'NULL' by default.
#' @param parallelize  is a logical, called TRUE if using parallel processing (multi-threading) is desired; default is FALSE.
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

map_HLA_clusters <- function(reads.list, cluster_coordinates, CB_rev_com = FALSE, additional_col_id = NULL, parallelize = FALSE) {
  ## parallelize
  if (parallelize) {
    multi_thread <- parallel::detectCores()
    # message(cat("\nMulti-threading! Available cores: ", parallel::detectCores(), sep = ""))
  } else {
    multi_thread <- 1
  }  
  if (class(reads.list) %in% c("data.frame", "data.table")) {
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
  return(reads.list)
}
