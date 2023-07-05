HLA_Matrix <- function(cts, seu, rp = character(), dn = character(), Ct = 0) {
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
  ## check format of 'rp' and 'dn'
  if (any(grepl(special, rp))) {
    rp <- gsub(special, "-", rp)
  } else if (all(grepl("-", rp))) {
    rp <- rp
  } else {
    stop("Incorrect format for recipient-defined and/or donor-defined HLA alleles. \nMake sure gene and allele are separated by standard nomenclature asterisk (or other special character).")
  }
  if (any(grepl(special, dn))) {
    dn <- gsub(special, "-", dn)
  } else if (all(grepl("-", dn))) {
    dn <- dn
  } else {
    stop("Incorrect format for recipient-defined and/or donor-defined HLA alleles. \nMake sure gene and allele are separated by standard nomenclature asterisk (or other special character).")
  }
  ## check class of 'Ct'
  if (class(Ct) != "integer"){
    if (class(Ct) != "numeric"){
      stop("Count threshold 'Ct' must be numeric or integer")
    }
  }
  ## check relevant col names exist
  if (!all(c("CB", "UMI", "gene", "samp") %in% colnames(cts))){
    stop("scrHLAtag output 'cts' dataframe must at least contain the columns 'CB', 'UMI', 'gene' (with HLA alleles), and 'samp' (matching the sample names in the corresponding Seurat object)")
  }
  
  ## see if more than 1 allele are present per umi at a time
  # count all the problematic CB:UMIs for which a molecular swap is suspected
  cts$mol_swap <- NA
  cts$mol_swap <- as.factor(cts$mol_swap)
  cts$class_swap <- NA
  cts$class_swap <- as.factor(cts$class_swap)
  # split,   this is computationally heavy (about 10min for 10M rows)
  cts.split <- with(cts, split(cts, list(cbumi=cbumi))) 
  print("1/6 - Estimating UMI duplication rate", quote=F)
  n_umi<-pbmclapply(cts.split, nrow, mc.cores = parallel::detectCores()) %>% unlist()
  umi_counts<- data.frame(n_umi)
  umi_counts$dummy <- 1
  umi_counts <- umi_counts[order(umi_counts$n_umi),]
  row.names(umi_counts)<- NULL
  g <- ggplot(umi_counts, aes(x= as.numeric(row.names(umi_counts)), y=n_umi, group= dummy))+
    #geom_smooth(size=2, method = "gam")+
    geom_line()+
    scale_y_log10(name = "PCR copies per UMI")+
    scale_x_continuous(name = "Rank", n.breaks = 8) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  print(g)
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
    cts.dedup.cb[[j]]$hla_conflict <- ifelse(any(cts.dedup.cb[[j]]$gene0 %in% rp) & any(cts.dedup.cb[[j]]$gene0 %in% dn),
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
  cts.dedup.cb <- pbmclapply(cts.dedup.cb, keep_two, recip = rp, donor = dn, mc.cores = parallel::detectCores())
  
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
