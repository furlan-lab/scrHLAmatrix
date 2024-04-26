#' Extracting the top HLA alleles from the scrHLAtag count files in a Pseudo-Bulk approach
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
#' Top_HLA_list_bulk(reads_1 = cts[["mRNA"]], reads_2 = cts[["gene"]], frac = 0.8, min_alleles_keep = 2, use_alt_align_ABC = TRUE)
#' @export

Top_HLA_list_bulk <- function(reads_1, reads_2 = NULL, frac = 0.75, min_alleles_keep = 5, min_reads_per_gene = 20, insert_pop_most_freq = TRUE, use_alt_align_ABC = FALSE){
  if (is.null(reads_2)) {
    reads_2 <- reads_1
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



#' Plotting the top HLA alleles extracted in a Pseudo-Bulk approach from the scrHLAtag count files 
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
#' Top_HLA_plot_bulk(reads_1 = cts[["mRNA"]], reads_2 = cts[["gene"]], use_alt_align_ABC = TRUE)
#' @export

Top_HLA_plot_bulk <- function(reads_1, reads_2 = NULL, top_hla = 10, min_reads_per_gene = 20, use_alt_align_ABC = FALSE, color_pal = NULL){
  if (is.null(reads_2)) {
    reads_2 <- reads_1
    warning("The molecule_info_gene.txt.gz count file does not seem to be included. The function will run but the argument 'use_alt_align_ABC' will be irrelevant.")
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
  print(pl2)
  return(pl2)
}
