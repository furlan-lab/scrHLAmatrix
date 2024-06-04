#' A comprehensive function to extract the top HLA alleles from the scrHLAtag count files (based on the most counts per Cell Barcode)
#' 
#' @param reads_1  is the primary scrHLAtag count file (1 of 2 files containing either the mRNA molecular info or the genomic (gene) molecular info). It includes columns for CB, UMI, and HLA alleles (https://github.com/furlan-lab/scrHLAtag).
#' @param reads_2  is the secondary scrHLAtag count file (the alternative file vs. the one designated in 'reads_1' argument). It includes columns for CB, UMI, and HLA alleles (https://github.com/furlan-lab/scrHLAtag). Default is NULL, in which case it will not be able to count alternative aligment.
#' @param allogeneic_entities  is the number 'k' of clusters to partition the UMAP space into, e.g. the number of entities or genotypes you 'think' there might be in your captured sample.
#' @param seu  is the Seurat object associated with the scrHLAtag count file (https://satijalab.org/seurat/index.html), and entered here if matching CBs in count file with Seurat colnames is desired.
#' @param CB_rev_com  is a logical, called TRUE if the need to obtained the reverse complement of Cell Barcodes (CBs) is desired; default is FALSE. 
#' @param hla_with_counts_above  is the number of total reads accross CBs at or above which an HLA allele is retained in the matrix.
#' @param CBs_with_counts_above  is the number of total reads accross HLA alleles at or above which a CB is retained in the matrix. Note: at present, the function will make sure that number of CBs is equal or more than available HLA alleles in the matrix.
#' @param match_CB_with_seu  is a logical, called TRUE if filtering CBs in the scrHLAtag count file with matching ones in the Seurat object is desired. 
#' @param top_by_read_frac.blk  in the Pseudo-Bulk algorithm, is the fraction (0 to 1) of total reads for a particular HLA gene, which incorporates the highest ranking alleles of that gene in terms of number of reads; default at 0.75
#' @param top_by_read_frac.cb  in the Single-Cell algorithm, is the fraction (0 to 1) of total reads for a particular HLA gene, which incorporates the highest ranking alleles of that gene in terms of number of reads; default at 0.85
#' @param allowed_alleles_per_cell  is a numeric (single or range) determining the minimum and maximum number of highest ranking allele genotypes per cell to keep if such number is beyond those limits when filtering by fraction 'frac'; default is c(1, 200), usefull in the early scrHLAtag iterations to give minimap2 lots of room to align; once you are ready to finalize the top HLA allele list, you can try c(1, 2) if you assume a cell can have a min of 1 allele (homozygous) and a max of 2 (heterozygous).
#' @param field_resolution  is a numeric, to select the HLA nomenclature level of Field resolution, where 1, 2, or 3 will take into consideration the first, the first two, or the first three field(s) of HLA designation; default is 3.
#' @param QC_mm2  is a logical, called TRUE if removing low quality reads based on minimap2 tags is desired.
#' @param s1_belowmax  is a proportion (0 to 1) of the maximum value (best quality) of the minimap2 's1' tag above which the quality of the read is acceptable; default at 0.75 of the max s1 score.
#' @param AS_belowmax  is a proportion (0 to 1) of the maximum value (best quality) of the minimap2 'AS' tag above which the quality of the read is acceptable; default at 0.85 of the max AS score.
#' @param NM_thresh  is the number of mismatches and gaps in the minimap2 alignment at or below which the quality of the read is acceptable; default is 15.
#' @param de_thresh  is the gap-compressed per-base sequence divergence in the minimap2 alignment at or below which the quality of the read is acceptable; the number is between 0 and 1, and default is 0.015.
#' @param parallelize  is a logical, called TRUE if using parallel processing (multi-threading) is desired; default is TRUE.
#' @param umap_spread  for uwot::umap(); effective scale of embedded points determining how clustered/clumped the embedded points are.
#' @param umap_min_dist  for uwot::umap(); effective minimum distance between embedded points. Smaller values will result in a more clustered/clumped embedding.
#' @param umap_repulsion  for uwot::umap(); weighting applied to negative samples in low dimensional embedding optimization.
#' @param ...  other arguments passed onto uwot::umap()
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
#' top_alleles <- Top_HLA_list(reads_1 = cts[["mRNA"]], reads_2 = cts[["gene"]], seu = your_Seurat_Obj)
#' # 
#' # Note: the function is optimized to choose whether to find top HLA alleles
#' # based on a Single-Cell approach (large number of uniquely mapped alleles) or 
#' # a Pseudo-Bulk approach (number of alleles is generally < 2000).
#' # 
#' # Note: if for a particular HLA, the alleles with the most counts are in a tie 
#' # between 3 or more alleles in a particular Cell Barcode, we cannot know which 
#' # ones are the top two alleles, so that CB is not counted. This is similar if  
#' # there were no counts for that allele (all zeros).
#' @export

Top_HLA_list <- function(reads_1, reads_2 = NULL, allogeneic_entities = 2, seu = NULL, CB_rev_com = FALSE,
                         hla_with_counts_above = 5, CBs_with_counts_above = 15, match_CB_with_seu = TRUE, 
                         QC_mm2 = TRUE, s1_belowmax = 0.8, AS_belowmax = 0.8, NM_thresh = 15, de_thresh = 0.01,
                         top_by_read_frac.blk = 0.75, top_by_read_frac.cb = 0.85,
                         allowed_alleles_per_cell = c(1, 200), field_resolution = 3, parallelize = TRUE, 
                         umap_spread = 3, umap_min_dist = 0.0001, umap_repulsion = 0.0001, ...) {
  s <- Sys.time()
  # creating HLA umap clusters
  HLA_umap_clusters <- HLA_clusters(reads = reads_1, 
                                    k = allogeneic_entities, 
                                    seu = seu, 
                                    CB_rev_com = CB_rev_com,  # TRUE for 3prime 10x on pacbio
                                    match_CB_with_seu = match_CB_with_seu,
                                    hla_with_counts_above = hla_with_counts_above, 
                                    CBs_with_counts_above = CBs_with_counts_above,
                                    QC_mm2 = QC_mm2, s1_belowmax = s1_belowmax, 
                                    AS_belowmax = AS_belowmax, NM_thresh = NM_thresh, de_thresh = de_thresh,
                                    return_heavy = TRUE, 
                                    spread = umap_spread, min_dist = umap_min_dist, repulsion_strength = umap_repulsion, ...)
  print(HLA_umap_clusters[[2]]+scale_color_manual(values=pals::glasbey())+theme_classic())
  if ((reads_1$gene %>% unique() %>% sort() %>% length()) > 2000) {
    message(cat("\nReads count file shows greater than 2000 mapped HLA alleles; extracting top alleles using the Pseudo-Bulk algorithm"))
    reads_1 <- map_HLA_clusters(reads.list = reads_1, HLA_umap_clusters[[1]], CB_rev_com = CB_rev_com)
    if (!is.null(reads_2)) {reads_2 <- map_HLA_clusters(reads.list = reads_2, HLA_umap_clusters[[1]], CB_rev_com = CB_rev_com)}
    cluster <- levels(reads_1$hla_clusters)
    alt_align <- rep(c(TRUE, FALSE), length.out = length(cluster))
    top_alleles_HLA <- c()
    for (i in seq_along(cluster)) {
      top_al <- Top_HLA_list_bulk(reads_1 = reads_1[reads_1$hla_clusters %in% cluster[i],], 
                                  reads_2 = reads_2[reads_2$hla_clusters %in% cluster[i],], 
                                  frac = top_by_read_frac.blk,
                                  insert_pop_most_freq = T,   # recommended to be TRUE for first iteration
                                  use_alt_align_ABC = alt_align[i])
      top_alleles_HLA <- c(top_alleles_HLA, top_al)
    }
    top_alleles_HLA <- top_alleles_HLA %>% unique() %>% sort()
  } else {
    message(cat("\nReads count file shows 2000 or fewer mapped HLA alleles; extracting top alleles using the Single-Cell algorithm"))
    reads_1 <- map_HLA_clusters(reads.list = HLA_umap_clusters[["reads"]], HLA_umap_clusters[[1]], CB_rev_com = F)
    top_alleles_HLA <- Top_HLA_list_byCB_preprocessed(reads = reads_1,
                                                      seu = seu,
                                                      match_CB_with_seu = match_CB_with_seu,
                                                      hla_with_counts_above = hla_with_counts_above,
                                                      CBs_with_counts_above = CBs_with_counts_above,
                                                      frac = top_by_read_frac.cb,
                                                      allowed_alleles_per_cell = allowed_alleles_per_cell,
                                                      field_resolution = field_resolution,
                                                      parallelize = parallelize)
  }
  e <- difftime(Sys.time(), s, units = "sec") %>% as.numeric() %>% abs()
  message(cat("\nDone!! (runtime: ", format(as.POSIXlt(e, origin = "1970-01-01", tz = "UTC"), "%H:%M:%S", tz = "UTC"), ")", sep = ""))
  return(top_alleles_HLA)
}