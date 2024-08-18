#' Loading the count file(s) as an output of scrHLAtag, including columns for CB, UMI, and HLA alleles (\url{https://github.com/furlan-lab/scrHLAtag}).
#' 
#' @param directories  a character list of one (or more) full-length directory(ies) where the the scrHLAtag count file(s) have been organized.
#' @param dir_names  a character list of one (or more) sample name(s) that \emph{strictly correspond(s)} to the directory(ies) where the the scrHLAtag count file(s) have been organized. Can remain \code{NULL} only if the \code{directories} are already named vectors. For 2 or more samples/directories, try to get the names as closely as possible to the sample ID names spposedly present in the corresponding Seurat object (usually concatenated to Seurat barcodes when experiments are merged).
#' @param seu  is the Seurat object associated with the scrHLAtag count file (\url{https://satijalab.org/seurat/index.html}).
#' @param scrHLAtag_outs  a character list of the proper scrHLAtag count file(s) to be loaded. For example, for each sample, the desired count files should be the \code{"molecule_info_mRNA.txt.gz"} and the \code{"molecule_info_gene.txt.gz"} files.
#' @param colnames  \code{NULL} or a character list of the colnames of the scrHLAtag count file(s), if they are different from default output (scrHLAtag, version 0.1.5).
#' @import stringr
#' @import magrittr
#' @import Seurat
#' @import dplyr
#' @return a list of data.frames/data.tables
#' @examples
#' dirs_path <- "path/to/scrHLAtag/out/files"
#' dirs<-list.dirs(path=dirs_path, full.names = T, recursive = F)
#' dirs<- lapply(dirs, list.dirs, recursive = F) %>% unlist
#' dirs<- lapply(dirs, dir, pattern = "unguided_hla_align_corrected", recursive = F, full.names = T) %>% unlist
#' dirnames <- c("AML_101_BM", "AML_101_34", "TN_BM", "TN_34") # this is how the samples were organized in the directories
#' ## Load the counts files
#' cts <- HLA_load(directories = dirs, dir_names = dirnames, seu = your_Seurat_obj)
#' @export

HLA_load <- function(directories, dir_names = NULL, seu, scrHLAtag_outs = c("molecule_info_mRNA.txt.gz", "molecule_info_gene.txt.gz"), colnames = NULL) {
  s <- Sys.time()
  if (!requireNamespace("stringdist", quietly = TRUE)) { stop("Package 'stringdist' needed for this function to work. Please install it.", call. = FALSE) }
  mol_info <- scrHLAtag_outs
  if (is.null(names(directories))) {
    if (is.null(dir_names)) { stop("'directories' need names. Please provide the names in 'dir_names' and make sure the naming order is strictly as organized in your directories.", call. = FALSE) }
    if (length(directories) != length(dir_names)) { stop("Number of 'directories' should match number of sample names in 'dir_names'. Also, make sure the naming order is strictly as organized in your directories.", call. = FALSE) }
    names(directories) <- dir_names
  }
  ## extract sample names from the Seurat colnames 'Cells(seu)'
  if (!("Seurat" %in% class(seu))) { stop("The function will extract sample ID(s) concatenated to Seurat barcodes. Object in 'seu' must be of class 'Seurat'", call. = FALSE) }
  v <- colnames(seu)
  prefx <- substr(v, 1, nchar(v) - sapply(base::gsub(".*[ATGC]+(.*)$", "\\1", v), nchar))
  prefx <- substr(prefx, 1, sapply(base::gregexpr("[^A-Za-z0-9]", prefx), function(pos) if (max(pos)>0) max(pos) else 0))
  sepr <- sapply(prefx, function(x) substr(x, nchar(x), nchar(x)))
  smpl <- substr(prefx, 1, nchar(prefx) - sapply(sepr, nchar))
  sufx <- base::gsub(".*[ATGC]+(.*)$", "\\1", v)
  prefx.sufx <- paste0(prefx, "*", sufx)
  rv <- unlist(lapply(v, function(x) intToUtf8(rev(utf8ToInt(x)))))
  cbrv <- gsub("^[^ATGC]*([ATGC]+).*", "\\1", rv)
  cbrv <- unlist(lapply(cbrv, function(x) intToUtf8(rev(utf8ToInt(x)))))
  smpls <- unique(smpl)
  prefxs <- unique(prefx)
  prefx.sufxs <- unique(prefx.sufx)
  if (length(smpls) < length(prefx.sufxs)) {
    message(cat("Sample IDs found in Seurat:", crayon::bgWhite(" ", prefx.sufxs, " ")))
    } else {
    message(cat("Sample IDs found in Seurat:", crayon::bgWhite(" ", smpls, " ")))
    }
  # length(unique(sepr))
  if (length(unique(sapply(cbrv, nchar))) > 1) { message(cat(crayon::red("The Barcode DNA sequences in the Seurat object have different lengths. Does the object contains different capture chemistries?"), sep = "")) }
  if (length(prefx.sufxs) > length(directories)) { stop("More samples found in Seurat than scrHLAtag selected directories. Make sure to select all experiements linked to the Seurate object, or subset your Seurat to the corresponding samples in your scrHLAtag experiments", call. = FALSE) }
  message(cat("Loading..."))
  cts <- lapply(mol_info, function(x) {
    dl<-lapply(prefx.sufxs, function(sample){
      prefx1 <- stringr::str_split_fixed(sample, "\\*", 2)[,1]
      samp1 <- ifelse(nchar(prefx1) > 0, substr(prefx1, 1, nchar(prefx1) - 1), "")
      sepr1 <- ifelse(nchar(prefx1) > 0, substr(prefx1, nchar(prefx1), nchar(prefx1)), "")
      sufx1 <- stringr::str_split_fixed(sample, "\\*", 2)[,2]
      d<-read.table(file.path(directories[which.min(stringdist::stringdist(samp1, names(directories), method = "lv"))], x), header = F, sep=" ", fill = T) 
      if (is.null(colnames)) {
        colnames(d)<-c("name","CB", "nb", "UMI", "gene", "query_len","start", "mapq", "cigar", "NM", "AS", "s1", "de", "seq")
      } else {
        colnames(d) <- colnames
      }
      d$samp <- samp1
      d$id_cb_separator <- sepr1
      d$id_cb_suffix <- sufx1
      return(d)
    }) #, mc.cores = parallel::detectCores())
    ctsu<-do.call(rbind,dl)
    if (nrow(which(is.na(ctsu), arr.ind=T)) > 0) {
      message(cat(nrow(which(is.na(ctsu), arr.ind=T)), " problematic row(s) detected; will be omitted.", sep = ""))
      print(which(is.na(ctsu), arr.ind=T))
      ctsu <- na.omit(ctsu)
      print(sapply(ctsu, class))
      if (is.null(colnames)) {
        # fix the classes in case there was a bug
        ctsu <- ctsu %>% mutate(across(c(13), as.numeric), across(c(3,6,7,8,10,11,12), as.integer), across(c(1,2,4,5,9,14,15), as.character))
      }
    }
    rm(dl)
    return(ctsu)
  }) #, mc.cores = parallel::detectCores())
  for (i in seq_along(mol_info)) {
    if (length(mol_info) == 2) {
      difchar <- which(strsplit(mol_info[1], "")[[1]] != strsplit(mol_info[2], "")[[1]]) 
      nms <- strsplit(mol_info[i], "")[[1]][difchar] %>% paste(., collapse = "")
      names(cts)[i] <- nms
    } else {
      nms <- mol_info[i]
      names(cts)[i] <- nms
    }
    message(cat("Unique alleles mapped through the ", nms, " reference: ", cts[[i]]$gene %>% unique() %>% sort() %>% length(), sep = ""))
  }
  e <- difftime(Sys.time(), s, units = "sec") %>% as.numeric() %>% abs()
  message(cat("Done!! (runtime: ", format(as.POSIXlt(e, origin = "1970-01-01", tz = "UTC"), "%H:%M:%S", tz = "UTC"), ")", sep = ""))
  return(cts)
}