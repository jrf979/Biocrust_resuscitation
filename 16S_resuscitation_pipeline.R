#title: "16S DADA2 Pipeline + 97% clustering + taxonomy assignment"
#author: "Raul Roman"
#date: "2025-10-08"


## ============================ 16S DADA2 pipeline ============================
rm(list=ls()); gc()
suppressPackageStartupMessages({
  library(dada2)
  library(BiocParallel)
  library(tidyverse)
})

set.seed(1)

## ---- Paths ----
WORKDIR <- "/storage/work/jxr6215/16S_DADA2_20250923"
PATH    <- file.path(WORKDIR, "sequences/cutadapt_trimmed")

dir.create(file.path(WORKDIR, "dada2_out"), showWarnings = FALSE, recursive = TRUE)
setwd(file.path(WORKDIR, "dada2_out"))

get_ncores <- function() {
  x <- Sys.getenv("SLURM_CPUS_PER_TASK")
  if (nzchar(x)) as.integer(x) else max(1L, parallel::detectCores() - 1L)
}
NCORES <- get_ncores()
register(MulticoreParam(workers = NCORES))
message("Using ", NCORES, " cores.")

## ===================== Input trimmed sequences from cutadapt ========================
F.sam <- sort(list.files(PATH, pattern = "_R1_001\\.trim\\.fastq\\.gz$", full.names = TRUE))
R.sam <- sort(list.files(PATH, pattern = "_R2_001\\.trim\\.fastq\\.gz$", full.names = TRUE))
stopifnot(length(F.sam) == length(R.sam), length(F.sam) > 0)

# Use unique filenames
bnF   <- basename(F.sam)
bnR   <- basename(R.sam)
sample.ids <- sub("_R1_001\\.trim\\.fastq\\.gz$", "", bnF)     # e.g., "100_S173_L001"
stopifnot(identical(sample.ids, sub("_R2_001\\.trim\\.fastq\\.gz$", "", bnR)))
stopifnot(anyDuplicated(sample.ids) == 0)

# Name vectors by sample ID
names(F.sam) <- names(R.sam) <- sample.ids

## ================== QC plots  ==========
suppressPackageStartupMessages({ library(ggplot2) })
plots_dir <- "00_qc_plots"; dir.create(plots_dir, showWarnings = FALSE)

idx <- seq_len(min(4, length(F.sam)))
QPforward <- plotQualityProfile(F.sam[idx])
QPreverse <- plotQualityProfile(R.sam[idx])
ggplot2::ggsave(file.path(plots_dir, "quality_forward.pdf"), QPforward, width = 10, height = 6)
ggplot2::ggsave(file.path(plots_dir, "quality_reverse.pdf"), QPreverse, width = 10, height = 6)
try({
  ggplot2::ggsave(file.path(plots_dir, "quality_forward.png"), QPforward, width = 10, height = 6, dpi = 150)
  ggplot2::ggsave(file.path(plots_dir, "quality_reverse.png"), QPreverse, width = 10, height = 6, dpi = 150)
}, silent = TRUE)

## =============== Filter and Trim =====================
filt_path <- "01_filtered"; dir.create(filt_path, showWarnings = FALSE)
filtFs <- file.path(filt_path, paste0(sample.ids, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.ids, "_R_filt.fastq.gz"))
names(filtFs) <- names(filtRs) <- sample.ids

out <- filterAndTrim(
  fwd = F.sam, filt = filtFs,
  rev = R.sam, filt.rev = filtRs,
  truncLen = c(230, 200),   # from visual inspection of QC plots
  trimLeft = c(0, 0),
  maxN = 0, maxEE = c(2, 5), truncQ = 2,
  rm.phix = TRUE, compress = TRUE,
  multithread = NCORES,
  matchIDs = TRUE, minLen = 20
)

# Ensure rownames are our sample IDs 
if (!identical(rownames(out), sample.ids)) rownames(out) <- sample.ids
write.csv(out, "01_filter_counts.csv")
saveRDS(list(filtFs=filtFs, filtRs=filtRs, out=out), "01_filter_env.rds")

## =============== Learn errors  =================
set.seed(2)
errF <- learnErrors(filtFs, multithread = FALSE, randomize = TRUE)
errR <- learnErrors(filtRs, multithread = FALSE, randomize = TRUE)
saveRDS(errF, "02_errF.rds")
saveRDS(errR, "02_errR.rds")

## =============== Dereplicate and denoise ========================
names(filtFs) <- names(filtRs) <- sample.ids

denoise_one <- function(sname) {
  drpF <- derepFastq(filtFs[[sname]], verbose = TRUE)
  drpR <- derepFastq(filtRs[[sname]], verbose = TRUE)
  ddF  <- dada(drpF, err = errF, multithread = FALSE, pool = FALSE)
  ddR  <- dada(drpR, err = errR, multithread = FALSE, pool = FALSE)
  mg   <- mergePairs(ddF, drpF, ddR, drpR, verbose = FALSE)
  list(dadaF = ddF, dadaR = ddR, merge = mg)
}

res_list <- bplapply(sample.ids, denoise_one, BPPARAM = MulticoreParam(NCORES))
saveRDS(res_list, "03_denoise_merge_list.rds")

mergers <- setNames(lapply(res_list, `[[`, "merge"), sample.ids)
seqtab  <- makeSequenceTable(mergers)
saveRDS(seqtab, "04_seqtab_raw.rds")

## =============== Check that primers are gone after cutadapt trimming step ===================
library(Biostrings)

F   <- DNAString("GTGYCAGCMGCCGCGGTAA")        # 515F-Y
R   <- DNAString("GGACTACNVGGGTWTCTAAT")       # 806R-B
Rrc <- reverseComplement(R)

asv <- DNAStringSet(colnames(seqtab))
hitsF  <- vmatchPattern(F, asv, fixed = FALSE)
hitsRr <- vmatchPattern(Rrc, asv, fixed = FALSE)

startsF <- vapply(hitsF, function(h) any(start(h) == 1), logical(1))
endsRr  <- mapply(function(h, w) any(end(h) == w), hitsRr, width(asv))

primer_hits <- sum(startsF | endsRr)
cat("Primer hits remaining after Cutadapt:", primer_hits, "\n")

## =============== Inspect distribution of sequence lengths ===================
length_distribution <- table(nchar(getSequences(seqtab))) 
print(length_distribution)# The tails (≤250 or ≥260, especially ≥270 or ≫300) are mostly chimeras/non-specifics/borderline merges and should be removed

# Filter length around 250:256
keep_lens <- 250:256
seqtab_len <- seqtab[, nchar(colnames(seqtab)) %in% keep_lens, drop = FALSE]
saveRDS(seqtab_len, "04_seqtab_lenfiltered.rds")

cat("Read retention (len-gate): ", sum(seqtab_len) / sum(seqtab), "\n")
cat("ASV  retention (len-gate): ", ncol(seqtab_len) / ncol(seqtab), "\n")

## =============== Chimera removal ===================
seqtab.nochim <- removeBimeraDenovo(
  seqtab_len,
  method = "consensus",
  multithread = NCORES,
  verbose = TRUE
)
saveRDS(seqtab.nochim, "05_seqtab_nochim.rds")

## =============== Track lost reads across the pipeline =================
getN <- function(x) sum(dada2::getUniques(x))

# Safety helper for 'out' columns
get_col <- function(m, nm) {
  if (nm %in% colnames(m)) m[, nm] else rep(NA_integer_, nrow(m))
}

out_counts <- as.matrix(out)
if (!identical(rownames(out_counts), sample.ids)) rownames(out_counts) <- sample.ids

track <- tibble(
  sample     = sample.ids,
  input      = get_col(out_counts, "reads.in"),
  filtered   = get_col(out_counts, "reads.out"),
  denoisedF  = sapply(res_list, function(z) getN(z$dadaF)),
  denoisedR  = sapply(res_list, function(z) getN(z$dadaR)),
  mergedSeqs = sapply(mergers, nrow),
  tabled     = rowSums(seqtab_len[sample.ids, , drop = FALSE]),
  nonchim    = rowSums(seqtab.nochim[sample.ids, , drop = FALSE])
)

write.csv(track, "06_track_reads.csv", row.names = FALSE)

## ========================= Export ASV table counts ===========================
# At this point, colnames(seqtab.nochim) are still the DNA sequences.
# Taxonomy is not assigned here, it is assigned later after 97% clustering.

asv_seqs <- colnames(seqtab.nochim)
asv_ids  <- paste0("ASV", seq_len(ncol(seqtab.nochim)))

colnames(seqtab.nochim) <- asv_ids

write.csv(seqtab.nochim, "ASV_table_counts.csv", quote = FALSE)
writeLines(paste0(">", asv_ids, "\n", asv_seqs), con = "ASV_sequences.fa")


## ======================= CLUSTER ASVs BY 97% SEQUENCE IDENTITY + Taxonomy =======================

library(Biostrings)
library(DECIPHER)
library(dada2)
library(data.table)
library(dplyr)

## =============== load data  =================
asv_fasta  <- "ASV_sequences.fa"
asv_counts <- "ASV_table_counts.csv"

SILVA_TRAIN <- "/storage/work/jxr6215/16S_DADA2_20250923/silva_nr99_v138.2_toGenus_trainset.fa.gz"

## =============== set parameters  =================
perc_identity <- 0.97
cutoff <- 1 - perc_identity

## =============== load counts  =================
cts <- read.csv(asv_counts, check.names = FALSE)
stopifnot(ncol(cts) >= 2)
sample_ids <- cts[[1]]
cts <- as.data.frame(cts[ , -1, drop = FALSE])
rownames(cts) <- sample_ids

## =============== load sequences  =================
dna <- readDNAStringSet(asv_fasta)
stopifnot(length(dna) > 0)

# Keep only ASVs present in counts; order to match counts
common_asvs <- intersect(colnames(cts), names(dna))
cts <- cts[, common_asvs, drop = FALSE]
dna <- dna[common_asvs]

# Check that there are no zero-length sequences
dna <- dna[width(dna) > 0]
stopifnot(identical(names(dna), colnames(cts)))

## =============== Orient, align + cluster  =================
dna <- OrientNucleotides(dna)

aln <- AlignSeqs(dna, processors = NCORES, verbose = TRUE)

d <- DistanceMatrix(
  aln,
  type = "dist",
  includeTerminalGaps = FALSE,
  penalizeGapLetterMatches = FALSE,
  correction = "none",
  processors = NCORES,
  verbose = TRUE
)

tree <- TreeLine(
  myDistMatrix = d,
  method = "UPGMA",
  type   = "both",# dendrogram + clusters
  cutoff = cutoff, # The 0.03 for 97% ID
  processors = NCORES,
  verbose = TRUE
)

clus_vec <- tree[[1]][["cluster"]]

# This is to make sure ASV names are attached
if (is.null(names(clus_vec)) || any(!nzchar(names(clus_vec)))) {
  labs <- attr(d, "Labels")
  if (is.null(labs)) labs <- names(dna)
  stopifnot(length(labs) == length(clus_vec))
  names(clus_vec) <- labs
}

cl_df <- data.frame(
  ASV        = as.character(names(clus_vec)),
  cluster_id = as.integer(clus_vec),
  row.names  = NULL,
  stringsAsFactors = FALSE
)

## =============== Select abundance centroids  =================

# Pick abundance centroid per cluster
asv_abund <- colSums(cts)
centroids <- cl_df |>
  dplyr::mutate(total = asv_abund[ASV]) |>
  dplyr::group_by(cluster_id) |>
  dplyr::slice_max(order_by = total, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::select(cluster_id, rep_ASV = ASV)

# Map ASV -> representative
map_df <- cl_df |>
  dplyr::left_join(centroids, by = "cluster_id") |>
  dplyr::select(ASV, cluster_id, rep_ASV)

# Collapse counts to representatives
dt <- data.table::as.data.table(cts, keep.rownames = "sample")
long <- data.table::melt(dt, id.vars = "sample", variable.name = "ASV", value.name = "count")[count > 0]
long <- long[map_df, on = "ASV"]
collapsed <- long[, .(count = sum(count)), by = .(sample, rep_ASV)]
wide <- data.table::dcast(collapsed, sample ~ rep_ASV, value.var = "count", fill = 0)

clustered_counts <- as.data.frame(wide)
rownames(clustered_counts) <- clustered_counts$sample
clustered_counts$sample <- NULL


## =============== Representative fasta  =================
rep_ids <- unique(centroids$rep_ASV)
rep_dna <- dna[rep_ids]
names(rep_dna) <- rep_ids

## =============== Export outputs  =================
write.csv(cbind(sample = rownames(clustered_counts), clustered_counts), "ASV_table_clustered.csv", row.names = FALSE)

write.csv(map_df,"ASV_cluster_map.csv",row.names = FALSE)

writeXStringSet(rep_dna,"ASV_representatives.fa",compress = FALSE)

# ==================== Taxonomy is assigned to cluster representatives ====================
taxa <- dada2::assignTaxonomy(
  seqs = as.character(rep_dna),
  refFasta = SILVA_TRAIN,
  multithread = NCORES,
  tryRC = TRUE
)

taxa_df <- as.data.frame(taxa, stringsAsFactors = FALSE)

taxa_df <- tibble::tibble(ASV = names(rep_dna)) %>%
  bind_cols(taxa_df) %>%
  relocate(ASV)

write.csv(
  taxa_df,
  "Tax_table_clustered.csv",
  row.names = FALSE
)
