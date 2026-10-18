#!/usr/bin/env Rscript
# FDR threshold for "high SNP enrichment" in functional ncRNAs
#
# Approach
#   Null distribution : enrichment values of pseudogenes (Gene_Type), i.e. the
#                       SNP enrichment expected for sequence without function.
#   Test              : for each gene, empirical upper-tail p-value
#                       p = (1 + #{null >= x}) / (N_null + 1)
#   Correction        : Benjamini-Hochberg across functional genes
#   Threshold         : smallest enrichment among functional genes with q <= FDR
#                       (genes with enrichment >= threshold = "high enrichment")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

# ---------------------------- configuration ---------------------------------
infile <- "../../results/combined_gene_data_with_snp_promoters.csv"
outdir <- "../../results"
fdr    <- 0.05

# Enrichment columns to test (one threshold is computed per SNP database)
enrich_cols <- c("enrichment_intergenic_1000genomes",
                 "enrichment_intergenic_gnomad",
                 "enrichment_intergenic_pangenome",
                 "enrichment_intergenic_dbsnp")

# Column used for the summary table's "Excess SNPs" column
primary_col <- "enrichment_intergenic_dbsnp"

# Values in Gene_Type -- ADJUST after checking the printed unique values.
# Matching is case-insensitive and by regex.
pseudo_pattern     <- "pseudo"
functional_pattern <- "^(?!.*pseudo)"   # anything that is not a pseudogene
# e.g. use "^functional$" if Gene_Type is strictly functional/pseudogene

# Functional groups: name, families (as they appear in Gene_group, lower case), polymerase
# NB: rnu5 is deliberately listed in both the major and minor spliceosome rows
# (same genes counted in both), as in the target table.
groups <- list(
  list(name = "Major spliceosome (Pol II)",  fam = c("rnu1","rnu2","rnu4","rnu5"),            pol = "Pol II"),
  list(name = "Major spliceosome (Pol III)", fam = c("rnu6"),                                 pol = "Pol III"),
  list(name = "Minor spliceosome (Pol II)",  fam = c("rnu11","rnu12","rnu4atac","rnu5"),      pol = "Pol II"),
  list(name = "Minor spliceosome (Pol III)", fam = c("rnu6atac"),                             pol = "Pol III"),
  list(name = "Translation",                 fam = c("trna"),                                 pol = "Pol III"),
  list(name = "RNP / RNA surveillance?",     fam = c("vtrna","rny"),                          pol = "Pol II, Pol III"),
  list(name = "Specialised RNA processing / targeting", fam = c("rn7sl","rnu7"),              pol = "Pol III, Pol II"),
  list(name = "Transcription regulation",    fam = c("rn7sk"),                                pol = "Pol III")
)

# Display labels for the "ncRNA families" column
fam_label <- c(rnu1="U1", rnu2="U2", rnu4="U4", rnu5="U5", rnu6="U6",
               rnu11="U11", rnu12="U12", rnu4atac="U4ATAC", rnu6atac="U6ATAC",
               rnu7="U7", trna="tRNA", vtrna="VTRNA", rny="Y RNA",
               rn7sk="7SK", rn7sl="7SL")

# ------------------------------ load data -----------------------------------
# The header contains duplicated column names; keep the first occurrence only
dat <- read_csv(infile, show_col_types = FALSE, name_repair = "minimal")
dat <- dat[, !duplicated(names(dat))]
dat$Gene_group <- tolower(trimws(dat$Gene_group))

cat("Unique Gene_Type values:\n");  print(table(dat$Gene_Type, useNA = "ifany"))
cat("\nUnique Gene_group values:\n"); print(table(dat$Gene_group, useNA = "ifany"))

dat <- dat %>%
  mutate(class = case_when(
    grepl(pseudo_pattern, Gene_Type, ignore.case = TRUE, perl = TRUE)     ~ "pseudogene",
    grepl(functional_pattern, Gene_Type, ignore.case = TRUE, perl = TRUE) ~ "functional",
    TRUE ~ NA_character_))

# ------------------------- FDR threshold function ---------------------------
fdr_threshold <- function(x, class, fdr = 0.05) {
  null <- x[class == "pseudogene" & !is.na(x)]
  func <- !is.na(x) & class == "functional"
  p <- rep(NA_real_, length(x))
  p[func] <- vapply(x[func], function(v) (1 + sum(null >= v)) / (length(null) + 1),
                    numeric(1))
  q <- rep(NA_real_, length(x))
  q[func] <- p.adjust(p[func], method = "BH")
  sig <- !is.na(q) & q <= fdr
  list(p = p, q = q, sig = sig,
       threshold = if (any(sig)) min(x[sig]) else NA_real_,
       n_null = length(null), n_func = sum(func))
}

# --------------------- run for every SNP database ---------------------------
thr_tab <- list()
for (col in enrich_cols) {
  if (!col %in% names(dat)) { warning("Missing column: ", col); next }
  r <- fdr_threshold(dat[[col]], dat$class, fdr)
  short <- sub("enrichment_intergenic_", "", col)
  dat[[paste0("p_", short)]]    <- r$p
  dat[[paste0("q_", short)]]    <- r$q
  dat[[paste0("high_", short)]] <- r$sig
  thr_tab[[col]] <- data.frame(
    snp_source = short, FDR = fdr,
    n_pseudogene_null = r$n_null, n_functional = r$n_func,
    n_high_enrichment = sum(r$sig), threshold = r$threshold)
}
thr_tab <- bind_rows(thr_tab)
cat("\nFDR thresholds for 'high SNP enrichment':\n"); print(thr_tab, row.names = FALSE)

# ------------------------- summary table by group ---------------------------
high_col <- paste0("high_", sub("enrichment_intergenic_", "", primary_col))

fmt <- function(total, parts) {
  if (length(parts) <= 1) as.character(total)
  else sprintf("%d (%s)", total, paste(parts, collapse = ", "))
}

rows <- lapply(groups, function(g) {
  d <- dat %>% filter(Gene_group %in% g$fam)
  per_fam <- function(cls, sig_only = FALSE) {
    vapply(g$fam, function(f) {
      s <- d %>% filter(Gene_group == f, class == cls)
      if (sig_only) sum(s[[high_col]], na.rm = TRUE) else nrow(s)
    }, numeric(1))
  }
  F <- per_fam("functional"); P <- per_fam("pseudogene"); E <- per_fam("functional", TRUE)
  data.frame(
    `Functional group` = g$name,
    `ncRNA families`   = paste(fam_label[g$fam], collapse = ", "),
    `RNA Polymerase`   = g$pol,
    `F(n)` = fmt(sum(F), F),
    `P(n)` = fmt(sum(P), P),
    `Functional ncRNA with Excess SNPs` = fmt(sum(E), E),
    check.names = FALSE)
})
summary_tab <- bind_rows(rows)
cat("\nSummary table (", primary_col, ", FDR = ", fdr, "):\n", sep = "")
print(summary_tab, row.names = FALSE)

# ------------------------------- output -------------------------------------
write_csv(thr_tab,     file.path(outdir, "snp_enrichment_fdr_thresholds.csv"))
write_csv(summary_tab, file.path(outdir, "snp_enrichment_group_summary.csv"))
write_csv(dat,         file.path(outdir, "gene_data_with_snp_enrichment_fdr.csv"))