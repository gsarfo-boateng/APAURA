#!/usr/bin/env Rscript

# APAURA plotting utility
# Generates publication-ready summary figures from PAS -> PAU -> APA -> RBP outputs.
# The script detects available files/columns and skips plots that are not applicable.

options(stringsAsFactors = FALSE)

required_packages <- c("ggplot2", "dplyr", "tidyr", "readr", "scales", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them with:\n  install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "), "))",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
})

usage <- function() {
  cat(paste0(
    "APAURA plotting script\n\n",
    "Usage:\n",
    "  Rscript plot_apaura.R --input RESULTS_DIR [options]\n\n",
    "Options:\n",
    "  --input DIR          APAURA results directory [required]\n",
    "  --rbp DIR            RBP results directory; defaults to INPUT/rbp or INPUT\n",
    "  --out DIR            Figure output directory [default: INPUT/figures]\n",
    "  --fdr FLOAT          FDR threshold [default: 0.05]\n",
    "  --delta-pau FLOAT    |delta PAU| threshold [default: 0.10]\n",
    "  --top INT            Number of top genes/RBPs to show [default: 20]\n",
    "  --width FLOAT        Default individual figure width in inches [default: 8]\n",
    "  --height FLOAT       Default individual figure height in inches [default: 6]\n",
    "  --help                Show this message\n\n",
    "Outputs:\n",
    "  Individual PDF and PNG figures\n",
    "  APAURA_figures.pdf (multi-page figure book)\n",
    "  figure_manifest.tsv\n"
  ))
}

parse_args <- function(x) {
  out <- list(
    input = NULL,
    rbp = NULL,
    out = NULL,
    fdr = 0.05,
    delta_pau = 0.10,
    top = 20L,
    width = 8,
    height = 6,
    help = FALSE
  )
  i <- 1L
  while (i <= length(x)) {
    key <- x[[i]]
    if (key == "--help" || key == "-h") {
      out$help <- TRUE
      i <- i + 1L
      next
    }
    if (!grepl("^--", key)) stop("Unknown argument: ", key, call. = FALSE)
    if (i == length(x)) stop("Missing value for ", key, call. = FALSE)
    val <- x[[i + 1L]]
    if (key == "--input") out$input <- val
    else if (key == "--rbp") out$rbp <- val
    else if (key == "--out") out$out <- val
    else if (key == "--fdr") out$fdr <- as.numeric(val)
    else if (key == "--delta-pau") out$delta_pau <- as.numeric(val)
    else if (key == "--top") out$top <- as.integer(val)
    else if (key == "--width") out$width <- as.numeric(val)
    else if (key == "--height") out$height <- as.numeric(val)
    else stop("Unknown argument: ", key, call. = FALSE)
    i <- i + 2L
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (args$help) {
  usage()
  quit(save = "no", status = 0)
}
if (is.null(args$input)) {
  usage()
  stop("\n--input is required.", call. = FALSE)
}
if (!dir.exists(args$input)) stop("Input directory does not exist: ", args$input, call. = FALSE)
if (is.null(args$out)) args$out <- file.path(args$input, "figures")
dir.create(args$out, recursive = TRUE, showWarnings = FALSE)

if (is.null(args$rbp)) {
  rbp_candidate <- file.path(args$input, "rbp")
  args$rbp <- if (dir.exists(rbp_candidate)) rbp_candidate else args$input
}

clean_names <- function(x) {
  x <- sub("^#", "", x)
  x <- tolower(x)
  gsub("[^a-z0-9]+", "_", x)
}

read_table_safe <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  d <- suppressMessages(readr::read_tsv(path, show_col_types = FALSE, progress = FALSE, comment = ""))
  names(d) <- clean_names(names(d))
  d
}

find_file <- function(name, dirs) {
  candidates <- unique(file.path(dirs[!is.na(dirs) & nzchar(dirs)], name))
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) NULL else hit[[1]]
}

read_bed_safe <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  d <- suppressMessages(readr::read_tsv(path, show_col_types = FALSE, progress = FALSE, comment = ""))
  names(d) <- clean_names(names(d))
  d
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

num <- function(x) suppressWarnings(as.numeric(x))
finite_only <- function(x) x[is.finite(x)]

has_cols <- function(d, cols) !is.null(d) && all(cols %in% names(d))

safe_log10 <- function(x) {
  x <- num(x)
  x[!is.finite(x) | x <= 0] <- .Machine$double.xmin
  -log10(x)
}

pretty_label <- function(x) {
  x <- gsub("_", " ", as.character(x))
  x <- gsub("\\bTutr\\b", "TUTR", tools::toTitleCase(tolower(x)))
  x <- gsub("\\bAle\\b", "ALE", x)
  x <- gsub("\\bIpa\\b", "IPA", x)
  x
}

theme_apaura <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

theme_set(theme_apaura())

plots <- list()
manifest <- tibble::tibble(
  figure = character(),
  pdf = character(),
  png = character(),
  description = character()
)

add_plot <- function(name, p, description, width = args$width, height = args$height) {
  pdf_file <- file.path(args$out, paste0(name, ".pdf"))
  png_file <- file.path(args$out, paste0(name, ".png"))
  ggsave(pdf_file, p, width = width, height = height, units = "in")
  ggsave(png_file, p, width = width, height = height, units = "in", dpi = 300)
  plots[[name]] <<- p
  manifest <<- dplyr::bind_rows(
    manifest,
    tibble::tibble(figure = name, pdf = pdf_file, png = png_file, description = description)
  )
  message("[figure] ", name)
}

message("APAURA plotting")
message("  input: ", normalizePath(args$input, mustWork = FALSE))
message("  rbp:   ", normalizePath(args$rbp, mustWork = FALSE))
message("  out:   ", normalizePath(args$out, mustWork = FALSE))

pas_path <- find_file("pas.bed", c(args$input))
pas_tsv_path <- find_file("pas.tsv", c(args$input))
pas_counts_path <- find_file("pas_counts.tsv", c(args$input))
pas_coverage_path <- find_file("pas_coverage.tsv", c(args$input))
pau_path <- find_file("pau.tsv", c(args$input))
apa_sites_path <- find_file("apa_sites.tsv", c(args$input))
apa_path <- find_file("apa.tsv", c(args$input))
apa_arch_path <- find_file("apa_architecture.tsv", c(args$input))
rbp_path <- find_file("rbp.tsv", c(args$rbp, args$input))
rbp_hits_path <- find_file("rbp_hits.tsv", c(args$rbp, args$input))
rbp_profile_path <- find_file("rbp_profile.tsv", c(args$rbp, args$input))

pas <- if (!is.null(pas_path)) read_bed_safe(pas_path) else read_table_safe(pas_tsv_path)
pas_counts <- read_table_safe(pas_counts_path)
pas_coverage <- read_table_safe(pas_coverage_path)
pau <- read_table_safe(pau_path)
apa_sites <- read_table_safe(apa_sites_path)
apa <- read_table_safe(apa_path)
apa_arch <- read_table_safe(apa_arch_path)
rbp <- read_table_safe(rbp_path)
rbp_hits <- read_table_safe(rbp_hits_path)
rbp_profile <- read_table_safe(rbp_profile_path)

# -----------------------------------------------------------------------------
# PAS figures
# -----------------------------------------------------------------------------

if (has_cols(pas, "region")) {
  d <- pas %>% filter(!is.na(region), nzchar(region)) %>% count(region, sort = TRUE)
  if (nrow(d) > 0) {
    p <- ggplot(d, aes(x = reorder(pretty_label(region), n), y = n)) +
      geom_col() + coord_flip() +
      labs(title = "Poly(A) sites by genomic region", x = NULL, y = "PAS count")
    add_plot("01_pas_region_counts", p, "Number of retained PAS by annotated genomic region")
  }
}

if (has_cols(pas, "architecture_role")) {
  d <- pas %>% filter(!is.na(architecture_role), nzchar(architecture_role)) %>% count(architecture_role, sort = TRUE)
  if (nrow(d) > 0) {
    p <- ggplot(d, aes(x = reorder(pretty_label(architecture_role), n), y = n)) +
      geom_col() + coord_flip() +
      labs(title = "Poly(A) sites by architecture role", x = NULL, y = "PAS count")
    add_plot("02_pas_architecture_role_counts", p, "PAS counts stratified by TUTR/ALE/other architecture role")
  }
}

if (has_cols(pas, c("gene_id", "site_id"))) {
  d <- pas %>% filter(!is.na(gene_id), nzchar(gene_id)) %>% count(gene_id, name = "n_pas")
  if (nrow(d) > 0) {
    p <- ggplot(d, aes(x = n_pas)) +
      geom_histogram(binwidth = 1, boundary = 0.5) +
      scale_x_continuous(breaks = pretty_breaks()) +
      labs(title = "Poly(A) sites per gene", x = "Number of PAS", y = "Genes")
    add_plot("03_pas_per_gene", p, "Distribution of retained PAS counts per gene")
  }
}

if (has_cols(pas_counts, c("site_id", "gene_id"))) {
  meta_cols <- intersect(c("site_id", "chrom", "pos0", "strand", "gene_id", "gene_name", "region"), names(pas_counts))
  sample_cols <- setdiff(names(pas_counts), meta_cols)
  if (length(sample_cols) > 0) {
    d <- pas_counts %>%
      pivot_longer(cols = all_of(sample_cols), names_to = "sample", values_to = "count") %>%
      mutate(count = num(count)) %>%
      group_by(sample) %>%
      summarise(total_pas_reads = sum(count, na.rm = TRUE), .groups = "drop")
    p <- ggplot(d, aes(x = reorder(sample, total_pas_reads), y = total_pas_reads)) +
      geom_col() + coord_flip() +
      scale_y_continuous(labels = comma) +
      labs(title = "PAS-supporting reads by sample", x = NULL, y = "Total PAS-supporting reads")
    add_plot("04_pas_reads_by_sample", p, "Total direct/long-read PAS-supporting molecules per sample")
  }
}

if (has_cols(pas_coverage, c("sample", "raw_mean_depth", "fitted_mean_depth"))) {
  d <- pas_coverage %>%
    transmute(sample, raw = num(raw_mean_depth), fitted = num(fitted_mean_depth)) %>%
    filter(is.finite(raw), is.finite(fitted))
  if (nrow(d) > 0) {
    p <- ggplot(d, aes(x = raw, y = fitted)) +
      geom_point(alpha = 0.55) +
      geom_abline(slope = 1, intercept = 0, linetype = 2) +
      facet_wrap(~sample, scales = "free") +
      labs(title = "Short-read PAS coverage fit", x = "Raw mean depth", y = "Fitted mean depth")
    add_plot("05_pas_coverage_fit", p, "Observed versus isotonic-fitted segment coverage for short-read analysis", width = 10, height = 7)
  }
}

# -----------------------------------------------------------------------------
# PAU figures
# -----------------------------------------------------------------------------

if (has_cols(pau, c("sample", "condition", "pau"))) {
  d <- pau %>% mutate(pau = num(pau)) %>% filter(is.finite(pau))
  if (nrow(d) > 0) {
    p <- ggplot(d, aes(x = condition, y = pau)) +
      geom_violin(trim = TRUE) +
      geom_boxplot(width = 0.16, outlier.shape = NA) +
      labs(title = "Poly(A)-site usage by condition", x = NULL, y = "PAU")
    add_plot("06_pau_distribution", p, "Distribution of per-sample poly(A)-site usage by condition")
  }
}

pau_matrix <- NULL
if (has_cols(pau, c("sample", "site_id", "pau"))) {
  wide <- pau %>%
    transmute(sample, site_id, pau = num(pau)) %>%
    filter(is.finite(pau)) %>%
    group_by(site_id, sample) %>%
    summarise(pau = mean(pau), .groups = "drop") %>%
    pivot_wider(names_from = sample, values_from = pau)
  if (nrow(wide) >= 2 && ncol(wide) >= 3) {
    pau_matrix <- as.matrix(wide[, -1, drop = FALSE])
    rownames(pau_matrix) <- wide$site_id
    storage.mode(pau_matrix) <- "numeric"
  }
}

if (!is.null(pau_matrix) && ncol(pau_matrix) >= 2) {
  cm <- suppressWarnings(cor(pau_matrix, use = "pairwise.complete.obs"))
  cd <- as.data.frame(as.table(cm), stringsAsFactors = FALSE)
  names(cd) <- c("sample1", "sample2", "correlation")
  p <- ggplot(cd, aes(x = sample1, y = sample2, fill = correlation)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", correlation)), size = 3) +
    scale_fill_gradient2(midpoint = 0, limits = c(-1, 1), oob = squish) +
    coord_equal() +
    labs(title = "Sample correlation from PAU", x = NULL, y = NULL, fill = "Pearson r") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  add_plot("07_pau_sample_correlation", p, "Pairwise sample correlation calculated from PAU profiles", width = 7, height = 6.5)
}

if (!is.null(pau_matrix) && ncol(pau_matrix) >= 3 && nrow(pau_matrix) >= 3) {
  keep <- apply(pau_matrix, 1, function(z) sum(is.finite(z)) >= 2 && sd(z, na.rm = TRUE) > 0)
  pm <- pau_matrix[keep, , drop = FALSE]
  if (nrow(pm) >= 2) {
    # Replace rare missing values with the site mean for PCA only.
    for (i in seq_len(nrow(pm))) {
      miss <- !is.finite(pm[i, ])
      if (any(miss)) pm[i, miss] <- mean(pm[i, !miss], na.rm = TRUE)
    }
    pc <- prcomp(t(pm), center = TRUE, scale. = FALSE)
    pct <- 100 * pc$sdev^2 / sum(pc$sdev^2)
    pd <- as.data.frame(pc$x[, 1:min(2, ncol(pc$x)), drop = FALSE])
    pd$sample <- rownames(pd)
    if (has_cols(pau, c("sample", "condition"))) {
      sm <- pau %>% distinct(sample, condition)
      pd <- left_join(pd, sm, by = "sample")
    } else {
      pd$condition <- "sample"
    }
    if (all(c("PC1", "PC2") %in% names(pd))) {
      p <- ggplot(pd, aes(x = PC1, y = PC2, shape = condition)) +
        geom_point(size = 3) +
        geom_text(aes(label = sample), nudge_y = 0.03 * diff(range(pd$PC2)), check_overlap = TRUE, size = 3) +
        labs(
          title = "PCA of poly(A)-site usage",
          x = sprintf("PC1 (%.1f%%)", pct[1]),
          y = sprintf("PC2 (%.1f%%)", pct[2]),
          shape = "Condition"
        )
      add_plot("08_pau_pca", p, "Principal-component analysis of per-sample PAU profiles")
    }
  }
}

# -----------------------------------------------------------------------------
# APA site-level figures
# -----------------------------------------------------------------------------

if (has_cols(apa_sites, c("delta_pau", "fdr"))) {
  d <- apa_sites %>%
    mutate(
      delta_pau = num(delta_pau),
      fdr = num(fdr),
      neglog10_fdr = safe_log10(fdr),
      significant = is.finite(fdr) & fdr <= args$fdr & abs(delta_pau) >= args$delta_pau
    ) %>%
    filter(is.finite(delta_pau), is.finite(neglog10_fdr))
  if (nrow(d) > 0) {
    p <- ggplot(d, aes(x = delta_pau, y = neglog10_fdr, shape = significant)) +
      geom_point(alpha = 0.7) +
      geom_vline(xintercept = c(-args$delta_pau, args$delta_pau), linetype = 2) +
      geom_hline(yintercept = -log10(args$fdr), linetype = 2) +
      labs(
        title = "Differential poly(A)-site usage",
        x = expression(Delta*"PAU (case - control)"),
        y = expression(-log[10]("FDR")),
        shape = "Significant"
      )
    add_plot("09_apa_site_volcano", p, "Site-level delta PAU versus FDR")
  }
}

if (has_cols(apa_sites, c("control_pau", "case_pau"))) {
  d <- apa_sites %>%
    mutate(control_pau = num(control_pau), case_pau = num(case_pau)) %>%
    filter(is.finite(control_pau), is.finite(case_pau))
  if (nrow(d) > 0) {
    p <- ggplot(d, aes(x = control_pau, y = case_pau)) +
      geom_point(alpha = 0.65) +
      geom_abline(slope = 1, intercept = 0, linetype = 2) +
      coord_equal() +
      labs(title = "Condition-level PAS usage", x = "Control PAU", y = "Case PAU")
    add_plot("10_apa_site_pau_control_vs_case", p, "Control versus case condition-level PAU at each PAS")
  }
}

if (has_cols(apa_sites, c("delta_pau", "architecture_role"))) {
  d <- apa_sites %>%
    mutate(delta_pau = num(delta_pau), architecture_role = pretty_label(architecture_role)) %>%
    filter(is.finite(delta_pau), !is.na(architecture_role), nzchar(architecture_role))
  if (nrow(d) > 0) {
    p <- ggplot(d, aes(x = architecture_role, y = delta_pau)) +
      geom_violin(trim = TRUE) +
      geom_boxplot(width = 0.14, outlier.shape = NA) +
      geom_hline(yintercept = 0, linetype = 2) +
      labs(title = "PAS usage shifts by architecture role", x = NULL, y = expression(Delta*"PAU")) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1))
    add_plot("11_apa_delta_pau_by_role", p, "Distribution of delta PAU by PAS architecture role")
  }
}

# -----------------------------------------------------------------------------
# APA gene/group-level figures
# -----------------------------------------------------------------------------

if (has_cols(apa, "apa_architecture")) {
  d <- apa %>% filter(!is.na(apa_architecture), nzchar(apa_architecture)) %>% count(apa_architecture, sort = TRUE)
  if (nrow(d) > 0) {
    p <- ggplot(d, aes(x = reorder(pretty_label(apa_architecture), n), y = n)) +
      geom_col() + coord_flip() +
      labs(title = "APA architecture", x = NULL, y = "Genes / APA groups")
    add_plot("12_apa_architecture_counts", p, "Counts of TUTR, ALE, mixed, and none architectures")
  }
}

if (has_cols(apa, "apa_direction")) {
  d <- apa %>% filter(!is.na(apa_direction), nzchar(apa_direction)) %>% count(apa_direction, sort = TRUE)
  if (nrow(d) > 0) {
    p <- ggplot(d, aes(x = reorder(pretty_label(apa_direction), n), y = n)) +
      geom_col() + coord_flip() +
      labs(title = "APA change direction", x = NULL, y = "Genes / APA groups")
    add_plot("13_apa_direction_counts", p, "Counts of shortening, lengthening, ALE switches, IPA gain/loss, and other calls")
  }
}

if (has_cols(apa, c("apa_architecture", "apa_direction"))) {
  d <- apa %>%
    filter(!is.na(apa_architecture), !is.na(apa_direction)) %>%
    count(apa_architecture, apa_direction)
  if (nrow(d) > 0) {
    p <- ggplot(d, aes(x = pretty_label(apa_architecture), y = n, fill = pretty_label(apa_direction))) +
      geom_col() +
      labs(title = "APA architecture and direction", x = "APA architecture", y = "Genes / APA groups", fill = "Direction")
    add_plot("14_apa_architecture_by_direction", p, "APA direction composition within each architecture", width = 9, height = 6)
  }
}

if (has_cols(apa, c("control_wpi", "case_wpi"))) {
  d <- apa %>%
    mutate(control_wpi = num(control_wpi), case_wpi = num(case_wpi)) %>%
    filter(is.finite(control_wpi), is.finite(case_wpi))
  if (nrow(d) > 0) {
    if (!"apa_architecture" %in% names(d)) d$apa_architecture <- "APA"
    p <- ggplot(d, aes(x = control_wpi, y = case_wpi, shape = pretty_label(apa_architecture))) +
      geom_point(alpha = 0.7) +
      geom_abline(slope = 1, intercept = 0, linetype = 2) +
      coord_equal() +
      labs(title = "Gene-level weighted PAS position", x = "Control WPI", y = "Case WPI", shape = "Architecture")
    add_plot("15_apa_wpi_control_vs_case", p, "Control versus case weighted poly(A)-position index")
  }
}

if (has_cols(apa, c("delta_wpi", "apa_architecture"))) {
  d <- apa %>%
    mutate(delta_wpi = num(delta_wpi), apa_architecture = pretty_label(apa_architecture)) %>%
    filter(is.finite(delta_wpi))
  if (nrow(d) > 0) {
    p <- ggplot(d, aes(x = apa_architecture, y = delta_wpi)) +
      geom_violin(trim = TRUE) + geom_boxplot(width = 0.14, outlier.shape = NA) +
      geom_hline(yintercept = 0, linetype = 2) +
      labs(title = "Global proximal-distal APA shift", x = NULL, y = expression(Delta*"WPI"))
    add_plot("16_apa_delta_wpi_by_architecture", p, "Distribution of gene/group delta WPI by APA architecture")
  }
}

if (has_cols(apa, c("control_tutr_wpi", "case_tutr_wpi"))) {
  d <- apa %>%
    mutate(control_tutr_wpi = num(control_tutr_wpi), case_tutr_wpi = num(case_tutr_wpi)) %>%
    filter(is.finite(control_tutr_wpi), is.finite(case_tutr_wpi))
  if (nrow(d) > 0) {
    p <- ggplot(d, aes(x = control_tutr_wpi, y = case_tutr_wpi)) +
      geom_point(alpha = 0.75) + geom_abline(slope = 1, intercept = 0, linetype = 2) +
      coord_equal() + labs(title = "TUTR positional shift", x = "Control TUTR WPI", y = "Case TUTR WPI")
    add_plot("17_tutr_wpi_control_vs_case", p, "Within-terminal-exon proximal/distal movement for TUTR-capable genes")
  }
}

if (has_cols(apa, c("control_ale_wpi", "case_ale_wpi"))) {
  d <- apa %>%
    mutate(control_ale_wpi = num(control_ale_wpi), case_ale_wpi = num(case_ale_wpi)) %>%
    filter(is.finite(control_ale_wpi), is.finite(case_ale_wpi))
  if (nrow(d) > 0) {
    p <- ggplot(d, aes(x = control_ale_wpi, y = case_ale_wpi)) +
      geom_point(alpha = 0.75) + geom_abline(slope = 1, intercept = 0, linetype = 2) +
      coord_equal() + labs(title = "Alternative last-exon branch shift", x = "Control ALE WPI", y = "Case ALE WPI")
    add_plot("18_ale_wpi_control_vs_case", p, "Movement among alternative terminal-exon branches")
  }
}

if (has_cols(apa, c("control_ipa_usage", "case_ipa_usage"))) {
  d <- apa %>%
    mutate(control_ipa_usage = num(control_ipa_usage), case_ipa_usage = num(case_ipa_usage)) %>%
    filter(is.finite(control_ipa_usage), is.finite(case_ipa_usage))
  if (nrow(d) > 0 && any(d$control_ipa_usage > 0 | d$case_ipa_usage > 0)) {
    p <- ggplot(d, aes(x = control_ipa_usage, y = case_ipa_usage)) +
      geom_point(alpha = 0.75) + geom_abline(slope = 1, intercept = 0, linetype = 2) +
      coord_equal() + labs(title = "Intronic polyadenylation usage", x = "Control IPA usage", y = "Case IPA usage")
    add_plot("19_ipa_usage_control_vs_case", p, "Condition-level intronic PAS usage")
  }
}

if (has_cols(apa, c("termination_shift_index", "termination_jsd"))) {
  d <- apa %>%
    mutate(termination_shift_index = num(termination_shift_index), termination_jsd = num(termination_jsd)) %>%
    filter(is.finite(termination_shift_index), is.finite(termination_jsd))
  if (nrow(d) > 0) {
    if (!"apa_direction" %in% names(d)) d$apa_direction <- "APA"
    p <- ggplot(d, aes(x = termination_shift_index, y = termination_jsd, shape = pretty_label(apa_direction))) +
      geom_point(alpha = 0.75) +
      labs(title = "Termination-compartment remodeling", x = "Termination shift index", y = "Termination JSD", shape = "Direction")
    add_plot("20_termination_shift", p, "Total-variation and Jensen-Shannon measures of termination-compartment change")
  }
}

if (has_cols(apa, c("tutr_direction", "ale_direction", "ipa_direction"))) {
  id_col <- if ("gene_name" %in% names(apa)) "gene_name" else if ("gene_id" %in% names(apa)) "gene_id" else NULL
  if (!is.null(id_col)) {
    d0 <- apa
    if ("apa_change" %in% names(d0)) d0 <- d0 %>% filter(apa_change == "yes")
    d <- d0 %>%
      select(all_of(id_col), tutr_direction, ale_direction, ipa_direction) %>%
      pivot_longer(cols = c(tutr_direction, ale_direction, ipa_direction), names_to = "mechanism", values_to = "direction") %>%
      filter(!direction %in% c("not_applicable", "stable", NA_character_))
    if (nrow(d) > 0) {
      genes <- unique(d[[id_col]])
      if (length(genes) > min(args$top, 40L)) genes <- genes[seq_len(min(args$top, 40L))]
      d <- d %>% filter(.data[[id_col]] %in% genes) %>% mutate(mechanism = pretty_label(mechanism), direction = pretty_label(direction))
      p <- ggplot(d, aes(x = mechanism, y = .data[[id_col]], fill = direction)) +
        geom_tile() +
        labs(title = "Mechanism-specific APA calls", x = NULL, y = NULL, fill = "Direction")
      add_plot("21_apa_mechanism_matrix", p, "TUTR, ALE, and IPA direction calls for changed genes", width = 8, height = max(5, 0.28 * length(genes) + 2))
    }
  }
}

# Top gene effects: architecture-aware effect metric.
if (!is.null(apa) && ("gene_id" %in% names(apa) || "gene_name" %in% names(apa))) {
  d <- apa
  d$gene_label <- if ("gene_name" %in% names(d)) ifelse(is.na(d$gene_name) | d$gene_name == "", d$gene_id, d$gene_name) else d$gene_id
  getn <- function(col) if (col %in% names(d)) num(d[[col]]) else rep(NA_real_, nrow(d))
  d$effect <- getn("delta_wpi")
  if ("apa_direction" %in% names(d)) {
    tw <- getn("delta_tutr_wpi")
    aw <- getn("delta_ale_wpi")
    iu <- getn("delta_ipa_usage")
    d$effect[d$apa_direction %in% c("shortening", "lengthening") & is.finite(tw)] <- tw[d$apa_direction %in% c("shortening", "lengthening") & is.finite(tw)]
    d$effect[d$apa_direction %in% c("proximal_switch", "distal_switch", "switch") & is.finite(aw)] <- aw[d$apa_direction %in% c("proximal_switch", "distal_switch", "switch") & is.finite(aw)]
    d$effect[d$apa_direction %in% c("gain", "loss") & is.finite(iu)] <- iu[d$apa_direction %in% c("gain", "loss") & is.finite(iu)]
  }
  fdr_col <- if ("min_site_fdr" %in% names(d)) "min_site_fdr" else if ("min_fdr" %in% names(d)) "min_fdr" else NULL
  if (!is.null(fdr_col)) d$summary_fdr <- num(d[[fdr_col]]) else d$summary_fdr <- NA_real_
  d <- d %>% filter(is.finite(effect))
  if (nrow(d) > 0) {
    d <- d %>% arrange(desc(abs(effect))) %>% slice_head(n = args$top)
    d$gene_label <- factor(d$gene_label, levels = rev(d$gene_label))
    if (!"apa_direction" %in% names(d)) d$apa_direction <- "APA"
    p <- ggplot(d, aes(x = effect, y = gene_label, shape = pretty_label(apa_direction))) +
      geom_segment(aes(x = 0, xend = effect, yend = gene_label), linewidth = 0.5) +
      geom_point(size = 2.5) +
      geom_vline(xintercept = 0, linetype = 2) +
      labs(title = paste0("Top ", min(args$top, nrow(d)), " APA effect sizes"), x = "Architecture-aware signed effect", y = NULL, shape = "Direction")
    add_plot("22_top_apa_effects", p, "Largest architecture-aware gene/group APA effects", width = 9, height = max(6, 0.3 * nrow(d) + 2))
  }
}

# -----------------------------------------------------------------------------
# RBP figures
# -----------------------------------------------------------------------------

if (has_cols(rbp, c("rbp_name", "odds_ratio", "fisher_fdr"))) {
  d <- rbp %>%
    mutate(
      odds_ratio = num(odds_ratio),
      fisher_fdr = num(fisher_fdr),
      neglog10_fdr = safe_log10(fisher_fdr),
      log2_or = log2(pmax(odds_ratio, .Machine$double.xmin))
    ) %>%
    filter(is.finite(log2_or), is.finite(neglog10_fdr))
  if (nrow(d) > 0) {
    dtop <- d %>% arrange(fisher_fdr, desc(odds_ratio)) %>% slice_head(n = args$top)
    dtop$rbp_label <- factor(dtop$rbp_name, levels = rev(unique(dtop$rbp_name)))
    if (!"group" %in% names(dtop)) dtop$group <- "RBP"
    if (!"target_type" %in% names(dtop)) dtop$target_type <- "site"
    p <- ggplot(dtop, aes(x = log2_or, y = rbp_label, size = neglog10_fdr, shape = target_type)) +
      geom_point(alpha = 0.8) +
      geom_vline(xintercept = 0, linetype = 2) +
      facet_wrap(~group, scales = "free_y") +
      labs(title = "RBP motif enrichment", x = expression(log[2]("odds ratio")), y = NULL, size = expression(-log[10]("FDR")), shape = "Target")
    add_plot("23_rbp_enrichment", p, "Top RBP motif enrichments by APA group and target type", width = 10, height = max(6, 0.28 * nrow(dtop) + 2))
  }
}

if (has_cols(rbp, c("rbp_name", "foreground_hit_rate", "background_hit_rate"))) {
  d <- rbp %>%
    mutate(
      foreground_hit_rate = num(foreground_hit_rate),
      background_hit_rate = num(background_hit_rate)
    )
  d$fdr_sort <- if ("fisher_fdr" %in% names(d)) num(d$fisher_fdr) else rep(1, nrow(d))
  d <- d %>%
    filter(is.finite(foreground_hit_rate), is.finite(background_hit_rate)) %>%
    arrange(fdr_sort) %>%
    slice_head(n = args$top)
  if (nrow(d) > 0) {
    if (!"group" %in% names(d)) d$group <- "RBP"
    d$label <- paste(d$rbp_name, d$group, sep = " | ")
    dl <- d %>% select(label, foreground_hit_rate, background_hit_rate) %>%
      pivot_longer(cols = c(foreground_hit_rate, background_hit_rate), names_to = "set", values_to = "hit_rate") %>%
      mutate(set = ifelse(set == "foreground_hit_rate", "foreground", "background"))
    dl$label <- factor(dl$label, levels = rev(unique(d$label)))
    p <- ggplot(dl, aes(x = hit_rate, y = label, shape = set)) +
      geom_point(size = 2.5, position = position_dodge(width = 0.4)) +
      scale_x_continuous(labels = percent_format(accuracy = 1)) +
      labs(title = "RBP motif prevalence", x = "Sites with motif hit", y = NULL, shape = "Set")
    add_plot("24_rbp_foreground_background", p, "Foreground versus matched-background motif-hit prevalence", width = 10, height = max(6, 0.28 * length(unique(dl$label)) + 2))
  }
}

if (has_cols(rbp_profile, c("rbp_name", "relative_start", "hit_count"))) {
  top_rbps <- NULL
  if (has_cols(rbp, c("rbp_name", "fisher_fdr"))) {
    top_rbps <- rbp %>% mutate(fisher_fdr = num(fisher_fdr)) %>% arrange(fisher_fdr) %>% pull(rbp_name) %>% unique() %>% head(min(args$top, 8L))
  }
  d <- rbp_profile %>%
    mutate(relative_start = num(relative_start), hit_count = num(hit_count)) %>%
    filter(is.finite(relative_start), is.finite(hit_count))
  if (!is.null(top_rbps)) d <- d %>% filter(rbp_name %in% top_rbps)
  if (nrow(d) > 0) {
    if (!"group" %in% names(d)) d$group <- "RBP"
    if (!"target_type" %in% names(d)) d$target_type <- "site"
    p <- ggplot(d, aes(x = relative_start, y = hit_count, group = rbp_name, linetype = rbp_name)) +
      geom_line(linewidth = 0.7) +
      geom_vline(xintercept = 0, linetype = 2) +
      facet_grid(group ~ target_type, scales = "free_y") +
      labs(title = "RBP positional motif profile", x = "Position relative to PAS / splice acceptor (nt)", y = "Motif hit count", linetype = "RBP")
    add_plot("25_rbp_positional_profile", p, "Position of enriched RBP motif hits relative to regulated PAS or ALE acceptor", width = 11, height = 7)
  }
}

if (has_cols(rbp, c("rbp_name", "group", "fisher_fdr"))) {
  d <- rbp %>%
    mutate(fisher_fdr = num(fisher_fdr), score = safe_log10(fisher_fdr)) %>%
    filter(is.finite(score)) %>%
    group_by(group, rbp_name) %>%
    summarise(score = max(score, na.rm = TRUE), .groups = "drop")
  if (nrow(d) > 0) {
    keep_rbp <- d %>% group_by(rbp_name) %>% summarise(best = max(score), .groups = "drop") %>% arrange(desc(best)) %>% slice_head(n = args$top) %>% pull(rbp_name)
    dh <- d %>% filter(rbp_name %in% keep_rbp)
    p <- ggplot(dh, aes(x = group, y = rbp_name, fill = score)) +
      geom_tile() +
      labs(title = "RBP enrichment across APA mechanisms", x = NULL, y = NULL, fill = expression(-log[10]("FDR"))) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1))
    add_plot("26_rbp_group_heatmap", p, "Heatmap of RBP motif-enrichment significance across APA mechanism groups", width = 9, height = max(6, 0.28 * length(keep_rbp) + 2))
  }
}

# -----------------------------------------------------------------------------
# Write figure book and manifest
# -----------------------------------------------------------------------------

if (length(plots) > 0) {
  master_pdf <- file.path(args$out, "APAURA_figures.pdf")
  pdf(master_pdf, width = args$width, height = args$height, onefile = TRUE)
  for (nm in names(plots)) print(plots[[nm]])
  dev.off()
  message("[figure book] ", master_pdf)
} else {
  warning("No applicable figures were generated. Check that --input points to a APAURA results directory.")
}

readr::write_tsv(manifest, file.path(args$out, "figure_manifest.tsv"))

run_info <- c(
  paste0("input=", normalizePath(args$input, mustWork = FALSE)),
  paste0("rbp=", normalizePath(args$rbp, mustWork = FALSE)),
  paste0("out=", normalizePath(args$out, mustWork = FALSE)),
  paste0("fdr=", args$fdr),
  paste0("delta_pau=", args$delta_pau),
  paste0("top=", args$top),
  paste0("figures_generated=", length(plots)),
  paste0("R_version=", R.version.string)
)
writeLines(run_info, file.path(args$out, "plot_run_info.txt"))

message("Done. Figures generated: ", length(plots))
