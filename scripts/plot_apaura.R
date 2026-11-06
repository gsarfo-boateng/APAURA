#!/usr/bin/env Rscript

# APAURA plotting utility
# Generates publication-ready summary figures from PAS -> PAU -> APA -> RBP outputs.
# The script detects available files/columns and skips plots that are not applicable.

options(stringsAsFactors = FALSE)

required_packages <- c("ggplot2", "dplyr", "tidyr", "readr", "scales", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop(
  "Missing R packages: ", paste(missing_packages, collapse = ", "),
  "\nInstall them with:\n  install.packages(c(", paste(sprintf('"%s"', missing_packages), collapse = ", "), "))", call. = FALSE)

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr); library(scales)
})

usage <- function() cat(paste0(
  "APAURA plotting script\n\nUsage:\n  Rscript plot_apaura.R --input RESULTS_DIR [options]\n\n",
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
  "Outputs:\n  Individual PDF and PNG figures\n  APAURA_figures.pdf (multi-page figure book)\n  figure_manifest.tsv\n"))

parse_args <- function(x) {
  out <- list(input = NULL, rbp = NULL, out = NULL, fdr = .05, delta_pau = .10,
              top = 20L, width = 8, height = 6, help = FALSE)
  value_keys <- c("--input"="input", "--rbp"="rbp", "--out"="out", "--fdr"="fdr",
                  "--delta-pau"="delta_pau", "--top"="top", "--width"="width", "--height"="height")
  i <- 1L
  while (i <= length(x)) {
    key <- x[[i]]
    if (key %in% c("--help", "-h")) { out$help <- TRUE; i <- i + 1L; next }
    if (!key %in% names(value_keys)) stop("Unknown argument: ", key, call. = FALSE)
    if (i == length(x)) stop("Missing value for ", key, call. = FALSE)
    nm <- unname(value_keys[[key]]); val <- x[[i + 1L]]
    out[[nm]] <- if (nm %in% c("fdr", "delta_pau", "width", "height")) as.numeric(val) else if (nm == "top") as.integer(val) else val
    i <- i + 2L
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (args$help) { usage(); quit(save = "no", status = 0) }
if (is.null(args$input)) { usage(); stop("\n--input is required.", call. = FALSE) }
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
  roots <- dirs[!is.na(dirs) & nzchar(dirs)]
  candidates <- unique(c(
    file.path(roots, name),
    file.path(roots, "details", name)))
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
      legend.title = element_text(face = "bold"))
}

theme_set(theme_apaura())

plots <- list()
manifest <- tibble::tibble(
  figure = character(),
  pdf = character(),
  png = character(),
  description = character())

add_plot <- function(name, p, description, width = args$width, height = args$height) {
  pdf_file <- file.path(args$out, paste0(name, ".pdf"))
  png_file <- file.path(args$out, paste0(name, ".png"))
  ggsave(pdf_file, p, width = width, height = height, units = "in")
  ggsave(png_file, p, width = width, height = height, units = "in", dpi = 300)
  plots[[name]] <<- p
  manifest <<- dplyr::bind_rows(
    manifest,
    tibble::tibble(figure = name, pdf = pdf_file, png = png_file, description = description))
  message("[figure] ", name)
}

count_plot <- function(d, col, name, title, ylab, desc, width = args$width, height = args$height) {
  if (!has_cols(d, col)) return(invisible(NULL))
  z <- d %>% transmute(value = .data[[col]]) %>% filter(!is.na(value), nzchar(value)) %>% count(value, sort = TRUE)
  if (!nrow(z)) return(invisible(NULL))
  p <- ggplot(z, aes(x = reorder(pretty_label(value), n), y = n)) + geom_col() + coord_flip() + labs(title = title, x = NULL, y = ylab)
  add_plot(name, p, desc, width, height)
}

xy_plot <- function(d, x, y, name, title, xlab, ylab, desc, shape_col = NULL, alpha = .75,
                    width = args$width, height = args$height, require_nonzero = FALSE) {
  if (!has_cols(d, c(x, y))) return(invisible(NULL))
  z <- d %>% transmute(x = num(.data[[x]]), y = num(.data[[y]]),
                       shape = if (is.null(shape_col) || !shape_col %in% names(d)) "" else as.character(.data[[shape_col]])) %>%
    filter(is.finite(x), is.finite(y))
  if (!nrow(z) || (require_nonzero && !any(z$x > 0 | z$y > 0))) return(invisible(NULL))
  p <- ggplot(z, aes(x, y)) + geom_point(alpha = alpha) + geom_abline(slope = 1, intercept = 0, linetype = 2) + coord_equal()
  if (!is.null(shape_col)) p <- p + aes(shape = pretty_label(shape))
  p <- p + labs(title = title, x = xlab, y = ylab, shape = if (is.null(shape_col)) NULL else "Architecture")
  add_plot(name, p, desc, width, height)
}

violin_plot <- function(d, x, y, name, title, ylab, desc, hline = NULL, rotate = FALSE) {
  if (!has_cols(d, c(x, y))) return(invisible(NULL))
  z <- d %>% transmute(x = pretty_label(.data[[x]]), y = num(.data[[y]])) %>% filter(!is.na(x), nzchar(x), is.finite(y))
  if (!nrow(z)) return(invisible(NULL))
  p <- ggplot(z, aes(x, y)) + geom_violin(trim = TRUE) + geom_boxplot(width = .14, outlier.shape = NA) + labs(title = title, x = NULL, y = ylab)
  if (!is.null(hline)) p <- p + geom_hline(yintercept = hline, linetype = 2)
  if (rotate) p <- p + theme(axis.text.x = element_text(angle = 35, hjust = 1))
  add_plot(name, p, desc)
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

# PAS figures

count_plot(pas, "region", "01_pas_region_counts", "Poly(A) sites by genomic region", "PAS count",
           "Number of retained PAS by annotated genomic region")
count_plot(pas, "architecture_role", "02_pas_architecture_role_counts", "Poly(A) sites by architecture role", "PAS count",
           "PAS counts stratified by TUTR/ALE/other architecture role")

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

# PAU figures

violin_plot(pau, "condition", "pau", "06_pau_distribution", "Poly(A)-site usage by condition", "PAU",
            "Distribution of per-sample poly(A)-site usage by condition")

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
          shape = "Condition")
      add_plot("08_pau_pca", p, "Principal-component analysis of per-sample PAU profiles")
    }
  }
}

# APA site-level figures

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
        shape = "Significant")
    add_plot("09_apa_site_volcano", p, "Site-level delta PAU versus FDR")
  }
}

xy_plot(apa_sites, "control_pau", "case_pau", "10_apa_site_pau_control_vs_case", "Condition-level PAS usage",
        "Control PAU", "Case PAU", "Control versus case condition-level PAU at each PAS", alpha = .65)
violin_plot(apa_sites, "architecture_role", "delta_pau", "11_apa_delta_pau_by_role", "PAS usage shifts by architecture role",
            expression(Delta*"PAU"), "Distribution of delta PAU by PAS architecture role", hline = 0, rotate = TRUE)

# APA gene/group-level figures# APA gene/group-level figures

count_plot(apa, "apa_architecture", "12_apa_architecture_counts", "APA architecture", "Genes / APA groups",
           "Counts of TUTR, ALE, mixed, and none architectures")
count_plot(apa, "apa_direction", "13_apa_direction_counts", "APA change direction", "Genes / APA groups",
           "Counts of shortening, lengthening, ALE switches, IPA gain/loss, and other calls")

if (has_cols(apa, c("apa_architecture", "apa_direction"))) {
  d <- apa %>% filter(!is.na(apa_architecture), !is.na(apa_direction)) %>% count(apa_architecture, apa_direction)
  if (nrow(d)) {
    p <- ggplot(d, aes(x = pretty_label(apa_architecture), y = n, fill = pretty_label(apa_direction))) + geom_col() +
      labs(title = "APA architecture and direction", x = "APA architecture", y = "Genes / APA groups", fill = "Direction")
    add_plot("14_apa_architecture_by_direction", p, "APA direction composition within each architecture", 9, 6)
  }
}

xy_plot(apa, "control_wpi", "case_wpi", "15_apa_wpi_control_vs_case", "Gene-level weighted PAS position",
        "Control WPI", "Case WPI", "Control versus case weighted poly(A)-position index", shape_col = "apa_architecture", alpha = .7)
violin_plot(apa, "apa_architecture", "delta_wpi", "16_apa_delta_wpi_by_architecture", "Global proximal-distal APA shift",
            expression(Delta*"WPI"), "Distribution of gene/group delta WPI by APA architecture", hline = 0)
xy_plot(apa, "control_tutr_wpi", "case_tutr_wpi", "17_tutr_wpi_control_vs_case", "TUTR positional shift",
        "Control TUTR WPI", "Case TUTR WPI", "Within-terminal-exon proximal/distal movement for TUTR-capable genes")
xy_plot(apa, "control_ale_wpi", "case_ale_wpi", "18_ale_wpi_control_vs_case", "Alternative last-exon branch shift",
        "Control ALE WPI", "Case ALE WPI", "Movement among alternative terminal-exon branches")
xy_plot(apa, "control_ipa_usage", "case_ipa_usage", "19_ipa_usage_control_vs_case", "Intronic polyadenylation usage",
        "Control IPA usage", "Case IPA usage", "Condition-level intronic PAS usage", require_nonzero = TRUE)

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

# RBP figures

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
      background_hit_rate = num(background_hit_rate))
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

# Write figure book and manifest

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
  paste0("R_version=", R.version.string))
writeLines(run_info, file.path(args$out, "plot_run_info.txt"))

message("Done. Figures generated: ", length(plots))
