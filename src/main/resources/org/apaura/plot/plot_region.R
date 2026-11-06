#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(patchwork)
  library(rtracklayer)
  library(GenomicRanges)
})

# CLI

parse_cli <- function(x){
  out <- list(); i <- 1
  while(i <= length(x)){
    a <- x[[i]]
    if(!startsWith(a,"--")) stop("Unexpected argument: ",a)
    k <- gsub("-","_",sub("^--","",a))
    if(i==length(x) || startsWith(x[[i+1]],"--")){
      out[[k]] <- TRUE; i <- i+1
    } else {
      out[[k]] <- x[[i+1]]; i <- i+2
    }
  }
  out
}

A <- parse_cli(commandArgs(trailingOnly=TRUE))
getv <- function(x,d=NULL) if(!is.null(A[[x]])) A[[x]] else d
flag <- function(x) isTRUE(A[[x]])
num <- function(x,d) as.numeric(getv(x,d))
int <- function(x,d) as.integer(getv(x,d))
empty <- function(x) is.null(x) || !length(x) || is.na(x[1]) || !nzchar(as.character(x[1]))

mode <- getv("mode")
results_dir <- getv("results_dir")
samples_file <- getv("samples")
gtf_file <- getv("gtf")
control <- getv("control")
case <- getv("case")
gene <- getv("gene")
# BEGIN APAURA_SC_GROUP_ARGUMENT
group_filter <- getv("group", "")
# END APAURA_SC_GROUP_ARGUMENT
region_arg <- getv("region")
out <- getv("out","APAURA_region")

flank <- int("flank",1000)
bin_size <- int("bin_size",10)
min_mapq <- int("min_mapq",10)
merge_method <- tolower(getv("merge","mean"))
library_mode <- tolower(getv("library","unstranded"))

show_replicates <- flag("show_replicates")
show_overview <- flag("show_overview")
verbose <- flag("verbose")
standalone <- flag("standalone")
show_ipa <- !flag("no_ipa")
keep_duplicates <- flag("keep_duplicates")
genomic_orientation <- flag("genomic_orientation")
show_transcripts <- flag("show_transcripts")
max_transcripts <- int("max_transcripts",20)

fdr_cut <- num("fdr",0.05)
apa_effect <- num("apa_effect",0.10)
ipa_effect <- num("ipa_effect",0.10)

format_out <- tolower(getv("format","all"))
plot_width <- num("width",12)
plot_height <- num("height",9)
SAMTOOLS <- getv("samtools","samtools")

if(empty(mode) || !mode %in% c("short-read","run-direct","run-isoquant"))
  stop("--mode must be short-read, run-direct, or run-isoquant")

if(any(vapply(list(results_dir,samples_file,gtf_file,control,case),empty,logical(1))))
  stop("--results-dir, --samples, --gtf, --control and --case are required")

if(empty(gene)==empty(region_arg))
  stop("Specify exactly one of --gene or --region")

if(!merge_method %in% c("mean","median"))
  stop("--merge must be mean or median")

if(!library_mode %in% c("unstranded","forward","reverse"))
  stop("--library must be unstranded, forward or reverse")

if(!format_out %in% c("all","pdf","png","svg"))
  stop("--format must be all, pdf, png or svg")

clean_ens <- function(x) sub("\\..*$","",x)
read_tab <- function(x) read_tsv(x,show_col_types=FALSE,progress=FALSE)

pretty_direction <- function(x){
  x <- as.character(x)
  x[x=="no_directional_shift"] <- "redistribution"
  x[x=="switch"] <- "redistribution"
  x[x=="proximal_switch"] <- "proximal switch"
  x[x=="distal_switch"] <- "distal switch"
  x
}

result_file <- function(name){
  primary <- file.path(results_dir,name)
  if(file.exists(primary)) return(primary)
  detail <- file.path(results_dir,"details",name)
  if(file.exists(detail)) return(detail)
  primary
}

pick_col <- function(df,x,required=TRUE){
  hit <- x[x %in% names(df)][1]
  if(is.na(hit) && required) stop("Expected one of: ",paste(x,collapse=", "))
  hit
}

parse_region <- function(x){
  m <- str_match(x,"^([^:]+):([0-9,]+)-([0-9,]+)$")
  if(any(is.na(m))) stop("Region format must be chr:start-end")
  tibble(chrom=m[2], start=as.integer(gsub(",","",m[3])), end=as.integer(gsub(",","",m[4])))
}

# Samples + annotation + locus

samples <- read_tab(samples_file)

sc <- pick_col(samples,c("sample","sample_id","name"))
cc <- pick_col(samples,c("condition","group"))
names(samples)[names(samples)==sc] <- "sample"
names(samples)[names(samples)==cc] <- "condition"

samples <- samples %>% filter(condition %in% c(control,case))

if(!all(c(control,case) %in% samples$condition))
  stop("Control/case conditions not found in sample sheet")

message("[APAURA plot] loading GTF")
gtf <- rtracklayer::import(gtf_file)

if(!empty(gene)){
  gx <- gtf[gtf$type %in% c("gene","exon")]
  gid <- as.character(gx$gene_id)
  gnm <- as.character(gx$gene_name)

  hit <- (!is.na(gnm) & gnm==gene) | (!is.na(gid) & clean_ens(gid)==clean_ens(gene))

  gm <- gx[hit]
  if(!length(gm)) stop("Gene not found in GTF: ",gene)

  chr <- unique(as.character(seqnames(gm)))
  if(length(chr)!=1) stop("Gene maps to multiple chromosomes/contigs")

  region <- tibble(chrom=chr, start=max(1,min(start(gm))-flank), end=max(end(gm))+flank)

  gn <- unique(na.omit(gnm[hit]))
  gene_label <- if(length(gn)) gn[1] else gene

  ss <- unique(as.character(strand(gm)))
  ss <- ss[ss %in% c("+","-")]
  strand0 <- if(length(ss)) ss[1] else "+"
} else {
  region <- parse_region(region_arg)
  gene_label <- region_arg
  strand0 <- NA_character_
}

REGION <- paste0(region$chrom,":",region$start,"-",region$end)

# Never reverse genomic x coordinates.
# Minus-strand short-read coverage is only mirrored vertically.
MIRROR_COV <- mode=="short-read" && !empty(gene) && identical(strand0,"-") && !genomic_orientation

message(
  "[APAURA plot] mode=",mode, " target=",gene_label,
  ifelse(is.na(strand0),"",paste0(" strand=",strand0)),
  " library=",library_mode,
  ifelse(MIRROR_COV," coverage=mirrored"," coverage=upright"))

xscale <- function(labels=TRUE){
  scale_x_continuous(limits=c(region$start,region$end),
    labels=if(labels) scales::comma else NULL,
    breaks=scales::breaks_pretty(n=5), expand=c(0,0))
}

empty_panel <- function(title){
  ggplot()+xscale(FALSE)+theme_void()+ggtitle(title)
}

# APA sites

apa_file <- file.path(
  results_dir,
  "apa_sites.tsv"
)

# BEGIN APAURA_SC_SITE_FILE
sc_site_file <- file.path(
  results_dir,
  "apaura_site_usage.tsv"
)

sc_plot_region <- (
  !file.exists(apa_file)
  && file.exists(sc_site_file)
)

if (sc_plot_region) {
  apa_file <- sc_site_file
}
# END APAURA_SC_SITE_FILE
if(!file.exists(apa_file)) stop("Missing: ",apa_file)

apa <- read_tab(apa_file)

# BEGIN APAURA_SC_SITE_ADAPTER

# SC event-level metadata.
#
# IMPORTANT:
# apaura_site_usage.tsv contains PAS usage/effect information.
# Its PAS rows do NOT currently provide a site-specific inferential FDR.
#
# For replicated SC analyses, statistical support therefore refers
# to the gene-level/event-level APA test, not to each individual PAS.
sc_event_architecture <- NA_character_
sc_event_classification <- NA_character_
sc_event_effect <- NA_real_
sc_event_q <- NA_real_
sc_analysis_mode <- NA_character_

sc_n_control <- NA_integer_
sc_n_case <- NA_integer_

sc_unreplicated <- FALSE
sc_event_supported <- FALSE

if (sc_plot_region) {

  required_sc_columns <- c(
    "group",
    "gene_id",
    "gene_symbol",
    "chrom",
    "strand",
    "site_id",
    "position",
    "control_pau",
    "treated_pau",
    "delta_pau"
  )

  missing_sc_columns <- setdiff(
    required_sc_columns,
    names(apa)
  )

  if (length(missing_sc_columns) > 0) {

    stop(
      paste0(
        "apaura_site_usage.tsv is missing required columns: ",
        paste(
          missing_sc_columns,
          collapse = ", "
        )
      )
    )
  }

  # Exact group selection

  available_groups <- unique(
    as.character(
      apa$group
    )
  )

  available_groups <- available_groups[
    !is.na(available_groups)
    & nzchar(available_groups)
  ]

  if (nzchar(group_filter)) {

    apa <- apa[
      as.character(apa$group) == group_filter,
      ,
      drop = FALSE
    ]

    if (nrow(apa) == 0) {

      stop(
        paste0(
          "No rows in apaura_site_usage.tsv for --group ",
          group_filter
        )
      )
    }

  } else if (length(available_groups) > 1) {

    stop(
      paste0(
        "Single-cell APA results contain multiple groups. ",
        "Supply --group. Available groups: ",
        paste(
          available_groups,
          collapse = ", "
        )
      )
    )
  }

  # Canonical APAURA plotting columns

  apa$gene_name <- as.character(
    apa$gene_symbol
  )

  apa$pas_position_1based <- suppressWarnings(
    as.integer(
      apa$position
    )
  )

  apa$control_pau <- suppressWarnings(
    as.numeric(
      apa$control_pau
    )
  )

  apa$case_pau <- suppressWarnings(
    as.numeric(
      apa$treated_pau
    )
  )

  apa$delta_pau <- suppressWarnings(
    as.numeric(
      apa$delta_pau
    )
  )

  # Biological-replicate counts from site table

  first_finite_int <- function(x) {

    z <- suppressWarnings(
      as.integer(x)
    )

    z <- z[
      is.finite(z)
    ]

    if (!length(z))
      return(NA_integer_)

    z[1]
  }

  if ("n_control_samples" %in% names(apa)) {

    sc_n_control <- first_finite_int(
      apa$n_control_samples
    )
  }

  if ("n_treated_samples" %in% names(apa)) {

    sc_n_case <- first_finite_int(
      apa$n_treated_samples
    )
  }

  # Read gene-level/event-level SC result

  sc_event_candidates <- c(
    file.path(
      results_dir,
      "all_results.with_architecture.tsv"
    ),
    file.path(
      results_dir,
      "unified_results.with_architecture.tsv"
    ),
    file.path(
      results_dir,
      "observed_unreplicated_events.with_architecture.tsv"
    )
  )

  sc_event_candidates <- sc_event_candidates[
    file.exists(sc_event_candidates)
  ]

  if (length(sc_event_candidates)) {

    sc_event_file <- sc_event_candidates[1]

    sc_event <- read_tab(
      sc_event_file
    )

    if (
      "group" %in% names(sc_event)
      && nzchar(group_filter)
    ) {

      sc_event <- sc_event[
        as.character(sc_event$group)
        == group_filter,
        ,
        drop = FALSE
      ]
    }

    # Support either gene_symbol or gene_name.
    if (!"gene_symbol" %in% names(sc_event)) {

      if ("gene_name" %in% names(sc_event)) {

        sc_event$gene_symbol <- as.character(
          sc_event$gene_name
        )

      } else {

        sc_event$gene_symbol <- ""
      }
    }

    if (!"gene_id" %in% names(sc_event)) {

      sc_event$gene_id <- ""
    }

    if (!empty(gene)) {

      sc_event <- sc_event %>%
        filter(
          gene_symbol == gene
          | clean_ens(gene_id)
            == clean_ens(gene)
        )
    }

    if (nrow(sc_event)) {

      ev <- sc_event[1, , drop = FALSE]

      if ("apa_architecture" %in% names(ev)) {

        sc_event_architecture <- toupper(
          as.character(
            ev$apa_architecture[1]
          )
        )
      }

      if ("classification" %in% names(ev)) {

        sc_event_classification <- as.character(
          ev$classification[1]
        )
      }

      effect_col <- c(
        "delta_utr_index",
        "delta_tutr_wpi",
        "delta_wpi",
        "delta_dui"
      )

      effect_col <- effect_col[
        effect_col %in% names(ev)
      ]

      if (length(effect_col)) {

        sc_event_effect <- suppressWarnings(
          as.numeric(
            ev[[effect_col[1]]][1]
          )
        )
      }

      if ("qvalue" %in% names(ev)) {

        sc_event_q <- suppressWarnings(
          as.numeric(
            ev$qvalue[1]
          )
        )
      }

      if ("analysis_mode" %in% names(ev)) {

        sc_analysis_mode <- as.character(
          ev$analysis_mode[1]
        )
      }

      if ("n_control_samples" %in% names(ev)) {

        z <- suppressWarnings(
          as.integer(
            ev$n_control_samples[1]
          )
        )

        if (is.finite(z))
          sc_n_control <- z
      }

      if ("n_treated_samples" %in% names(ev)) {

        z <- suppressWarnings(
          as.integer(
            ev$n_treated_samples[1]
          )
        )

        if (is.finite(z))
          sc_n_case <- z
      }
    }
  }

  # Determine replication state
  #
  # Cells are never treated as biological replicates.
  # These counts refer to biological samples/pseudobulks.

  sc_unreplicated <- (
    is.na(sc_n_control)
    || is.na(sc_n_case)
    || sc_n_control < 2
    || sc_n_case < 2
    || (
      !is.na(sc_analysis_mode)
      && grepl(
        "UNREPLICATED",
        sc_analysis_mode,
        ignore.case = TRUE
      )
    )
  )

  # Event-level statistical support
  #
  # Only valid when biological replication exists.
  # This is NOT a per-PAS FDR.

  sc_event_supported <- (
    !sc_unreplicated
    && is.finite(sc_event_q)
    && sc_event_q <= fdr_cut
    && (
      !is.finite(sc_event_effect)
      || abs(sc_event_effect) >= apa_effect
    )
  )

  # PAS rows themselves do not currently carry a site-specific FDR.
  # Keep the canonical site FDR explicitly NA.
  apa$fdr <- NA_real_

  if (sc_unreplicated) {

    apa$status <- "observed_unreplicated"

  } else {

    apa$status <- "replicated_sc"
  }

  message(
    "[APAURA plot] single-cell result group: ",
    ifelse(
      nzchar(group_filter),
      group_filter,
      available_groups[1]
    )
  )

  message(
    "[APAURA plot] single-cell replication: control=",
    sc_n_control,
    " case=",
    sc_n_case,
    " mode=",
    ifelse(
      sc_unreplicated,
      "unreplicated descriptive",
      "replicated inferential"
    )
  )

  if (!is.na(sc_event_architecture)) {

    message(
      "[APAURA plot] single-cell event: architecture=",
      sc_event_architecture,
      " classification=",
      sc_event_classification,
      " deltaUTR=",
      signif(
        sc_event_effect,
        4
      ),
      if (
        sc_unreplicated
        || !is.finite(sc_event_q)
      ) {
        ""
      } else {
        paste0(
          " event_q=",
          signif(
            sc_event_q,
            4
          )
        )
      }
    )
  }
}

# END APAURA_SC_SITE_ADAPTER
if(!"pas_position_1based" %in% names(apa)){
  pc <- pick_col(apa,c("position_1based","pas_position","position","pos1"),FALSE)

  if(!is.na(pc)){
    apa$pas_position_1based <- as.numeric(apa[[pc]])
  } else if("pas_start0" %in% names(apa)){
    apa$pas_position_1based <- as.numeric(apa$pas_start0)+1
  } else if("start0" %in% names(apa)){
    apa$pas_position_1based <- as.numeric(apa$start0)+1
  } else {
    stop("Cannot identify PAS coordinate in apa_sites.tsv")
  }
}

if(!"delta_pau" %in% names(apa)){
  dc <- pick_col(apa,c("delta_PAU","deltaPAU","delta_usage"),FALSE)
  if(is.na(dc)) stop("Cannot identify delta PAU column")
  apa$delta_pau <- as.numeric(apa[[dc]])
}

if(!"fdr" %in% names(apa)){
  fc <- pick_col(apa,c("FDR","padj","qvalue"),FALSE)
  apa$fdr <- if(is.na(fc)) NA_real_ else as.numeric(apa[[fc]])
}

if(!"status" %in% names(apa)) apa$status <- "ok"
if(!"gene_id" %in% names(apa)) apa$gene_id <- ""
if(!"gene_name" %in% names(apa)) apa$gene_name <- ""

if(!"chrom" %in% names(apa)){
  cc <- pick_col(apa,c("chr","seqname","contig"))
  apa$chrom <- apa[[cc]]
}

apa <- apa %>%
  filter(chrom==region$chrom, between(pas_position_1based,region$start,region$end))

if(!empty(gene)){
  apa <- apa %>%
    filter(gene_name==gene | clean_ens(gene_id)==clean_ens(gene))
}

apa <- apa %>%
  mutate(
    direction=case_when(
      delta_pau > 0 ~ "Higher PAU in case",
      delta_pau < 0 ~ "Lower PAU in case",
      TRUE ~ "No change"
    )
  )

if (sc_plot_region) {

  # SC PAS rows do not have a per-site inferential FDR.
  #
  # Unreplicated:
  #   show all observed PAS changes at full opacity.
  #
  # Replicated:
  #   opacity reflects gene/event-level FDR support.
  #   It must not be interpreted as PAS-specific significance.
  apa$significant <- if (
    sc_unreplicated
  ) {
    TRUE
  } else {
    sc_event_supported
  }

} else {

  # Preserve existing bulk/non-SC behavior exactly.
  apa <- apa %>%
    mutate(
      significant=
        status=="ok"
        & !is.na(fdr)
        & fdr<=fdr_cut
        & abs(delta_pau)>=apa_effect
    )
}

if (
  sc_plot_region
  && "site_id" %in% names(apa)
) {

  sid <- as.character(
    apa$site_id
  )

  apa$site_label <- ifelse(
    grepl(
      "_PAS[0-9]+$",
      sid
    ),
    sub(
      "^.*_(PAS[0-9]+)$",
      "\\1",
      sid
    ),
    sid
  )

} else {

  apa$site_label <- ""
}

message("[APAURA plot] terminal PAS in region: ",nrow(apa))

# IPA

ipa <- tibble()
ipa_file <- result_file("ipa_sites.tsv")

if(show_ipa && file.exists(ipa_file)){
  ipa <- read_tab(ipa_file)

  if(!"representative_pas_1based" %in% names(ipa))
    stop("ipa_sites.tsv lacks representative_pas_1based")

  if(!"status" %in% names(ipa)) ipa$status <- "ok"
  if(!"gene_id" %in% names(ipa)) ipa$gene_id <- ""
  if(!"gene_name" %in% names(ipa)) ipa$gene_name <- ""

  ipa <- ipa %>%
    filter(chrom==region$chrom, between(representative_pas_1based,region$start,region$end))

  if(!empty(gene)){
    ipa <- ipa %>%
      filter(gene_name==gene | clean_ens(gene_id)==clean_ens(gene))
  }

  ipa <- ipa %>%
    mutate(significant=status=="ok" & !is.na(fdr) & fdr<=fdr_cut & abs(delta_ipa_usage)>=ipa_effect,
      direction=case_when(delta_ipa_usage>0 ~ "IPA gain in case", delta_ipa_usage<0 ~ "IPA loss in case", TRUE ~ "No change"))
}

message("[APAURA plot] IPA events in region: ",nrow(ipa))

# Colours + PAS guides

COND_COL <- c(setNames("#4C78A8",control), setNames("#E45756",case))

EVENT_COL <- c(
  "Higher PAU in case"="#D55E00",
  "Lower PAU in case"="#0072B2",
  "IPA gain in case"="#CC6677",
  "IPA loss in case"="#4477AA",
  "No change"="#999999")

apa_guides <- apa %>% distinct(pas_position_1based)

ipa_guides <- if(nrow(ipa)){
  ipa %>% distinct(representative_pas_1based)
} else {
  tibble()
}

# Evidence

if(mode=="short-read"){

  pc <- pick_col(samples,c("bam","file","path","bam/file/path"))
  names(samples)[names(samples)==pc] <- "bam"

  sample_dir <- dirname(normalizePath(samples_file,mustWork=TRUE))

  samples <- samples %>%
    mutate(bam=ifelse(grepl("^/",bam),bam,file.path(sample_dir,bam)))

  if(Sys.which(SAMTOOLS)=="" && !file.exists(SAMTOOLS))
    stop("samtools not found: ",SAMTOOLS)

  # unmapped + secondary + QC-fail + supplementary
  BASE_EXCLUDE <- if(keep_duplicates) 2820L else 3844L
  FORWARD_EXCLUDE <- bitwOr(BASE_EXCLUDE,16L)

  # Paired-end transcript strand assignment

  strand_specs <- function(s){

    if(library_mode=="unstranded") return(NULL)

    if(library_mode=="reverse"){
      if(s=="+"){
        return(tibble(req=c(81L,129L),exc=c(BASE_EXCLUDE,FORWARD_EXCLUDE)))
      } else {
        return(tibble(req=c(65L,145L),exc=c(FORWARD_EXCLUDE,BASE_EXCLUDE)))
      }
    }

    if(s=="+"){
      tibble(req=c(65L,145L),exc=c(FORWARD_EXCLUDE,BASE_EXCLUDE))
    } else {
      tibble(req=c(81L,129L),exc=c(BASE_EXCLUDE,FORWARD_EXCLUDE))
    }
  }

  sam_count <- function(bam,req=NULL,exc=BASE_EXCLUDE){
    a <- c("view","-c","-q",min_mapq)
    if(!is.null(req)) a <- c(a,"-f",req)
    a <- c(a,"-F",exc,bam)

    x <- system2(SAMTOOLS,a,stdout=TRUE)
    n <- suppressWarnings(as.numeric(tail(x,1)))

    if(!length(n) || !is.finite(n))
      stop("Could not count alignments: ",bam)

    n
  }

  # Same denominator as strand-specific bigWig:
  # total assigned = plus transcript + minus transcript
  # Persistent normalization cache

  CACHE_FILE <- file.path(results_dir,"plot_normalization.tsv")

  empty_cache <- function(){
    tibble(
      sample=character(),
      bam=character(),
      library=character(),
      min_mapq=integer(),
      keep_duplicates=logical(),
      bam_size=double(),
      bam_mtime=double(),
      total_assigned=double()
    )
  }

  norm_cache <- if(file.exists(CACHE_FILE)){
    suppressWarnings(read_tab(CACHE_FILE))
  } else {
    empty_cache()
  }

  get_norm_total <- function(sample,bam){

    bp <- normalizePath(bam,mustWork=TRUE)
    fi <- file.info(bp)
    size0 <- as.numeric(fi$size)
    mtime0 <- as.numeric(fi$mtime)

    hit <- which(
      norm_cache$sample==sample &
      norm_cache$bam==bp &
      norm_cache$library==library_mode &
      norm_cache$min_mapq==min_mapq &
      norm_cache$keep_duplicates==keep_duplicates &
      norm_cache$bam_size==size0 &
      abs(norm_cache$bam_mtime-mtime0)<1
    )

    if(length(hit)){
      total <- norm_cache$total_assigned[hit[1]]
      if(verbose)
        message("[APAURA plot] normalization cache hit: ",sample,
                " | total=",round(total))
      return(total)
    }

    if(verbose) message("[APAURA plot] normalization cache miss: ",sample, " | scanning BAM once")

    total <- if(library_mode=="unstranded"){
      sam_count(bam)
    } else {
      sam_count(bam,req=1L,exc=BASE_EXCLUDE)
    }

    if(!is.finite(total) || total<=0)
      stop("No usable alignments for normalization: ",bam)

    row <- tibble(
      sample=sample,
      bam=bp,
      library=library_mode,
      min_mapq=min_mapq,
      keep_duplicates=keep_duplicates, bam_size=size0, bam_mtime=mtime0, total_assigned=total)

    same <- norm_cache$sample==sample &
            norm_cache$bam==bp &
            norm_cache$library==library_mode &
            norm_cache$min_mapq==min_mapq &
            norm_cache$keep_duplicates==keep_duplicates

    norm_cache <<- bind_rows(norm_cache[!same,],row)

    tmp <- paste0(CACHE_FILE,".tmp.",Sys.getpid())
    write_tsv(norm_cache,tmp)

    if(!file.rename(tmp,CACHE_FILE)){
      unlink(tmp)
      warning("Could not update normalization cache: ",CACHE_FILE)
    }

    total
  }

  # Build regional BAM for selected transcript strand

  make_region_bam <- function(bam,out_bam){

    if(library_mode=="unstranded"){
      s <- system2(SAMTOOLS, c("view","-b", "-q",as.character(min_mapq), "-F",as.character(BASE_EXCLUDE), bam,REGION ), stdout=out_bam)

      if(s!=0) stop("samtools view failed: ",bam)
      return(invisible(TRUE))
    }

    specs <- if(!empty(gene)){
      strand_specs(strand0)
    } else {
      bind_rows(strand_specs("+"),strand_specs("-"))
    }

    parts <- paste0(out_bam,".part",seq_len(nrow(specs)),".bam")
    on.exit(unlink(parts),add=TRUE)

    for(i in seq_len(nrow(specs))){
      s <- system2(SAMTOOLS, c("view","-b", "-q",as.character(min_mapq), "-f",as.character(specs$req[i]),  "-F",as.character(specs$exc[i]), bam,REGION), stdout=parts[i])

      if(s!=0) stop("strand-specific samtools view failed: ",bam)
    }

    s <- system2(SAMTOOLS, c("merge","-f",out_bam,parts))

    if(s!=0) stop("samtools merge failed: ",bam)
  }

  # Coverage

  get_cov <- function(sample,condition,bam){

    total <- get_norm_total(sample,bam)

    message("[APAURA plot] ",sample," | ",condition, " | library=",library_mode,
      ifelse(!empty(gene) && library_mode!="unstranded",
        paste0(" | transcript_strand=",strand0), "" ),
      " | normalization_total=",round(total))

    tmp_bam <- tempfile(fileext=".bam")
    tmp_dep <- tempfile(fileext=".depth")
    tmp_bai <- paste0(tmp_bam,".bai")

    on.exit(unlink(c(tmp_bam,tmp_bai,tmp_dep)),add=TRUE)

    make_region_bam(bam,tmp_bam)

    if(system2(SAMTOOLS,c("index",tmp_bam))!=0)
      stop("samtools index failed: ",sample)

    if(system2(SAMTOOLS,
      c("depth","-a","-r",REGION,tmp_bam), stdout=tmp_dep)!=0) stop("samtools depth failed: ",sample)

    x <- read_tsv(tmp_dep, col_names=c("chrom","pos","depth"), col_types="cid", show_col_types=FALSE, progress=FALSE)

    if(!nrow(x))
      stop("No coverage returned for ",sample)

    x %>%
      mutate(sample=sample, condition=condition, norm_cov=depth/total*1e6,
        bin=region$start+((pos-region$start)%/%bin_size)*bin_size) %>%
      group_by(sample,condition,bin) %>%
      summarise(norm_cov=mean(norm_cov,na.rm=TRUE), .groups="drop") %>%
      mutate(pos=bin+(bin_size-1)/2)
  }

  cov_df <- pmap_dfr(samples %>% select(sample,condition,bam), get_cov)

  rep_n <- samples %>% count(condition,name="n")
  cov_df <- cov_df %>%
    left_join(rep_n,by="condition") %>%
    mutate(label=paste0(condition," (n=",n,")"))

  merged_df <- cov_df %>%
    group_by(condition,pos) %>%
    summarise(norm_cov=if(merge_method=="mean"){
        mean(norm_cov,na.rm=TRUE)
      } else {
        median(norm_cov,na.rm=TRUE)
      },
      .groups="drop") %>%
    left_join(rep_n,by="condition") %>%
    mutate(label=paste0(condition," (n=",n,")"))

  mult <- if(MIRROR_COV) -1 else 1

  cov_df <- cov_df %>% mutate(plot_cov=norm_cov*mult)
  merged_df <- merged_df %>% mutate(plot_cov=norm_cov*mult)

  order0 <- c(control,case)

  label_levels <- paste0(order0," (n=", rep_n$n[match(order0,rep_n$condition)],
    ")")

  cov_df <- cov_df %>%
    mutate(condition=factor(condition,levels=order0),
      label=factor(label,levels=label_levels))

  merged_df <- merged_df %>%
    mutate(condition=factor(condition,levels=order0),
      label=factor(label,levels=label_levels))

  # Coverage plot

  # Whole-gene overview + 3' terminal-region coverage zoom

  tested_terminal_pas <- numeric()

  if(nrow(apa) && "pas_position_1based" %in% names(apa)){
    tested_terminal_pas <- suppressWarnings(
      as.numeric(apa$pas_position_1based)
    )
    tested_terminal_pas <- tested_terminal_pas[
      is.finite(tested_terminal_pas)
    ]
  }

  terminal_exon_df <- tibble()

  if(!is.null(gene) && length(gene) && !is.na(gene[1]) &&
     nzchar(as.character(gene[1]))){

    gx_terminal <- gtf[
      gtf$type=="exon" &
      as.character(seqnames(gtf))==region$chrom
    ]

    if(length(gx_terminal)){
      gid0 <- as.character(gx_terminal$gene_id)
      gnm0 <- as.character(gx_terminal$gene_name)

      gx_terminal <- gx_terminal[
        (!is.na(gnm0) & gnm0==gene) |
        (!is.na(gid0) & clean_ens(gid0)==clean_ens(gene))
      ]
    }

    if(length(gx_terminal) &&
       "transcript_id" %in% names(mcols(gx_terminal))){

      tx_exons0 <- tibble(
        transcript=as.character(gx_terminal$transcript_id),
        start=start(gx_terminal),
        end=end(gx_terminal)
      ) %>%
        filter(!is.na(transcript), transcript!="")

      if(nrow(tx_exons0)){
        if(identical(strand0,"-")){
          terminal_exon_df <- tx_exons0 %>%
            group_by(transcript) %>%
            slice_min(order_by=start,n=1,with_ties=FALSE) %>%
            ungroup()
        } else {
          terminal_exon_df <- tx_exons0 %>%
            group_by(transcript) %>%
            slice_max(order_by=end,n=1,with_ties=FALSE) %>%
            ungroup()
        }
      }
    }
  }

  # Prefer the exact terminal PAS represented in the ΔPAU panel.
  # Pure-ALE loci may not have such PAS, so fall back to the
  # annotated alternative last exons.
  if(length(tested_terminal_pas)){
    feature_min <- min(tested_terminal_pas)
    feature_max <- max(tested_terminal_pas)
    terminal_zoom_basis <- "tested terminal PAS"
  } else if(nrow(terminal_exon_df)){
    feature_min <- min(terminal_exon_df$start)
    feature_max <- max(terminal_exon_df$end)
    terminal_zoom_basis <- "annotated alternative last exons"
  } else {
    feature_min <- region$start
    feature_max <- region$end
    terminal_zoom_basis <- "requested region"
  }

  feature_span <- max(1,feature_max-feature_min)

  # Give a little sequence context without letting a tiny terminal
  # exon become a single-pixel plot.
  terminal_flank <- as.integer(
    max(500,min(1000,round(feature_span*.20)))
  )

  terminal_zoom_start <- max(
    region$start,
    floor(feature_min-terminal_flank)
  )

  terminal_zoom_end <- min(
    region$end,
    ceiling(feature_max+terminal_flank)
  )

  min_zoom_width <- 2000

  if((terminal_zoom_end-terminal_zoom_start) < min_zoom_width){
    centre0 <- (terminal_zoom_start+terminal_zoom_end)/2

    terminal_zoom_start <- max(
      region$start,
      floor(centre0-min_zoom_width/2)
    )

    terminal_zoom_end <- min(
      region$end,
      ceiling(centre0+min_zoom_width/2)
    )
  }

  terminal_merged <- merged_df %>%
    filter(
      pos>=terminal_zoom_start,
      pos<=terminal_zoom_end
    )

  terminal_reps <- cov_df %>%
    filter(
      pos>=terminal_zoom_start,
      pos<=terminal_zoom_end
    )

  if(!nrow(terminal_merged)){
    terminal_merged <- merged_df
    terminal_reps <- cov_df
    terminal_zoom_start <- region$start
    terminal_zoom_end <- region$end
    terminal_zoom_basis <- "requested region"
  }

  coverage_panel <- function(
    merged_data,
    replicate_data,
    xmin0,
    xmax0,
    title0,
    subtitle0,
    show_x=FALSE,
    show_zoom=FALSE
  ){

    q <- ggplot(
      merged_data,
      aes(pos,plot_cov)
    )

    if(show_zoom &&
       is.finite(terminal_zoom_start) &&
       is.finite(terminal_zoom_end)){

      q <- q +
        annotate(
          "rect",
          xmin=terminal_zoom_start,
          xmax=terminal_zoom_end,
          ymin=-Inf,
          ymax=Inf,
          fill="grey45",
          alpha=.08
        )
    }

    q <- q +
      geom_hline(
        yintercept=0,
        linewidth=.30,
        colour="grey35"
      ) +
      geom_area(
        aes(fill=condition),
        alpha=.16,
        colour=NA
      )

    if(show_replicates && nrow(replicate_data)){
      q <- q +
        geom_line(
          data=replicate_data,
          aes(
            pos,
            plot_cov,
            colour=condition,
            group=interaction(condition,sample)
          ),
          linewidth=.18,
          alpha=.16,
          inherit.aes=FALSE
        )
    }

    q <- q +
      geom_line(
        aes(colour=condition),
        linewidth=.65
      )

    # Mark tested PAS only on the zoomed panel so coverage changes
    # can be compared directly with the ΔPAU panel.
    if(!show_zoom && length(tested_terminal_pas)){
      pas_lines0 <- tibble(
        pas_position_1based=tested_terminal_pas
      )

      q <- q +
        geom_vline(
          data=pas_lines0,
          aes(xintercept=pas_position_1based),
          linewidth=.30,
          linetype=3,
          colour="grey68"
        )
    }

    q <- q +
      facet_grid(label~.,switch="y") +
      scale_colour_manual(values=COND_COL) +
      scale_fill_manual(values=COND_COL) +
      scale_y_continuous(
        labels=function(x)
          scales::number(abs(x),accuracy=.1)
      ) +
      scale_x_continuous(
        limits=c(xmin0,xmax0),
        labels=scales::comma,
        expand=c(0,0)
      ) +
      theme_classic(base_size=10) +
      theme(
        legend.position="none",
        legend.title=element_blank(),
        strip.background=element_blank(),
        strip.placement="outside",
        strip.text.y.left=element_text(
          angle=0,
          face="bold"
        ),
        axis.text.x=if(show_x)
          element_text(size=8)
        else
          element_blank(),
        axis.ticks.x=if(show_x)
          element_line()
        else
          element_blank(),
        plot.title=element_text(
          face="bold",
          size=11
        ),
        plot.subtitle=element_text(
          size=8.5,
          colour="grey30"
        )
      ) +
      labs(
        x=if(show_x)
          paste0(region$chrom," genomic position")
        else
          NULL,
        y=paste0(
          ifelse(
            merge_method=="mean",
            "Mean",
            "Median"
          ),
          " CPM coverage"
        ),
        title=title0,
        subtitle=subtitle0,
        colour=NULL,
        fill=NULL
      )

    q
  }

  three_prime_text <- if(
    identical(strand0,"-")
  ){
    "transcript 3′ direction: left"
  } else if(
    identical(strand0,"+")
  ){
    "transcript 3′ direction: right"
  } else {
    "transcript strand unavailable"
  }

  p_overview <- coverage_panel(
    merged_df,
    cov_df,
    region$start,
    region$end,
    paste0(
      gene_label,
      ifelse(
        is.na(strand0),
        "",
        paste0(" | ",strand0," strand")
      )
    ),
    paste0(
      "Whole-gene coverage overview • shaded area = terminal zoom • ",
      control," vs ",case
    ),
    show_x=FALSE,
    show_zoom=TRUE
  )

  p_evidence <- coverage_panel(
    terminal_merged,
    terminal_reps,
    terminal_zoom_start,
    terminal_zoom_end,
    "3′ terminal-region coverage",
    paste0(
      "Zoom based on ",
      terminal_zoom_basis,
      " • shared CPM scale across conditions • ",
      three_prime_text
    ),
    show_x=TRUE,
    show_zoom=FALSE
  )

# Direct RNA / IsoQuant

} else {

  if(!all(c("control_pau","case_pau") %in% names(apa)))
    stop("apa_sites.tsv must contain control_pau and case_pau")

  evidence_df <- bind_rows(
    apa %>% transmute(condition=control,pos=pas_position_1based,usage=control_pau),
    apa %>% transmute(condition=case,pos=pas_position_1based,usage=case_pau))

  evidence_df$condition <- factor(evidence_df$condition, levels=c(control,case))

  p_evidence <- ggplot(evidence_df,aes(pos,usage,colour=condition))+
    geom_segment(aes(xend=pos,y=0,yend=usage),linewidth=.65)+
    geom_point(size=3)+
    facet_grid(condition~.,switch="y")+
    scale_colour_manual(values=COND_COL)+
    xscale(FALSE)+
    theme_classic(base_size=11)+
    theme(legend.position="none",
      strip.background=element_blank(),
      strip.placement="outside",
      strip.text.y.left=element_text(angle=0,face="bold"),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())+
    labs(x=NULL, y="Mean PAU", title=paste0(gene_label,
        ifelse(is.na(strand0),"",paste0(" | ",strand0," strand"))),
      subtitle=ifelse(mode=="run-direct", "Direct-RNA 3′-end/PAS usage", "IsoQuant corrected-read PAS usage"))
}

# APA panel

summary_architecture <- NA_character_
summary_tutr_call <- "not significant"

apa_panel_title <- "Terminal poly(A) site usage"
apa_panel_subtitle <- paste0(
  "ΔPAU = ",
  case,
  " − ",
  control
)

if (sc_plot_region) {

  if (!is.na(sc_event_architecture)) {

    summary_architecture <-
      sc_event_architecture
  }

  cls <- toupper(
    ifelse(
      is.na(sc_event_classification),
      "",
      sc_event_classification
    )
  )

  sc_direction_label <- case_when(
    grepl("LENGTHEN", cls) ~ "lengthening",
    grepl("SHORTEN", cls) ~ "shortening",
    grepl("COMPLEX", cls) ~ "complex remodeling",
    TRUE ~ "remodeling"
  )

  arch_label <- ifelse(
    is.na(summary_architecture)
    || !nzchar(summary_architecture),
    "APA",
    summary_architecture
  )

  if (sc_unreplicated) {

    apa_panel_title <- paste0(
      "Terminal poly(A) site usage — ",
      arch_label,
      " • observed ",
      sc_direction_label
    )

  } else {

    support_label <- if (
      sc_event_supported
    ) {
      "FDR-supported "
    } else {
      "not FDR-supported "
    }

    apa_panel_title <- paste0(
      "Terminal poly(A) site usage — ",
      arch_label,
      " • ",
      support_label,
      sc_direction_label
    )
  }

  effect_label <- if (
    is.finite(sc_event_effect)
  ) {
    sprintf(
      "ΔUTR index = %+.3f",
      sc_event_effect
    )
  } else {
    "ΔUTR index = NA"
  }

  max_dpau <- suppressWarnings(
    max(
      abs(apa$delta_pau),
      na.rm = TRUE
    )
  )

  dpau_label <- if (
    is.finite(max_dpau)
  ) {
    sprintf(
      "max |ΔPAU| = %.3f",
      max_dpau
    )
  } else {
    "max |ΔPAU| = NA"
  }

  group_label <- if (
    nzchar(group_filter)
  ) {
    group_filter
  } else {
    "single-cell group"
  }

  if (sc_unreplicated) {

    apa_panel_subtitle <- paste0(
      effect_label,
      " • ",
      dpau_label,
      " • ΔPAU = ",
      case,
      " − ",
      control,
      "
",
      group_label,
      " • unreplicated descriptive comparison"
    )

  } else {

    q_label <- if (
      is.finite(sc_event_q)
    ) {
      format(
        sc_event_q,
        digits = 3,
        scientific = TRUE
      )
    } else {
      "NA"
    }

    apa_panel_subtitle <- paste0(
      effect_label,
      " • event q = ",
      q_label,
      " • ",
      dpau_label,
      "
",
      group_label,
      " • n=",
      sc_n_control,
      " vs ",
      sc_n_case,
      " biological samples • ΔPAU = ",
      case,
      " − ",
      control
    )
  }

} else if (!empty(gene)) {

  # Existing bulk/non-SC summary behavior.
  summary_file0 <- result_file(
    "apa.tsv"
  )

  if (file.exists(summary_file0)) {

    summary0 <- read_tab(
      summary_file0
    )

    if (!"gene_id" %in% names(summary0))
      summary0$gene_id <- ""

    if (!"gene_name" %in% names(summary0))
      summary0$gene_name <- ""

    summary0 <- summary0 %>%
      filter(
        gene_name == gene
        | clean_ens(gene_id)
          == clean_ens(gene)
      )

    if (
      nrow(summary0)
      && "apa_architecture" %in% names(summary0)
    ) {

      summary_architecture <- as.character(
        summary0$apa_architecture[1]
      )
    }

    if (
      nrow(summary0)
      && all(
        c(
          "tutr_change",
          "tutr_direction"
        ) %in% names(summary0)
      )
    ) {

      if (
        as.character(
          summary0$tutr_change[1]
        ) == "yes"
      ) {

        summary_tutr_call <- pretty_direction(
          summary0$tutr_direction[1]
        )

      } else if (
        !is.na(summary_architecture)
        && !summary_architecture
           %in% c(
             "TUTR",
             "mixed"
           )
      ) {

        summary_tutr_call <- "not applicable"
      }
    }
  }

  apa_panel_title <- paste0(
    "Terminal poly(A) site usage — TUTR: ",
    summary_tutr_call
  )
}
if(nrow(apa)){

  p_apa <- ggplot(
    apa,
    aes(
      pas_position_1based,
      delta_pau
    )
  )+
    geom_hline(
      yintercept=0,
      linewidth=.35,
      colour="grey25"
    )+
    geom_hline(
      yintercept=c(
        -apa_effect,
        apa_effect
      ),
      linetype=3,
      linewidth=.3,
      colour="grey70"
    )+
    geom_segment(
      aes(
        xend=pas_position_1based,
        y=0,
        yend=delta_pau,
        colour=direction,
        alpha=significant
      ),
      linewidth=.7
    )+
    geom_point(
      aes(
        colour=direction,
        alpha=significant
      ),
      size=3
    )+
    geom_text(
      aes(
        label=site_label
      ),
      vjust=ifelse(
        apa$delta_pau >= 0,
        -0.7,
        1.5
      ),
      size=2.6,
      check_overlap=TRUE,
      show.legend=FALSE
    )+
    scale_colour_manual(
      values=EVENT_COL
    )+
    scale_alpha_manual(
      values=c(
        `FALSE`=.42,
        `TRUE`=1
      ),
      guide="none"
    )+
    xscale(FALSE)+
    theme_classic(
      base_size=10
    )+
    theme(
      legend.position="right",
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank()
    )+
    labs(
      x=NULL,
      y=expression(
        Delta*"PAU"
      ),
      colour=NULL,
      title=apa_panel_title,
      subtitle=apa_panel_subtitle
    )

} else {
  p_apa <- if(identical(summary_architecture,"ALE")){
    empty_panel("TUTR usage: not applicable")
  } else {
    empty_panel("TUTR usage: no tested sites in region")
  }
}

# IPA panel

if(show_ipa){

  if(nrow(ipa)){

    p_ipa <- ggplot(ipa, aes(representative_pas_1based,delta_ipa_usage))+
      geom_hline(yintercept=0,linewidth=.35,colour="grey25")+
      geom_hline(yintercept=c(-ipa_effect,ipa_effect), linetype=3, linewidth=.3, colour="grey70")+
      geom_segment(aes(xend=representative_pas_1based, y=0, yend=delta_ipa_usage, colour=direction, alpha=significant), linewidth=.7)+
      geom_point(aes(colour=direction,alpha=significant), size=3,shape=17)+
      scale_colour_manual(values=EVENT_COL)+
      scale_alpha_manual(values=c(`FALSE`=.42,`TRUE`=1), guide="none")+
      xscale(FALSE)+
      theme_classic(base_size=10)+
      theme(legend.position="right", axis.text.x=element_blank(), axis.ticks.x=element_blank())+
      labs(x=NULL, y=expression(Delta*"IPA"), colour=NULL, title="Intronic polyadenylation — IPA usage",
        subtitle=paste0("ΔIPA = ",case," − ",control))

  } else {
    p_ipa <- empty_panel("Intronic polyadenylation: no tested event in region")
  }
}

# Gene model
# APA architecture / gene model

p_ale <- NULL
show_ale_panel <- FALSE

exons <- gtf[
  gtf$type=="exon" &
  as.character(seqnames(gtf))==region$chrom &
  start(gtf)<=region$end &
  end(gtf)>=region$start
]

if(!empty(gene) && length(exons)){
  gid <- as.character(exons$gene_id)
  gnm <- as.character(exons$gene_name)

  exons <- exons[
    gnm==gene |
    clean_ens(gid)==clean_ens(gene)
  ]
}

gene_height <- 1.4

apa_summary <- tibble()
apa_summary_file <- result_file("apa.tsv")

if(!empty(gene) && file.exists(apa_summary_file)){
  apa_summary <- read_tab(apa_summary_file)

  if(!"gene_id" %in% names(apa_summary)) apa_summary$gene_id <- ""
  if(!"gene_name" %in% names(apa_summary)) apa_summary$gene_name <- ""

  apa_summary <- apa_summary %>%
    filter(
      gene_name==gene |
      clean_ens(gene_id)==clean_ens(gene)
    )
}

pick_architecture <- function(x){
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]

  for(v in c("mixed","ALE","TUTR","none"))
    if(v %in% x) return(v)

  if(length(x)) x[1] else NA_character_
}

reported_architecture <- NA_character_
architecture_detail <- ""

if(nrow(apa_summary)){
  if("apa_architecture" %in% names(apa_summary))
    reported_architecture <- pick_architecture(apa_summary$apa_architecture)

  if("mechanism_directions" %in% names(apa_summary)){
    z <- unique(as.character(apa_summary$mechanism_directions))
    z <- z[!is.na(z) & nzchar(z) & z!="not_applicable"]
    if(length(z)){
      architecture_detail <- z[1]
      architecture_detail <- str_replace_all(architecture_detail,"no_directional_shift","redistribution")
      architecture_detail <- str_replace_all(architecture_detail,"proximal_switch","proximal switch")
      architecture_detail <- str_replace_all(architecture_detail,"distal_switch","distal switch")
      architecture_detail <- str_replace_all(architecture_detail,"(?<=ALE=)switch","redistribution")
    }
  }

  if(!nzchar(architecture_detail)){
    pieces <- character()

    if("tutr_direction" %in% names(apa_summary)){
      z <- unique(as.character(apa_summary$tutr_direction))
      z <- z[!is.na(z) & !z %in% c("","stable","not_applicable")]
      if(length(z)) pieces <- c(pieces,paste0("TUTR=",pretty_direction(z[1])))
    }

    if("ale_direction" %in% names(apa_summary)){
      z <- unique(as.character(apa_summary$ale_direction))
      z <- z[!is.na(z) & !z %in% c("","stable","not_applicable")]
      if(length(z)) pieces <- c(pieces,paste0("ALE=",pretty_direction(z[1])))
    }

    if(length(pieces)) architecture_detail <- paste(pieces,collapse="; ")
  }

  if(!nzchar(architecture_detail) &&
     "apa_direction" %in% names(apa_summary)){
    z <- unique(as.character(apa_summary$apa_direction))
    z <- z[!is.na(z) & !z %in% c("","stable","not_applicable")]
    if(length(z)) architecture_detail <- pretty_direction(z[1])
  }
}

cluster_anchors <- function(x,tol=5L){
  if(!length(x)) return(integer())

  o <- order(x)
  out <- integer(length(x))
  cluster <- 1L
  out[o[1]] <- cluster
  last <- x[o[1]]

  if(length(o)>1){
    for(j in o[-1]){
      if(abs(x[j]-last)>tol) cluster <- cluster+1L
      out[j] <- cluster
      last <- x[j]
    }
  }

  out
}

if(length(exons) && show_transcripts &&
   "transcript_id" %in% names(mcols(exons))){

  tx <- tibble(
    transcript=as.character(exons$transcript_id),
    start=start(exons),
    end=end(exons)
  ) %>%
    filter(!is.na(transcript),transcript!="")

  keep <- tx %>%
    count(transcript,name="n") %>%
    arrange(desc(n)) %>%
    slice_head(n=max_transcripts) %>%
    pull(transcript)

  tx <- tx %>%
    filter(transcript %in% keep) %>%
    mutate(
      transcript=factor(transcript,levels=rev(keep)),
      y=as.numeric(transcript)
    )

  lines <- tx %>%
    group_by(transcript,y) %>%
    summarise(
      start=min(start),
      end=max(end),
      .groups="drop"
    )

  p_gene <- ggplot()+
    geom_segment(
      data=lines,
      aes(x=start,xend=end,y=y,yend=y),
      linewidth=.3,
      colour="grey50"
    )+
    geom_rect(
      data=tx,
      aes(
        xmin=start,xmax=end,
        ymin=y-.28,ymax=y+.28
      ),
      fill="grey15"
    )+
    scale_y_continuous(
      breaks=seq_along(rev(keep)),
      labels=rev(keep)
    )+
    xscale(TRUE)+
    theme_classic(base_size=9)+
    theme(
      axis.title.y=element_blank(),
      axis.ticks.y=element_blank(),
      axis.line.y=element_blank()
    )+
    labs(
      x=paste0(region$chrom," genomic position"),
      title="Annotated transcripts"
    )

  gene_height <- min(5,max(1.5,length(keep)*.22+.5))

} else if(length(exons) && !empty(gene) &&
          "transcript_id" %in% names(mcols(exons))){

  tx <- tibble(
    transcript=as.character(exons$transcript_id),
    start=start(exons),
    end=end(exons)
  ) %>%
    filter(!is.na(transcript),transcript!="")

  terminal_tx <- if(identical(strand0,"-")){
    tx %>%
      group_by(transcript) %>%
      slice_min(order_by=start,n=1,with_ties=FALSE) %>%
      ungroup()
  } else {
    tx %>%
      group_by(transcript) %>%
      slice_max(order_by=end,n=1,with_ties=FALSE) %>%
      ungroup()
  }

  cds <- gtf[
    gtf$type=="CDS" &
    as.character(seqnames(gtf))==region$chrom
  ]

  if(length(cds)){
    gid <- as.character(cds$gene_id)
    gnm <- as.character(cds$gene_name)

    cds <- cds[
      gnm==gene |
      clean_ens(gid)==clean_ens(gene)
    ]
  }

  if(length(cds) &&
     "transcript_id" %in% names(mcols(cds))){

    cds_df <- tibble(
      transcript=as.character(cds$transcript_id),
      start=start(cds),
      end=end(cds)
    ) %>%
      filter(!is.na(transcript),transcript!="") %>%
      group_by(transcript) %>%
      summarise(
        cds_start=min(start),
        cds_end=max(end),
        .groups="drop"
      )

    terminal_tx <- terminal_tx %>%
      left_join(cds_df,by="transcript")

  } else {
    terminal_tx$cds_start <- NA_real_
    terminal_tx$cds_end <- NA_real_
  }

  if(identical(strand0,"-")){
    terminal_tx <- terminal_tx %>%
      mutate(
        anchor=end,
        utr_start=start,
        utr_end=ifelse(
          is.na(cds_start),
          end,
          pmin(end,cds_start-1)
        ),
        has_utr=utr_start<=utr_end
      )
  } else {
    terminal_tx <- terminal_tx %>%
      mutate(
        anchor=start,
        utr_start=ifelse(
          is.na(cds_end),
          start,
          pmax(start,cds_end+1)
        ),
        utr_end=end,
        has_utr=utr_start<=utr_end
      )
  }

  terminal_tx$branch <- cluster_anchors(terminal_tx$anchor,5L)

  branches <- terminal_tx %>%
    group_by(branch) %>%
    summarise(
      anchor=round(median(anchor)),
      exon_start=min(start),
      exon_end=max(end),
      utr_start=if(any(has_utr)) min(utr_start[has_utr]) else NA_real_,
      utr_end=if(any(has_utr)) max(utr_end[has_utr]) else NA_real_,
      has_utr=any(has_utr),
      n_transcripts=n_distinct(transcript),
      .groups="drop"
    )

  # Prefer APAURA's engine-derived ALE branch catalog when available. It is
  # built by the same TerminalFamilyReconciler used for apa_architecture.tsv,
  # so the plot shows every structural terminal branch (including ALE-only
  # branches that have no tandem PAS differential row).
  ale_ann <- tibble()
  ale_ann_file <- result_file("ale_annotation.tsv")
  if(file.exists(ale_ann_file)){
    ale_ann <- read_tab(ale_ann_file)
    if(!"gene_id" %in% names(ale_ann)) ale_ann$gene_id <- ""
    if(!"gene_name" %in% names(ale_ann)) ale_ann$gene_name <- ""
    ale_ann <- ale_ann %>%
      filter(
        gene_name==gene |
        clean_ens(gene_id)==clean_ens(gene)
      )
  }

  ale_diff <- tibble()
  ale_diff_file <- result_file("ale_branches.tsv")
  if(file.exists(ale_diff_file)){
    ale_diff <- read_tab(ale_diff_file)
    if(!"gene_id" %in% names(ale_diff)) ale_diff$gene_id <- ""
    if(!"gene_name" %in% names(ale_diff)) ale_diff$gene_name <- ""
    ale_diff <- ale_diff %>%
      filter(
        gene_name==gene |
        clean_ens(gene_id)==clean_ens(gene)
      )
  }

  annotation_pas <- tibble()

  if(nrow(ale_ann)){
    branches <- map_dfr(seq_len(nrow(ale_ann)),function(i){
      a <- ale_ann[i,]
      candidates <- terminal_tx %>%
        filter(abs(anchor-a$anchor0)<=5)

      if(nrow(candidates)){
        candidates <- candidates %>%
          arrange(abs(anchor-a$anchor0))

        exon_start0 <- min(candidates$start)
        exon_end0 <- max(candidates$end)
        has_utr0 <- any(candidates$has_utr)
        utr_start0 <- if(has_utr0) min(candidates$utr_start[candidates$has_utr]) else NA_real_
        utr_end0 <- if(has_utr0) max(candidates$utr_end[candidates$has_utr]) else NA_real_
      } else {
        pos <- suppressWarnings(as.numeric(
          unlist(strsplit(as.character(a$annotated_pas_positions_1based),",",fixed=TRUE))
        ))
        pos <- pos[is.finite(pos)]
        exon_start0 <- min(c(a$anchor0,pos),na.rm=TRUE)
        exon_end0 <- max(c(a$anchor0,pos),na.rm=TRUE)
        has_utr0 <- FALSE
        utr_start0 <- NA_real_
        utr_end0 <- NA_real_
      }

      tibble(
        branch=as.integer(a$branch_rank),
        branch_id=as.character(a$branch_id),
        branch_rank=as.integer(a$branch_rank),
        anchor=as.numeric(a$anchor0),
        exon_start=exon_start0,
        exon_end=exon_end0,
        utr_start=utr_start0,
        utr_end=utr_end0,
        has_utr=has_utr0,
        n_transcripts=NA_integer_,
        n_pas=as.integer(a$n_annotated_pas),
        annotated_pas_positions_1based=as.character(a$annotated_pas_positions_1based)
      )
    })

    annotation_pas <- map_dfr(seq_len(nrow(branches)),function(i){
      z <- suppressWarnings(as.numeric(
        unlist(strsplit(branches$annotated_pas_positions_1based[i],",",fixed=TRUE))
      ))
      z <- z[is.finite(z)]
      if(!length(z)) return(tibble())
      tibble(
        branch=branches$branch[i],
        branch_id=branches$branch_id[i],
        pas_position_1based=z
      )
    })

  } else {
    # Backward-compatible fallback for older result directories.
    pas_map0 <- map_dfr(seq_len(nrow(branches)),function(i){
      b <- branches[i,]
      hit <- apa %>%
        filter(
          pas_position_1based >= b$exon_start-24,
          pas_position_1based <= b$exon_end+24
        )
      if(!nrow(hit)) return(tibble())
      hit %>%
        transmute(
          branch=b$branch,
          pas_position_1based
        )
    })

    pas_counts <- if(nrow(pas_map0)){
      pas_map0 %>%
        distinct(branch,pas_position_1based) %>%
        count(branch,name="n_pas")
    } else {
      tibble(branch=integer(),n_pas=integer())
    }

    branches <- branches %>%
      left_join(pas_counts,by="branch") %>%
      mutate(
        n_pas=coalesce(n_pas,0L),
        branch_id=paste0("branch_",branch),
        branch_rank=branch
      )
    annotation_pas <- if(nrow(pas_map0)){
      pas_map0 %>%
        left_join(branches %>% select(branch,branch_id),by="branch")
    } else {
      tibble(
        branch=integer(),
        branch_id=character(),
        pas_position_1based=double()
      )
    }
  }

  # Map differential tandem-PAS results to the structural branches by genomic
  # PAS proximity. This keeps structural PAS (grey) separate from tested ΔPAU.
  pas_map <- map_dfr(seq_len(nrow(branches)),function(i){
    b <- branches[i,]

    if(nrow(annotation_pas)){
      branch_positions <- annotation_pas %>%
        filter(branch==b$branch) %>%
        pull(pas_position_1based)
      hit <- apa %>%
        filter(vapply(
          pas_position_1based,
          function(x) any(abs(x-branch_positions)<=24),
          logical(1)
        ))
    } else {
      hit <- apa %>%
        filter(
          pas_position_1based >= b$exon_start-24,
          pas_position_1based <= b$exon_end+24
        )
    }

    if(!nrow(hit)) return(tibble())

    hit %>%
      transmute(
        branch=b$branch,
        branch_id=b$branch_id,
        pas_position_1based,
        delta_pau,
        direction,
        significant
      )
  })

  has_tutr <- any(branches$n_pas>=2)
  has_ale <- nrow(branches)>=2

  inferred_architecture <- case_when(
    has_tutr && has_ale ~ "mixed",
    has_tutr ~ "TUTR",
    has_ale ~ "ALE",
    TRUE ~ "none"
  )

  architecture0 <- if(
    !is.na(reported_architecture) &&
    nzchar(reported_architecture)
  ){
    reported_architecture
  } else {
    inferred_architecture
  }

  branch_view <- branches

  if("branch_rank" %in% names(branch_view)){
    branch_view <- branch_view %>% arrange(branch_rank)
  } else if(identical(strand0,"-")){
    branch_view <- branch_view %>% arrange(desc(anchor))
  } else {
    branch_view <- branch_view %>% arrange(anchor)
  }

  if(nrow(branch_view)>max_transcripts){
    message(
      "[APAURA plot] terminal branches truncated: ",
      nrow(branch_view)," -> ",max_transcripts
    )
    branch_view <- branch_view %>% slice_head(n=max_transcripts)
  }

  if(nrow(ale_diff)){
    branch_view <- branch_view %>%
      left_join(
        ale_diff %>%
          select(branch_id,delta_usage,fdr,status),
        by="branch_id"
      )
  } else {
    branch_view <- branch_view %>%
      mutate(delta_usage=NA_real_,fdr=NA_real_,status=NA_character_)
  }

  branch_view <- branch_view %>%
    mutate(
      display_rank=row_number(),
      y=display_rank+1,
      branch_label=paste0(
        "Last exon ",display_rank,
        " • ",n_pas," PAS"
      )
    )

  # Dedicated ALE panel. Branch/BU remain internal engine terminology; the
  # user-facing plot describes these as alternative last exons and Δusage.
  last_exon_view <- branch_view %>%
    filter(is.finite(delta_usage)) %>%
    mutate(
      usage_direction=case_when(
        delta_usage>0 ~ "Higher last-exon usage in case",
        delta_usage<0 ~ "Lower last-exon usage in case",
        TRUE ~ "No change"
      ),
      significant=coalesce(status=="ok" & is.finite(fdr) & fdr<=fdr_cut,FALSE),
      last_exon_label=paste0("LE",display_rank)
    )

  if(nrow(last_exon_view) && architecture0 %in% c("ALE","mixed")){
    ale_call <- "not significant"
    if(nrow(apa_summary) && all(c("ale_change","ale_direction") %in% names(apa_summary))){
      if(as.character(apa_summary$ale_change[1])=="yes")
        ale_call <- pretty_direction(apa_summary$ale_direction[1])
    }

    LAST_EXON_COL <- c(
      "Higher last-exon usage in case"="#D55E00",
      "Lower last-exon usage in case"="#0072B2",
      "No change"="#999999"
    )

    last_exon_labels <- last_exon_view %>%
      arrange(desc(abs(delta_usage))) %>%
      slice_head(n=4) %>%
      filter(abs(delta_usage)>=0.04 | significant)

    p_ale <- ggplot(last_exon_view,aes(anchor,delta_usage))+
      geom_hline(yintercept=0,linewidth=.35,colour="grey25")+
      geom_segment(
        aes(xend=anchor,y=0,yend=delta_usage,colour=usage_direction,alpha=significant),
        linewidth=.75
      )+
      geom_point(aes(colour=usage_direction,alpha=significant),size=3)+
      geom_text(
        data=last_exon_labels,
        aes(label=last_exon_label),
        vjust=ifelse(last_exon_labels$delta_usage>=0,-.8,1.4),
        size=2.7,
        colour="grey25",
        check_overlap=TRUE
      )+
      scale_colour_manual(values=LAST_EXON_COL)+
      scale_alpha_manual(values=c(`FALSE`=.48,`TRUE`=1),guide="none")+
      xscale(FALSE)+
      theme_classic(base_size=10)+
      theme(
        legend.position="right",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
      )+
      labs(
        x=NULL,
        y=expression(Delta*"usage"),
        colour=NULL,
        title=paste0("Alternative last-exon usage — ALE: ",ale_call),
        subtitle=paste0("Δusage = ",case," − ",control)
      )

    show_ale_panel <- TRUE
  }

  pas_view <- if(nrow(pas_map)){
    pas_map %>%
      inner_join(
        branch_view %>% select(branch,y),
        by="branch"
      )
  } else {
    tibble()
  }

  annotation_pas_view <- if(nrow(annotation_pas)){
    annotation_pas %>%
      inner_join(
        branch_view %>% select(branch,y),
        by="branch"
      )
  } else {
    tibble()
  }

  er <- GenomicRanges::reduce(granges(exons))
  edf <- tibble(
    start=start(er),
    end=end(er)
  )

  gs <- min(edf$start)
  ge <- max(edf$end)

  arrow_from <- ifelse(
    !is.na(strand0) && strand0=="-",
    ge,gs
  )
  arrow_to <- ifelse(
    !is.na(strand0) && strand0=="-",
    gs,ge
  )

  y_breaks <- c(1,branch_view$y)
  y_labels <- c("Gene overview",branch_view$branch_label)

  # APAURA standalone: annotation-only last-exon labels.
  # Do not report "0 PAS" when no APAURA PAS results were supplied.
  if(standalone){
    y_labels <- c(
      "Gene overview",
      paste0("Last exon ", seq_len(nrow(branch_view)))
    )
  }

  subtitle0 <- "Blue = terminal region; open circles = annotated PAS; filled circles = tested PAS"

  p_gene <- ggplot()+
    geom_segment(
      aes(
        x=arrow_from,xend=arrow_to,
        y=1,yend=1
      ),
      linewidth=.35,
      colour="grey45",
      arrow=grid::arrow(
        length=grid::unit(.07,"inches"),
        type="closed"
      )
    )+
    geom_rect(
      data=edf,
      aes(
        xmin=start,xmax=end,
        ymin=.82,ymax=1.18
      ),
      fill="grey20"
    )+
    geom_segment(
      data=branch_view,
      aes(
        x=exon_start,xend=exon_end,
        y=y,yend=y
      ),
      linewidth=.35,
      colour="grey55"
    )+
    geom_rect(
      data=branch_view,
      aes(
        xmin=exon_start,xmax=exon_end,
        ymin=y-.22,ymax=y+.22
      ),
      fill="grey25"
    )

  if(any(branch_view$has_utr)){
    p_gene <- p_gene+
      geom_rect(
        data=branch_view %>% filter(has_utr),
        aes(
          xmin=utr_start,xmax=utr_end,
          ymin=y-.27,ymax=y+.27
        ),
        fill="#4C78A8",
        alpha=.72
      )
  }

  if(nrow(annotation_pas_view)){
    p_gene <- p_gene+
      geom_point(
        data=annotation_pas_view,
        aes(x=pas_position_1based,y=y),
        inherit.aes=FALSE,
        size=2.1,
        shape=21,
        fill="white",
        colour="grey40",
        stroke=.55
      )
  }

  if(nrow(pas_view)){
    p_gene <- p_gene+
      geom_segment(
        data=pas_view,
        aes(
          x=pas_position_1based,
          xend=pas_position_1based,
          y=y-.31,
          yend=y+.31,
          colour=direction
        ),
        linewidth=.55
      )+
      geom_point(
        data=pas_view,
        aes(
          x=pas_position_1based,
          y=y,
          colour=direction,
          alpha=significant
        ),
        size=2.5
      )+
      scale_colour_manual(
        values=EVENT_COL,
        guide="none"
      )+
      scale_alpha_manual(
        values=c(`FALSE`=.45,`TRUE`=1),
        guide="none"
      )
  }

  p_gene <- p_gene+
    scale_y_continuous(
      breaks=y_breaks,
      labels=y_labels,
      expand=expansion(mult=c(.08,.08))
    )+
    xscale(TRUE)+
    theme_classic(base_size=9)+
    theme(
      axis.title.y=element_blank(),
      axis.ticks.y=element_blank(),
      axis.line.y=element_blank(),
      panel.grid=element_blank()
    )+
    labs(
      x=paste0(region$chrom," genomic position"),
      title=paste0("3′ terminal APA architecture | ",architecture0),
      subtitle=subtitle0
    )

  gene_height <- min(
    5,
    max(2,1.1+nrow(branch_view)*.48)
  )

} else if(length(exons)){

  er <- GenomicRanges::reduce(granges(exons))

  edf <- tibble(
    start=start(er),
    end=end(er)
  )

  gs <- min(edf$start)
  ge <- max(edf$end)

  arrow_from <- ifelse(!is.na(strand0) && strand0=="-",ge,gs)
  arrow_to <- ifelse(!is.na(strand0) && strand0=="-",gs,ge)

  p_gene <- ggplot()+
    geom_segment(
      aes(
        x=arrow_from,xend=arrow_to,
        y=.5,yend=.5
      ),
      linewidth=.45,
      arrow=grid::arrow(
        length=grid::unit(.08,"inches"),
        type="closed"
      )
    )+
    geom_rect(
      data=edf,
      aes(
        xmin=start,xmax=end,
        ymin=.22,ymax=.78
      ),
      fill="grey15"
    )+
    xscale(TRUE)+
    coord_cartesian(ylim=c(0,1),clip="off")+
    theme_classic(base_size=10)+
    theme(
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      axis.title.y=element_blank()
    )+
    labs(
      x=paste0(region$chrom," genomic position"),
      title=paste0(
        "Gene model",
        ifelse(
          is.na(strand0),
          "",
          paste0(" | ",strand0," strand")
        )
      )
    )

} else {
  p_gene <- empty_panel("Gene model")
}

# Combine

# Keep all terminal APA panels aligned to the coverage zoom.
# This is deliberately applied after the individual ggplots are built.
if(mode=="short-read" &&
   exists("terminal_zoom_start") &&
   exists("terminal_zoom_end")){

  terminal_xlim <- c(
    terminal_zoom_start,
    terminal_zoom_end
  )

  if(exists("p_apa") && !is.null(p_apa)){
    p_apa <- p_apa +
      coord_cartesian(
        xlim=terminal_xlim,
        clip="off"
      )
  }

  if(exists("p_ipa") && !is.null(p_ipa)){
    p_ipa <- p_ipa +
      coord_cartesian(
        xlim=terminal_xlim,
        clip="off"
      )
  }

  if(exists("p_ale") && !is.null(p_ale)){
    p_ale <- p_ale +
      coord_cartesian(
        xlim=terminal_xlim,
        clip="off"
      )
  }

  if(exists("p_gene") && !is.null(p_gene)){
    p_gene <- p_gene +
      coord_cartesian(
        xlim=terminal_xlim,
        clip="off"
      )
  }
}

# Standalone plots show GTF-derived terminal structure only.
# Do not present annotation-derived branches as a tested APAURA
# TUTR/ALE/mixed architecture call.
if(standalone && exists("p_gene") && !is.null(p_gene)){
  p_gene <- p_gene +
    labs(
      title=if(empty(gene)){
        "Annotated terminal-exon structure"
      } else {
        "Annotated 3′ terminal-exon architecture"
      },
      subtitle=paste0(
        "Blue = annotated terminal region",
        " • differential APA/PAS statistics not available",
        " — no APAURA results supplied"
      )
    )
}

show_tutr_panel <- nrow(apa)>0
show_ipa_panel <- show_ipa && nrow(ipa)>0

if(mode=="short-read" &&
   show_overview &&
   exists("p_overview") &&
   !is.null(p_overview)){
  plot_parts <- list(
    p_overview,
    p_evidence
  )
  plot_heights <- c(
    1.65,
    4.25
  )
} else {
  plot_parts <- list(p_evidence)
  plot_heights <- c(5)
}

if(show_tutr_panel){
  plot_parts <- c(plot_parts,list(p_apa))
  plot_heights <- c(plot_heights,1.55)
}

if(show_ipa_panel){
  plot_parts <- c(plot_parts,list(p_ipa))
  plot_heights <- c(plot_heights,1.55)
}

if(show_ale_panel){
  plot_parts <- c(plot_parts,list(p_ale))
  plot_heights <- c(plot_heights,1.70)
}

plot_parts <- c(plot_parts,list(p_gene))
plot_heights <- c(plot_heights,gene_height)

final_plot <- wrap_plots(
  plotlist=plot_parts,
  ncol=1,
  heights=plot_heights
)

plot_height_final <- if(
  mode=="short-read" &&
  show_overview &&
  exists("p_overview") &&
  !is.null(p_overview)
){
  plot_height + 1.35
} else {
  plot_height
}

# Save

dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)

save_one <- function(ext){
  f <- paste0(out,".",ext)

  if(ext=="pdf"){
    ggsave(f,final_plot,width=plot_width,height=plot_height_final)
  } else if(ext=="png"){
    ggsave(f,final_plot,width=plot_width,height=plot_height_final,dpi=300)
  } else if(ext=="svg"){
    ggsave(f,final_plot,width=plot_width,height=plot_height_final,device="svg")
  }

  message("[APAURA plot] wrote ",f)
}

if(format_out=="all"){
  save_one("pdf")
  save_one("png")
  save_one("svg")
} else {
  save_one(format_out)
}

write_tsv(apa,paste0(out,".apa_sites.tsv"))

if(show_ipa)
  write_tsv(ipa,paste0(out,".ipa_sites.tsv"))

message("[APAURA plot] finished")
