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

# ============================================================
# CLI
# ============================================================

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
region_arg <- getv("region")
out <- getv("out","APAURA_region")

flank <- int("flank",1000)
bin_size <- int("bin_size",10)
min_mapq <- int("min_mapq",10)
merge_method <- tolower(getv("merge","mean"))
library_mode <- tolower(getv("library","unstranded"))

show_replicates <- flag("show_replicates")
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

# ============================================================
# APA sites
# ============================================================

apa_file <- file.path(results_dir,"apa_sites.tsv")
if(!file.exists(apa_file)) stop("Missing: ",apa_file)

apa <- read_tab(apa_file)

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
  mutate(significant=status=="ok" & !is.na(fdr) & fdr<=fdr_cut & abs(delta_pau)>=apa_effect,
    direction=case_when(delta_pau>0 ~ "Higher PAU in case", delta_pau<0 ~ "Lower PAU in case", TRUE ~ "No change"))

message("[APAURA plot] terminal PAS in region: ",nrow(apa))

# IPA

ipa <- tibble()
ipa_file <- file.path(results_dir,"ipa_sites.tsv")

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

# ============================================================
# Evidence
# ============================================================

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
  # ============================================================
  # Persistent normalization cache
  # ============================================================

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
      message("[APAURA plot] normalization cache hit: ",sample,
              " | total=",round(total))
      return(total)
    }

    message("[APAURA plot] normalization cache miss: ",sample, " | scanning BAM once")

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

  # ----------------------------------------------------------
  # Coverage plot
  # ----------------------------------------------------------

 p_evidence <- ggplot(merged_df,aes(pos,plot_cov))+
   geom_hline(yintercept=0,linewidth=.3,colour="grey35")+
   geom_area(aes(fill=condition),alpha=.16,colour=NA)

 if(show_replicates){
   p_evidence <- p_evidence+
     geom_line(data=cov_df, aes(pos,plot_cov, colour=condition,
         group=interaction(condition,sample)), linewidth=.18, alpha=.15, inherit.aes=FALSE)
 }

 p_evidence <- p_evidence+
   geom_line(aes(colour=condition), linewidth=.65)+
   facet_grid(label~.,switch="y")+
   scale_colour_manual(values=COND_COL)+
   scale_fill_manual(values=COND_COL)+
   scale_y_continuous(labels=function(x) scales::number(abs(x),accuracy=.1))+
   xscale(FALSE)+
   theme_classic(base_size=11)+
   theme(legend.position="none",
     legend.title=element_blank(),
     strip.background=element_blank(),
     strip.placement="outside",
     strip.text.y.left=element_text(angle=0,face="bold"),
     axis.text.x=element_blank(),
     axis.ticks.x=element_blank())+
   labs(x=NULL, y=paste0(ifelse(merge_method=="mean","Mean","Median"), " CPM coverage"),
     title=paste0(gene_label, ifelse(is.na(strand0),"",paste0(" | ",strand0," strand"))),
     subtitle=paste0(control," vs ",case,
       ifelse(library_mode=="unstranded", "", paste0(" • ",library_mode,"-stranded"))))

# ============================================================
# Direct RNA / IsoQuant
# ============================================================

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

# ============================================================
# APA panel
# ============================================================

if(nrow(apa)){

  p_apa <- ggplot(apa,aes(pas_position_1based,delta_pau))+
    geom_hline(yintercept=0,linewidth=.35,colour="grey25")+
    geom_hline(yintercept=c(-apa_effect,apa_effect), linetype=3, linewidth=.3, colour="grey70")+
    geom_segment(aes(xend=pas_position_1based, y=0, yend=delta_pau, colour=direction, alpha=significant), linewidth=.7)+
    geom_point(aes(colour=direction,alpha=significant), size=3)+
    scale_colour_manual(values=EVENT_COL)+
    scale_alpha_manual(values=c(`FALSE`=.42,`TRUE`=1), guide="none")+
    xscale(FALSE)+
    theme_classic(base_size=10)+
    theme(legend.position="right",
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())+
    labs(x=NULL, y=expression(Delta*"PAU"),
      colour=NULL, title="Terminal poly(A) site usage • case − control")

} else {
  p_apa <- empty_panel("Terminal PAS: none in region")
}

# ============================================================
# IPA panel
# ============================================================

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
      labs(x=NULL, y=expression(Delta*"IPA"), colour=NULL, title="Intronic polyadenylation usage • case − control")

  } else {
    p_ipa <- empty_panel("Intronic PAS: none in region")
  }
}

# Gene model
exons <- gtf[gtf$type=="exon" & as.character(seqnames(gtf))==region$chrom & start(gtf)<=region$end & end(gtf)>=region$start]

if(!empty(gene) && length(exons)){
  gid <- as.character(exons$gene_id)
  gnm <- as.character(exons$gene_name)

  exons <- exons[gnm==gene | clean_ens(gid)==clean_ens(gene)]
}

gene_height <- 1.4

if(length(exons) && show_transcripts &&
   "transcript_id" %in% names(mcols(exons))){

  tx <- tibble(transcript=as.character(exons$transcript_id),
    start=start(exons), end=end(exons)) %>%
    filter(!is.na(transcript),transcript!="")

  keep <- tx %>%
    count(transcript,name="n") %>%
    arrange(desc(n)) %>%
    slice_head(n=max_transcripts) %>%
    pull(transcript)

  tx <- tx %>%
    filter(transcript %in% keep) %>%
    mutate(transcript=factor(transcript,levels=rev(keep)), y=as.numeric(transcript))

  lines <- tx %>%
    group_by(transcript,y) %>%
    summarise(start=min(start), end=max(end), .groups="drop")

  p_gene <- ggplot()+
    geom_segment(data=lines, aes(x=start,xend=end,y=y,yend=y), linewidth=.3, colour="grey50")+
    geom_rect(data=tx, aes(xmin=start,xmax=end, ymin=y-.28,ymax=y+.28), fill="grey15")+
    scale_y_continuous(breaks=seq_along(rev(keep)), labels=rev(keep))+
    xscale(TRUE)+
    theme_classic(base_size=9)+
    theme(axis.title.y=element_blank(),
      axis.ticks.y=element_blank(),
      axis.line.y=element_blank())+
    labs(x=paste0(region$chrom," genomic position"),
      title="Annotated transcripts")

  gene_height <- min(5,max(1.5,length(keep)*.22+.5))

} else if(length(exons)){

  er <- GenomicRanges::reduce(granges(exons))

  edf <- tibble(start=start(er), end=end(er))

  gs <- min(edf$start)
  ge <- max(edf$end)

  arrow_from <- ifelse(!is.na(strand0) && strand0=="-",ge,gs)
  arrow_to <- ifelse(!is.na(strand0) && strand0=="-",gs,ge)

  p_gene <- ggplot()+
    geom_segment(aes(x=arrow_from,xend=arrow_to, y=.5,yend=.5),
      linewidth=.45,
      arrow=grid::arrow(length=grid::unit(.08,"inches"), type="closed"))+
    geom_rect(data=edf, aes(xmin=start,xmax=end, ymin=.22,ymax=.78), fill="grey15")+
    xscale(TRUE)+
    coord_cartesian(ylim=c(0,1),clip="off")+
    theme_classic(base_size=10)+
    theme(axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      axis.title.y=element_blank())+
    labs(x=paste0(region$chrom," genomic position"),
      title=paste0("Gene model",
        ifelse(is.na(strand0), "", paste0(" | ",strand0," strand | 5′→3′ arrow"))))

} else {
  p_gene <- empty_panel("Gene model")
}

# Combine

if(show_ipa){
  final_plot <- p_evidence/p_apa/p_ipa/p_gene +
    plot_layout(heights=c(5,1.6,1.6,gene_height))
} else {
  final_plot <- p_evidence/p_apa/p_gene +
    plot_layout(heights=c(5,1.7,gene_height))
}

# Save

dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)

save_one <- function(ext){
  f <- paste0(out,".",ext)

  if(ext=="pdf"){
    ggsave(f,final_plot,width=plot_width,height=plot_height)
  } else if(ext=="png"){
    ggsave(f,final_plot,width=plot_width,height=plot_height,dpi=300)
  } else if(ext=="svg"){
    ggsave(f,final_plot,width=plot_width,height=plot_height,device="svg")
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