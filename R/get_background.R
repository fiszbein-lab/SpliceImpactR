#' Parse a genomic interval string into components
#'
#' Converts a coordinate string of the form \code{"chr:start-end"} into a
#' \link[data.table]{data.table} with separate columns for chromosome, start,
#' and stop positions.
#'
#' @param x Character vector of genomic interval strings.
#'
#' @return A \link[data.table]{data.table} with columns \code{chr}, \code{start},
#'   and \code{stop}.
#' @keywords internal
#' @noRd
#'
#' @importFrom data.table data.table tstrsplit
parse_one_interval <- function(x) {
  x <- trimws(as.character(x))
  # expected "chrX:start-end"
  p <- data.table::tstrsplit(x, ":", fixed = FALSE)
  chr <- p[[1]]
  ab  <- data.table::tstrsplit(p[[2]], "-", fixed = TRUE)
  data.table::data.table(
    chr   = chr,
    start = suppressWarnings(as.integer(ab[[1]])),
    stop  = suppressWarnings(as.integer(ab[[2]]))
  )
}

#' Match input exon coordinates to annotated exons
#'
#' Matches exon intervals from a user-supplied table (e.g., background or
#' candidate exons) to reference exon annotations using genomic overlap.
#' The function reports best matches per query exon, optionally preferring
#' annotated exons in protein-coding transcripts.
#'
#' @param exon_df A \code{data.frame} or data.table with
#'   at least the columns \code{gene}, \code{exon}, and \code{ID}, where
#'   \code{exon} is formatted as \code{"chr:start-end"}.
#' @param annotations A reference annotation table as produced by
#'   \code{build_from_annotations()}.
#' @param minOverlap Numeric value in [0,1]; minimum fraction of the query
#'   interval that must overlap a reference exon. Default: 0.8.
#' @return A data.table with the subset of exons successfully
#'   matched to annotations, containing columns:
#'   \code{gene_id}, \code{exon}, \code{transcript_id}, and \code{exon_id}.
#'
#' @details
#' When multiple annotated exons overlap a query, the best match is chosen
#' according to (1) matching classification (first/internal/last),
#' (2) presence of a protein-linked transcript, (3) highest query coverage,
#' and (4) intersection width.
#'
#'
#' @importFrom data.table as.data.table setorder setcolorder fifelse
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges width pintersect
#' @keywords internal
#' @noRd
match_exon_table <- function(exon_df,
                             annotations,
                             minOverlap = 0.8) {

  X  <- data.table::as.data.table(exon_df)
  need <- c("gene","exon","ID")
  miss <- setdiff(need, names(X)); if (length(miss)) stop("Missing columns: ", paste(miss, collapse=", "))

  # normalize
  X[, gene_id := sub("\\.\\d+$","", as.character(gene))]
  X[, `:=`(ID = as.character(ID),
           exon = as.character(exon))]
  coords <- parse_one_interval(X$exon)
  X[, c("chr","start","stop") := coords]

  # build once
  AA <- build_from_annotations(annotations)
  EX <- AA$exons
  TX <- AA$transcripts
  protein_tx <- AA$protein_tx
  tx_rownum  <- setNames(TX$rownum, TX$transcript_id)

  # subset exons to genes & chromosomes we actually need (huge speedup)
  want_genes <- unique(X$gene_id)
  want_chr   <- unique(X$chr)
  EX_sub <- EX[gene_id %chin% want_genes & chr %chin% want_chr]

  if (!nrow(EX_sub)) {
    return(X[, .(row = .I, gene, exon, ID,
                 hit = FALSE,
                 transcript_id = NA_character_, exon_id = NA_character_,
                 classification = NA_character_,
                 recip_q = NA_real_, recip_s = NA_real_, int_width = NA_integer_,
                 protein_linked = NA)])
  }

  # make GRanges (query has no strand; gene join makes strand moot)
  gr_q <- GenomicRanges::GRanges(seqnames = X$chr,
                                 ranges = IRanges::IRanges(X$start, X$stop))
  gr_s <- GenomicRanges::GRanges(seqnames = EX_sub$chr,
                                 ranges = IRanges::IRanges(EX_sub$start, EX_sub$stop))
  H    <- GenomicRanges::findOverlaps(gr_q, gr_s, ignore.strand = TRUE)
  if (!length(H)) {
    return(X[, .(row = .I, gene, exon, ID,
                 hit = FALSE,
                 transcript_id = NA_character_, exon_id = NA_character_,
                 classification = NA_character_,
                 recip_q = NA_real_, recip_s = NA_real_, int_width = NA_integer_,
                 protein_linked = NA)])
  }

  # compute overlaps
  q  <- gr_q[queryHits(H)]
  s  <- gr_s[subjectHits(H)]
  iw <- as.numeric(IRanges::width(pintersect(q, s)))
  rq <- iw / as.numeric(IRanges::width(q))   # coverage of the query (our exon)
  rs <- iw / as.numeric(IRanges::width(s))   # coverage of the annotated exon

  keep <- which(rq >= minOverlap)
  if (!length(keep)) {
    return(X[, .(row = .I, gene, exon, ID,
                 hit = FALSE,
                 transcript_id = NA_character_, exon_id = NA_character_,
                 classification = NA_character_,
                 recip_q = NA_real_, recip_s = NA_real_, int_width = NA_integer_,
                 protein_linked = NA)])
  }

  hits <- data.table(
    row            = queryHits(H)[keep],
    gene_id_hit    = EX_sub$gene_id[subjectHits(H)[keep]],
    transcript_id  = EX_sub$transcript_id[subjectHits(H)[keep]],
    exon_id        = EX_sub$exon_id[subjectHits(H)[keep]],
    classification = EX_sub$classification[subjectHits(H)[keep]],
    recip_q        = rq[keep],
    recip_s        = rs[keep],
    int_width      = iw[keep]
  )

  # restrict to same gene (rigid input)
  hits <- hits[gene_id_hit %chin% X$gene_id[row]]
  if (!nrow(hits)) {
    return(X[, .(row = .I, gene, exon, ID,
                 hit = FALSE,
                 transcript_id = NA_character_, exon_id = NA_character_,
                 classification = NA_character_,
                 recip_q = NA_real_, recip_s = NA_real_, int_width = NA_integer_,
                 protein_linked = NA)])
  }

  hits[, protein_linked := transcript_id %chin% protein_tx]

  # prefer class that matches provided ID (first/internal/last)
  hits[, class_match := (classification == tolower(X$ID[row])) |
         (classification == X$ID[row])]

  # choose best per input row (vectorized by ordering + first)
  setorder(hits, row, -class_match, -protein_linked, -recip_q, -int_width)
  best <- hits[, .SD[1L], by = row]

  # assemble output
  out <- best[X[, .(row = .I, gene, exon, ID)], on = "row"]
  out[, `:=`(
    hit = !is.na(transcript_id),
    transcript_row = ifelse(is.na(transcript_id), NA_integer_, tx_rownum[transcript_id])
  )]

  # tidy columns
  out[, gene_id := sub("\\.\\d+$","", as.character(gene))]
  setcolorder(out, c("row","gene_id","exon","ID",
                     "hit","transcript_id","transcript_row","exon_id","classification",
                     "recip_q","recip_s","int_width","protein_linked"))
  background_out <- out[class_match == TRUE & hit == TRUE & gene_id_hit == gene_id, .(gene_id, exon, transcript_id, exon_id)]
  return(background_out)
}

#' Read and filter exon background files
#'
#' Internal helper that loads exon-level summary files from specified directories,
#' filters low-coverage or unannotated exons, and standardizes exon classification
#' labels for downstream background construction.
#'
#' @param paths_df Data frame containing at least a column \code{path} with
#'   directories to read. An optional column \code{sample_name} may also be supplied.
#' @param keep_annotated_first_last Logical; if \code{TRUE}, retain only annotated
#'   first/last exons with valid coverage flags (\code{nFE}, \code{nLE}). Default: \code{FALSE}.
#'
#' @return Adata.table (or tibble-like) with columns
#'   \code{gene}, \code{exon}, and standardized \code{ID}.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom data.table rbindlist
#' @importFrom dplyr select
read_background <- function(paths_df, keep_annotated_first_last=FALSE) {
  stopifnot(is.data.frame(paths_df))
  req <- c("path")
  miss <- setdiff(req, names(paths_df))
  if (length(miss)) stop("Data frame must include: ", paste(req, collapse=", "))

  # fill sample_name if missing
  if (!"sample_name" %in% names(paths_df)) {
    paths_df$sample_name <- basename(normalizePath(paths_df$path, mustWork=FALSE))
  }

  parts <- lapply(seq_len(nrow(paths_df)), function(i){
    files <- .read_exon_files(paste0(paths_df$path[i], basename(paths_df$path[i]), '.'), columns = c("gene", "exon", "ID", "nFE", "nLE", "nUP", "nDOWN"))
    if (!length(files)) stop("No .exon files under: ", paths_df$path[i])
    files
  })
  exon_files <- rbindlist(parts)
  if (keep_annotated_first_last == TRUE) {
    exon_files <- exon_files[!((ID == "first" | ID == "FirstInternal_high" | ID == "FirstInternal_medium") & nFE == 0),]
    exon_files <- exon_files[!((ID == "last" | ID == "InternalLast_high" | ID == "InternalLast_medium") & nLE == 0),]
  }

  exon_files <- exon_files[nUP >= 10 | nDOWN >= 10] %>% dplyr::select(gene, exon, ID)
  map <- c(
    "first"          = "first",
    "FirstInternal_high"    = "hybrid",
    "FirstInternal_medium" = "hybrd",
    "internal" = "internal",
    "last"          = "last",
    "InternalLast_high"    = "hybrid",
    "InternalLast_medium" = "hybrd"
  )
  exon_files$ID <- map[exon_files$ID]

  exon_files <- exon_files[!duplicated(exon_files)]

  return(exon_files)

}


#' Build domain difference maps for background transcript pairs
#'
#' Internal helper that annotates background transcript pairs with protein domain
#' sets derived from annotated protein features. Computes domains unique to each
#' transcript and summary counts of gained/lost domains.
#'
#' @param background A data.table containing paired transcript IDs in columns
#'   \code{transcript_id_1} and \code{transcript_id_2}.
#' @param protein_features A data.frame or data.table of protein domain features,
#'   typically from \code{get_comprehensive_annotations}
#'
#' @return A data.table with domain lists and summary metrics
#'   per transcript pair, including counts of unique and shared domains.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom data.table as.data.table setnames
#' @importFrom BiocParallel bpparam bplapply
get_domain_background <- function(background, protein_features) {

  Pf <- data.table::as.data.table(protein_features)
  # Pf <- Pf[, domain_id := paste0(database, ";", gsub("[|]", " ", gsub(";", " ", name)))][!is.na(domain_id) & nzchar(domain_id)]

  Pf[, domain_id := paste0(
    database, ";",
    gsub("[|;]", " ", clean_name), " ",
    sub(".*;", "", name)
  )]
  Pf <- Pf[!is.na(domain_id) & nzchar(domain_id)]

  # keep one row per (tx, protein, domain)
  D <- unique(Pf[,
                 .(transcript_id = as.character(ensembl_transcript_id),
                   protein_id    = as.character(ensembl_peptide_id),
                   domain_id     = as.character(domain_id)
                 )
  ])



  doms <- D[, .(domains = list(unique(domain_id))), by = transcript_id]

  bg <- doms[background, on = .(transcript_id = transcript_id_1)]
  data.table::setnames(bg, "domains", "d1")
  bg <- doms[bg, on = .(transcript_id = transcript_id_2)]
  data.table::setnames(bg, "domains", "d2")

  bg[, `:=`(
    d1 = lapply(d1, function(x) if (is.null(x)) character() else x),
    d2 = lapply(d2, function(x) if (is.null(x)) character() else x)
  )]
  bg <- bg[lengths(d1) > 0 & lengths(d2) > 0]
  bg <- bg[vapply(Map(setequal, d1, d2), isFALSE, logical(1))]
  bg[, c("d1_pruned", "d2_pruned") := {
    tmp <- Map(function(x, y) {
      common <- intersect(x, y)
      list(setdiff(x, common), setdiff(y, common))
    }, d1, d2)
    list(lapply(tmp, `[[`, 1), lapply(tmp, `[[`, 2))
  }]
  bg <- bg[lengths(d1_pruned) > 0 & lengths(d2_pruned) > 0]

  bg[, `:=`(
    d1_parsed = lapply(d1_pruned, .parse_domain_coords),
    d2_parsed = lapply(d2_pruned, .parse_domain_coords)
  )]

  param <- BiocParallel::bpparam()
  res <- BiocParallel::bplapply(seq_len(nrow(bg)), function(i) list(
    d1 = .diff_domains(bg$d1_parsed[[i]], bg$d2_parsed[[i]])$domain,
    d2 = .diff_domains(bg$d2_parsed[[i]], bg$d1_parsed[[i]])$domain
  ), BPPARAM = param)
  bg[, domains_1 := lapply(res, `[[`, "d1")]
  bg[, domains_2 := lapply(res, `[[`, "d2")]

  bg[, `:=`(n_domains_1 = lengths(d1),
            n_domains_2 = lengths(d2),
            d1 = NULL,
            d2 = NULL,
            d1_pruned = NULL,
            d2_pruned = NULL,
            d1_parsed = NULL,
            d2_parsed = NULL,
            sd_n_domains_1 = lengths(domains_1),
            sd_n_domains_2 = lengths(domains_2),
            total_sd_n_domains = lengths(domains_1)+lengths(domains_2))]
  bg[, total_sd_domains := Map(union, domains_1, domains_2)]

  bg[]
}

#' Generate within-gene transcript pairs for background analysis
#'
#' Internal helper that constructs all unique transcript pairs per gene from a
#' background data table. Used to define background isoform pairs for domain or
#' feature comparison analyses.
#'
#' @param background A \link[data.table]{data.table} or \code{data.frame}
#'   containing at least the columns \code{gene_id} and \code{transcript_id},
#'   \code{matched_background} from \code{get_background}
#'
#' @return A \link[data.table]{data.table} with columns:
#'   \describe{
#'     \item{gene_id}{Ensembl gene identifier.}
#'     \item{transcript_id_1, transcript_id_2}{Paired transcript identifiers.}
#'     \item{n_pairs_in_gene}{Number of total transcript pairs within the gene.}
#'   }
#'
#' @details
#' The function removes missing or duplicated transcript IDs, performs an
#' in-gene Cartesian self-join, and keeps only ordered transcript pairs
#' (\code{transcript_id_1 < transcript_id_2}) to avoid duplication.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom data.table as.data.table setkey
make_transcript_pairs <- function(background) {
  B <- data.table::as.data.table(background)

  # 1) unique transcripts per gene (drop exon-level duplication, NAs)
  U <- unique(B[!is.na(transcript_id), .(gene_id, transcript_id)])

  # 2) self cross-join by gene, then keep ordered i<j to remove redundant pairs
  data.table::setkey(U, gene_id, transcript_id)
  P <- U[U, on = "gene_id", allow.cartesian = TRUE][
    transcript_id < i.transcript_id,
    .(gene_id,
      transcript_id_1 = transcript_id,
      transcript_id_2 = i.transcript_id)
  ]

  # (optional) attach pair counts per gene
  P[, n_pairs_in_gene := .N, by = gene_id][]
}

#' Build a per-transcript length index from annotations
#'
#' Internal helper that aggregates exon-level annotation data into
#' transcript-level lengths (CDS, total exon, and amino-acid length).
#' Optionally attaches peptide identifiers if present.
#'
#' @param annotations A \link[data.table]{data.table} annotations from get_annotations
#'
#' @return A \link[data.table]{data.table} with one row per transcript, including:
#'   \itemize{
#'     \item \code{transcript_id}
#'     \item \code{cds_bp} – total CDS bases
#'     \item \code{exon_bp_all} – total exon span bases
#'     \item \code{prot_aa} – translated amino-acid length
#'     \item \code{peptide_id} – optional peptide identifier
#'   }
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom data.table as.data.table setnames
build_tx_length_index <- function(annotations) {
  A <- data.table::as.data.table(annotations)

  # Keep exon rows only; de-dup by (tx, exon) before summing
  EX <- A[type %chin% "exon",
          .(transcript_id = as.character(transcript_id),
            exon_id       = as.character(exon_id),
            cds_len       = as.numeric(cds_len),
            feature_len   = as.numeric(feature_length),
            type          = as.character(type)
          )
  ]

  # Some GTFs can repeat exon rows; keep one per (tx, exon)
  EX <- unique(EX, by = c("transcript_id","exon_id"))

  # Aggregate per transcript
  TX <- EX[type == "exon", .(
    cds_bp      = sum(cds_len,     na.rm = TRUE),
    exon_bp_all = sum(feature_len, na.rm = TRUE)
  ), by = transcript_id]

  # Convert cds_bp to amino-acid length (floor; stop not subtracted)
  TX[, prot_aa := ifelse(is.finite(cds_bp), floor(cds_bp / 3), NA_integer_)]

  # Optional: attach peptide IDs if available anywhere in annotations
  pep_col <- intersect(c("ensembl_peptide_id","protein_id","peptide_id"), names(A))
  if (length(pep_col)) {
    PEP <- A[!is.na(get(pep_col[1])) & nzchar(as.character(get(pep_col[1]))),
             .(pep = as.character(get(pep_col[1]))), by = transcript_id]
    # One peptide per transcript: take the first
    PEP <- PEP[!duplicated(transcript_id)]
    TX <- PEP[TX, on = "transcript_id"]
    data.table::setnames(TX, "pep", "peptide_id")
  } else {
    TX[, peptide_id := NA_character_]
  }

  TX[]
}

#' Annotate transcript pairs with length differences
#'
#' Internal helper that attaches CDS, exon, and amino-acid lengths to transcript
#' pairs (generated by \code{make_transcript_pairs()}) and computes signed and
#' absolute differences between each pair.
#'
#' @param background A data.frame or \link[data.table]{data.table} containing
#'   transcript pairs, as produced by \code{make_transcript_pairs()}.
#' @param annotations GTF-like annotations containing exon-level features, used
#'   by \code{build_tx_length_index()} to derive transcript lengths. From
#'   \code{get_annotation} annotations
#'
#' @return A \link[data.table]{data.table} with per-pair length statistics,
#'   including:
#'   \itemize{
#'     \item \code{cds_bp_*}, \code{exon_bp_all_*}, \code{prot_aa_*}
#'     \item \code{d_*} (signed differences)
#'     \item \code{d_*_abs} (absolute differences)
#'   }
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom data.table setkey setnames setcolorder
make_transcript_pairs_with_lengths <- function(background, annotations) {
  pairs <- make_transcript_pairs(background)
  TX    <- build_tx_length_index(annotations)

  data.table::setkey(TX, transcript_id)

  out <- pairs[
    TX, on = c(transcript_id_1 = "transcript_id")
  ][
    TX, on = c(transcript_id_2 = "transcript_id"),
    nomatch = 0L
  ]

  # tidy names
  data.table::setnames(out,
           old = c("cds_bp","exon_bp_all","prot_aa","peptide_id"),
           new = c("cds_bp_1","exon_bp_all_1","prot_aa_1","peptide_id_1"))
  data.table::setnames(out,
           old = c("i.cds_bp","i.exon_bp_all","i.prot_aa","i.peptide_id"),
           new = c("cds_bp_2","exon_bp_all_2","prot_aa_2","peptide_id_2"),
           skip_absent = TRUE)  # second join’s cols already have these names

  # Compute signed and absolute changes
  out[, `:=`(
    d_prot_aa      = prot_aa_1      - prot_aa_2,
    d_prot_aa_abs  = abs(prot_aa_1  - prot_aa_2),
    d_cds_bp       = cds_bp_1       - cds_bp_2,
    d_cds_bp_abs   = abs(cds_bp_1   - cds_bp_2),
    d_exon_bp_all  = exon_bp_all_1  - exon_bp_all_2,
    d_exon_bp_abs  = abs(exon_bp_all_1 - exon_bp_all_2)
  )]

  # Nice ordering
  data.table::setcolorder(out, c(
    "gene_id","n_pairs_in_gene",
    "transcript_id_1","peptide_id_1","prot_aa_1","cds_bp_1","exon_bp_all_1",
    "transcript_id_2","peptide_id_2","prot_aa_2","cds_bp_2","exon_bp_all_2",
    "d_prot_aa","d_prot_aa_abs","d_cds_bp","d_cds_bp_abs","d_exon_bp_all","d_exon_bp_abs"
  ))
  out[]
}

#' Attach gene IDs to transcripts
#'
#' Internal helper that maps \code{transcript_id} values to their
#' corresponding \code{gene_id}s using an annotation table.
#'
#' @param transcripts A \link[data.table]{data.table} or \code{data.frame}
#'   containing a column \code{transcript_id}.
#' @param annotations A \link[data.table]{data.table} with at least
#'   \code{transcript_id} and \code{gene_id} columns.
#'
#' @return The input data with an added \code{gene_id_raw} column.
#' @importFrom data.table setDT
#' @keywords internal
#' @noRd
#'
#' @importFrom data.table as.data.table
get_genes_from_transcripts <- function(matched_background, A) {
  A <- data.table::as.data.table(A)
  MAP <- A[!is.na(transcript_id) & !is.na(gene_id),
           .(transcript_id = as.character(transcript_id),
             gene_id   = as.character(gene_id))]
  MAP <- unique(MAP, by = "transcript_id")
  data.table::setDT(matched_background)
  matched_background[, transcript_id := sub("\\.\\d+$", "", transcript_id)]
  matched_genes_background <- MAP[matched_background, on  = 'transcript_id']
  matched_genes_background[]
}


#' Build a transcript-pair background with domain and length annotations
#'
#' Constructs a background dataset of transcript pairs suitable for
#' domain-level or exon-level enrichment analyses. Depending on the
#' \code{source} parameter, the function can derive the background from
#' HIT index results, annotated transcripts, or a user-provided list
#' of transcript IDs.
#'
#'
#' @param source Character string specifying the source of the background.
#'   One of:
#'   \itemize{
#'     \item \code{"hit_index"} use HIT index output directories.
#'     \item \code{"annotated"} use all transcripts from the annotation.
#'     \item \code{"user-given"} use a user-supplied list of transcript IDs.
#'   }
#' @param input Source-specific input:
#'   \itemize{
#'     \item For \code{"hit_index"}: a data.frame of paths with a \code{path} column.
#'     \item For \code{"annotated"}: ignored.
#'     \item For \code{"user-given"}: a character vector or data.frame containing \code{transcript_id}.
#'   }
#' @param annotations A data.table annotations from \code{get_annotations}
#' @param protein_features A data.frame of protein domain or motif features
#'   (e.g. InterPro, Pfam) with at least
#'   \code{ensembl_transcript_id}, \code{ensembl_peptide_id},
#'   \code{database}, and \code{name} columns from
#'   \code{get_comprehensive_annotations}
#' @param keep_annotated_first_last Logical; passed to
#'   read_background to control filtering of annotated
#'   first/last exons. Defaults to \code{TRUE}.
#' @param minOverlap Numeric (0-1); minimum fraction overlap required when
#'   matching exons between HIT index data and annotations. Defaults to 0.8.
#'
#' @return A \link[data.table]{data.table} in which each row represents
#'   a transcript pair annotated with gene ID, CDS/exon length differences,
#'   and protein domain differences. Columns include:
#'   \code{gene_id}, \code{transcript_id_1}, \code{transcript_id_2},
#'   \code{prot_aa_1}, \code{prot_aa_2}, \code{domains_1}, \code{domains_2},
#'   and related summary metrics.
#'
#'
#' @examples
#'
#' annots <- get_annotation(load = "test")
#' interpro_features <- get_protein_features(c("interpro"), annots$annotations, timeout = 600, test = TRUE)
#' protein_feature_total <- get_comprehensive_annotations(list(interpro_features))
#'
#' # Build background from HIT index paths
#' sample_frame <- data.frame(path = c(check_extdata_dir('rawData/control_S5/'),
#'                                     check_extdata_dir('rawData/control_S6/'),
#'                                     check_extdata_dir('rawData/control_S7/'),
#'                                     check_extdata_dir('rawData/control_S8/'),
#'                                     check_extdata_dir('rawData/case_S1/'),
#'                                     check_extdata_dir('rawData/case_S2/'),
#'                                     check_extdata_dir('rawData/case_S3/'),
#'                                     check_extdata_dir('rawData/case_S4/')),
#'                            sample_name  = c("S5", "S6", "S7", "S8", "S1", "S2", "S3", "S4"),
#'                            condition    = c("control", "control", "control", "control", "case",  "case",  "case",  "case"),
#'                            stringsAsFactors = FALSE)
#' bg <- get_background(source = "hit_index",
#'                      input = sample_frame,
#'                      annotations = annots$annotations,
#'                      protein_features = protein_feature_total)
#'
#' # Or from user-supplied transcript IDs
#' tx_ids <- c("ENST00000466994","ENST00000484435")
#' bg_user <- get_background(source = "user-given",
#'                           input = tx_ids,
#'                           annotations = annots$annotations,
#'                           protein_features = protein_feature_total)
#'
#'
#' @export
get_background <- function(source = c("hit_index", "annotated", "user-given"),
                           input,
                           annotations,
                           protein_features,
                           keep_annotated_first_last = TRUE,
                           minOverlap = 0.8) {
  if (source == "hit_index") {
    background_init <- read_background(paths_df = input, keep_annotated_first_last)
    matched_background <- match_exon_table(exon_df = background_init, annotations = annotations, minOverlap)
  } else if (source == "annotated") {
    matched_background <- annotations[, .(transcript_id, gene_id)]
  } else if (source == "user-given") {
    if (sum(input %in% annotations$transcript_id) == 0) {
      stop("No user-given transcripts found in annotations", call. = FALSE)
    }
    matched_background <- get_genes_from_transcripts(data.table(transcript_id = input), annotations)
  }
  matched_paired_background <- make_transcript_pairs_with_lengths(matched_background, annotations)
  background_domains <- get_domain_background(matched_paired_background, protein_features)
  return(background_domains)
}
