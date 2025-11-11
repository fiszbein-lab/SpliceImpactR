#' Build GENCODE download URLs (internal helper)
#'
#' Internal function to construct URLs for downloading GENCODE
#' annotation and sequence files for either human or mouse.
#' This is used by higher-level SpliceImpactR functions.
#'
#' @param species Character string, either `"human"` or `"mouse"`.
#' @param release Integer or string specifying the GENCODE release number.
#'
#' @return A named list containing URLs and the release tag.
#' @keywords internal
.si_gencode_urls <- function(species = c("human","mouse"),
                             release) {
  species <- match.arg(species)

  if (species == "human") {
    stopifnot(is.numeric(release) || grepl("^[0-9]+$", release))
    rel  <- as.integer(release)
    base <- sprintf("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_%d", rel)
    tag  <- sprintf("v%d", rel)
  } else {
    if (is.numeric(release) || grepl("^[0-9]+$", release)) release <- paste0("M", as.integer(release))
    stopifnot(grepl("^M[0-9]+$", release))
    base <- sprintf("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_%s", release)
    tag  <- sprintf("v%s", release)
  }

  list(
    gtf  = sprintf("%s/gencode.%s.annotation.gtf.gz",     base, tag),
    txfa = sprintf("%s/gencode.%s.pc_transcripts.fa.gz",  base, tag),
    aafa = sprintf("%s/gencode.%s.pc_translations.fa.gz", base, tag),
    tag  = tag
  )
}

#' Prepare GENCODE annotation and sequence assets (internal helper)
#'
#' Internal function to download, import, and cache GENCODE annotation
#' and sequence files (GTF, transcript FASTA, and protein FASTA) for a
#' given species and release. Used internally by higher-level
#' SpliceImpactR functions to ensure annotation assets are available.
#'
#' @param base_dir Character string giving the base directory where
#'   downloaded or cached files should be stored.
#' @param species Character string, either `"human"` or `"mouse"`.
#' @param release Integer or string specifying the GENCODE release
#'   number (e.g., `45` for human or `"M35"` for mouse).
#' @param mode Character string, one of `"download"` or
#'   `"import_then_cache"`. Determines whether to only download the
#'   compressed files or to import the GTF into R and cache as RDS.
#' @param use_rds_cache Logical; if `TRUE`, loads cached `.rds` GTF file
#'   if available.
#'
#' @return A list containing:
#'   \describe{
#'     \item{`paths`}{File paths for the downloaded or cached assets.}
#'     \item{`gtf_df`}{Imported GTF data frame if loaded or created.}
#'     \item{`meta`}{Metadata list with species, release, and tag.}
#'   }
#'
#' @importFrom utils download.file
#' @keywords internal
.si_prepare_assets <- function(base_dir,
                              species = c("human","mouse"),
                              release,
                              mode = c("download","import_then_cache"),
                              use_rds_cache = TRUE) {
  if (is.null(base_dir)) {
    # follow Bioconductor convention: user cache for persistent assets
    base_dir <- tools::R_user_dir("SpliceImpactR", "cache")
    message("[INFO] Using default cache directory: ", base_dir)
  }

  species <- match.arg(species); mode <- match.arg(mode)
  urls <- .si_gencode_urls(species, release)
  tag  <- urls$tag

  dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
  paths <- list(
    gtf_gz  = file.path(base_dir, sprintf("gencode.%s.annotation.gtf.gz",     tag)),
    gtf_rds = file.path(base_dir, sprintf("gencode.%s.annotation.gtf.rds",    tag)),
    txfa_gz = file.path(base_dir, sprintf("gencode.%s.pc_transcripts.fa.gz",  tag)),
    aafa_gz = file.path(base_dir, sprintf("gencode.%s.pc_translations.fa.gz", tag))
  )

  .dl <- function(url, dest) {
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    if (nzchar(Sys.which("wget")))      system2("wget", c("-q","-c","-O", shQuote(dest), shQuote(url)))
    else if (nzchar(Sys.which("curl"))) system2("curl", c("-L","-C","-","-o", shQuote(dest), shQuote(url)))
    else                                utils::download.file(url, destfile = dest, mode = "wb", quiet = TRUE)
    dest
  }

  # GTF: load RDS if present; otherwise download or import+cache
  gtf_df <- NULL
  if (use_rds_cache && file.exists(paths$gtf_rds)) {
    gtf_df <- readRDS(paths$gtf_rds)
  } else if (!file.exists(paths$gtf_gz)) {
    if (mode == "download") {
      .dl(urls$gtf, paths$gtf_gz)
    } else {
      message("Importing GTF from URL via rtracklayer (first run)")
      gtf_df <- rtracklayer::readGFF(urls$gtf)
      saveRDS(gtf_df, paths$gtf_rds)
      .dl(urls$gtf, paths$gtf_gz)
    }
  }

  if (!file.exists(paths$aafa_gz)) .dl(urls$aafa, paths$aafa_gz)
  if (!file.exists(paths$txfa_gz)) .dl(urls$txfa, paths$txfa_gz)

  list(paths = paths, gtf_df = gtf_df,
       meta = list(species = species, release = release, tag = tag))
}

#' Null-coalescing operator (internal)
#'
#' Internal infix operator returning its left-hand side if not `NULL`,
#' otherwise its right-hand side. Commonly used to provide defaults
#' for optional arguments.
#'
#' @param a Left-hand side value.
#' @param b Right-hand side value (default returned if `a` is `NULL`).
#' @noRd
#' @return `a` if not `NULL`, otherwise `b`.
#' @keywords internal
`%||%` <- function(a,b) {
  if (!is.null(a)) a else b
}

#' Strip version suffix from identifiers (internal)
#'
#' Removes trailing version components (e.g., ".1", ".2") from
#' Ensembl or similar identifiers.
#'
#' @param x Character vector of identifiers.
#' @return A character vector with version suffixes removed.
#'
#' @keywords internal
strip_ver <- function(x) {
  sub("\\..*$", "", x)
}


#' Extract an attribute value from GTF/GFF attribute fields (internal)
#'
#' Parses a specific key-value pair from the `attributes` column of a
#' GTF/GFF line, returning the value associated with the requested key.
#'
#' @param x Character vector of attribute strings.
#' @param key Character string naming the attribute key to extract
#'   (e.g., `"gene_id"` or `"transcript_id"`).
#' @param strip Logical; if `TRUE`, removes version suffixes using
#'   [strip_ver()].
#'
#' @return Character vector of extracted values, `NA` where the key
#'   was not found.
#' @keywords internal
attr_get <- function(x, key, strip = FALSE) {
  m <- regexec(paste0(key, ' "([^"]+)"'), x, perl = TRUE)
  val <- regmatches(x, m)
  out <- vapply(
    val,
    function(z) if (length(z)>=2) z[2] else NA_character_,
    character(1))
  if (strip) out <- strip_ver(out)
  out
}


#' Extract transcript identifier from protein FASTA header
#'
#' @param h Character string giving a single FASTA header line
#'   (including the leading `>` symbol if present).
#'
#' @return A character string containing the transcript identifier
#'   (without version suffix). Returns an empty string if no match
#'   is found.
#'
#' @details
#' GENCODE protein FASTA headers typically contain a transcript
#' reference such as:
#'
#' ```
#' >ENSP00000369497.3|ENST00000355832.3|ENSG00000141510.15|...
#' ```
#'
#' or may encode it in tagged fields like `transcript:ENST...` or
#' `transcript_id=ENSMUST...`. This function uses regular-expression
#' pattern matching to locate the transcript identifier, strips the
#' version suffix (e.g. `.3`), and returns the cleaned transcript ID.
#'
#' @keywords internal
prot_hdr_to_enst <- function(h) {
  for (p in c("transcript:([^| ]+)", "transcript_id:([^| ]+)", "transcript_id=([^| ]+)", "tr:([^| ]+)")) {
    if (grepl(p, h, perl = TRUE)) return(strip_ver(sub(paste0(".*\\b", p, ".*"), "\\1", h, perl = TRUE)))
  }
  f <- strsplit(h, "\\|", perl = TRUE)[[1]]
  cand <- grep("^(ENST|ENSMUST)", f, value = TRUE)  # human or mouse
  if (length(cand)) return(strip_ver(cand[1]))
  ""
}

#' Load a GTF file into a long-form data.table (internal)
#'
#' Internal helper to read a GTF or GFF file (optionally gzipped) and
#' return it as a tidy \pkg{data.table}. Ensures consistent column
#' naming, fills in missing identifiers from the `attributes` column,
#' and optionally adds a unique row identifier.
#'
#' @param gtf_path_or_df Character string giving the path or URL to a
#'   GTF/GFF file (optionally prefixed with `file://`), or an existing
#'   \code{data.frame}/\code{data.table} containing GTF-like data.
#' @param add_row_uid Logical; if `TRUE`, adds a unique `row_uid` column
#'   for reference.
#'
#' @return A \code{data.table} with standardized GTF fields and a
#'   consistent set of identifier columns.
#'
#' @details
#' The function uses \pkg{rtracklayer} to import GTF/GFF files, which
#' returns an S4 \code{DataFrame}. It is then coerced to a standard
#' \code{data.table}. When certain ID columns are missing, they are
#' extracted from the `attributes` column using [attr_get()].
#'
#' @examples
#' \dontrun{
#' gtf_dt <- load_gtf_long("gencode.v45.annotation.gtf.gz")
#' data.table::head(gtf_dt)
#' }
#'
#' @importFrom rtracklayer readGFF
#' @importFrom data.table as.data.table :=
#' @keywords internal
load_gtf_long <- function(gtf_path_or_df, add_row_uid = TRUE) {
  # --- decide how to read ----------------------------------------------------
  if (is.character(gtf_path_or_df) && length(gtf_path_or_df) == 1L) {
    p <- sub("^file://", "", gtf_path_or_df)    # allow file:// URLs
    gtf <- rtracklayer::readGFF(p)

  } else if (inherits(gtf_path_or_df, c("data.table","data.frame"))) {
    gtf <- gtf_path_or_df

  } else {
    stop("`gtf_path_or_df` must be a single character path/URL or a data.frame/data.table.")
  }

  # rtracklayer returns S4 DataFrame; coerce to data.frame -> data.table
  if (inherits(gtf, "DataFrame")) gtf <- as.data.frame(gtf)
  gtf <- data.table::as.data.table(gtf)

  # --- keep only known columns (whatever exists) -----------------------------
  keep <- intersect(
    c("seqid","source","type","start","end","score","strand","phase","attributes",
      "gene_id","gene_type","gene_name","transcript_id","transcript_type","transcript_name",
      "exon_id","exon_number","protein_id","tag","level","transcript_support_level"),
    names(gtf)
  )
  gtf <- gtf[, ..keep]

  # Ensure attributes is character for parsing helpers
  if ("attributes" %in% names(gtf) && !is.character(gtf$attributes)) {
    gtf[, attributes := as.character(attributes)]
  }

  # --- ensure IDs (fallback via attributes) ----------------------------------
  if (!"gene_id"       %in% names(gtf)) gtf[, gene_id       := attr_get(attributes, "gene_id",       TRUE)] else gtf[, gene_id       := strip_ver(gene_id)]
  if (!"transcript_id" %in% names(gtf)) gtf[, transcript_id := attr_get(attributes, "transcript_id", TRUE)] else gtf[, transcript_id := strip_ver(transcript_id)]
  if (!"exon_id"       %in% names(gtf)) gtf[, exon_id       := attr_get(attributes, "exon_id",       TRUE)] else gtf[, exon_id       := strip_ver(exon_id)]
  if (!"protein_id"    %in% names(gtf)) gtf[, protein_id    := attr_get(attributes, "protein_id",    TRUE)] else gtf[, protein_id    := strip_ver(protein_id)]
  if (!"exon_number"   %in% names(gtf)) gtf[, exon_number   := suppressWarnings(as.integer(attr_get(attributes, "exon_number")))]

  # --- row uid + column order ------------------------------------------------
  if (isTRUE(add_row_uid)) {
    gtf[, row_uid := .I]
    data.table::setcolorder(gtf, c("row_uid", setdiff(names(gtf), "row_uid")))
  }

  colnames(gtf)[which(colnames(gtf) == 'seqid')] <- 'chr'
  gtf[]
}


#' Build a compact transcript protein sequence map (internal)
#'
#' Internal helper that integrates information from a GTF, transcript FASTA,
#' and protein FASTA to produce a compact mapping of genes, transcripts,
#' and proteins, optionally including nucleotide and amino-acid sequences.
#'
#' @param gtf_df A \code{data.frame} or \code{data.table} produced by
#'   [load_gtf_long()], containing at least \code{gene_id},
#'   \code{transcript_id}, \code{protein_id}, and optionally
#'   \code{row_uid}.
#' @param txfa_path Character string giving the path to the GENCODE
#'   transcript FASTA file (e.g. \code{gencode.v45.pc_transcripts.fa.gz}).
#' @param aafa_path Character string giving the path to the GENCODE
#'   protein FASTA file (e.g. \code{gencode.v45.pc_translations.fa.gz}).
#' @param take_tx_cds_slice Logical; if \code{TRUE}, restricts transcript
#'   sequences to their coding region based on the \code{CDS:start-end}
#'   annotation in the FASTA header.
#' @param keep_sequences Logical; if \code{TRUE}, include sequence strings
#'   in the output table, otherwise store \code{NA_character_}.
#' @param add_row_uids Logical; if \code{TRUE}, attaches corresponding
#'   \code{row_uid}s for genes and transcripts from \code{gtf_df}.
#'
#' @return A \code{data.table} containing, for each transcript:
#'   \describe{
#'     \item{row_uid}{Row index within the output mapping.}
#'     \item{gene_id, transcript_id, protein_id}{Identifiers from GENCODE.}
#'     \item{transcript_seq, protein_seq}{Optional sequence strings.}
#'     \item{gene_row_uid, transcript_row_uid}{Original UID references (if added).}
#'   }
#'
#' @details
#' The function combines identifiers from GTF features, transcript FASTA
#' headers, and protein FASTA headers to construct a unified mapping.
#' It uses \pkg{Biostrings} to load FASTA data and
#' \pkg{data.table} for efficient joins.
#' @importFrom data.table as.data.table data.table rbindlist setcolorder :=
#' @importFrom Biostrings readAAStringSet DNAStringSet subseq
#' @importFrom stats na.omit
#' @keywords internal
load_seq_map <- function(gtf_df, txfa_path, aafa_path,
                         take_tx_cds_slice = TRUE,
                         keep_sequences    = TRUE,
                         add_row_uids      = TRUE) {
  gtf_dt <- data.table::as.data.table(gtf_df)

  if (add_row_uids && !"row_uid" %in% names(gtf_dt)) {
    stop("gtf_df must include a 'row_uid' column to attach gene/transcript UIDs. ",
         "Run load_gtf_long(..., add_row_uid = TRUE) first.")
  }

  # Map from GTF CDS rows (best source for protein_id)
  cds_map <- unique(na.omit(gtf_dt[type == "CDS", .(gene_id, transcript_id, protein_id)]))
  tx_gene <- unique(na.omit(gtf_dt[type == "transcript", .(transcript_id, gene_id)]))

  AA <- Biostrings::readAAStringSet(aafa_path)
  hdrs <- names(AA)
  protein_id <- strip_ver(sub("^([^|]+).*", "\\1", hdrs))
  transcript_id <- vapply(hdrs, prot_hdr_to_enst, character(1))
  prot_dt <- data.table::data.table(
    protein_id    = protein_id,
    transcript_id = transcript_id,
    protein_seq   = if (keep_sequences) as.character(AA) else NA_character_
  )
  prot_map <- unique(prot_dt[, .(transcript_id, protein_id)])

  # ---- Merge maps and attach gene_id ----
  tx_prot  <- unique(data.table::rbindlist(list(cds_map[, .(transcript_id, protein_id)],
                                                prot_map), use.names = TRUE, fill = TRUE))
  map_all  <- merge(tx_prot, tx_gene, by = "transcript_id", all.x = TRUE)
  data.table::setcolorder(map_all, c("gene_id","transcript_id","protein_id"))

  # ---- Optional: attach transcript_row_uid and gene_row_uid from gtf_df ----
  if (add_row_uids) {
    tmap <- gtf_dt[type == "transcript" & !is.na(transcript_id),
                   .(transcript_row_uid = min(row_uid)), by = transcript_id]
    gmap <- gtf_dt[type == "gene" & !is.na(gene_id),
                   .(gene_row_uid = min(row_uid)), by = gene_id]
    map_all <- merge(map_all, tmap, by = "transcript_id", all.x = TRUE)
    map_all <- merge(map_all, gmap, by = "gene_id",       all.x = TRUE)
  }

  # ---- Protein sequences (optional) ----
  protein_seq_df <- if (keep_sequences) unique(prot_dt[!is.na(protein_seq),
                                                       .(protein_id, protein_seq)]) else data.table::data.table()
  need_enst <- unique(map_all$transcript_id)
  transcript_seq_df <- data.table::data.table()

  if (length(need_enst)) {
    # Load transcript FASTA
    TX <- Biostrings::readDNAStringSet(txfa_path)

    hdrs <- names(TX)
    enst_ids <- strip_ver(sub("^([^|]+).*", "\\1", hdrs))

    # Keep only transcripts we need
    keep_idx <- enst_ids %in% need_enst
    if (any(keep_idx)) {
      TX <- TX[keep_idx]
      enst_ids <- enst_ids[keep_idx]

      # Optional CDS trimming (header may contain CDS:start-end)
      if (isTRUE(take_tx_cds_slice)) {
        cds_bounds <- regexec("\\bCDS:([0-9]+)-([0-9]+)", hdrs[keep_idx], perl = TRUE)
        cds_matches <- regmatches(hdrs[keep_idx], cds_bounds)

        TX <- Biostrings::DNAStringSet(mapply(function(seq, match) {
          if (length(match) >= 3L) {
            start <- as.integer(match[2])
            end   <- as.integer(match[3])
            if (!is.na(start) && !is.na(end) &&
                start >= 1L && end <= nchar(as.character(seq)) && end >= start) {
              return(Biostrings::subseq(seq, start = start, end = end))
            }
          }
          seq
        }, seq = TX, match = cds_matches, SIMPLIFY = FALSE))
      }

      # Build data.table of transcript sequences
      transcript_seq_df <- data.table::data.table(
        transcript_id  = enst_ids,
        transcript_seq = if (isTRUE(keep_sequences))
          as.character(TX)
        else
          NA_character_
      )
    }
  }

  # ---- Final sequence-map table ----
  seq_map_df <- merge(
    map_all,
    transcript_seq_df,
    by = "transcript_id",
    all.x = TRUE
  )


  if (nrow(protein_seq_df)) seq_map_df <- merge(seq_map_df, protein_seq_df, by = "protein_id", all.x = TRUE)

  # stable uid for this table (independent of gtf_df row_uids)
  seq_map_df[, row_uid := .I]

  # nice column order (row_uids first if present)
  front <- c("row_uid",
             intersect(c("gene_row_uid","transcript_row_uid"), names(seq_map_df)),
             "gene_id","transcript_id","protein_id","transcript_seq","protein_seq")
  rest  <- setdiff(names(seq_map_df), front)
  data.table::setcolorder(seq_map_df, c(front, rest))
  seq_map_df[, `:=`(row_uid = transcript_row_uid, transcript_row_uid = NULL)][]
  seq_map_df[]
}

#' Restrict GTF entries by gene type (internal)
#'
#' Internal helper to subset a GTF \code{data.frame} or \code{data.table}
#' to include only selected gene biotypes (e.g. protein_coding, lncRNA).
#'
#' @param gtf_df A \code{data.frame} or \code{data.table} containing a
#'   \code{gene_type} column.
#' @param restrictions Character vector of allowed \code{gene_type}
#'   values. Defaults to common gene biotypes.
#'
#' @return A subset of \code{gtf_df} containing only rows whose
#'   \code{gene_type} matches one of the specified \code{restrictions}.
#' @keywords internal
restrict_gtf_genetype <- function(gtf_df,
                                  restrictions = c('protein_coding', 'lncRNA',
                                                  'snoRNA', 'snRNA', 'rRNA',
                                                  'miRNA', 'MT_tRNA', 'MT_rRNA')) {
  return(gtf_df[gtf_df$gene_type %in% restrictions,])
}

#' Restrict GTF entries by feature type (internal)
#'
#' Internal helper to subset a GTF \code{data.frame} or \code{data.table}
#' to retain only \code{gene}, \code{exon}, and \code{transcript} rows.
#'
#' @param gtf_df A \code{data.frame} or \code{data.table} containing a
#'   \code{type} column (e.g. output of [load_gtf_long()]).
#'
#' @return A subset of \code{gtf_df} including only rows where
#'   \code{type} is one of \code{"gene"}, \code{"exon"}, or
#'   \code{"transcript"}.
#'
#' @keywords internal
restrict_gtf_rowtype <- function(gtf_df) {
  return(gtf_df[gtf_df$type %in% c('gene', 'exon', 'transcript'),])
}




#' Add per-exon coding and feature coordinates (internal)
#'
#' Internal helper to annotate each exon in a GTF \code{data.table} with
#' strand-aware coding information. For every exon, the function records
#' whether it overlaps coding sequence (CDS), untranslated region (UTR),
#' start codon, or stop codon, together with absolute, genomic, and
#' transcript-relative coordinates.
#'
#' @param gtf_df A \code{data.frame} or \code{data.table} containing GTF
#'   annotations, typically produced by [load_gtf_long()]. Must include
#'   columns \code{type}, \code{start}, \code{end}, \code{strand},
#'   \code{transcript_id}, \code{exon_id}, and a unique \code{row_uid}.
#'
#' @return A \code{data.table} identical to the input \code{gtf_df} but
#'   with additional per-exon fields:
#'   \itemize{
#'     \item \code{cds_has}, \code{utr_has}, etc logical indicators
#'     \item Absolute (within-exon) start/stop coordinates (\code{_abs_*})
#'     \item Genomic coordinates (\code{_gen_*})
#'     \item Transcript-relative coordinates (\code{_rel_*})
#'     \item Feature lengths (\code{_len})
#'   }
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Orders exons per transcript, respecting strand.
#'   \item Computes exon-wise coordinates for CDS, UTR, start/stop codons.
#'   \item Merges these annotations back into the exon-level GTF table.
#' }
#' It relies on \pkg{data.table} joins and by-group operations for speed.
#' @importFrom data.table setDF := setnames
#' @keywords internal
add_exon_coding_information <- function(gtf_df) {
  data.table::setDT(gtf_df)

  # ---------- 0) Strand-aware absolute exon order per transcript ----------
  ex <- gtf_df[type == "exon",
               .(row_uid, transcript_id, exon_id, chr, strand,
                 exon_start = start, exon_end = end)]

  # order: + strand by start asc; - strand by start desc
  ex[order(transcript_id,
           ifelse(strand=="+", exon_start, -exon_start),
           ifelse(strand=="+", exon_end,   -exon_end)),
     absolute_exon_position := seq_len(.N), by = transcript_id]

  # ---------- helper to pull feature rows keyed by exon ----------
  feat_by_exon <- function(t) {
    dt <- gtf_df[type == t, .(transcript_id, exon_id, fstart = start, fend = end)]
    dt[!is.na(exon_id)]
  }

  # ---------- core annotator: adds abs/gen/rel + per-exon length for one feature kind ----------
  annotate_feature <- function(ex, feat_dt, prefix) {
    # empty (no rows of this feature)
    if (!nrow(feat_dt)) {
      out <- ex[, .(row_uid, transcript_id, absolute_exon_position)]
      out[, (paste0(prefix, "_has")) := FALSE]
      for (nm in c("_abs_start","_abs_stop","_gen_start","_gen_stop","_rel_start","_rel_stop","_len"))
        out[, (paste0(prefix, nm)) := NA_integer_]
      # lengths for empty = 0 not NA (so they sum nicely)
      out[, (paste0(prefix, "_len")) := 0L]
      return(out[])
    }

    # join by (transcript_id, exon_id), allow multiple fragments
    dt <- merge(ex, feat_dt, by = c("transcript_id","exon_id"), allow.cartesian = TRUE)

    # intersect with exon bounds (defensive)
    dt[, `:=`(
      gstart   = pmax(exon_start, fstart),
      gend     = pmin(exon_end,   fend)
    )]
    dt <- dt[gend >= gstart]

    # exon-relative 1-based coords + fragment length
    dt[, `:=`(
      rstart   = gstart - exon_start + 1L,
      rend     = gend   - exon_start + 1L,
      frag_len = gend - gstart + 1L
    )]

    # per exon: outer span (abs), genomic min/max, and TRUE feature length (sum of fragments)
    exon_agg <- dt[, .(
      abs_start = min(rstart),           # exon-relative outer start
      abs_stop  = max(rend),             # exon-relative outer stop
      gen_start = min(gstart),           # genomic outer start
      gen_stop  = max(gend),             # genomic outer stop
      feat_len  = sum(frag_len)          # true length inside this exon
    ), by = .(transcript_id, row_uid)]

    # ensure every exon is present; fill NAs
    out <- merge(ex[, .(row_uid, transcript_id, absolute_exon_position)], exon_agg,
                 by = c("transcript_id","row_uid"), all.x = TRUE)
    out[, has := !is.na(abs_start)]
    for (c in c("abs_start","abs_stop","gen_start","gen_stop")) out[is.na(get(c)), (c) := NA_integer_]
    out[is.na(feat_len), feat_len := 0L]

    # transcript-relative stacked coords:
    # cumulative sum of prior feature lengths in transcript order (feature-only counting)
    out[order(transcript_id, absolute_exon_position),
        rel_offset := cumsum(data.table::shift(feat_len, fill = 0L)), by = transcript_id]

    # relative coords count ONLY feature bases; first base across transcript is 1
    out[has == TRUE, `:=`(
      rel_start = rel_offset + 1L,
      rel_stop  = rel_offset + feat_len
    )]
    out[has == FALSE, `:=`(rel_start = NA_integer_, rel_stop = NA_integer_)]

    # finalize names
    nm <- function(s) paste0(prefix, s)
    out <- out[, .(row_uid,
                   has_col = has,
                   abs_start, abs_stop, gen_start, gen_stop,
                   rel_start, rel_stop,
                   feat_len)]
    data.table::setnames(out,
             c("has_col","abs_start","abs_stop","gen_start","gen_stop","rel_start","rel_stop","feat_len"),
             c(nm("_has"), nm("_abs_start"), nm("_abs_stop"),
               nm("_gen_start"), nm("_gen_stop"),
               nm("_rel_start"), nm("_rel_stop"),
               nm("_len")))
    out[]
  }

  # ---------- 1) build feature tables ----------
  cds_dt <- feat_by_exon("CDS")
  utr_dt <- feat_by_exon("UTR")
  sc_dt  <- feat_by_exon("start_codon")
  tc_dt  <- feat_by_exon("stop_codon")

  # ---------- 2) annotate each feature kind ----------
  cds_anno <- annotate_feature(ex, cds_dt, "cds")
  utr_anno <- annotate_feature(ex, utr_dt, "utr")
  sc_anno  <- annotate_feature(ex, sc_dt,  "start_codon")
  tc_anno  <- annotate_feature(ex, tc_dt,  "stop_codon")

  # ---------- 3) merge back to exon rows in gtf_df ----------
  ex_anno <- Reduce(function(a,b) merge(a,b, by = "row_uid", all = TRUE),
                    list(cds_anno, utr_anno, sc_anno, tc_anno))

  setkey(gtf_df, row_uid)
  setkey(ex_anno, row_uid)
  gtf_df[ex_anno, `:=`(
    cds_has            = i.cds_has,
    cds_abs_start      = i.cds_abs_start,
    cds_abs_stop       = i.cds_abs_stop,
    cds_gen_start      = i.cds_gen_start,
    cds_gen_stop       = i.cds_gen_stop,
    cds_rel_start      = i.cds_rel_start,
    cds_rel_stop       = i.cds_rel_stop,
    cds_len            = i.cds_len,        # <-- NEW

    utr_has            = i.utr_has,
    utr_abs_start      = i.utr_abs_start,
    utr_abs_stop       = i.utr_abs_stop,
    utr_gen_start      = i.utr_gen_start,
    utr_gen_stop       = i.utr_gen_stop,
    utr_rel_start      = i.utr_rel_start,
    utr_rel_stop       = i.utr_rel_stop,
    utr_len            = i.utr_len,        # <-- NEW

    start_codon_has        = i.start_codon_has,
    start_codon_abs_start  = i.start_codon_abs_start,
    start_codon_abs_stop   = i.start_codon_abs_stop,
    start_codon_gen_start  = i.start_codon_gen_start,
    start_codon_gen_stop   = i.start_codon_gen_stop,
    start_codon_rel_start  = i.start_codon_rel_start,
    start_codon_rel_stop   = i.start_codon_rel_stop,

    stop_codon_has         = i.stop_codon_has,
    stop_codon_abs_start   = i.stop_codon_abs_start,
    stop_codon_abs_stop    = i.stop_codon_abs_stop,
    stop_codon_gen_start   = i.stop_codon_gen_start,
    stop_codon_gen_stop    = i.stop_codon_gen_stop,
    stop_codon_rel_start   = i.stop_codon_rel_start,
    stop_codon_rel_stop    = i.stop_codon_rel_stop
  )]

  return(gtf_df[])
}


#' Add exon order and positional classification (internal)
#'
#' Internal helper that annotates exon rows within a GTF table with
#' strand-aware order and coding position information.  Adds both absolute
#' (by transcript) and coding-region-specific ordering, along with
#' classification labels such as \code{"first"}, \code{"internal"},
#' \code{"last"}, or \code{"single_exon"}.
#'
#' @param gtf_df A \code{data.frame} or \code{data.table} containing GTF
#'   annotations (usually from [add_exon_coding_information()]).  Must
#'   include columns \code{type}, \code{start}, \code{end},
#'   \code{strand}, \code{transcript_id}, \code{exon_number},
#'   \code{cds_has}, and \code{cds_rel_start}.
#'
#' @return A \code{data.table} identical to \code{gtf_df} but with added
#'   columns:
#'   \itemize{
#'     \item \code{absolute_exon_position} exon index by transcript
#'     \item \code{coding_exon_position} exon index among coding exons
#'     \item \code{absolute_exon_class} positional label for all exons
#'     \item \code{coding_exon_class} positional label for coding exons
#'   }
#'
#' @details
#' Exons are ordered by genomic position within each transcript, taking
#' strand into account.  If \code{exon_number} is provided, it is used to
#' set the absolute order; otherwise order is inferred from coordinates.
#'
#' @importFrom data.table := frank fifelse
#' @keywords internal
add_exon_order_information <- function(gtf_df) {
  data.table::setDT(gtf_df)

  # Work on exon rows only
  ex <- gtf_df[type == "exon",
               .(row_uid, transcript_id, exon_id, strand, start, end,
                 exon_number = get("exon_number"),
                 cds_has = get("cds_has"),
                 cds_rel_start = get("cds_rel_start"))]

  # 1) absolute_exon_position
  if (!all(is.na(ex$exon_number))) {
    ex[, absolute_exon_position := as.integer(exon_number)]
    # fill any remaining NAs from strand-aware genomic order
    if (anyNA(ex$absolute_exon_position)) {
      ex[is.na(absolute_exon_position)][
        order(transcript_id,
              ifelse(strand=="+", start, -start),
              ifelse(strand=="+", end,   -end)),
        absolute_exon_position := seq_len(.N), by = transcript_id]
    }
  } else {
    ex[order(transcript_id,
             ifelse(strand=="+", start, -start),
             ifelse(strand=="+", end,   -end)),
       absolute_exon_position := seq_len(.N), by = transcript_id]
  }

  # ensure logical cds_has (default FALSE if missing)
  ex[is.na(cds_has), cds_has := FALSE]

  # 2) coding_exon_position using cds_rel_start (ascending) within transcript
  # initialize to -1 for all exons
  ex[, coding_exon_position := -1L]

  # rank only coding exons that have a non-NA cds_rel_start
  ex[cds_has == TRUE & !is.na(cds_rel_start),
     coding_exon_position := data.table::frank(cds_rel_start, ties.method = "first"),
     by = transcript_id]

  # 3) absolute_exon_class with single-exon handling
  ex[, n_abs := .N, by = transcript_id]
  ex[, absolute_exon_class :=
       data.table::fifelse(n_abs == 1L, "single_exon",
               fifelse(absolute_exon_position == 1L, "first",
                       fifelse(absolute_exon_position == n_abs, "last", "internal")))]

  # 4) coding_exon_class with single-coding-exon handling
  ex[, n_coding := sum(cds_has == TRUE & !is.na(cds_rel_start)), by = transcript_id]
  ex[cds_has == FALSE | is.na(cds_rel_start), coding_exon_class := "noncoding"]
  ex[cds_has == TRUE  & !is.na(cds_rel_start) & n_coding == 1L, coding_exon_class := "single_exon"]
  ex[cds_has == TRUE  & !is.na(cds_rel_start) & n_coding > 1L,
     coding_exon_class :=
       data.table::fifelse(coding_exon_position == 1L,         "first",
               fifelse(coding_exon_position == n_coding,   "last",  "internal"))]

  # write back by row_uid
  setkey(gtf_df, row_uid)
  setkey(ex,     row_uid)
  gtf_df[ex, `:=`(
    absolute_exon_position = i.absolute_exon_position,
    coding_exon_position   = i.coding_exon_position,
    absolute_exon_class    = i.absolute_exon_class,
    coding_exon_class      = i.coding_exon_class
  )]

  gtf_df[]
}

#' Add exon count per transcript (internal)
#'
#' Internal helper that computes the number of exon entries per transcript
#' in a GTF table and attaches it as a new column.
#'
#' @param gtf_df A \code{data.frame} or \code{data.table} containing at
#'   least \code{type == "exon"} and \code{transcript_id} columns.
#' @param col Character string giving the name of the output column for
#'   the exon count (default \code{"n_exons"}).
#'
#' @return A \code{data.table} identical to \code{gtf_df} but with one
#'   additional column containing the number of exons per transcript.
#'
#' @importFrom data.table setnames
#' @keywords internal
add_exon_count_per_transcript <- function(gtf_df, col = "n_exons") {
  cnt <- gtf_df[type == "exon" & !is.na(transcript_id),
                .(N = .N), by = transcript_id]
  data.table::setnames(cnt, "N", col)

  res <- cnt[gtf_df, on = "transcript_id"]

  res[]
}

#' Add feature length column (internal)
#'
#' Internal helper that computes the genomic length of each feature in a
#' GTF table and appends it as \code{feature_length}.
#'
#' @param gtf_df A \code{data.frame} or \code{data.table} containing
#'   \code{start} and \code{end} columns.
#'
#' @return A \code{data.table} identical to the input but with an added
#'   \code{feature_length} column (\code{end - start + 1}).
#'
#' @importFrom data.table :=
#' @keywords internal
add_feature_length <- function(gtf_df) {
  return(gtf_df[, feature_length := abs(gtf_df$start - gtf_df$end)+1])
}

#' Identify potential hybrid exons (internal)
#'
#' Internal helper to find exons that overlap between internal and
#' terminal (first/last) exons of different transcripts within the same
#' gene. Used to detect possible hybrid exon configurations.
#'
#' @param gtf_df A \code{data.frame} or \code{data.table} containing GTF
#'   annotations with \code{gene_id}, \code{transcript_id},
#'   \code{exon_id}, \code{chr}, \code{start}, \code{end}, and
#'   \code{absolute_exon_class}.
#'
#' @return A list of two \code{data.table}s:
#'   \itemize{
#'     \item \code{first_hybrids} overlaps between first and internal exons
#'     \item \code{last_hybrids} overlaps between last and internal exons
#'   }
#'   Each includes transcript IDs, exon IDs, and mapped transcript
#'   \code{row_uid}s.
#'
#' @details
#' Uses \pkg{data.table::foverlaps()} to identify exons that share
#' genomic coordinates but belong to different transcripts within the
#' same gene.
#'
#' @importFrom data.table foverlaps setkey as.data.table
#' @keywords internal
identify_hybrid_exons_split <- function(gtf_df) {
  DT <- data.table::as.data.table(gtf_df)

  tx_uid <- DT[type == "transcript", .(transcript_id, transcript_uid = row_uid)]
  data.table::setkey(tx_uid, transcript_id)

  ex <- DT[type == "exon" & !is.na(gene_id) & !is.na(transcript_id) & !is.na(exon_id),
           .(gene_id, chr, start, end, transcript_id, exon_id, absolute_exon_class)]

  ex_internal <- copy(ex[absolute_exon_class == "internal"])
  ex_first    <- copy(ex[absolute_exon_class == "first"])
  ex_last     <- copy(ex[absolute_exon_class == "last"])

  # helper: overlap terminal (x) vs internal (i) within same gene+seqid
  overlap_pairs <- function(terminal_dt, internal_dt) {
    if (!nrow(terminal_dt) || !nrow(internal_dt)) {
      return(data.table(
        gene_id = character(), transcript_id_terminal = character(),
        transcript_id_internal = character(),
        exon_id_terminal = character(), exon_id_internal = character()
      ))
    }
    data.table::setkey(terminal_dt, gene_id, chr, start, end)
    data.table::setkey(internal_dt, gene_id, chr, start, end)

    ov <- data.table::foverlaps(terminal_dt, internal_dt,
                    by.x = c("gene_id","chr","start","end"),
                    by.y = c("gene_id","chr","start","end"),
                    type = "any", nomatch = 0L)

    # x.* are bare; i.* are from internal table
    out <- ov[
      transcript_id != i.transcript_id,
      .(gene_id,
        transcript_id_terminal = transcript_id,
        transcript_id_internal = i.transcript_id,
        exon_id_terminal       = exon_id,
        exon_id_internal       = i.exon_id)
    ]
    unique(out)
  }

  # compute separately
  first_hybrids <- overlap_pairs(ex_first, ex_internal)
  last_hybrids  <- overlap_pairs(ex_last,  ex_internal)

  add_uids <- function(df, tx_uid) {
    df[tx_uid, on = .(transcript_id_terminal = transcript_id),
       transcript_uid_terminal := transcript_uid]
    df[tx_uid, on = .(transcript_id_internal = transcript_id),
       transcript_uid_internal := transcript_uid]
    return(df)
  }
  return(list(
    first_hybrids = add_uids(first_hybrids, tx_uid),
    last_hybrids  = add_uids(last_hybrids, tx_uid)
  ))
}



#' Resolve and create an output directory
#'
#' Determines where to write package output files, supporting both
#' session-temporary and persistent cache locations. By default,
#' temporary directories are created under [tempdir()], while persistent
#' ones follow Bioconductors cache convention via
#' [tools::R_user_dir()].
#'
#' @param dir Optional user-supplied path. If provided, this path is
#'   created (recursively) if necessary and returned as an absolute,
#'   normalized directory path.
#' @param persist One of \code{"temp"} or \code{"cache"}. When
#'   \code{"temp"}, a session-temporary directory under [tempdir()] is
#'   used; when \code{"cache"}, a persistent user cache directory
#'   following [tools::R_user_dir()] is used.
#' @param pkg Package name used to determine the cache location when
#'   \code{persist = "cache"}. Defaults to \code{"SpliceImpactR"}.
#'
#' @return A character scalar giving an absolute, writable directory path.
#'   The directory is created if it does not exist.
#'
#'   @importFrom tools R_user_dir
#' @keywords internal
resolve_output_dir <- function(dir = NULL,
                               persist = c("temp","cache")[1],
                               pkg = "SpliceImpactR") {
  persist <- match.arg(persist)

  if (!is.null(dir)) {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    if (!dir.exists(dir)) stop("Could not create dir: ", dir)
    return(normalizePath(dir, winslash = "/", mustWork = TRUE))
  }

  if (persist == "temp") {
    d <- tempfile(pattern = paste0(pkg, "_"), tmpdir = tempdir())
  } else { # persist == "cache"
    d <- tools::R_user_dir(pkg, "cache")
  }

  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  normalizePath(d, winslash = "/", mustWork = TRUE)
}

#' Retrieve transcript and protein sequences (internal)
#'
#' Internal helper that calls [load_seq_map()] to obtain transcript
#' and protein sequences for a given GTF annotation and corresponding
#' FASTA files. Used within higher-level SpliceImpactR workflows.
#'
#' @param gtf_df A \code{data.frame} or \code{data.table} containing
#'   GTF annotation data, typically from [load_gtf_long()].
#' @param transcript_path Path to the transcript FASTA file
#'   (e.g. \code{gencode.v45.pc_transcripts.fa.gz}).
#' @param translation_path Path to the protein translation FASTA file
#'   (e.g. \code{gencode.v45.pc_translations.fa.gz}).
#'
#' @return A \code{data.table} with gene, transcript, and protein IDs
#'   and their corresponding nucleotide and amino acid sequences.
#'
#' @keywords internal
get_sequences <- function(gtf_df, transcript_path, translation_path) {
  seq_map_df <- load_seq_map(
    gtf_df,
    txfa_path = transcript_path,
    aafa_path = translation_path,
    take_tx_cds_slice = TRUE,
    keep_sequences    = TRUE
  )
  return(seq_map_df)
}

#' Add frames to exons dependent on cds location
#'
#' Internal helper
#'
#' @param gtf_df A \code{data.frame} or \code{data.table} containing
#'   GTF annotation data, typically from [load_gtf_long()].
#'
#' @return A \code{data.table} with added start_frame and stop_frame
#'
#' @keywords internal
add_exon_frames <- function(gtf_df) {
  data.table::setDT(gtf_df)

  # compute frames only for exon rows that have CDS
  idx <- (gtf_df$type == "exon") & gtf_df$cds_has
  # guard NAs as well
  idx <- idx & !is.na(gtf_df$cds_rel_start) & !is.na(gtf_df$cds_rel_stop)

  # frame is (position - 1) %% 3, with the first coding nt of a transcript = frame 0
  gtf_df[idx, `:=`(
    start_frame = as.integer((cds_rel_start - 1L) %% 3L),
    stop_frame  = as.integer((cds_rel_stop  - 1L) %% 3L)
  )]

  return(gtf_df[])
}


#' Load and cache GENCODE annotations, sequences, and hybrid exon annotations
#'
#' This function loads GENCODE gene models (GTF), processes exon annotations,
#' extracts transcript and protein sequences, identifies hybrid exons, and
#' optionally caches the processed objects for future fast access.
#'
#' Four loading modes are supported:
#' \describe{
#'   \item{`test`}{Load small internal test data shipped with the package.}
#'   \item{`cached`}{Load previously cached GTF, transcript sequences, protein sequences, and hybrid tables from `base_dir`.}
#'   \item{`path`}{Read local GTF and FASTA files, process, then cache results in `base_dir`.}
#'   \item{`link`}{Download GENCODE files from URLs, process, then cache results in `base_dir`.}
#' }
#'
#' Cached files created in `base_dir`:
#' \preformatted{
#' {species}_gencode_v{release}.gtf.rds
#' {species}_gencode_v{release}_transcripts.rds
#' {species}_gencode_v{release}_proteins.rds
#' {species}_gencode_v{release}_hybrids.rds
#' }
#'
#' @param load Character string specifying load mode:
#'   one of `"link"`, `"path"`, `"cached"`, `"test"`.
#' @param base_dir Directory to cache/load processed objects.
#' @param species Species label used in filenames (default `"human"`).
#' @param release GENCODE release version (default `45`).
#' @param gtf_path Path to a GTF file when `load = "path"`.
#' @param transcript_path Path to transcript FASTA (.fa/.fa.gz) when `load = "path"`.
#' @param translation_path Path to protein FASTA (.fa/.fa.gz) when `load = "path"`.
#' @param filter_tsl Transcript support levels to retain (default `c("1","2","3")`).
#'   Transcripts outside this set are dropped unless the row is a gene record.
#'
#' @return A list with:
#' \describe{
#'   \item{`annotations`}{Processed long-format GTF as `data.table`}
#'   \item{`sequences`}{List with elements `transcripts` and `proteins` (or `NULL` if not loaded)}
#'   \item{`hybrids`}{Hybrid exon annotation list}
#' }
#'
#'
#' @examples
#' # Load bundled test data
#' ann <- get_annotation(load = "test")
#'
#' #
#' #  # Load from local files and cache results
#' #  ann <- get_annotation(
#' #    load = "path",
#' #    base_dir = "/project/annotation_cache/",
#' #    gtf_path = "gencode.v45.annotation.gtf.gz",
#' #    transcript_path = "gencode.v45.pc_transcripts.fa.gz",
#' #    translation_path = "gencode.v45.pc_translations.fa.gz"
#' #  )
#'
#' #  # Download files, process, and cache
#' #  ann <- get_annotation(
#' #    load = "link",
#' #    base_dir = "/project/annotation_cache/"
#' #  )
#'
#' # Load from cached RDS (fast)
#' #  ann <- get_annotation(
#' #    load = "cached",
#' #    base_dir = "/project/annotation_cache/"
#' #  )
#'
#' @importFrom magrittr %>%
#' @export
get_annotation <- function(
    load = c("link", "path", "cached", "test"),
    base_dir = NULL,
    species = c("human", 'mouse')[1],
    release = 45,
    gtf_path = NULL,
    transcript_path = NULL,
    translation_path = NULL,
    filter_tsl = c('1','2','3', '4', '5')[c(1, 2, 3)]
) {
  load <- match.arg(load)
  base_dir <- resolve_output_dir(base_dir)

  ### ---- Define cache files ----
  rds_gtf  <- file.path(base_dir, sprintf("%s_gencode_v%s.gtf.rds", species, release))
  rds_seq   <- file.path(base_dir, sprintf("%s_gencode_v%s_sequences.rds", species, release))
  rds_hyb  <- file.path(base_dir, sprintf("%s_gencode_v%s_hybrids.rds", species, release))

  message("[INFO] get_annotation mode: ", load)

  ### ---- TEST MODE ----
  if (load == "test") {
    message("[INFO] Loading bundled test annotation data")
    return(list(
      annotations = fread(get_example_data("human_test_gencode_v45_annotations.csv")),
      sequences = fread(get_example_data("human_test_gencode_v45_sequences.csv")),
      hybrids = list(
        first_hybrids = fread(get_example_data("human_test_gencode_v45_first_hybrids.csv")),
        last_hybrids = fread(get_example_data("human_test_gencode_v45_last_hybrids.csv"))
      )
    ))
  }

  ### ---- CACHED MODE ----
  if (load == "cached") {
    message("[FAST] Loading cached annotation objects from: ", base_dir)

    req <- c(rds_gtf, rds_seq, rds_hyb)
    if (!all(file.exists(req))) {
      missing <- req[!file.exists(req)]
      msg <- paste0(
        "Missing cached file(s):\n",
        paste(missing, collapse = "\n"),
        "\nRun with load='path' or load='link' first to generate caches."
      )
      stop(msg, call. = FALSE)
    }

    return(list(
      annotations = readRDS(rds_gtf),
      sequences = readRDS(rds_seq),
      hybrids = readRDS(rds_hyb)
    ))
  }

  ### ---- LINK / PATH: Acquire GTF/FASTA paths ----
  if (load == "link") {
    message("[STEP] Downloading GENCODE assets (first run only)")
    assets <- .si_prepare_assets(base_dir, species, release, mode = "download")
    gtf_file  <- assets$paths$gtf_gz
    tx_fa     <- transcript_path  %||% assets$paths$txfa_gz
    aa_fa     <- translation_path %||% assets$paths$aafa_gz
  } else if (load == "path") {
    if (is.null(gtf_path))
      stop("load='path' requires gtf_path=")

    gtf_file <- gtf_path
    tx_fa <- transcript_path
    aa_fa <- translation_path

  }

  message("[STEP] Reading GTF: ", gtf_file)
  gtf_dt <- load_gtf_long(gtf_file)

  ### ---- Annotation pipeline ----
  message("[STEP] Annotating GTF")
  gtf_df <- gtf_dt %>%
    restrict_gtf_genetype %>%
    add_exon_coding_information %>%
    add_exon_order_information %>%
    restrict_gtf_rowtype %>%
    add_exon_frames %>%
    add_feature_length

  gtf_df <- gtf_df[
    transcript_support_level %in% filter_tsl | type == 'gene'
  ][
    !(tag %in% c("cds_start_NF", "cds_end_NF")) | type == 'gene'
  ]

  ### ---- Hybrid exons ----
  message("[STEP] Identifying hybrid exons")
  hybrids <- identify_hybrid_exons_split(gtf_df)

  ### ---- Sequence loading ----
  message("[STEP] Loading transcript + protein sequences")
  # tx_local <- if (grepl("^https?://", tx_fa)) cached_resource(tx_fa) else tx_fa
  # aa_local <- if (grepl("^https?://", aa_fa)) cached_resource(aa_fa) else aa_fa
  tx_local <- tx_fa
  aa_local <- aa_fa

  if (!file.exists(tx_local)) stop("Transcript FASTA missing: ", tx_local)
  if (!file.exists(aa_local)) stop("Protein FASTA missing: ", aa_local)

  seq_map <- get_sequences(
    gtf_df,
    transcript_path = tx_local,
    translation_path = aa_local
  )


  ### ---- Save cache ----
  message("[CACHE] Saving cleaned annotations and sequences ->", base_dir)
  saveRDS(gtf_df,  rds_gtf)
  saveRDS(seq_map, rds_seq)
  saveRDS(hybrids,    rds_hyb)

  ### ---- Return ----
  list(
    annotations = gtf_df,
    sequences = seq_map,
    hybrids = hybrids
  )
}

#' Helper to fetch test files for vignettes / tests
#' @param filename file to probe
#' @keywords internal
#' @return proper path to example data
get_example_data <- function(filename) {
  # Locate extdata folder inside the installed package
  data_dir <- system.file("extdata", package = "SpliceImpactR")
  if (data_dir == "") stop("extdata directory not found in SpliceImpactR")

  # Build full file path
  fpath <- file.path(data_dir, filename)
  if (!file.exists(fpath)) stop("File not found in extdata: ", filename)

  return(fpath)
}


#' Helper to fetch test dirs for vignettes / tests
#' @keywords internal
#' @param dir_name Name (or relative path) of the directory inside `inst/extdata/`.
#' @param package Package name (defaults to `"SpliceImpactR"`).
#' @param error Logical; if TRUE (default), stops with an informative message
#' @examples
#' check_extdata_dir("rawData")
#'
#' @return proper directory for data
#'
#' @export
check_extdata_dir <- function(dir_name,
                               package = "SpliceImpactR",
                               error = TRUE) {
  base_dir <- system.file("extdata", package = package)
  if (base_dir == "") {
    msg <- sprintf("No extdata folder found for package '%s'.", package)
    if (error) stop(msg, call. = FALSE) else return(FALSE)
  }

  full_dir <- file.path(base_dir, dir_name)
  if (!dir.exists(full_dir)) {
    msg <- sprintf("Directory not found in extdata: %s", dir_name)
    if (error) stop(msg, call. = FALSE) else return(FALSE)
  }

  return(invisible(full_dir))
}














