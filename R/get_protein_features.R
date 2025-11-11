#' Bin elements under a cumulative cap (internal)
#'
#' Internal helper that groups elements sequentially into bins such that
#' the cumulative value in each bin does not exceed a specified cap.
#' Often used to split items for batch processing or job chunking based
#' on approximate size or cost.
#'
#' @param df A \code{data.frame} containing at least two columns: one
#'   with element names and one with numeric values.
#' @param name_col Character string giving the column name for element
#'   identifiers.
#' @param value_col Character string giving the column name for numeric
#'   values used in cumulative binning.
#' @param cap Numeric scalar giving the maximum cumulative value allowed
#'   per bin (default \code{3100}).
#'
#' @return A list where each element is a character vector of names that
#'   belong to one bin.  The bins are created sequentially in the order
#'   of the input rows.
#'
#' @details
#' The function walks through rows in order, adding elements to the
#' current bin until the running total exceeds \code{cap}, then starts a
#' new bin. The algorithm is greedy: it does not reorder or rebalance
#' after bin formation.
#'
#'
#' @keywords internal
bin_under_cap <- function(df, name_col, value_col, cap = 3100) {
  stopifnot(is.data.frame(df),
            is.character(name_col),
            is.character(value_col))
  nms  <- as.character(df[[name_col]])
  vals <- as.numeric(df[[value_col]])

  bins  <- list()
  totals <- numeric(0)
  cur_names <- character(0)
  cur_sum <- 0

  for (i in seq_along(vals)) {
    v  <- vals[i]
    nm <- nms[i]

    if (cur_sum + v <= cap || cur_sum == 0) {
      cur_names <- c(cur_names, nm)
      cur_sum   <- cur_sum + v
    } else {
      bins[[length(bins) + 1L]] <- cur_names
      totals <- c(totals, cur_sum)
      cur_names <- nm
      cur_sum   <- v
    }
  }

  if (length(cur_names)) {
    bins[[length(bins) + 1L]] <- cur_names
    totals <- c(totals, cur_sum)
  }

  return(bins)
}


#' Split GTF transcripts into balanced chromosome groups (internal)
#'
#' Internal helper that groups protein-coding transcripts by chromosome,
#' ensuring that each group remains below a cumulative size threshold.
#' Useful for dividing GTF processing or annotation tasks into balanced
#' batches.
#'
#' @param gtf_df A \code{data.frame} or \code{data.table} containing
#'   GTF annotations, typically from [load_gtf_long()].
#' @param max_group_size Numeric scalar giving the maximum total number
#'   of transcripts allowed per group (passed to \code{bin_under_cap()}).
#'
#' @return A list of character vectors, where each element corresponds
#'   to a bin of chromosome names grouped under the cumulative cap.
#'
#' @details
#' The function considers only transcripts where
#' \code{type == "transcript"}, \code{transcript_type == "protein_coding"}
#'
#' Chromosome names are simplified by removing a leading \code{"chr"}
#' prefix before grouping.
#'
#' @keywords internal
split_into_bits <- function(gtf_df, max_group_size) {
  chr_counts <- table(gsub("chr", "",
                           gtf_df[type == 'transcript' & transcript_type == 'protein_coding', chr]))
  return(bin_under_cap(data.frame(index = names(chr_counts),
                                  vals = as.integer(chr_counts)),
                       "index", "vals"))
}

#' Retrieve protein feature annotations from Ensembl BioMart (internal)
#'
#' Internal helper to query Ensembl BioMart for per-transcript protein
#' feature annotations such as InterPro or Pfam domains. Designed for
#' batched retrieval using [split_into_bits()] to avoid large single
#' queries.
#'
#' @param protein_features Character vector specifying which feature
#'   types to request (e.g. \code{"interpro"}, \code{"pfam"}).
#' @param gtf_df A \code{data.frame} or \code{data.table} containing GTF
#'   annotations; used to derive chromosome groups for batching.
#' @param max_accession_size Integer specifying the maximum total number
#'   of transcripts per query batch (default \code{3500}).
#' @param species_dataset Character string giving the Ensembl BioMart
#'   dataset (default \code{"hsapiens_gene_ensembl"}). For mouse, use
#'   \code{"mmusculus_gene_ensembl"}.
#' @param release release version from ensemble associated with the gencode
#' version provided in the get_annotation (This can be found in gencode, for
#' example in humans: https://www.gencodegenes.org/human/releases.html)
#'
#' @return A \code{data.table} containing protein feature annotations
#'   including transcript and peptide IDs, feature start/end positions,
#'   and database-specific identifiers (e.g., InterPro accession).
#'
#' @details
#' The function queries the Ensembl BioMart service using
#' \pkg{biomaRt::getBM()} with filters \code{"chromosome_name"} and
#' \code{"transcript_biotype"}, restricted to
#' \code{"protein_coding"} transcripts. Queries are executed in chunks
#' per chromosome group to avoid API timeouts.
#'
#'
#' @keywords internal
get_biomart_protein_features <- function(protein_features = c("interpro"), gtf_df, max_accession_size = 3500, species_dataset = "hsapiens_gene_ensembl", release = 109) {
  options(biomaRt.cache = FALSE)
  if (!(exists('mart'))) {
    mart <- biomaRt::useEnsembl(
      biomart = "genes",
      dataset = species_dataset,
      version = release
    )
  }
  atts <- c("ensembl_transcript_id", "ensembl_peptide_id",
            if ("interpro" %in% protein_features) c("interpro", "interpro_short_description", "interpro_description", "interpro_start", "interpro_end"),
            if (length(protein_features[protein_features != "interpro"]) > 0 && length(protein_features) > 0)
              c(t(outer(protein_features[protein_features != "interpro"], c("", "_start", "_end"), paste0))))

  message(paste0("[PROCESSING] Accessing biomaRt for protein features, retrieving: ", paste0(atts, collapse = ", ")))

  access_groups <- split_into_bits(gtf_df, max_accession_size)
  biomart_list <- lapply(seq_along(access_groups), function(x) {
    bm0 <- biomaRt::getBM(attributes = atts,
                          mart = mart,
                          values = list(chromosome_name = access_groups[[x]],
                                        transcript_biotype = "protein_coding"),
                          filters = c('chromosome_name', "transcript_biotype"))
    out_mes <- paste0("[PROCESSING] Protein feature chunk ", x, " access complete")
    message(out_mes)
    return(bm0)
  })

  res <- do.call(rbind, biomart_list)
  return(setDT(res))
}

#' Convert wide BioMart feature table to long format (internal)
#'
#' Internal helper that reshapes a wide BioMart results table containing
#' per-feature columns (e.g. InterPro, Pfam, TMHMM) into a unified long
#' format suitable for downstream annotation or visualization.
#'
#' @param ipr A \code{data.frame} or \code{data.table} returned from
#'   [biomaRt::getBM()] containing per-transcript protein feature data.
#' @param features Character vector of feature prefixes to include
#'   (default: \code{c("mobidblite","seg","ncoils","tmhmm","signalp")}).
#' @param include_interpro Logical; whether to also include InterPro
#'   annotations if present (default \code{TRUE}).
#'
#' @return A \code{data.table} in long format with columns:
#'   \code{ensembl_transcript_id}, \code{ensembl_peptide_id},
#'   \code{database}, \code{feature_id}, \code{name}, \code{alt_name},
#'   \code{start}, \code{stop}, and \code{method = "biomaRt"}.
#'
#' @details
#' For each feature type \code{x}, the function expects columns
#' \code{x}, \code{x_start}, and \code{x_end}. These are combined into a
#' single long table. Non-InterPro features use the feature ID for all
#' name fields. InterPro entries optionally use the description fields
#' \code{interpro_description} and
#' \code{interpro_short_description} if available.
#'
#' @importFrom data.table as.data.table melt := rbindlist
#' @importFrom magrittr %>%
#'
#' @keywords internal
to_long_features <- function(ipr,
                             features = c("mobidblite","seg","ncoils","tmhmm","signalp"),
                             include_interpro = TRUE) {
  ipr <- data.table::as.data.table(ipr)

  id_cols <- c("ensembl_transcript_id", "ensembl_peptide_id")

  ## --- pattern features: x, x_start, x_end ---
  # keep only features that actually exist in the table
  feats <- features[
    paste0(features)       %chin% names(ipr) &
      paste0(features,"_start") %chin% names(ipr) &
      paste0(features,"_end")   %chin% names(ipr)
  ]
  long_pat <- data.table()

  if (length(feats)) {
    melted <- data.table::melt(
      ipr,
      id.vars   = id_cols,
      measure   = list(
        feature_id = feats,
        start      = paste0(feats, "_start"),
        stop       = paste0(feats, "_end")
      ),
      variable.name = "database",
    )

    # map 1..k → feature names, drop empties/NAs
    melted[, database := feats[database]]
    long_pat <- melted[
      !is.na(feature_id) & nzchar(as.character(feature_id)) &
        !is.na(start) & !is.na(stop),
      .(ensembl_transcript_id, ensembl_peptide_id,
        database,
        feature_id = paste0(as.character(feature_id)),#, ";", ensembl_peptide_id, ";", start, "-", stop),
        name = paste0(as.character(feature_id)),#, ";", ensembl_peptide_id, ";", start, "-", stop),
        alt_name = paste0(as.character(feature_id)),#, ";", ensembl_peptide_id, ";", start, "-", stop),
        start = as.integer(start),
        stop = as.integer(stop))
    ]
  }

  ## --- InterPro (special: extra name fields) ---
  long_ipr <- data.table()
  if (include_interpro &&
      all(c("interpro","interpro_start","interpro_end") %chin% names(ipr))) {
    nm  <- intersect("interpro_description", names(ipr))
    alt <- intersect("interpro_short_description", names(ipr))
    nm  <- if (length(nm)) nm else "interpro"
    alt <- if (length(alt)) alt else NA_character_

    long_ipr <- ipr[
      !is.na(interpro) & nzchar(interpro) &
        !is.na(interpro_start) & !is.na(interpro_end),
      .(ensembl_transcript_id, ensembl_peptide_id,
        database   = "interpro",
        feature_id = interpro,
        name       = get(nm),
        alt_name   = if (is.character(alt)) get(alt) else NA_character_,
        start      = as.integer(interpro_start),
        stop       = as.integer(interpro_end))
    ]
  }

  # bind and order
  out <- data.table::rbindlist(list(long_ipr, long_pat), use.names = TRUE, fill = TRUE)
  out[, method := "biomaRt"]
  return(unique(out)[order(ensembl_transcript_id, database, start, stop)])
}

#' Standardize user-supplied protein feature annotations
#'
#' Internal helper that converts a user-provided feature table into the
#' standardized long-format schema used throughout SpliceImpactR.
#' Ensures consistent column names, data types, and coordinate logic.
#'
#' @param x A \code{data.frame} or \code{data.table} supplied by the
#'   user containing at least \code{name}, \code{start}, and \code{stop},
#'   and optionally one or both of \code{ensembl_transcript_id} or
#'   \code{ensembl_peptide_id}.
#' @param default_database Character scalar giving the default value for
#'   the \code{database} column when none is provided
#'   (default \code{"user"}).
#'
#' @return A \code{data.table} with standardized columns:
#'   \code{ensembl_transcript_id}, \code{ensembl_peptide_id},
#'   \code{database}, \code{feature_id}, \code{name}, \code{alt_name},
#'   \code{start}, \code{stop}, and \code{method = "manual"}.
#'
#' @details
#' The function verifies the presence of required coordinate columns and
#' at least one Ensembl identifier, fills in any missing optional
#' columns with default values, coerces all types to their expected
#' formats, removes invalid coordinates (\code{stop < start}), and
#' reorders columns into the canonical schema used by
#' [to_long_features()].
#'
#'
#' @importFrom data.table := setcolorder as.data.table
#' @keywords internal
add_user_features <- function(x, default_database = "user") {
  DT <- data.table::as.data.table(x)

  # --- required columns ---
  has_enst <- "ensembl_transcript_id" %in% names(DT)
  has_ensp <- "ensembl_peptide_id"    %in% names(DT)
  if (!has_enst && !has_ensp) {
    stop("Provide at least one of: ensembl_transcript_id or ensembl_peptide_id.")
  }
  req <- c("name","start","stop")
  miss <- setdiff(req, names(DT))
  if (length(miss)) stop("Missing required column(s): ", paste(miss, collapse=", "))

  # --- coerce types & fill optionals ---
  # if absent, create columns so downstream is uniform
  if (!has_enst) DT[, ensembl_transcript_id := NA_character_]
  if (!has_ensp) DT[, ensembl_peptide_id    := NA_character_]
  if (!("database"   %in% names(DT))) DT[, database   := default_database]
  if (!("alt_name"   %in% names(DT))) DT[, alt_name   := NA_character_]
  if (!("feature_id" %in% names(DT))) DT[, feature_id := name]

  # type safety
  DT[, `:=`(
    ensembl_transcript_id = as.character(ensembl_transcript_id),
    ensembl_peptide_id    = as.character(ensembl_peptide_id),
    database              = as.character(database),
    feature_id            = as.character(feature_id),
    name                  = as.character(name),
    alt_name              = as.character(alt_name),
    start                 = as.integer(start),
    stop                  = as.integer(stop),
    method                = "manual"
  )]

  # basic sanity filters
  bad <- DT[is.na(start) | is.na(stop) | stop < start]
  if (nrow(bad)) {
    warning("Dropped ", nrow(bad), " row(s) with invalid start/stop.")
    DT <- DT[!(is.na(start) | is.na(stop) | stop < start)]
  }

  # standard column order, de-dup
  data.table::setcolorder(DT, c("ensembl_transcript_id","ensembl_peptide_id",
                                "database","feature_id","name","alt_name",
                                "start","stop"))
  return(unique(DT))
}


#' Read an object from a supported file format (.rds, .csv, .tsv)
#'
#' Internal helper to load serialized or tabular data in a uniform way.
#'
#' @param path Character scalar; path to a file ending in `.rds`, `.csv`, or `.tsv`.
#'
#' @return An R object (data.table, data.frame, or arbitrary R object from `.rds`).
#' @keywords internal
#' @importFrom tools file_ext
#' @importFrom data.table fread
.read_any <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    readRDS(path)
  } else if (ext %in% c("csv","tsv")) {
    data.table::fread(path)
  } else {
    stop("Unsupported load file extension: ", ext, " (use .rds or .csv/.tsv)")
  }
}

#' Write an object to a supported file format (.rds, .csv, .tsv)
#'
#' Internal helper to save serialized or tabular data in a uniform way.
#'
#' @param dt Object to save (typically a data.frame or data.table).
#' @param path Character scalar; output file path ending in `.rds`, `.csv`, or `.tsv`.
#'
#' @return Invisibly returns the output file path.
#' @keywords internal
#' @importFrom tools file_ext
#' @importFrom data.table fwrite
.write_any <- function(dt, path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    saveRDS(dt, path)
  } else if (ext %in% c("csv","tsv")) {
    data.table::fwrite(dt, path)
  } else {
    stop("Unsupported save file extension: ", ext, " (use .rds or .csv/.tsv)")
  }
}


#' @title External function to fetch protein features from biomaRt
#' @description Here we also remove any duplicate and overlapping domains
#' We also add the genomic location to the name of the protein feature for
#' downstream safeguarding and precision. This is to prevent different occurences
#' of the same domain being called as the same in domain identification and
#' enrichment.
#' @param biomaRt_databases choose what biomaRt attribute to access, defaulting
#' to interpro, mobidblite, seg, ncoils, tmhmm, signalp
#' @param gtf_df annotations from get_annotation()
#' @param load_path path to load prior protein features from
#' @param save_path path to save prior protein features from
#' @param timeout ability to extend timeout if biomaRt is not cooperating
#' @param species Character string giving the Ensembl BioMart
#' dataset (default \code{"human"}). For mouse, use
#' \code{"mouse"}.
#' @param release release version from ensemble associated with the gencode
#' version provided in the get_annotation (This can be found in gencode, for
#' example in human: https://www.gencodegenes.org/human/releases.html)
#' example in mouse: https://www.gencodegenes.org/mouse/releases.html)
#'
#' @param test Logical; bool for whether to load from reduced test set
#'
#' @importFrom data.table rbindlist rleid
#' @examples
#' annotation_df <- get_annotation(load = "test")
#' interpro_features <- get_protein_features(c("interpro"), annotation_df$annotations, timeout = 600, test = TRUE)
#' @return
#' A `data.table` with one row per protein feature and transcript coupling
#'
#' @export
get_protein_features <- function(biomaRt_databases = c("interpro", "mobidblite", "seg", "ncoils", "tmhmm", "signalp"),
                                 gtf_df,
                                 load_path = NULL,
                                 save_path = NULL,
                                 timeout = 600,
                                 species = c('human', 'mouse')[1],
                                 release = c(109, 115)[1],
                                 test = FALSE) {
  if (species == 'human') {
    species <- 'hsapiens_gene_ensembl'
  } else if (species == 'mouse') {
    species <- 'mmusculus_gene_ensembl'
  }
  if (test == TRUE) {
    return(data.table::rbindlist(lapply(biomaRt_databases[biomaRt_databases %in% c("interpro", "signalp")], function(x) {
      fread(get_example_data(paste0("test_", x, ".csv")))
    })))
  }
  options(timeout = timeout)
  if (!is.null(load_path) && file.exists(load_path)) {
    message("[LOADING] Protein features loaded from: ", load_path)
    pf <- .read_any(load_path)
    return(as.data.table(pf))
  }
  pf <- get_biomart_protein_features(protein_features = biomaRt_databases,
                                     gtf_df = gtf_df,
                                     species_dataset = species,
                                     release = release) %>% to_long_features
  message("[PROCESSING] Deduping and locating genomic coordinates")
  ## Dedup overlapping identical domains
  pf <- pf[ensembl_transcript_id %in% gtf_df$transcript_id]

  pf_i <- pf

  pf <- as.data.table(pf)[, .(ensembl_transcript_id, feature_id, start, stop)]

  pf[, `:=`(start = as.integer(start), stop = as.integer(stop))]

  pf <- unique(pf)

  setorder(pf, ensembl_transcript_id, feature_id, start, stop)

  pf[, k := rleid(ensembl_transcript_id, feature_id)]

  pf[, rmax := cummax(stop), by = k]
  pf[, grp  := cumsum(start > shift(rmax, fill = -1)), by = k]

  # aggregate per (key, grp) in one shot
  ans <- pf[, .(start = min(start), stop = max(stop)),
            by = .(ensembl_transcript_id, feature_id, grp)]
  ans[, grp := NULL][]

  pf <- merge(ans,
        unique(pf_i[, .(ensembl_transcript_id, ensembl_peptide_id, database, feature_id, name, alt_name, method)]),
        by = c("ensembl_transcript_id", "feature_id"))

  ## get a unique genomic coord for domain
  cds_map <- as.data.table(gtf_df)[
    type == "exon" & cds_has == TRUE,
    .(
      ensembl_transcript_id = transcript_id,
      chr,
      strand,
      cds_rel_start,
      cds_rel_stop,
      cds_gen_start,
      cds_gen_stop
    )
  ][order(ensembl_transcript_id, cds_rel_start)]
  setkey(cds_map, ensembl_transcript_id, cds_rel_start, cds_rel_stop)

  cds_len_nt <- cds_map[, .(cds_len_nt = max(cds_rel_stop)), by = ensembl_transcript_id]
  cds_len_aa <- cds_len_nt[, .(ensembl_transcript_id, cds_len_aa = floor(cds_len_nt / 3L))]
  pf <- merge(pf, cds_len_aa, by = "ensembl_transcript_id", all.x = TRUE)
  pf[, stop := pmin(stop, cds_len_aa, na.rm = TRUE)]
  pf[, start := pmin(start, stop, na.rm = TRUE)]  # just in case start > stop

  # Drop any zero-length or NA domains that got truncated to nothing
  pf <- pf[!is.na(start) & !is.na(stop)]

  out <- as.data.table(pf)
  out[, cds_nt_start := start * 3L - 2L]
  out[, cds_nt_end   := stop  * 3L]

  strand_info <- unique(gtf_df[
    type == "transcript" & !is.na(strand),
    .(ensembl_transcript_id = transcript_id, chr, strand)
  ])
  out <- merge(out, strand_info, by = "ensembl_transcript_id", all.x = TRUE)
  out <- out[!is.na(strand)]

  map_to_genomic <- function(features, cds_map) {
    features[, c("genomic_start", "genomic_end") := {

      cds_sub <- cds_map[.BY$ensembl_transcript_id]
      if (nrow(cds_sub) == 0L) {
        list(NA_integer_, NA_integer_)
      } else {

      i1 <- pmin(findInterval(cds_nt_start, cds_sub$cds_rel_stop+1)+1, nrow(cds_sub))
      i2 <- pmin(findInterval(cds_nt_end,   cds_sub$cds_rel_stop+1)+1, nrow(cds_sub))

      if (.BY$strand == "+") {
        gstart <- cds_sub$cds_gen_start[i1] + (cds_nt_start - cds_sub$cds_rel_start[i1])
        gend   <- cds_sub$cds_gen_start[i2] + (cds_nt_end   - cds_sub$cds_rel_start[i2])
      } else {
        gstart <- cds_sub$cds_gen_stop[i1] - (cds_sub$cds_rel_stop[i1] - cds_nt_start)
        gend   <- cds_sub$cds_gen_stop[i2] - (cds_sub$cds_rel_stop[i2] - cds_nt_end)
      }
      list(gstart, gend)
      }
    }, by = .(ensembl_transcript_id, strand)]

    features[]
  }

  mapped <- map_to_genomic(out, cds_map)
  mapped <- mapped[!is.na(genomic_start) & !is.na(genomic_end)]
  mapped <- mapped[
    , .(
      chr = first(chr),
      strand = first(strand),
      genomic_start = if (all(is.na(genomic_start))) NA_real_ else min(genomic_start, na.rm = TRUE),
      genomic_end   = if (all(is.na(genomic_end))) NA_real_ else max(genomic_end,   na.rm = TRUE),
      feature_id = first(feature_id),
      clean_name = first(name),
      alt_name = first(alt_name),
      database = first(database),
      ensembl_peptide_id = first(ensembl_peptide_id),
      method = first(method)
    ),
    by = .(ensembl_transcript_id, start, stop)
  ]

  pf <- mapped[, name := paste0(clean_name, ";", chr, ":",
                                format(as.numeric(genomic_start), scientific = FALSE, trim = TRUE), "-",
                                format(as.numeric(genomic_end), scientific = FALSE, trim = TRUE)
  )]
  pf <- pf[, `:=` (genomic_start = NULL, genomic_end = NULL)]


  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    message("[SAVING] Protein features saved to: ", save_path)
    .write_any(pf, save_path)
  }
  return(pf)
}


#' Incorporate user-supplied protein features
#'
#' Converts a manual feature table into the standardized long format
#' and optionally merges it with biomaRt-derived features.
#'
#' @param manual_features Data.frame or data.table with at least
#'   \code{name}, \code{start}, \code{stop} amino acid, and one of
#'   \code{ensembl_transcript_id} or \code{ensembl_peptide_id}.
#' @param gtf_df get_annotation annotation output
#' @param biomaRt_features Optional \code{data.table} of features from
#'   [get_protein_features()] to merge with.
#' @param load_path Optional path to load precomputed manual features.
#' @param save_path Optional path to save the combined feature table.
#'
#' @return A \code{data.table} of manual (and optionally combined)
#'   protein features.
#'
#' @examples
#' annotation_df <- get_annotation(load = "test")
#' user_df <- data.frame(
#'   ensembl_transcript_id = c(
#'     "ENST00000511072","ENST00000374900","ENST00000373020","ENST00000456328",
#'     "ENST00000367770","ENST00000331789","ENST00000335137","ENST00000361567",
#'     NA,                    "ENST00000380152"
#'   ),
#'   ensembl_peptide_id = c(
#'     "ENSP00000426975", NA,                   "ENSP00000362048","ENSP00000407743",
#'     "ENSP00000356802","ENSP00000326734", NA,                  "ENSP00000354587",
#'     "ENSP00000364035", NA
#'   ),
#'   name = c(
#'     "Low complexity","Transmembrane helix","Coiled-coil","Signal peptide",
#'     "Transmembrane helix","Low complexity","Coiled-coil","Transmembrane helix",
#'     "Signal peptide","Low complexity"
#'   ),
#'   start = c(80L, 201L, 35L, 1L, 410L, 150L, 220L, 30L, 1L, 300L),
#'   stop  = c(120L,223L, 80L, 20L, 430L, 190L, 260L, 55L, 24L, 360L),
#'   database   = c("seg","tmhmm","ncoils","signalp","tmhmm","seg","ncoils","tmhmm","signalp", NA),
#'   alt_name   = c(NA,"TMhelix",NA,"SignalP-noTM", "TMhelix", NA, NA, "TMhelix", "SignalP-TAT", NA),
#'   feature_id = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
#' )
#' user_features <- get_manual_features(user_df, annotation_df$annotation)
#'
#' @export
#' @importFrom data.table first
#' @seealso [add_user_features()]
get_manual_features <- function(manual_features,
                                gtf_df,
                                biomaRt_features = NULL,
                                load_path = NULL,
                                save_path = NULL) {
  if (!is.null(load_path) && file.exists(load_path)) {
    message("[LOADING] Manual protein features loaded from: ", load_path)
    manual_features_loaded <- .read_any(load_path)
    return(as.data.table(manual_features_loaded))
  }
  manual_features_loaded <- add_user_features(manual_features)

  cds_map <- as.data.table(gtf_df)[
    type == "exon" & cds_has == TRUE,
    .(
      ensembl_transcript_id = transcript_id,
      chr,
      strand,
      cds_rel_start,
      cds_rel_stop,
      cds_gen_start,
      cds_gen_stop
    )
  ][order(ensembl_transcript_id, cds_rel_start)]
  setkey(cds_map, ensembl_transcript_id, cds_rel_start, cds_rel_stop)

  cds_len_nt <- cds_map[, .(cds_len_nt = max(cds_rel_stop)), by = ensembl_transcript_id]
  cds_len_aa <- cds_len_nt[, .(ensembl_transcript_id, cds_len_aa = floor(cds_len_nt / 3L))]
  manual_features_loaded <- merge(manual_features_loaded, cds_len_aa, by = "ensembl_transcript_id", all.x = TRUE)
  manual_features_loaded[, stop := pmin(stop, cds_len_aa, na.rm = TRUE)]
  manual_features_loaded[, start := pmin(start, stop, na.rm = TRUE)]

  # Drop any zero-length or NA domains that got truncated to nothing
  manual_features_loaded <- manual_features_loaded[!is.na(start) & !is.na(stop)]

  out <- as.data.table(manual_features_loaded)
  out[, cds_nt_start := start * 3L - 2L]
  out[, cds_nt_end   := stop  * 3L]

  strand_info <- unique(gtf_df[
    type == "transcript" & !is.na(strand),
    .(ensembl_transcript_id = transcript_id, chr, strand)
  ])
  out <- merge(out, strand_info, by = "ensembl_transcript_id", all.x = TRUE)
  out <- out[!is.na(strand)]

  map_to_genomic <- function(features, cds_map) {
    features[, c("genomic_start", "genomic_end") := {

      cds_sub <- cds_map[.BY$ensembl_transcript_id]
      if (nrow(cds_sub) == 0L) {
        list(NA_integer_, NA_integer_)
      } else {

      i1 <- pmin(findInterval(cds_nt_start, cds_sub$cds_rel_stop+1)+1, nrow(cds_sub))
      i2 <- pmin(findInterval(cds_nt_end,   cds_sub$cds_rel_stop+1)+1, nrow(cds_sub))

      if (.BY$strand == "+") {
        gstart <- cds_sub$cds_gen_start[i1] + (cds_nt_start - cds_sub$cds_rel_start[i1])
        gend   <- cds_sub$cds_gen_start[i2] + (cds_nt_end   - cds_sub$cds_rel_start[i2])
      } else {
        gstart <- cds_sub$cds_gen_stop[i1] - (cds_sub$cds_rel_stop[i1] - cds_nt_start)
        gend   <- cds_sub$cds_gen_stop[i2] - (cds_sub$cds_rel_stop[i2] - cds_nt_end)
      }
      list(gstart, gend)
      }
    }, by = .(ensembl_transcript_id, strand)]

    features[]
  }
  mapped <- map_to_genomic(out, cds_map)
  mapped$na_genomic <- is.na(mapped$genomic_end) | is.na(mapped$genomic_start)
  badRows <- paste0(sum(mapped$na_genomic), " domain(s) not matched to genomic coords")
  message(badRows)
  mapped <- mapped[!(na_genomic)]
  mapped <- mapped[!is.na(genomic_start) & !is.na(genomic_end)]
  mapped <- mapped[
    , .(
      chr = first(chr),
      strand = first(strand),
      genomic_start = if (all(is.na(genomic_start))) NA_real_ else min(genomic_start, na.rm = TRUE),
      genomic_end   = if (all(is.na(genomic_end))) NA_real_ else max(genomic_end,   na.rm = TRUE),
      # feature_id = first(feature_id),
      clean_name = first(name),
      alt_name = first(alt_name),
      database = first(database),
      ensembl_peptide_id = first(ensembl_peptide_id),
      method = first(method),
      name = first(name),
      na_genomic = first(na_genomic)
    ),
    by = .(ensembl_transcript_id, feature_id, start, stop)
  ]


  manual_features_loaded <- mapped[, name := paste0(clean_name, ";", chr, ":",
                                                    format(as.numeric(genomic_start), scientific = FALSE, trim = TRUE), "-",
                                                    format(as.numeric(genomic_end), scientific = FALSE, trim = TRUE)
  )]

  manual_features_loaded <- manual_features_loaded[, `:=` (na_genomic = NULL,
                                                         genomic_start = NULL,
                                                         genomic_end = NULL)]

  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    message("[SAVING] Protein features saved to: ", save_path)
    .write_any(manual_features_loaded, save_path)
  }

  if (is.null(biomaRt_features)) {
    return(manual_features_loaded)
  } else {
    return(rbind(manual_features_loaded, biomaRt_features))
  }
}

#' Combine multiple sources of protein feature annotations
#'
#' Aggregates all available feature tables (e.g., biomaRt + manual)
#' into a single unified long-format annotation table.
#'
#' @param protein_feature_list List of \code{data.table} objects,
#'   typically from [get_protein_features()] or [get_manual_features()].
#' @param load_path_list Optional vector of file paths to load each
#'   feature source from disk (instead of providing in memory).
#' @param save_path Optional path to cache the combined annotations.
#'
#' @examples
#' annotation_df <- get_annotation(load = "test")
#' user_df <- data.frame(
#'   ensembl_transcript_id = c(
#'     "ENST00000511072","ENST00000374900","ENST00000373020","ENST00000456328",
#'     "ENST00000367770","ENST00000331789","ENST00000335137","ENST00000361567",
#'     NA,                    "ENST00000380152"
#'   ),
#'   ensembl_peptide_id = c(
#'     "ENSP00000426975", NA,                   "ENSP00000362048","ENSP00000407743",
#'     "ENSP00000356802","ENSP00000326734", NA,                  "ENSP00000354587",
#'     "ENSP00000364035", NA
#'   ),
#'   name = c(
#'     "Low complexity","Transmembrane helix","Coiled-coil","Signal peptide",
#'     "Transmembrane helix","Low complexity","Coiled-coil","Transmembrane helix",
#'     "Signal peptide","Low complexity"
#'   ),
#'   start = c(80L, 201L, 35L, 1L, 410L, 150L, 220L, 30L, 1L, 300L),
#'   stop  = c(120L,223L, 80L, 20L, 430L, 190L, 260L, 55L, 24L, 360L),
#'   database   = c("seg","tmhmm","ncoils","signalp","tmhmm","seg","ncoils","tmhmm","signalp", NA),
#'   alt_name   = c(NA,"TMhelix",NA,"SignalP-noTM", "TMhelix", NA, NA, "TMhelix", "SignalP-TAT", NA),
#'   feature_id = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
#' )
#' user_features <- get_manual_features(user_df, annotation_df$annotation)
#' interpro_features <- get_protein_features(c("interpro"), annotation_df$annotations, timeout = 600, test = TRUE)
#' protein_feature_total <- get_comprehensive_annotations(list(user_features, interpro_features))
#'
#' @return A combined \code{data.table} containing all unique feature rows.
#' @importFrom data.table data.table
#' @export
#'
get_comprehensive_annotations <- function(protein_feature_list, load_path_list=NULL, save_path=NULL) {
  if (!is.null(load_path_list) && lapply(load_path_list, file.exists)) {
    message("[LOADING] All protein features loaded from: ", paste0(load_path_list, collapse = ', '))
    protein_feature_list <- lapply(load_path_list, .read_any)
    full_features_loaded <- do.call(rbind, protein_feature_list)
    return(as.data.table(full_features_loaded))
  }
  full_features_loaded <- do.call(rbind, protein_feature_list)
  return(do.call(rbind, protein_feature_list))
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    message("[SAVING] All protein features saved to: ", save_path)
    .write_any(full_features_loaded, save_path)
  }
}


#' Map protein features to coding exons
#'
#' Overlaps amino acid-based protein features (e.g., InterPro, TMHMM)
#' with exon coding spans defined by transcript-relative CDS coordinates.
#'
#' @param gtf_dt A \code{data.frame} or \code{data.table} containing
#'   transcript and exon annotation, including \code{type},
#'   \code{transcript_id}, \code{cds_rel_start}, and \code{cds_rel_stop}
#'   columns (see [add_exon_coding_information()]).
#' @param feat A \code{data.frame} or \code{data.table} of long-format
#'   protein features (from [get_protein_features()] or
#'   [get_comprehensive_annotations()]) containing columns
#'   \code{ensembl_transcript_id}, \code{start}, \code{stop},
#'   \code{database}, \code{feature_id}, \code{name}, and
#'   \code{alt_name}.
#' @param inclusive Logical; whether to round both feature and exon
#'   boundaries upward when converting from nucleotide to amino acid
#'   coordinates (default \code{TRUE}). If \code{FALSE}, downstream
#'   exons own partial codons to avoid double counting.
#'
#' @return A \code{data.table} of overlapping feature-exon pairs with
#'   amino acid coordinates (\code{overlap_aa_start}, \code{overlap_aa_end},
#'   \code{overlap_aa_len}) and associated exon and feature metadata.
#'
#' @examples
#' annotation_df <- get_annotation(load = 'test')
#' interpro_features <- get_protein_features(c("interpro"), annotations$annotations, timeout = 600, test = TRUE)
#' protein_feature_total <- get_comprehensive_annotations(list(interpro_features))
#'
#' exon_features <- get_exon_features(annotation_df$annotations, protein_feature_total)
#'
#' @importFrom data.table as.data.table first setkey setnames
#'
#' @export
get_exon_features <- function(gtf_dt, feat, inclusive = TRUE) {
  gtf <- data.table::as.data.table(gtf_dt)
  Fe   <- data.table::as.data.table(feat)

  # minimal checks
  stopifnot(all(c("type","transcript_id","cds_rel_start","cds_rel_stop") %in% names(gtf)))
  stopifnot(all(c("ensembl_transcript_id","start","stop","database","feature_id","name","alt_name") %in% names(Fe)))

  # nt -> aa helpers
  if (inclusive) {
    nt_to_aa_start <- function(nt) ifelse(is.na(nt), NA_integer_, ceiling(nt / 3))
    nt_to_aa_end   <- function(nt) ifelse(is.na(nt), NA_integer_, ceiling(nt / 3))
  } else {
    # “downstream owns partial codon” (no double count)
    nt_to_aa_start <- function(nt) ifelse(is.na(nt), NA_integer_, ceiling(nt / 3))
    nt_to_aa_end   <- function(nt) ifelse(is.na(nt), NA_integer_, floor(nt / 3))
  }

  # build exon AA spans from cds_rel_* on exon rows
  ex <- gtf[type == "exon" & !is.na(transcript_id)]
  if (!nrow(ex)) return(data.table())

  ex[, `:=`(
    exon_aa_start = nt_to_aa_start(cds_rel_start),
    exon_aa_end   = nt_to_aa_end(cds_rel_stop)
  )]
  ex <- ex[!is.na(exon_aa_start) & !is.na(exon_aa_end) & exon_aa_end >= exon_aa_start]

  # group by exon identity (prefer exon_id; else exon_number)
  grp <- if ("exon_id" %in% names(ex) && any(!is.na(ex$exon_id)))
    c("transcript_id","exon_id") else c("transcript_id","exon_number")

  exon_aa <- ex[, .(
    gene_id       = data.table::first(na.omit(gene_id)),
    exon_id       = data.table::first(na.omit(exon_id)),
    exon_number   = suppressWarnings(as.integer(
      data.table::first(na.omit(exon_number)))),
    strand        = data.table::first(strand),
    exon_aa_start = min(exon_aa_start),
    exon_aa_end   = max(exon_aa_end)
  ), by = grp]
  data.table::setnames(exon_aa, "transcript_id", "ensembl_transcript_id")

  # overlap features (AA) with exon AA spans
  keep_feat <- Fe[!is.na(start) & !is.na(stop) & start <= stop]
  if (!nrow(keep_feat) || !nrow(exon_aa)) return(data.table())
  data.table::setkey(keep_feat, ensembl_transcript_id)
  data.table::setkey(exon_aa,   ensembl_transcript_id)

  j <- exon_aa[keep_feat, allow.cartesian = TRUE, nomatch = 0L]
  if (!nrow(j)) return(data.table())

  j[, `:=`(
    overlap_start = pmax(exon_aa_start, start),
    overlap_end   = pmin(exon_aa_end,   stop)
  )]
  j <- j[overlap_start <= overlap_end]

  j[, .(
    gene_id,
    ensembl_transcript_id,
    ensembl_peptide_id,
    exon_id,
    exon_number,
    strand,
    database, feature_id, name, alt_name,
    prot_start = start, prot_stop = stop,
    exon_aa_start, exon_aa_end,
    overlap_aa_start = overlap_start,
    overlap_aa_end   = overlap_end,
    overlap_aa_len   = overlap_end - overlap_start + 1L
  )][order(ensembl_transcript_id, exon_number, database, prot_start, prot_stop)]
}


















