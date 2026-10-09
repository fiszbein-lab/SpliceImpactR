# SpliceImpactR 1.1.4

## Bug fixes

* BioMart queries for `release` go to `https://e<release>.ensembl.org`, which
  Ensembl forwards to that release's archive, instead of to the archive host
  itself. With biomaRt before 2.70, new installations could no longer connect
  to archive hosts ("Unable to contact any Ensembl mirror"): biomaRt first
  downloads Ensembl's list of archives, which is no longer available. Release
  116 now works with any biomaRt version. The data come from the same
  archives, so results and cached features are unchanged.
* With biomaRt before 2.70, connection errors for an archive host passed as
  `ensembl_host` explain this and suggest the e<release> address.

# SpliceImpactR 1.1.3

## Breaking changes and changed results

* `get_protein_features()` queries Ensembl 111 by default (was 109), the
  Ensembl release of the default GENCODE v45 annotation from
  `get_annotation()`. Features from full runs change, and cached features are
  fetched again once.
* `ensembl_mirror` is deprecated and ignored: Ensembl retired its BioMart
  mirrors.

## Bug fixes

* BioMart queries go to the Ensembl archive of `release` by its explicit host
  (releases 105 to 116). Ensembl 116 (June 2026) is the last release with
  BioMart, and www.ensembl.org no longer serves it. The new `ensembl_host`
  argument selects another host. Release 116 needs biomaRt 2.70 or later
  (earlier versions send its queries to www.ensembl.org), and connection errors
  name the host.

# SpliceImpactR 1.1.2

## Breaking changes and changed results

* `get_splicing_impact()` selects transcript pairs with the ORF-aware matcher
  (`matching = "orf"`) by default; `matching = "legacy"` restores the previous
  matcher. The ORF matcher checks each form's exact event structure, never
  pairs a transcript with itself, labels pairs without a structural match as
  `approximate`, and returns unresolved comparisons in `matching`. Selected
  transcripts and pair counts change, and matching is slower.

## New features

* `get_ranked_pairs(return_class = "S4")` stores the selected pairs, form rows
  and diagnostics in a `SpliceImpactResult`, where `get_splicing_impact()`
  stores them, so the stepwise S4 workflow uses the ORF matcher. The default
  still returns the list. `get_pairs()` S4 output sets
  `metadata$matching = "legacy"` and drops earlier ORF diagnostics.

## Documentation

* The README and vignettes describe the ORF matcher as the default, including
  in the stepwise S4 workflow. The README also covers one-based rMATS
  coordinates, the domain-enrichment population, gene universes, new output
  columns and input requirements. The matching vignette notes that diagnostics
  are empty tables, not `NULL`, when no forms are significant.

# SpliceImpactR 1.1.1

## Breaking changes and changed results

* rMATS coordinates are converted to one-based, closed intervals: starts move
  by +1 and retained introns exclude the flanking exon ends. Differential
  inclusion tables saved with earlier versions miss exact matches by 1 nt and
  must be regenerated; `get_rmats_post_di()` returns package-format tables
  unchanged.
* A3SS and A5SS forms include the partner (flanking) exon in `inc`, in genomic
  order. `inc` strings and event counts change, and events that differ only in
  their partner exon are no longer merged.
* `get_rmats_post_di()` numbers events across files and stops when one event
  type appears under more than one grp1/grp2 comparison. The new `case_group`
  argument selects the case group (default 1, rMATS `--b1`, as before), and
  each import reports which group it used.
* `get_background()` defaults to `source = "annotated"` (was `"hit_index"`).
  Passing `input` without `source` warns that it is ignored; `"hit_index"` or
  `"user-given"` without `input` is an error. The `gene_universe` and
  `feature_gene_universe` attributes give the eligible gene universes.
* Domain enrichment counts only pairs with a tested domain change, in both
  foreground and background, and counts repeated transcript pairs once.
  Foreground pairs and their domains must be in the background, and the odds
  ratio compares the foreground with the rest of the background.
  `domain_col_fg` and `domain_col_bg` are honoured, and `delim` defaults to
  `"\\|"`. Enrichment results change.
* The legacy matcher ranks candidates by structural fit (exon class, overlap,
  width) before protein-coding status, with transcript support level as the
  last tie-break, and candidate exons must belong to the event's gene. Selected
  transcripts can change.
* Forms without a matched transcript no longer produce pairs; they were paired
  with an `NA` transcript.
* `summary_classification` has a new `"NMD"` class, with top precedence, when
  either pair member is annotated as nonsense-mediated decay.
* `import_di_table()` stops on strands other than `+`/`-` and uses `event_id`,
  `form` and `padj` columns when present.
* `get_rmats_hit()` stops on unknown event types and returns only the requested
  HIT index types.
* `get_user_data_post_di()` requires `event_id` and stops on events without
  INC and EXC forms or a SITE form, or that mix them. Without `event_id` it
  used to give each row its own ID, so forms could never pair.
* Protein-feature cache keys hash the full input content, so features are
  fetched again once; old cache files are left on disk.
* Requires R (>= 4.5.0).

## New features

* Opt-in ORF-aware transcript-pair matching: `get_ranked_pairs()` and
  `get_splicing_impact(matching = "orf")`, with `matching_max_candidates` and
  `matching_fallback`. Pairs carry `matching_tier` and ranking diagnostics, and
  unresolved comparisons are returned with reasons. The legacy matcher remains
  the default. The new `transcript_pair_ranking` vignette describes the
  protocol.
* `get_annotation(filter_tsl = NULL)` keeps transcripts of every support
  level, including missing ones, in a separate processed cache. The default
  remains TSL 1-3.
* `keep_sig_pairs()` adds `site_significant`, and pairs carry
  `site_significant_case`, `site_significant_control` and
  `n_event_comparisons` for events with several sites.
* Matched and paired tables carry each transcript's `transcript_type`.
* `get_matched_events_chunked()` gains `verbose`; chunk progress is reported
  with `message()`.
* `get_manual_features()` places rows given only `ensembl_peptide_id` and
  reports rows it cannot place.
* Empty analyses return zero-row tables with their usual columns, in table and
  S4 modes, and plots return placeholders.

## Bug fixes

* `get_pairs(source = "paired")` failed whenever an event had two forms; it now
  takes positive delta PSI as the case, like `"multi"`.
* PPI switches use the side of each interaction the event gene is on.
  Domain-motif changes with the motif on gene A were never detected, and
  domain-domain changes could be called from a partner's domain.
* HIT index import reads plain and `.gz` files once (PSI was halved), finds
  `.gz` PSI and `.exon.gz` files, and uses anchored file patterns. rMATS
  `.txt.gz` files are read, and one-base segments are kept.
* `get_rmats_hit()` binds rMATS and HIT index tables by column name.
* User-given size factors (`'user-given'`) failed; sample folders without a
  trailing slash now work.
* Frame checks follow transcript order upstream of the event, A5SS frame checks
  use the variable exon, and event lengths sum the selected transcript's exons.
* Manual features on minus-strand transcripts map to the right genomic
  positions, and BioMart features at identical coordinates are kept when their
  database or ID differ.
* `plot_alignment_summary(output_file = )` saved an undefined object. Feature
  plots draw every interval of multi-interval `inc`/`exc`, and
  `plot_enriched_domains_counts()` counts unique events.
* `get_splicing_impact()` works from `res` alone and no longer reloads raw data
  when given `sample_frame` and `res`. S4 metadata is replaced, not appended,
  on re-runs, and `has_*` flags reflect content.
* Processed annotation caches are stored inside the cache folder. BioMart
  queries use bounded transcript-ID batches and restore `options()`.
* Inputs are no longer modified in place.

## Other changes

* Imports R.utils (>= 2.13.0), which `data.table::fread()` needs for `.gz`
  input.
* Removed unused internal helpers.

# SpliceImpactR 0.99.4

* Version bump for Bioconductor

# SpliceImpactR 0.99.2

# SpliceImpactR 0.99.1

* Version bump for Bioconductor

# SpliceImpactR 0.99.0

* Added a `NEWS.md` file to track changes to the package.
