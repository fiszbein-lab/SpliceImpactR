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
