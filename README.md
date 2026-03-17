# SpliceImpactR <img src="./inst/screenshot1.png" alt="Fiszbein Lab Logo" width="110" align="right"/> <img src="./inst/screenshot2.png" width="200" align="right"/>

SpliceImpactR is an R package designed for studying the impact of alternative splicing on protein structure and function. 
It provides tools for analyzing RNA-seq data to identify differentially included splicing events and predict their consequences 
on the resulting protein products. SpliceImpactR output involves identifying key changes in proteins at various levels: primary sequence, 
domain content, and transcript-transcript interactions.

The suite of funcitons is designed to anaylyze the consequences of AFE, ALE, SE, MXE, A5SS, RI, and A3SS, along with hybrid exons. 
SpliceImpactR is built to take output from any source with custom functions to take data processed by from the [HIT Index](https://github.com/thepailab/HITindex) and [rMATS](https://github.com/Xinglab/rmats-turbo). 

The package is built to work with human and mouse data, primarily from GENCODE and biomaRt. We also allow for user-defined events and protein features for flexibility of use.

HIT Index data outputs, such as .exon files, are also incorporated into part of the process

## Features
Identification of alternative splicing events from RNA-seq data.
Analysis of the potential impact of splicing events on protein structure.
Functional annotation of spliced isoforms to predict their biological impact.
Integration with existing bioinformatics tools and databases for comprehensive analysis.
Holistic analysis of how the use of different RNA processing events differs.

## Installation
Choose one installation path below.

BiocManager
```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("SpliceImpactR")
```

devtools (GitHub)
```r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("fiszbein-lab/SpliceImpactR")
```

### Load package
```r
library(SpliceImpactR)
```




# Usage
## External tools and acronyms
- `rMATS` (replicate Multivariate Analysis of Transcript Splicing) provides
  event-level splicing tables.
- `HITindex` quantifies first/last exon usage and labels hybrid exons.
- `PPIDM` Protein-Protein Interaction Domain Miner for domain-domain
  interaction derived from PPI and 3did's domain-domain interactions.
- `ELM` Eukaryotic Linear Motif Database for short linear motif occurrences
  and domain-motif interactions.
- `BiomaRt` accesses data from InterPro, PFAM, SignalP, TMHMM, CDD,
  and Mobidb-lite.

## Workflow map
The standard analysis path is:
1. Load reference resources (annotations + protein features).
2. Build sample manifest (`sample_frame`).
3. (Optional) run the quick-start wrapper for end-to-end execution.
4. Read splicing events for stepwise analysis.
5. Run QC summaries and plots.
6. Run differential inclusion (`get_differential_inclusion`) and significance filtering (`keep_sig_pairs`).
7. Match significant events to annotation and sequences (`get_matched_events_chunked`, `attach_sequences`).
8. Build case/control transcript pairs (`get_pairs`).
9. Compute sequence/frame consequences (`compare_sequence_frame`).
10. Call domain changes (`get_domains`) and optional enrichment (`enrich_*`).
11. Infer PPI rewiring (`get_ppi_switches`) and summarize (`integrated_event_summary`, `plot_*`).
12. Build enrichment foreground gene sets from table or S4 inputs.
13. Use S4 container workflows and accessors.
14. Run optional custom-input entry points.

At a high level, each stage narrows from many sample-level rows to event-level effects:
`data` (sample/form rows) -> `res` (DI rows) -> `hits_final` (paired event impact rows).

## Load reference resources
### Load GENCODE annotations
__SpliceImpactR__ requires the accession of various genome annotations, accessed through biomaRt and directly through gencode, 
here we access the gencode files. SpliceImpactR is built to work with either human or mouse data. 
We will initially load a test set
```r
annotation_df <- get_annotation(load = "test")
```
If we were looking to load the full annotations, we'd run the following
```r
annotation_df <- get_annotation(load = "link", species = 'human', release = 45, base_dir = "./")
```
Or to load from cached annotations previously loaded:
```r
annotation_df <- get_annotation(load="cached", base_dir="./path/")
```

### Load protein features (BioMart / ELM / manual)
`get_protein_features()` supports:
- `interpro`: integrated domain/family/superfamily signatures.
- `pfam`: protein domain HMM families.
- `cdd`: NCBI Conserved Domain Database annotations.
- `gene3d`: CATH/Gene3D structural domains.
- `signalp`: signal peptide predictions.
- `tmhmm`: transmembrane helices.
- `ncoils`: coiled-coil predictions.
- `seg`: low-complexity regions.
- `mobidblite`: intrinsically disordered regions.
- `elm`: short linear motifs (SLiMs; loaded from ELM, not BioMart).

For BioMart-backed features, attributes must follow:
`{feature}`, `{feature}_start`, `{feature}_end` (for example
`pfam`, `pfam_start`, `pfam_end`).
Additional custom sources can be added with `get_manual_features()`
(described in more detail in the custom input section below).

We're loading test data here, but set test = FALSE to get the full set.
```r
interpro_features <- get_protein_features(c("interpro"), annotation_df$annotations, test = TRUE)
signalp_features <- get_protein_features(c("signalp"), annotation_df$annotations, test = TRUE)
elm_features <- get_protein_features(c("elm"), annotation_df$annotations, test = TRUE)
```

We can also load user-defined protein features by transcript/protein ensembl ids and the location of the protein feature within 
```r
user_df <- data.frame(
 ensembl_transcript_id = c(
   "ENST00000511072","ENST00000374900","ENST00000373020","ENST00000456328",
   "ENST00000367770","ENST00000331789","ENST00000335137","ENST00000361567",
   NA,                    "ENST00000380152"
 ),
 ensembl_peptide_id = c(
   "ENSP00000426975", NA,                   "ENSP00000362048","ENSP00000407743",
   "ENSP00000356802","ENSP00000326734", NA,                  "ENSP00000354587",
   "ENSP00000364035", NA
 ),
 name = c(
   "Low complexity","Transmembrane helix","Coiled-coil","Signal peptide",
   "Transmembrane helix","Low complexity","Coiled-coil","Transmembrane helix",
   "Signal peptide","Low complexity"
 ),
 start = c(80L, 201L, 35L, 1L, 410L, 150L, 220L, 30L, 1L, 300L),
 stop  = c(120L,223L, 80L, 20L, 430L, 190L, 260L, 55L, 24L, 360L),
 database   = c("seg","tmhmm","ncoils","signalp","tmhmm","seg","ncoils","tmhmm","signalp", NA),
 alt_name   = c(NA,"TMhelix",NA,"SignalP-noTM", "TMhelix", NA, NA, "TMhelix", "SignalP-TAT", NA),
 feature_id = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
)
user_features <- get_manual_features(user_df, gtf_df = annotation_df$annotations)
```

We use this function to combine multiple protein features and the user-defined features. 
If no user_features are added, remove user_features from get_comprehensive_annotations()
We also get the exon-level protein features from the prior overall features.
```r
protein_feature_total <- get_comprehensive_annotations(list(signalp_features, interpro_features, user_features))
exon_features <- get_exon_features(annotation_df$annotations, protein_feature_total)
```

## Build sample manifest
For the standard workflow, each row in `sample_frame` represents one sample.
The sample directory should contain outputs from `rMATS` and/or `HITindex`.
`rMATS` readers look for `{AS}.MATS.JC/JCEC.txt`, and HITindex readers look
for `.AFEPSI`, `.ALEPSI`, and `.exon` files.

Required sample manifest columns:
- `path`: per-sample directory path.
- `sample_name`: unique sample identifier.
- `condition`: case/control label.

Manifest expectations:
- One row per sample.
- Paths point to readable sample directories.
- Exactly two condition labels for the default DI workflow (`case`, `control`).
- Replicates per condition are recommended for stable DI estimation.

For the sake of this intro, we use toy versions (limited to a handful of genes).
The data files should be organized as such for each sample:
```r
print(list.files(file.path(system.file("extdata", package = "SpliceImpactR"), "rawData/control_S5/")))
```

Standard combined-directory example:
```r
sample_frame <- data.frame(path = c(file.path(system.file("extdata", package = "SpliceImpactR"), "rawData/control_S5/"),
                                    file.path(system.file("extdata", package = "SpliceImpactR"), "rawData/control_S6/"),
                                    file.path(system.file("extdata", package = "SpliceImpactR"), "rawData/control_S7/"),
                                    file.path(system.file("extdata", package = "SpliceImpactR"), "rawData/control_S8/"),
                                    file.path(system.file("extdata", package = "SpliceImpactR"), "rawData/case_S1/"),
                                    file.path(system.file("extdata", package = "SpliceImpactR"), "rawData/case_S2/"),
                                    file.path(system.file("extdata", package = "SpliceImpactR"), "rawData/case_S3/"),
                                    file.path(system.file("extdata", package = "SpliceImpactR"), "rawData/case_S4/")),
                           sample_name  = c("S5", "S6", "S7", "S8", "S1", "S2", "S3", "S4"),
                           condition    = c("control", "control", "control", "control", "case",  "case",  "case",  "case"),
                           stringsAsFactors = FALSE)
```

## Quick start wrapper (`get_splicing_impact`)
If your inputs already follow the standard layout, run the full workflow with
the wrapper:

```r
# data.table-first return
out <- get_splicing_impact(
  sample_frame = sample_frame,
  source_data = "both",  # "hitindex" | "rmats" | "both"
  annotation_df = annotation_df,
  protein_feature_total = protein_feature_total,
  return_class = "data.table"
)

data <- out$data
res <- out$res
hits_final <- out$hits_final

# S4-first return
obj <- get_splicing_impact(
  sample_frame = sample_frame,
  source_data = "both",
  annotation_df = annotation_df,
  protein_feature_total = protein_feature_total,
  return_class = "S4"
)
```

## Read splicing events (rmats + hit index example)
If you want stepwise control, load splicing events directly:

```r
data <- get_rmats_hit(
  sample_frame,
  event_types = c("ALE", "AFE", "MXE", "SE", "A3SS", "A5SS", "RI")
)

DT <- data.table::as.data.table(data)
DT[, .(
  n_rows = .N,
  n_events = data.table::uniqueN(event_id),
  n_genes = data.table::uniqueN(gene_id)
)]
```

Alternative when `rMATS` and `HITindex` are in separate directory trees:
```r
sample_frame_rmats <- sample_frame
sample_frame_hit <- sample_frame

# In real usage, set different path columns:
# sample_frame_rmats$path <- c("/rmats/S5/", "/rmats/S6/", ...)
# sample_frame_hit$path <- c("/hitindex/S5/", "/hitindex/S6/", ...)

rmats_only <- get_rmats(
  load_rmats(
    sample_frame_rmats,
    use = "JCEC",
    event_types = c("MXE", "SE", "A3SS", "A5SS", "RI")
  )
)

hit_only <- get_hitindex(
  sample_frame_hit,
  keep_annotated_first_last = TRUE
)

shared_cols <- intersect(names(rmats_only), names(hit_only))
data <- data.table::rbindlist(
  list(
    rmats_only[, ..shared_cols],
    hit_only[, ..shared_cols]
  ),
  use.names = TRUE,
  fill = TRUE
)
```

## Sample-Level QC and Exploration
### Compare HITindex Between Conditions
This summary compares event-level mean HIT values across conditions.

```r
hit_compare <- compare_hit_index(
  sample_frame,
  condition_map = c(control = "control", test = "case")
)

hit_compare$plot
head(hit_compare$results[order(fdr)], 6)
```

### PSI Distribution Overview
This panel summarizes event count depth-normalized metrics and PSI eCDF.

```r
overview_plot <- overview_spicing_comparison(
  events = data,
  sample_df = sample_frame,
  depth_norm = "exon_files",
  event_type = "AFE"
)
overview_plot
```


## Differential Inclusion
We then perform differential inclusion analysis. 
This uses a quasibinomial glm and subsequent F test to identify significant changes in PSI across condition. 
The default here is 10 minimum read count, 
at least present (nonzero) in half of the samples within either of the conditions. 
This step does various filtering and with `verbose = TRUE` prints wrapped,
readable progress lines (no horizontal scrolling needed).

We filter for fdr < 0.05 and delta_psi > 0.1 and output a volcano plot
If using real data, this may take a long time. To speed it up, 
set `parallel_glm = TRUE` and pass a BiocParallel backend via `BPPARAM`
(for example `BiocParallel::MulticoreParam(workers = 4)` on Linux/macOS or
`BiocParallel::SnowParam(workers = 4)` on Windows).
If needed, load the package explicitly with `library(BiocParallel)`.
Using keep_sig_pairs we filter for significance/threshold and plot with 
plot_di_volcano
```r
res <- get_differential_inclusion(
  data,
  min_total_reads = 10,
  parallel_glm = TRUE,
  BPPARAM = BiocParallel::SerialParam()
)
res_di <- keep_sig_pairs(res)
volcano_plot <- plot_di_volcano_dt(res)
volcano_plot
```


## Matching and pairing
Then we match the significant output to annotation. Here, we attach associated transcript and protein sequences and then extract pairs of 'swapping' events.
This matching is done through a strict hierarchy:

1. Prefilter by `chr`, `strand`, and `gene_id` to keep only compatible
   annotation intervals.
2. Match inclusion (`inc`) coordinates to exons and require all inclusion
   parts to be covered for a transcript candidate.
3. Remove candidates that overlap exclusion (`exc`) coordinates.
4. Prioritize transcript/exon choices by event-type-consistent exon class
   (`first`, `internal`, `last`), then reciprocal overlap and intersection
   width; protein-linked transcripts are preferred when available.
5. Build case/control pairs in `get_pairs(source = "multi")` by joining all
   positive `delta_psi` rows (case) with all negative rows (control) for each
   `event_id`, then ordering by strongest `|delta_psi|`.

```r
matched <- get_matched_events_chunked(res_di, annotation_df$annotations, chunk_size = 2000)
hits_sequences <- attach_sequences(matched, annotation_df$sequences)
pairs <- get_pairs(hits_sequences, source="multi")
```

We can also perform analysis looking at how events impact proximal/distal use of terminal exons
```r
proximal_output <- get_proximal_shift_from_hits(pairs)
```

## Inspect PSI for a single event
Use `probe_individual_event()` to visualize PSI distributions for a specific event across samples. For terminal exon events
(`AFE`/`ALE`), PSI is separated by the `inc` entry to highlight proximal vs distal choices.

Identify an event of interest from the differential inclusion results
```r
event_to_probe <- res$event_id[1]
```
Plot PSI by sample/condition; missing combinations are filled with zeros by default
```r
probe <- probe_individual_event(data, event = event_to_probe)
probe
```

If you want to start from transcript pairs directly (instead of DI/matching), see the
`Advanced / custom input workflows` section at the end.

## Primary sequence comparisons
Here, we compare sequence using protein-coding status, sequence alignment percent identity, protein length, and whether frame shifts / rescues are produced.
And make a summary plot
```r
seq_compare <-compare_sequence_frame(pairs, annotation_df$annotations)
alignment_summary <- plot_alignment_summary(seq_compare)
alignment_summary
```

Key labels used in sequence/frame outputs:
- `protein_coding`: both isoforms have protein IDs.
- `onePC`: only one isoform has a protein ID.
- `noPC`: neither isoform has a protein ID.
- `Match`: protein sequences are identical.
- `FrameShift`: reading frame is disrupted between isoforms.

We can also perform analysis looking at how events impact protein length
```r
length_output <- plot_length_comparison(seq_compare)
length_output
```

## Get background
We next must get a background set for domain enrichment analysis. We can do this through 
all annotated transcripts, a given set of possible transcripts, or the hit-index's .exon files
```r
bg <- get_background(source = "annotated",
                     annotations = annotation_df$annotations,
                     protein_features = protein_feature_total)
```


## Get domain changes
Here we identify when the alternative RNA processing event drives a change in protein features, then identify enriched domains using the backgrond set
First get the domains that change across pairs
```r
hits_domain <- get_domains(seq_compare, exon_features)
```

Then we can probe for any enriched domains that are changing and plot
```r
enriched_domains <- enrich_domains_hypergeo(hits_domain, bg, db_filter = 'interpro')
domain_plot <- plot_enriched_domains_counts(enriched_domains, top_n = 20)
domain_plot
```

And we're able to search for A) specific events enrichment (AFE, ALE, etc)
or by database (Interpro, SignalP, etc)
```r
enriched_domains <- enrich_by_event(hits_domain, bg, events = 'AFE', db_filter = 'interpro')
enriched_domains <- enrich_by_db(hits_domain, bg, dbs = 'interpro')
```

## Isoform-Isoform interaction network (Only available for human data currently)
For PPI analysis, we first obtain protein-protein interactions from a pre-derived network built from biogrid, ppidm, and elm
And use these to show when a change in domain / SLiM changes ppi
```r
ppi <- get_ppi_interactions()             
hits_final <- get_ppi_switches(hits_domain, ppi, protein_feature_total)
ppi_plot <- plot_ppi_summary(hits_final)
ppi_plot
```

## Gene Enrichment Foreground Probes (table and S4)
`get_gene_enrichment()` supports both table and S4 input.

- For `mode = "di"`, pass `res`.
- For `mode = "ppi"` / `mode = "domain"`, pass `hits`.
- S4 examples are shown in the S4 section so setup is centralized.

```r
fg_di <- get_gene_enrichment(
  mode = "di",
  res = res,
  padj_threshold = 0.05,
  delta_psi_threshold = 0.1
)
fg_domain <- get_gene_enrichment(mode = "domain", hits = hits_domain)
fg_ppi <- get_gene_enrichment(mode = "ppi", hits = hits_final)
```

If DI enrichment returns all `NA` statistics, the usual issue is sparse
Ensembl-to-Entrez mapping relative to the selected background. In practice,
increase foreground size (relax DI cutoffs), broaden background, and/or lower
`min_size` in `get_enrichment()`.

```r
enrichment_di <- get_enrichment(
  foreground = fg_di,
  background = bg$gene_id,
  species = "human",
  gene_id_type = "ensembl",
  sources = "GO:BP",
  min_size = 5
)
```

## Visualize specific transcript changes
You can visualize one paired event in transcript-centric and protein-centric
views.

```r
viz_pair <- hits_final[
  !is.na(transcript_id_case) &
    !is.na(transcript_id_control) &
    transcript_id_case != "" &
    transcript_id_control != ""
][1]

tx_pair <- c(viz_pair$transcript_id_case, viz_pair$transcript_id_control)

transcript_centric <- plot_two_transcripts_with_domains_unified(
  transcripts = tx_pair,
  gtf_df = annotation_df$annotations,
  protein_features = protein_feature_total,
  feature_db = c("interpro"),
  combine_domains = TRUE,
  view = "transcript"
)

protein_centric <- plot_two_transcripts_with_domains_unified(
  transcripts = tx_pair,
  gtf_df = annotation_df$annotations,
  protein_features = protein_feature_total,
  feature_db = c("interpro"),
  combine_domains = TRUE,
  view = "protein"
)
```

## Integrative visualization
Use `integrated_event_summary()` for a multi-panel overview across DI, sequence,
domain, and PPI layers.

Panel guide:
- Top-left: event classification composition by event type.
- Top-right: alignment score distributions by event type.
- Middle-left: domain-change prevalence (`Any`, case-only, control-only, both).
- Middle-center/right: PPI rewiring prevalence and gain distributions.
- Bottom-left: relative retention from DI input to final integrated hits.
- Bottom-right: gene-level event-type coordination (Jaccard heatmap).

```r
int_summary <- integrated_event_summary(hits_final, res)
int_summary$plot
int_summary$summaries$relative_use
```

### Understanding output columns (`hits_final`, `data`, `res`)
The pipeline returns three core tables in `data.table` mode:
- `data`: sample-level event measurements before differential modeling.
- `res`: differential inclusion results per tested event/site.
- `hits_final`: paired case/control isoform effects with sequence, domain, and PPI annotations.

Suffix convention used throughout:
- `_case` = case-preferred isoform values.
- `_control` = control-preferred isoform values.

#### `hits_final` (integrated event-level output)
Use this table for biological interpretation and downstream plotting.

**1) Event and isoform identifiers**
- `event_id`: event identifier used across all outputs.
- `event_type`: splicing class (`SE`, `A3SS`, `A5SS`, `MXE`, `RI`, `AFE`, `ALE`, `HFE`, `HLE`).
- `gene_id`: Ensembl gene identifier.
- `chr`, `strand`: genomic chromosome and strand.
- `transcript_id_case`, `transcript_id_control`: paired transcript IDs.
- `protein_id_case`, `protein_id_control`: paired protein IDs (if protein-coding).
- `form_case`, `form_control`: row form labels used during pairing.
- `exons_case`, `exons_control`: event exon IDs used for case/control mapping.

**2) Event coordinates and differential statistics**
- `inc_case`, `inc_control`: inclusion coordinate strings for each isoform.
- `exc_case`, `exc_control`: exclusion coordinate strings for each isoform.
- `delta_psi_case`, `delta_psi_control`: signed PSI shift for each side of the pair.
- `p.value_case`, `p.value_control`: differential model p-values.
- `padj_case`, `padj_control`: multiple-testing-adjusted p-values.
- `n_samples_case`, `n_samples_control`: total samples used.
- `n_case_case`, `n_case_control`: case sample counts.
- `n_control_case`, `n_control_control`: control sample counts.

**3) Sequence content and coding context**
- `transcript_seq_case`, `transcript_seq_control`: transcript nucleotide sequences.
- `protein_seq_case`, `protein_seq_control`: translated protein sequences.
- `pc_class`: coding relationship class for the pair.
- Length metrics: `prot_len_*`, `tx_len_*`, `exon_cds_len_*`, `exon_len_*`, and associated `*_diff` / `*_diff_abs` columns.

**4) Alignment and frame classification**
- DNA alignment: `dna_pid`, `dna_score`, `dna_width`.
- Protein alignment: `prot_pid`, `prot_score`, `prot_width`.
- Frame diagnostics: `frame_call`, `rescue`, `frame_check_exon1`, `frame_check_exon2`.
- Final summary label: `summary_classification`.

**5) Domain-level change annotations**
- `domains_exons_case`, `domains_exons_control`: domains mapped on event exons.
- `case_only_domains`, `control_only_domains`: collapsed domain strings unique to each side.
- `case_only_domains_list`, `control_only_domains_list`, `either_domains_list`: list-columns of domain tokens.
- Counts: `case_only_n`, `control_only_n`, `diff_n`.

**6) Predicted interaction rewiring (PPI/DDI/DMI-aware)**
- Partners: `case_ppi`, `control_ppi` (list-columns).
- Counts: `n_case_ppi`, `n_control_ppi`, `n_ppi`.
- Feature drivers: `case_ppi_drivers`, `control_ppi_drivers` (merged PFAM/ELM tokens, prefixed as `pfam;...` or `elm;...`).

#### `data` (raw sample-level input table)
Use `data` to inspect per-sample evidence feeding differential inclusion.

Core columns:
- `event_id`, `event_type`, `form`, `gene_id`, `chr`, `strand`.
- `inc`, `exc`: coordinate strings for inclusion/exclusion forms.
- `inclusion_reads`, `exclusion_reads`: read support.
- `psi`: sample-level PSI value.
- `sample`, `condition`: sample metadata.
- `source_file`: source path used during import.

Often present depending on import path:
- HITindex metadata such as `HITindex`, `class`, `nFE`, `nLE`, `nUP`, `nDOWN`, `nTXPT`, `psi_original`, `total_reads`, `source`.

#### `res` (differential inclusion output)
Use `res` to rank significant events before downstream pairing/domain/PPI steps.

Core columns:
- `site_id`: tested site/event key used by the model.
- `event_id`, `event_type`, `gene_id`, `chr`, `strand`, `inc`, `exc`, `form`.
- `n_samples`, `n_control`, `n_case`: sample counts used.
- `mean_psi_ctrl`, `mean_psi_case`: group PSI means.
- `delta_psi`: case minus control PSI shift.
- `p.value`, `padj`: statistical significance.
- `cooks_max`: maximum Cook's distance seen for the fitted site.
- `n`, `n_used`: total rows and rows retained after model filtering.

## S4 Applications and Accessors
You can run the full pipeline with `get_splicing_impact()` and choose either compact `data.table` outputs or a single S4 object.
This is the canonical place to initialize `obj` for S4 workflows.

`SpliceImpactResult` is a custom S4 container that keeps all major pipeline parts synchronized:
- `raw_events` (`SummarizedExperiment`): sample-level table + ranges/assays.
- `di_events` / `res_di` (`GRanges`): differential inclusion rows.
- `matched` (`DFrame`): annotation-matched rows (and sequence-attached rows).
- `paired_hits` (`GRanges`) + `segments` (`GRangesList`): final case/control event impacts.
- `sample_frame` (`DFrame`): sample manifest metadata.

```r
# End-to-end (combined HITindex + rMATS)
out <- get_splicing_impact(
  sample_frame = sample_frame,
  source_data = "both",                 # "hitindex" | "rmats" | "both"
  event_types = c("ALE", "AFE", "MXE", "SE", "A3SS", "A5SS", "RI"),
  annotation_df = annotation_df,
  protein_feature_total = protein_feature_total,
  return_class = "data.table"
)

# Compact returns in data.table mode
data <- out$data
res <- out$res
hits_final <- out$hits_final
```

```r
# Return a single S4 object
obj <- get_splicing_impact(
  sample_frame = sample_frame,
  source_data = "both",
  annotation_df = annotation_df,
  protein_feature_total = protein_feature_total,
  return_class = "S4"
)

# Convert slots back to data.table
raw_dt <- as_dt_from_s4(obj, "raw_events")
di_dt <- as_dt_from_s4(obj, "di_events")
hits_dt <- as_dt_from_s4(obj, "paired_hits")
```

When to prefer S4:
- You want one object to pass through multiple steps with consistent state.
- You want slot-level validation and synchronized filtering via `filter_spliceimpact_hits()`.
- You want Bioconductor-native structures (`SummarizedExperiment`/`GRanges`) for downstream tooling.

For paired-hit summaries, use fast accessors:
- `get_hits_core()`: identifiers, event metadata, and core comparison fields.
- `get_hits_domain()`: domain gain/loss content and domain-count summaries.
- `get_hits_ppi()`: PPI partner switches and feature-driver columns.
- `get_hits_sequence()`: sequence/alignment/frame and length-delta fields.

```r
# Generic subset accessor
core_dt <- get_hits_final_view(obj, col_subset = "core")
dom_dt <- get_hits_final_view(obj, col_subset = "domain")
ppi_dt <- get_hits_final_view(obj, col_subset = "ppi")
seq_dt <- get_hits_final_view(obj, col_subset = "sequence")

# Tiny wrappers
core_dt <- get_hits_core(obj)
dom_dt <- get_hits_domain(obj)
ppi_dt <- get_hits_ppi(obj)
seq_dt <- get_hits_sequence(obj)
```

To inspect S4 schema and slot usage:

```r
spliceimpact_s4_schema()
spliceimpact_s4_guide()
```

### Filtering an S4 object by `hits_final` columns
Use `filter_spliceimpact_hits()` to subset by any `paired_hits` column and keep all event-linked slots synchronized
(`paired_hits`, `segments`, `res_di`, `di_events`, `matched`, `raw_events`).

```r
# keep one event
obj_one <- filter_spliceimpact_hits(obj, event_id == "A3SS:44")

# keep one gene with PPI change
obj_gene <- filter_spliceimpact_hits(obj, gene_id == "ENSG00000142599", n_ppi > 0)

# keep coding events with domain and frame criteria
obj_focus <- filter_spliceimpact_hits(
  obj,
  pc_class == "protein_coding",
  diff_n > 0,
  frame_call %in% c("Match", "PartialMatch")
)
```

Quick validation and extraction after filtering:

```r
# confirm retained events
as_dt_from_s4(obj_focus, "paired_hits")[, .N, by = event_type][order(-N)]

# get compact outputs
get_hits_core(obj_focus)
get_hits_domain(obj_focus)
get_hits_ppi(obj_focus)
```

Notes:
- Filter expressions are evaluated in `paired_hits` context.
- Multiple expressions are combined with `AND`.
- `sample_frame` is intentionally not filtered (sample metadata stays intact).

### Common S4 applications
Use the same S4 object directly in downstream helpers without converting to
tables.

```r
fg_di_s4 <- get_gene_enrichment(
  mode = "di",
  x = obj,
  padj_threshold = 0.05,
  delta_psi_threshold = 0.1
)

fg_ppi_s4 <- get_gene_enrichment(mode = "ppi", x = obj)

fg_ppi_focus <- get_gene_enrichment(mode = "ppi", x = obj_focus)
fg_di_focus <- get_gene_enrichment(
  mode = "di",
  x = obj_focus,
  padj_threshold = 0.05,
  delta_psi_threshold = 0.1
)

int_summary_focus <- integrated_event_summary(obj_focus, obj_focus)
int_summary_focus$plot
```

### S4-first main workflow (code only)
This mirrors the table workflow but keeps all updates inside one
`SpliceImpactResult` object. The same core functions accept S4 input and
return an updated S4 object when `return_class = "S4"` is set.

```r
obj_flow <- as_splice_impact_result(
  data = data,
  sample_frame = sample_frame
)

obj_flow <- get_differential_inclusion(
  obj_flow,
  min_total_reads = 10,
  parallel_glm = TRUE,
  BPPARAM = BiocParallel::SerialParam(),
  return_class = "S4"
)
obj_flow <- keep_sig_pairs(obj_flow, return_class = "S4")

obj_flow <- get_matched_events_chunked(
  obj_flow,
  annotation_df$annotations,
  chunk_size = 2000,
  return_class = "S4"
)
obj_flow <- attach_sequences(
  obj_flow,
  annotation_df$sequences,
  return_class = "S4"
)
obj_flow <- get_pairs(obj_flow, source = "multi", return_class = "S4")
obj_flow <- compare_sequence_frame(
  obj_flow,
  annotation_df$annotations,
  return_class = "S4"
)
obj_flow <- get_domains(obj_flow, exon_features, return_class = "S4")
obj_flow <- get_ppi_switches(
  obj_flow,
  ppi,
  protein_feature_total,
  return_class = "S4"
)

hits_core_flow <- get_hits_core(obj_flow)
hits_domain_flow <- get_hits_domain(obj_flow)
hits_ppi_flow <- get_hits_ppi(obj_flow)
hits_sequence_flow <- get_hits_sequence(obj_flow)
```

## Advanced / custom input workflows (optional)
Use these entry points when your data starts outside the default
`get_rmats_hit()` to `get_differential_inclusion()` flow.

### Add user-defined protein features (`get_manual_features`)
```r
ann_dt <- data.table::as.data.table(annotation_df$annotations)
coding_tx <- unique(ann_dt[type == "exon" & cds_has == TRUE, transcript_id])
n_manual <- min(3L, length(coding_tx))
stopifnot(n_manual >= 1L)

manual_df <- data.frame(
  ensembl_transcript_id = coding_tx[seq_len(n_manual)],
  ensembl_peptide_id = rep(NA_character_, n_manual),
  name = paste0("demo_feature_", seq_len(n_manual)),
  start = c(20L, 45L, 80L)[seq_len(n_manual)],
  stop = c(35L, 58L, 92L)[seq_len(n_manual)],
  database = rep("manual", n_manual),
  alt_name = rep(NA_character_, n_manual),
  feature_id = rep(NA_character_, n_manual),
  stringsAsFactors = FALSE
)

manual_features <- get_manual_features(
  manual_features = manual_df,
  gtf_df = annotation_df$annotations
)
```

### Bring your own pre-DI event table
```r
example_df <- data.frame(
  event_id = rep("A3SS:1", 8),
  event_type = "A3SS",
  form = rep(c("INC", "EXC"), each = 4),
  gene_id = "ENSG00000158286",
  chr = "chrX",
  strand = "-",
  inc = c(rep("149608626-149608834", 4), rep("149608626-149608829", 4)),
  exc = c(rep("", 4), rep("149608830-149608834", 4)),
  inclusion_reads = c(30, 32, 29, 31, 2, 3, 4, 3),
  exclusion_reads = c(1, 1, 2, 1, 28, 27, 26, 30),
  sample = c("S1", "S2", "S3", "S4", "S1", "S2", "S3", "S4"),
  condition = rep(c("case", "case", "control", "control"), 2),
  stringsAsFactors = FALSE
)
example_df$psi <- example_df$inclusion_reads /
  (example_df$inclusion_reads + example_df$exclusion_reads)

user_data <- get_user_data(example_df)
```

### Bring your own post-DI table
```r
example_user_data <- data.frame(
  event_id = rep("A3SS:1", 8),
  event_type = "A3SS",
  gene_id = "ENSG00000158286",
  chr = "chrX",
  strand = "-",
  form = rep(c("INC", "EXC"), each = 4),
  inc = c(rep("149608626-149608834", 4), rep("149608626-149608829", 4)),
  exc = c(rep("", 4), rep("149608830-149608834", 4)),
  stringsAsFactors = FALSE
)

user_res <- get_user_data_post_di(example_user_data)
```

### Import rMATS post-DI results directly
Multiple files:
```r
input <- data.frame(
  path = c("/path/A3SS.MATS.JC.txt", "/path/A5SS.MATS.JC.txt"),
  event_type = c("A3SS", "A5SS"),
  stringsAsFactors = FALSE
)

# res_rmats_di <- get_rmats_post_di(input)
```

Single preloaded rMATS table:
```r
rmats_df <- data.frame(
  ID = 1L,
  GeneID = "ENSG00000182871",
  geneSymbol = "COL18A1",
  chr = "chr21",
  strand = "+",
  longExonStart_0base = 45505834L,
  longExonEnd = 45505966L,
  shortES = 45505837L,
  shortEE = 45505966L,
  flankingES = 45505357L,
  flankingEE = 45505431L,
  ID.2 = 2L,
  IJC_SAMPLE_1 = "1,1,1",
  SJC_SAMPLE_1 = "1,1,1",
  IJC_SAMPLE_2 = "1,1,1",
  SJC_SAMPLE_2 = "1,1,1",
  IncFormLen = 52L,
  SkipFormLen = 49L,
  PValue = 0.6967562,
  FDR = 1,
  IncLevel1 = "0.0,0.0,0.0",
  IncLevel2 = "1.0,1.0,1.0",
  IncLevelDifference = 1.0,
  stringsAsFactors = FALSE
)

res_rmats_di <- get_rmats_post_di(rmats_df, event_type = "A3SS")
```

### Start from transcript pairs instead of event matching
```r
tx_ids <- unique(annotation_df$annotations$transcript_id)
tx_ids <- tx_ids[!is.na(tx_ids) & tx_ids != ""]
stopifnot(length(tx_ids) >= 4L)

transcript_pairs <- data.frame(
  transcript1 = tx_ids[1:2],
  transcript2 = tx_ids[3:4],
  stringsAsFactors = FALSE
)

user_matched <- compare_transcript_pairs(
  transcript_pairs = transcript_pairs,
  annotations = annotation_df$annotations
)
```

## Contributing
Contributions to SpliceImpactR are welcome, including bug reports, feature requests, and pull requests. Please see CONTRIBUTING.md for guidelines on how to contribute.

## Support
If you encounter any problems or have suggestions, please file an issue on the GitHub issue tracker. Or contact zachpw@bu.edu

##Citation
If you use SpliceImpactR in your research, please cite:

```bibtex
Zachary Peters Wakefield, Ana Fiszbein
SpliceImpactR maps alternative RNA processing events driving protein functional diversity
2025
https://www.biorxiv.org/content/10.1101/2025.06.20.660706v1
https://github.com/fiszbein-lab/SpliceImpactR
```
