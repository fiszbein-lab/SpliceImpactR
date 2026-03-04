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
You can install SpliceImpactR directly from GitHub using the devtools package. If you do not have devtools installed, you can install it first with:

```r
install.packages("devtools")
```
Then, to install the package (devtools installation recommnded):
```r
devtools::install_github("fiszbein-lab/SpliceImpactR")
```
Or 
``` r
BiocManager::install("fiszbein-lab/SpliceImpactR", version="devel")
```




# Usage
## Access gencode information
__SpliceImpactR__ requires the acceession of various genome annotations, accessed through biomaRt and directly through gencode, 
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

## Access biomaRt information
Here we obtain further annotations through biomaRt. The typical protein features accessed through biomaRt are: interpro, pfam, 
signalp, ncoils, mobidblite. We also access the eukaryotic linear motif (ELM) database of short linear motifs.
Any features added are able to be accessed if they have the three attributes (biomaRt::listAttributes(mart)): 
{feature}, {feature}_start, {feature}_end. If there are more protein features desired, you can manually access them and upload through 
get_manual_features() shown in the chunk below using peptide locations.

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

## Loading data (rmats + hit index example)
For the sake of this intro, we use toy versions (limited to a handful of genes)
The sample data frame must have a path column pointing to where the files (rMATS output and hit_index is contained). We must also have sample_name and condition columns
For the standard workflow, we require all output files within the same directory for each sample. rMATS analysis will look for the {AS}.MATS.JC/JCEC.txt and HIT Index will look for the .AFEPSI, .ALEPSI, and .exon files
The data files should be organized as such for each sample:
```r
print(list.files(check_extdata_dir('rawData/control_S5/')))
```

If the rmats and hit index output are in separate directories, you can use get_hitindex() and get_rmats(),
then rbind on the shared columns (to avoid reorganizing data files)
```r
sample_frame <- data.frame(path = c(check_extdata_dir('rawData/control_S5/'),
                                    check_extdata_dir('rawData/control_S6/'),
                                    check_extdata_dir('rawData/control_S7/'),
                                    check_extdata_dir('rawData/control_S8/'),
                                    check_extdata_dir('rawData/case_S1/'),
                                    check_extdata_dir('rawData/case_S2/'),
                                    check_extdata_dir('rawData/case_S3/'),
                                    check_extdata_dir('rawData/case_S4/')),
                           sample_name  = c("S5", "S6", "S7", "S8", "S1", "S2", "S3", "S4"),
                           condition    = c("control", "control", "control", "control", "case",  "case",  "case",  "case"),
                           stringsAsFactors = FALSE)

data <- get_rmats_hit(sample_frame, event_types = c("ALE", "AFE", "MXE", "SE", "A3SS", "A5SS", "RI"))

```

Now that the data is loaded, we can make some initial comparisons. Such as comparing the hit index values across condition,
if hit index is supplied. This identifies how classification/use of exons may change
Or an overview of how the conditions compare as well. 
Here, we probe whether there are changes in depth-normalized counts of AFE or the distribution of PSI across condition
```r
hit_compare <- compare_hit_index(sample_frame, condition_map = c(control = "control", test = "case"))

ov <- overview_splicing_comparison_fixed(data, 
                                         sample_frame, 
                                         depth_norm = 'exon_files', 
                                         event_type = "AFE")
```


## Differential Inclusion
We then perform differential inclusion analysis. 
This uses a quasibinomial glm and subsequent F test to identify significant changes in PSI across condition. 
The default here is 10 minimum read count, 
at least present (nonzero) in half of the samples within either of the conditions. 
This step does various filtering and with verbose = TRUE prints out how many events are filtered / kept.

We filter for fdr < 0.05 and delta_psi > 0.1 and output a volcano plot
If using real data, this may take a long time. To speed it up, 
add parallel_glm = TRUE and set cores_glm > 1
Using keep_sig_pairs we filter for significance/threshold and plot with 
plot_di_volcano
```r
res <- get_differential_inclusion(data, min_total_reads = 10)
res_di <- keep_sig_pairs(res)
volcano_plot <- plot_di_volcano_dt(res)
```

We can also load user-supplied data from any source through: get_user_data and get_user_data_post_di. Here we use minimal example data frames
pre-differential-inclusion:
```r
example_df <- data.frame(
  event_id = rep("A3SS:1", 8),
  event_type = "A3SS",
  form = rep(c("INC","EXC"), each = 4),
  gene_id = "ENSG00000158286",
  chr = "chrX",
  strand = "-",
  inc = c(rep("149608626-149608834",4), rep("149608626-149608829",4)),
  exc = c(rep("",4), rep("149608830-149608834",4)),
  inclusion_reads = c(30,32,29,31, 2,3,4,3),
  exclusion_reads = c(1,1,2,1, 28,27,26,30),
  sample = c("S1","S2","S3","S4","S1","S2","S3","S4"),
  condition = rep(c("case","case","control","control"), 2),
  stringsAsFactors = FALSE
)
example_df$psi <- example_df$inclusion_reads / example_df$exclusion_reads
user_data <- get_user_data(example_df)
```
post-differential-inclusion:
```r
example_user_data <- data.frame(
  event_id = rep("A3SS:1", 8),
  event_type = "A3SS",
  gene_id = "ENSG00000158286",
  chr = "chrX",
  strand = "-",
  form = rep(c("INC","EXC"), each = 4),
  inc = c(
    rep("149608626-149608834", 4),
    rep("149608626-149608829", 4)
  ),
  exc = c(
    rep("", 4),
    rep("149608830-149608834", 4)
  ),
  stringsAsFactors = FALSE
)

user_data <- get_user_data_post_di(example_user_data)
```

SpliceImpactR also handles data with differential inclusion values generated
by rMATS
Loading multiple splice types are loaded through:
```r
input <- data.frame(
  path = c('/path/A3SS.MATS.JC.txt', '/path2/A5SS.MATS.JC.txt'),
  grp1 = c("WT","WT"),
  grp2 = c("KO","KO"),
  event_type = c("A3SS", "A5SS")
)
# res_rmats_di <- get_rmats_post_di(input)
```

We can also load from a single rMATS table, preloaded:
```r
df <- data.frame(
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
user_res <- get_rmats_post_di(df, event_type = "A3SS")
```


## Matching and pairing
Then we match the significant output to annotation. Here, we attach associated transcript and protein sequences and then extract pairs of 'swapping' events.
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
```

If you already have sets of transcripts you want to compare, you can feed them into 
compare_transcript_pairs and get output at the matched stage of the workflow. This will
position you to compare all protein features downstream.
```r
annotation_df <- get_annotation(load = "test")
transcript_pairs <- data.frame(
    transcript1 = c("ENST00000337907", "ENST00000426559"),
    transcript2 = c("ENST00000400908", "ENST00000399728")
)
user_matched <- compare_transcript_pairs(transcript_pairs, annotation_df$annotations)
```

## Primary sequence comparisons
Here, we compare sequence using protein-coding status, sequence alignment percent identity, protein length, and whether frame shifts / rescues are produced.
And make a summary plot
```r
seq_compare <-compare_sequence_frame(pairs, annotation_df$annotations)
alignment_summary <- plot_alignment_summary(seq_compare)
```

We can also perform analysis looking at how events impact protein length
```r
length_output <- plot_length_comparison(seq_compare)
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
```

We can also probe for enrichment of relevant genes:
```r
enrichment_ppi <- get_enrichment(get_ppi_gene_enrichment(hits_final), bg$gene_id, species = 'human', 'ensembl', 'GO:BP')
enrichment_domain <- get_enrichment(get_domain_gene_for_enrichment(hits_domain), bg$gene_id, species = 'human', 'ensembl', 'GO:BP')
enrichment_di <- get_enrichment(get_di_gene_enrichment(res, .05, .1), bg$gene_id, species = 'human', 'ensembl', 'GO:BP')
```

## Visualize specific transcript changes
Here we can look at how protein features are changing across the matched 
transcripts / proteins in two views: transcript (genomic-oriented) and protein
(amino-acid-oriented). This orients everything reading from L-R, so - strand is
reversed before visualization
```
transcript_centric <- plot_two_transcripts_with_domains_unified(
  transcripts = c("ENST00000371514", "ENST00000430330"),
  gtf_df = annotation_df$annotations,
  protein_features = protein_feature_total,
  feature_db = c("interpro"),
  combine_domains = TRUE,
  view = "transcript"
)

protein_centric <- plot_two_transcripts_with_domains_unified(
  transcripts = c("ENST00000371514", "ENST00000430330"),
  gtf_df = annotation_df$annotations,
  protein_features = protein_feature_total,
  feature_db = c("interpro"),
  combine_domains = TRUE,
  view = "protein"
)
```

## Integrative Analysis
Here we identify some holistic patterns / integrative analysis using all event types
```
int_summary <- integrated_event_summary(hits_final, res)
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
- Feature drivers: `case_pfam_changed`, `control_pfam_changed`, `case_elm_changed`, `control_elm_changed`.

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

## S4 Workflow and Accessors
You can run the full pipeline with `get_splicing_impact()` and choose either compact `data.table` outputs or a single S4 object.

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

For paired-hit summaries, use fast accessors:

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

