# SpliceImpactR <img src="./inst/screenshot1.png" alt="Fiszbein Lab Logo" width="110" align="right"/> <img src="./inst/screenshot2.png" width="200" align="right"/>

SpliceImpactR is an R package designed for studying the impact of alternative splicing on protein structure and function. It provides tools for analyzing RNA-seq data to identify differentially included splicing events and predict their consequences on the resulting protein products. SpliceImpactR output involves identifying key changes in proteins at various levels: primary sequence, domain content, and transcript-transcript interactions.

The suite of funcitons is designed to anaylyze the consequences of AFE, ALE, SE, MXE, A5SS, RI, and A3SS, along with hybrid exons. SpliceImpactR is built to take output from any source with custom functions to take data processed by from the [HIT Index](https://github.com/thepailab/HITindex) and [rMATS](https://github.com/Xinglab/rmats-turbo). 

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
Then, to install the package:
```r
devtools::install("fiszbein-lab/SpliceImpactR")
```
Or 
``` r
BiocManager::install("fiszbein-lab/SpliceImpactR", version="devel")
```




# Usage
## Access gencode information
__SpliceImpactR__ requires the acceession of various genome annotations, accessed through biomaRt and directly through gencode, here we access the gencode files
```r
## Loading annotations (if they aren't previously cached takes a bit of time)
## We will initally load a test set
annotation_df <- get_annotation(load = "link")

## If we were looking to load the full annotations, we'd run the following (or load from paths of already downloaded gtf/fa files)
# annotation_df <- get_annotation(load = "link", species = 'human', version = 45, base_dir = "/.")

## After the initial lengthy loading of annotations, we could quickly load from cached rds files
# annotation_df <- get_annotation(load="cached", base_dir="/.")
```

## Access biomaRt information
Here we obtain further annotations through biomaRt
```r
## We're loading test data here, but set test = FALSE to get the full set
interpro_features <- get_protein_features(c("interpro"), annotation_df$annotations, timeout = 600, test = TRUE)
signalp_features <- get_protein_features(c("signalp"), annotation_df$annotations, timeout = 600, test = TRUE)


## When loading multiple features from biomaRt, we suggest loading in separate get_protein_features calls for each individual feature database
# interpro_features <- get_protein_features(c("signalp"), annotation_df$annotations, timeout = 600)
# interpro_features <- get_protein_features(c("interpro"), annotation_df$annotations, timeout = 600)

## We can also load user-defined protein features by transcript/protein ensembl ids and the location of the protein feature within 
# user_df <- data.table(
#  ensembl_transcript_id = c(
#    "ENST00000511072","ENST00000374900","ENST00000373020","ENST00000456328",
#    "ENST00000367770","ENST00000331789","ENST00000335137","ENST00000361567",
#    NA,                    "ENST00000380152"
#  ),
#  ensembl_peptide_id = c(
#    "ENSP00000426975", NA,                   "ENSP00000362048","ENSP00000407743",
#    "ENSP00000356802","ENSP00000326734", NA,                  "ENSP00000354587",
#    "ENSP00000364035", NA
#  ),
#  name = c(
#    "Low complexity","Transmembrane helix","Coiled-coil","Signal peptide",
#    "Transmembrane helix","Low complexity","Coiled-coil","Transmembrane helix",
#    "Signal peptide","Low complexity"
#  ),
#  start = c(80L, 201L, 35L, 1L, 410L, 150L, 220L, 30L, 1L, 300L),
#  stop  = c(120L,223L, 80L, 20L, 430L, 190L, 260L, 55L, 24L, 360L),
#  database   = c("seg","tmhmm","ncoils","signalp","tmhmm","seg","ncoils","tmhmm","signalp", NA),
#  alt_name   = c(NA,"TMhelix",NA,"SignalP-noTM", "TMhelix", NA, NA, "TMhelix", "SignalP-TAT", NA),
#  feature_id = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
# )
# user_features <- get_manual_features(user_df)

## We use this function to combine multiple protein features and the user-defined features
protein_feature_total <- get_comprehensive_annotations(list(signalp_features, interpro_features))


## Finally, we get the exon-level protein features from the prior overall features
exon_features <- get_exon_features(annotation_df$annotations, protein_feature_total)

```

## Loading data (rmats + hit index example)
For the sake of this intro, we use toy versions (limited to a handful of genes)
The sample data frame must have a path column pointing to where the files (rMATS output and hit_index is contained). We must also have sample_name and condition columns
```r
# For purposes of these examples, data directory is in extdata, however this should be replaced with the directory path to  data
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
data <- get_rmats_hit(sample_frame)

## We can plot an overview of how the conditions compare as well
ov <- overview_splicing_comparison_fixed(hit_index, sample_frame, 'exon_files')

```


## Differential Inclusion
First we perform differential inclusion analysis. This uses a quasibinomial glm and subsequent F test to identify significant changes in PSI across condition.
We filter for fdr < 0.05 and delta_psi > 0.1 and can output a volcano plot
```
res <- get_differential_inclusion(data)
res_di <- keep_sig_pairs(res)
volcano_plot <- plot_di_volcano_dt(res)
```

## Matching and pairing
Then we match the signficant output to annotation. Here, we attach associated transcript and protein sequences and then extract pairs of 'swapping' events.
```
matched <- get_matched_events_chunked(res_di, annotation_df$annotations, chunk_size = 2000)
x_seq <- attach_sequences(matched, annotation_df$sequences)
pairs <- get_pairs(x_seq, source="multi")

## We can also perform analysis looking at how events impact protein length and 
## proximal/distal use of terminal exons
proximal_output <- get_proximal_shift_from_hits(pairs)
length_output <- plot_length_comparison(seq_compare)

```

## Primary sequence comparisons
Here, we compare sequence using protein-coding status, sequence alignment percent identity, changing length, and whether frame shifts are produced.
And make a summary plot
```
seq_compare <-compare_sequence_frame(pairs, annotation_df$annotations)
alignment_summary <- plot_alignment_summary(seq_compare)

```

## Get background
We next must get a background set for domain enrichment analysis
```
bg <- get_background(source = "annotated",
                     annotations = annotation_df$annotations,
                     protein_features = protein_feature_total)
```

## Get domain changes
```
## First get the domains that change across pairs
hits_domain <- get_domains(seq_compare, exon_features)

## Then we can probe for any enriched domains that are changing and plot
enriched_domains <- enrich_domains_hypergeo(hits_domain, bg, db_filter = 'interpro')
domain_plot <- plot_enriched_domains_counts(enriched_domains, top_n = 20)

## And we're able to search for A) specific events enrichment (AFE, ALE, etc)
## or by database (Interpro, SignalP, etc)
# enriched_domains <- enrich_by_event(hits_domain, bg, events = 'AFE', db_filter = 'interpro')
# enriched_domains <- enrich_by_db(hits_domain, bg, dbs = 'interpro')

```

## Isoform-Isoform interaction network
For PPI analysis, we first obtain protein-protein interaction domain miner (ppidm)
And use those domain-domain interactions to predict ppis 
```
## First we load from ppidm
ppidm <- get_ppidm(test=TRUE)

## Then we identify the predicted ppi
ppi <- get_isoform_interactions(protein_feature_total, ppidm, save = FALSE, load_dir = NULL, init = TRUE)

## Now we can identify ppi switches and plot
hits_final <- get_ppi_switches(hits_domain, ppi)
ppi_plot <- plot_ppi_summary(hits_final)

```


## Integrative Analysis
Here we identify some holistic patterns / integrative analysis using all event types
```
int_summary <- integrated_event_summary(hits_domain, res)
```

## Contributing
Contributions to SpliceImpactR are welcome, including bug reports, feature requests, and pull requests. Please see CONTRIBUTING.md for guidelines on how to contribute.

## Support
If you encounter any problems or have suggestions, please file an issue on the GitHub issue tracker. Or contact zachpw@bu.edu

##Citation
If you use SpliceImpactR in your research, please cite:

```bibtex
Zachary Wakefield
SpliceImpactR: Analyzing the Impact of Alternative Splicing on Protein Structure and Function
2025
https://github.com/fiszbein-lab/SpliceImpactR
```

