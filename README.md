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
Then, to install the package:
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
```r
## Loading annotations (if they aren't previously cached takes a bit of time)
## We will initally load a test set
annotation_df <- get_annotation(load = "test")

## If we were looking to load the full annotations, we'd run the following (or load from paths of already downloaded gtf/fa files)
# annotation_df <- get_annotation(load = "link", species = 'human', release = 45, base_dir = "./")

## After the initial lengthy loading of annotations, we could quickly load from cached rds files
# annotation_df <- get_annotation(load="cached", base_dir="./path/")
```

## Access biomaRt information
Here we obtain further annotations through biomaRt. The typical protein features accessed through biomaRt are: interpro, pfam, 
signalp, ncoils, mobidblite. Any features added are able to be accessed if they have the three attributes (biomaRt::listAttributes(mart)): 
{feature}, {feature}_start, {feature}_end. If there are more protein features desired, you can manually access them and upload through 
get_manual_features() shown in the chunk below using peptide locations.
```r
## We're loading test data here, but set test = FALSE to get the full set
interpro_features <- get_protein_features(c("interpro"), annotation_df$annotations, timeout = 600, test = TRUE)
signalp_features <- get_protein_features(c("signalp"), annotation_df$annotations, timeout = 600, test = TRUE)


## When loading multiple features from biomaRt, we suggest loading in separate get_protein_features calls for each individual feature database
# interpro_features <- get_protein_features(c("signalp"), annotation_df$annotations, timeout = 600)
# interpro_features <- get_protein_features(c("interpro"), annotation_df$annotations, timeout = 600)

## We can also load user-defined protein features by transcript/protein ensembl ids and the location of the protein feature within 
# user_df <- data.frame(
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
# user_features <- get_manual_features(user_df, gtf_df = annotation_df$annotations)

## We use this function to combine multiple protein features and the user-defined features
protein_feature_total <- get_comprehensive_annotations(list(signalp_features, interpro_features))


## Finally, we get the exon-level protein features from the prior overall features
exon_features <- get_exon_features(annotation_df$annotations, protein_feature_total)

```

## Loading data (rmats + hit index example)
For the sake of this intro, we use toy versions (limited to a handful of genes)
The sample data frame must have a path column pointing to where the files (rMATS output and hit_index is contained). We must also have sample_name and condition columns
```r
# For purposes of these examples, data directory is in extdata to point to toy data, however this should be replaced with the directory path to  data

# For the standard workflow, we require all output files within the same directory for each sample. rMATS analysis will look for the {AS}.MATS.JC/JCEC.txt and HIT Index will look for the .AFEPSI, .ALEPSI, and .exon files
# The data files should be organized as such for each sample:
print(list.files(check_extdata_dir('rawData/control_S5/')))

# If the rmats and hit index output are in separate directories, you can use get_hitindex() and get_rmats() to avoid reorganizing data files
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
# Here we load both rmats and hit index data
data <- get_rmats_hit(sample_frame, event_types = c("ALE", "AFE", "MXE", "SE", "A3SS", "A5SS", "RI"))

## We can compare the HIT Index values across condition, to identify how classification/use of exons may change
hit_compare <- compare_hit_index(sample_frame, condition_map = c(control = "control", test = "case"))

## We can plot an overview of how the conditions compare as well. Here, we probe whether there are changes in depth-normalized counts of AFE or the distribution of PSI across condition
ov <- overview_splicing_comparison_fixed(data, 
                                         sample_frame, 
                                         depth_norm = 'exon_files', 
                                         event_type = "AFE")


```


## Differential Inclusion
First we perform differential inclusion analysis. This uses a quasibinomial glm and subsequent F test to identify significant changes in PSI across condition. The default here is 10 minimum read count, at least present (nonzero) in half of the samples within either of the conditions. This step does various filtering and with verbose = TRUE prints out how many events are filtered / kept.

We filter for fdr < 0.05 and delta_psi > 0.1 and output a volcano plot
```r
## If not using the toy data, we can add cores_glm = n, parallel_glm = TRUE to run the statistical analysis chunked in parallel
res <- get_differential_inclusion(data, min_total_reads = 10)
res_di <- keep_sig_pairs(res)
volcano_plot <- plot_di_volcano_dt(res)
```

```r
## We can also load user-supplied data from any source through: get_user_data and get_user_data_post_di. Here we use minimal example data frames
# pre-di
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

# post-di
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

## We are also able to load post differential inclusion rMATS data if this data
## was supplied
# input <- data.frame(
#   path = c('/path/A3SS.MATS.JC.txt', '/path2/A5SS.MATS.JC.txt'),
#   grp1 = c("WT","WT"),
#   grp2 = c("KO","KO"),
#   event_type = c("A3SS", "A5SS")
# )
# res <- get_rmats_post_di(meta)
```

## We can also load from after rmats differential inclusion
```r
# # Single rMATS table already loaded as df, see documentation for taking multiple as a data frame or more details on inputs
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

## We can also perform analysis looking at how events impact proximal/distal use of terminal exons
proximal_output <- get_proximal_shift_from_hits(pairs)
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

## We can also perform analysis looking at how events impact protein length
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
```r
## First get the domains that change across pairs
hits_domain <- get_domains(seq_compare, exon_features)

## Then we can probe for any enriched domains that are changing and plot
enriched_domains <- enrich_domains_hypergeo(hits_domain, bg, db_filter = 'interpro')
domain_plot <- plot_enriched_domains_counts(enriched_domains, top_n = 20)

## And we're able to search for A) specific events enrichment (AFE, ALE, etc)
## or by database (Interpro, SignalP, etc)
enriched_domains <- enrich_by_event(hits_domain, bg, events = 'AFE', db_filter = 'interpro')
enriched_domains <- enrich_by_db(hits_domain, bg, dbs = 'interpro')

```

## Isoform-Isoform interaction network
For PPI analysis, we first obtain protein-protein interaction domain miner domain-domain interactions (ppidm)
And use those ddis to predict ppis 
```r
## First we load from ppidm
ppidm <- get_ppidm(test=TRUE)

## Then we identify the predicted ppi, this can take some while. We should restrict to transcripts that we've found in the samples for computational efficiency
restrict_protein_features <- protein_feature_total[ensembl_transcript_id %in% 
                                    c(hits_domain$transcript_id_exc, hits_domain$transcript_id_inc)]
ppi <- get_isoform_interactions(restrict_protein_features, ppidm, save = FALSE, load_dir = NULL, init = TRUE)

## Now we can identify ppi switches and plot
hits_final <- get_ppi_switches(hits_domain, ppi)
ppi_plot <- plot_ppi_summary(hits_final)

## We can also probe for enrichment of relevant genes:
# enrichment_ppi <- get_enrichment(get_ppi_gene_enrichment(hits_final), bg$gene_id, species = 'human', 'ensembl', 'GO:BP')
# enrichment_domain <- get_enrichment(get_domain_gene_for_enrichment(hits_domain), bg$gene_id, species = 'human', 'ensembl', 'GO:BP')
# enrichment_di <- get_enrichment(get_di_gene_enrichment(res, .05, .1), bg$gene_id, species = 'human', 'ensembl', 'GO:BP')
```


## Integrative Analysis
Here we identify some holistic patterns / integrative analysis using all event types
```
int_summary <- integrated_event_summary(hits_final, res)
```

### Understanding the `hits_final` columns
The `hits_final` table contains the paired isoform metadata together with alignment, domain, and PPI annotations that flow through the workflow. Key columns include:

* **Event metadata** – `event_id` plus the paired exon IDs (`exons_inc`, `exons_exc`) and event type labels (`event_type_inc`, `event_type_exc`) used throughout downstream comparisons. Transcript identifiers for each isoform (`transcript_id_inc`, `transcript_id_exc`) anchor all subsequent annotations.
* **Sequence content and alignment metrics** – nucleotide and protein sequences for inclusion and exclusion isoforms (`transcript_seq_inc/exc`, `protein_seq_inc/exc`) together with coding status (`pc_class`), protein and transcript lengths, exon/CDS length differences, and alignment statistics (`dna_pid/score/width`, `prot_pid/score/width`).
* **Frame-shift classification** – frame comparison and rescue outcomes (`frame_call`, `rescue`) and the consolidated `summary_classification` label used in plots (e.g., `Match`, `FrameShift`, `Rescue`, or inherited protein-coding class).
* **Domain changes** – domains observed on the event exons for inclusion vs. exclusion isoforms (`domains_exons_inc`, `domains_exons_exc`); isoform-specific domain sets (`inc_only_domains`, `exc_only_domains`), list-columns holding the underlying identifiers, and counts of changed domains (`inc_only_n`, `exc_only_n`, `diff_n`).
* **Predicted PPI switches** – partners unique to the inclusion or exclusion isoform (`inc_ppi`, `exc_ppi`) and their counts (`n_inc_ppi`, `n_exc_ppi`, with `n_ppi` as the total number of altered interactions).

These columns provide the necessary context to trace how an alternative splicing event alters coding potential, protein domains, and predicted interaction partners.

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

