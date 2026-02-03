# Data within /extData/rawData was generated through manually generating a
# limited number of alternative RNA processing events across AFE ALE exon A3SS
# A5SS MXE RI SE that
# span events which will be significantly different across condition and events
# which will not be significantly different across conditions. Along with this,
# we generated .exon files to go with the AFE data.
#
# Having both significant events ensures the rest of the
# function of the package will be able to be minimally used.
#
# We also have limited datasets for annotations, protein features, ppidm setup
# for tests, examples, and vignettes for efficient runtime


## For PPI generation from DMI and DDI and Biogrid:
bg <- fread("/Users/zacharywakefield/Downloads/BIOGRID-ALL-5.0.253.mitab 2.txt", nThread = 8)

extract_genes <- function(x) {
  m <- regexec("entrez gene/locuslink:([^|]+)", x, perl = TRUE)
  r <- regmatches(x, m)
  vapply(r, function(z) if (length(z) >= 2) z[2] else NA_character_, character(1))
}

bg[, geneA := extract_genes(`Alt IDs Interactor A`)]
bg[, geneB := extract_genes(`Alt IDs Interactor B`)]
bg_hgnc <- bg[`Taxid Interactor A` == 'taxid:9606' & `Taxid Interactor B` == "taxid:9606", .(gA =geneA, gB = geneB)]
bg_min <- unique(bg_hgnc[, .(geneA = pmin(gA, gB), geneB = pmax(gA, gB))])

annotation_df <- get_annotation(load = "cached", species = 'human', release = 45, base_dir = "./")
annotation_dt <- as.data.table(annotation_df$annotations)

g2p <- unique(annotation_dt[
  !is.na(gene_name) & nzchar(gene_name) &
    !is.na(protein_id) & nzchar(protein_id),
  .(gene = gene_name, ensg = gene_id)
])

incA <- g2p[bg_min, on = .(gene = geneA)][, geneA := gene][, gene := NULL][, ensg_a := ensg][, ensg := NULL]
bg_ensg <- g2p[incA, on = .(gene = geneB)][, geneB := gene][, gene := NULL][, ensg_b := ensg][, ensg := NULL]

biogrid <- bg_ensg[!is.na(ensg_a) & !is.na(ensg_b)][, .(geneA = ensg_a, geneB = ensg_b, biogrid = TRUE)]




## 2
# elm_instances gives protein name, primary uniprot accession, start/end of SLiM -> can load like a protein feature instead of specialized thing
elm_taxon <- function(species) {
  species <- match.arg(species, c("human","mouse"))
  if (species == "human") "Homo sapiens" else "Mus musculus"
}

taxon <- elm_taxon('human')
url_instances <- paste0(
  "http://elm.eu.org/instances.tsv?q=*&taxon=",
  URLencode(taxon, reserved = TRUE),
  "&instance_logic=true%20positive"
)
elm_instances <- fread(url_instances, showProgress = FALSE, skip = 5, nThread = 4)

url_classes <- "http://elm.eu.org/elms/elms_index.tsv"
elm_classes <- fread(url_classes, showProgress = FALSE, skip = 5, nThread = 4)

url_interactions <- "http://elm.eu.org/interactions/as_tsv"
elm_interactions <- fread(url_interactions, showProgress = FALSE,  nThread = 4)
ei <- unique(elm_interactions[, .(feature_id = Elm, domain_id = Domain)])

pfam_ei_int <- ei[grepl("PF", domain_id)]
ip_ei_int <- ei[grepl("IPR", domain_id)]

library(PFAM.db)
x <- PFAM.db::PFAMINTERPRO2AC
mk <- mappedkeys(x)
xx <- as.list(x[mk])
pfConvert <- rbindlist(lapply(seq_len(nrow(ip_ei_int)), function(i) {
  y <- ip_ei_int$domain_id[i]
  slimIn <- ip_ei_int$feature_id[i]
  out <- xx[y]
  names(out) <- NULL
  if (is.null(unlist(out))) {
    data.frame(row = i, feature_id = NA, domain_id = NA, type = "NULL")
  } else if (length(unlist(out)) > 1) {
    data.frame(row = i, feature_id = slimIn, domain_id = unlist(out), type = "multi")
  } else {
    data.frame(row = i, feature_id = slimIn, domain_id = unlist(out), type = "single")
  }
}))[type != "NULL", .(feature_id, domain_id)]
pf_slim <- unique(rbind(pfConvert, pfam_ei_int))
lms <- get_linear_motifs(annotation_df$annotations, annotation_df$sequences, species = 'hsapiens_gene_ensembl', release = c(109, 115)[1])
gp <- unique(annotation_df$annotations[!is.na(transcript_id), .(gene_id, transcript_id)])
gene_slim <- unique(gp[lms, on = .(transcript_id = ensembl_transcript_id)][, .(gene_id, feature_id)])

if (!(exists('mart'))) {
  mart <- biomaRt::useEnsembl(
    biomart = "genes",
    dataset = 'hsapiens_gene_ensembl',
    version = 109
  )
}
atts <- c("ensembl_gene_id", "ensembl_peptide_id", "pfam")
gene_interpro <- data.table(biomaRt::getBM(attributes = atts,
                                           mart = mart),
                            values = list(gene_biotype = "protein_coding"),
                            filters = c("gene_biotype"))

gene_interpro <- data.table(gene_interpro[pfam != "",])

gene_interpro_cont <- unique(gene_interpro[, .(ensembl_gene_id, pfam)])
slim_merge <- gene_slim[pf_slim, on = "feature_id"][!is.na(gene_id)]

dmi_set <- unique(gene_interpro_cont[slim_merge, 
                                     on = .(pfam = domain_id), 
                                     nomatch = 0, 
                                     allow.cartesian = TRUE][
                                       !is.na(gene_id) & !is.na(ensembl_gene_id)
                                     ][, 
                                       .(geneA = ensembl_gene_id, 
                                         geneB = gene_id, 
                                         DMI = TRUE,
                                         DMI_A = pfam,
                                         DMI_B = feature_id)])

ppidm <- get_ppidm(download = TRUE)
ddi <- ppidm[ppidm$CLASS %in% c("Gold", "Silver") | IN_GOLDSTANDARD == "yes", .(D1, D2)]

ddi_int <- unique(gene_interpro_cont[gene_interpro_cont[ddi, on = .(pfam = D1), allow.cartesian = TRUE, nomatch = 0], 
                                     on = .(pfam = D2), allow.cartesian = TRUE, nomatch = 0][ensembl_gene_id != i.ensembl_gene_id, 
                                                                                             .(geneA = ensembl_gene_id,
                                                                                               geneB = i.ensembl_gene_id,
                                                                                               DDI = TRUE,
                                                                                               DDI_A = pfam,
                                                                                               DDI_B = i.pfam)])


ddi_fin <- ddi_int[geneA %in% annotation_df$annotations$gene_id & geneB %in% annotation_df$annotations$gene_id]


get_elm_dmi_ppi_gene <- function(elm_interactions,
                                 species = c("human","mouse")[1],
                                 release = c(109,115)[1]) {
  
  species <- match.arg(species, c("human","mouse"))
  dataset <- if (species == "human") "hsapiens_gene_ensembl" else "mmusculus_gene_ensembl"
  
  options(biomaRt.cache = FALSE)
  if (!(exists('mart'))) {
    mart <- biomaRt::useEnsembl(
      biomart = "genes",
      dataset = species,
      version = release
    )
  }
  EI <- as.data.table(elm_interactions)
  unip <- unique(c(
    EI$interactorElm,
    EI$interactorDomain
  ))
  unip <- unip[!is.na(unip) & nzchar(unip)]
  
  map <- data.table(biomaRt::getBM(
    attributes = c("uniprotswissprot", "ensembl_gene_id"),
    mart = mart,
    filters = "uniprotswissprot",
    values = unip
  ))
  
  map <- map[
    uniprotswissprot != "" &
      !is.na(ensembl_gene_id) & nzchar(ensembl_gene_id),
    .(uniprotswissprot, ensembl_gene_id)
  ]
  map <- unique(map)
  setkey(map, uniprotswissprot)
  
  x <- map[EI, on = .(uniprotswissprot = interactorElm),
           allow.cartesian = TRUE, nomatch = 0L]
  setnames(x, "ensembl_gene_id", "geneElm")
  x[, interactorElm := uniprotswissprot][, uniprotswissprot := NULL]
  
  # map domain interactor
  x <- map[x, on = .(uniprotswissprot = interactorDomain),
           allow.cartesian = TRUE, nomatch = 0L]
  setnames(x, "ensembl_gene_id", "geneDomain")
  x[, interactorDomain := uniprotswissprot][, uniprotswissprot := NULL]
  
  x <- x[
    !is.na(geneElm) & !is.na(geneDomain) &
      geneElm != geneDomain
  ]
  
  x[, `:=`(geneA = geneElm, geneB = geneDomain)]
  
  # ------------------------------------------------------------
  # 5) final gene-centric PPI table
  # ------------------------------------------------------------
  ppi_gene <- unique(x[, .(
    geneA,
    geneB,
    DMI = TRUE,
    DMI_A = Elm,
    DMI_B = Domain
  )])
  
  ppi_gene[]
}


dmi_i <- get_elm_dmi_ppi_gene(elm_interactions, "human", 109)

dmi <- unique(rbind(dmi_i, dmi_set))
bg_pairs <- paste(
  pmin(biogrid$geneA, biogrid$geneB),
  pmax(biogrid$geneA, biogrid$geneB),
  sep = "||"
)

# canonicalize PPI pairs
ddi_pairs <- paste(
  pmin(ddi_fin$geneA, ddi_fin$geneB),
  pmax(ddi_fin$geneA, ddi_fin$geneB),
  sep = "||"
)

dmi_pairs <- paste(
  pmin(dmi$geneA, dmi$geneB),
  pmax(dmi$geneA, dmi$geneB),
  sep = "||"
)

library(VennDiagram)
# bg_pairs and ppi_pairs should be character vectors of canonicalized pairs, e.g. "ENSP..||ENSP.."
sets <- list(
  BioGRID = unique(bg_pairs),
  DDI     = unique(ddi_pairs),
  DMI = unique(dmi_pairs)
)


grid::grid.newpage()
venn.diagram(
  x = sets,
  filename = NULL,
  category.names = names(sets),
  main = "Overlap of interaction pairs"
) |> grid::grid.draw()




canon_pairs <- function(dt, a = "geneA", b = "geneB") {
  dt <- as.data.table(dt)
  stopifnot(a %in% names(dt), b %in% names(dt))
  dt[, `:=`(
    g1 = get(a),
    g2 = get(b)
  )]
  dt[, `:=`(
    geneA = pmin(g1, g2),
    geneB = pmax(g1, g2)
  )]
  # dt[, c("g1","g2") := NULL]
  dt[]
}

# --- 1) canonicalize all three ---
bg  <- canon_pairs(copy(biogrid))
ddi <- canon_pairs(copy(ddi_fin))
dmi2 <- canon_pairs(copy(dmi))

ddi[, `:=`(
  DDI_A_aligned = fifelse(geneA == g1, DDI_A, DDI_B),
  DDI_B_aligned = fifelse(geneB == g2, DDI_B, DDI_A)
)]
ddi[, c("DDI_A","DDI_B", "g1", "g2") := NULL]
setnames(ddi, c("DDI_A_aligned","DDI_B_aligned"), c("DDI_A","DDI_B"))

dmi2[, `:=`(
  DMI_A_aligned = fifelse(geneA == g1, DMI_A, DMI_B),
  DMI_B_aligned = fifelse(geneB == g2, DMI_B, DMI_A)
)]
dmi2[, c("DMI_A","DMI_B", "g1", "g2") := NULL]
setnames(dmi2, c("DMI_A_aligned","DMI_B_aligned"), c("DMI_A","DMI_B"))

bg_u <- unique(bg[, .(geneA = pmin(g1, g2), geneB = pmax(g1, g2), biogrid)])

dmi_w <- dmi2[, .(
  DMI = TRUE,
  DMI_A = list(unique(DMI_A)),
  DMI_B = list(unique(DMI_B))
), by = .(geneA, geneB)]

# aggregate DDI evidence per pair
ddi_w <- ddi[, .(
  DDI = TRUE,
  DDI_A = list(unique(DDI_A)),
  DDI_B = list(unique(DDI_B))
), by = .(geneA, geneB)]

# merge (safe now: one row per pair in each table)
setkey(bg_u, geneA, geneB)
setkey(dmi_w, geneA, geneB)
setkey(ddi_w, geneA, geneB)

wide_list <- merge(bg_u[, .(geneA, geneB, biogrid)], dmi_w, by=c("geneA","geneB"), all=TRUE)
wide_list <- merge(wide_list, ddi_w, by=c("geneA","geneB"), all=TRUE)

# fill booleans
wide_list[, biogrid := fifelse(is.na(biogrid), FALSE, biogrid)]
wide_list[, DMI     := fifelse(is.na(DMI),     FALSE, DMI)]
wide_list[, DDI     := fifelse(is.na(DDI),     FALSE, DDI)]

setcolorder(wide_list, c("geneA","geneB","biogrid","DMI","DMI_A","DMI_B","DDI","DDI_A","DDI_B"))

ppi <- wide_list[biogrid == TRUE]
saveRDS(ppi, '/inst/extdata/ppi.RDS')
