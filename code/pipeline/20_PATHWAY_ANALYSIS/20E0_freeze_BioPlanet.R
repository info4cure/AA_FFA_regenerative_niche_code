
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
})

cat("============================================================\n")
cat("PIPELINE_V2 — 20E0 FREEZE BIOPLANET\n")
cat("============================================================\n")

BASE <- file.path(
  "PIPELINE_V2","results","disease_biology",
  "20_PATHWAY_ANALYSIS","20E_BIOPLANET"
)

SRC <- file.path(
  BASE,"source",
  "BioPlanet_pathway_2019-05-14.csv"
)

OUT_GS <- file.path(BASE,"gene_sets")
OUT_TAB <- file.path(BASE,"tables")

dir.create(OUT_GS,recursive=TRUE,showWarnings=FALSE)
dir.create(OUT_TAB,recursive=TRUE,showWarnings=FALSE)

stopifnot(file.exists(SRC))

# ============================================================
# 1. Load immutable source
# ============================================================

cat("\nSTEP 1/7 — loading official BioPlanet source\n")

x <- fread(
  SRC,
  colClasses="character",
  na.strings=c("","NA")
)

expected <- c(
  "PATHWAY_ID",
  "PATHWAY_NAME",
  "GENE_ID",
  "GENE_SYMBOL"
)

stopifnot(
  identical(names(x),expected)
)

cat("Raw rows:",nrow(x),"\n")
cat("Columns :",ncol(x),"\n")

# ============================================================
# 2. Raw QC
# ============================================================

cat("\nSTEP 2/7 — raw QC\n")

qc_raw <- tibble(
  metric=c(
    "raw_rows",
    "unique_pathway_ids",
    "unique_pathway_names",
    "unique_gene_ids",
    "unique_gene_symbols",
    "missing_pathway_id",
    "missing_pathway_name",
    "missing_gene_id",
    "missing_gene_symbol",
    "exact_duplicate_rows"
  ),
  value=c(
    nrow(x),
    n_distinct(x$PATHWAY_ID),
    n_distinct(x$PATHWAY_NAME),
    n_distinct(x$GENE_ID,na.rm=TRUE),
    n_distinct(x$GENE_SYMBOL,na.rm=TRUE),
    sum(is.na(x$PATHWAY_ID)),
    sum(is.na(x$PATHWAY_NAME)),
    sum(is.na(x$GENE_ID)),
    sum(is.na(x$GENE_SYMBOL)),
    nrow(x)-nrow(distinct(x))
  )
)

print(as.data.frame(qc_raw),row.names=FALSE)

write_csv(
  qc_raw,
  file.path(OUT_TAB,"20E0_raw_QC.csv")
)

# ============================================================
# 3. Pathway identity integrity
# ============================================================

cat("\nSTEP 3/7 — pathway identity integrity\n")

id_to_name <- x %>%
  distinct(PATHWAY_ID,PATHWAY_NAME) %>%
  count(PATHWAY_ID,name="n_names") %>%
  filter(n_names>1)

name_to_id <- x %>%
  distinct(PATHWAY_ID,PATHWAY_NAME) %>%
  count(PATHWAY_NAME,name="n_ids") %>%
  filter(n_ids>1)

cat(
  "Pathway IDs with >1 name:",
  nrow(id_to_name),"\n"
)

cat(
  "Pathway names with >1 ID:",
  nrow(name_to_id),"\n"
)

write_csv(
  id_to_name,
  file.path(
    OUT_TAB,
    "20E0_pathway_ID_name_conflicts.csv"
  )
)

write_csv(
  name_to_id,
  file.path(
    OUT_TAB,
    "20E0_pathway_name_ID_conflicts.csv"
  )
)

# ============================================================
# 4. Freeze clean pathway-gene universe
# ============================================================

cat("\nSTEP 4/7 — freezing pathway-gene universe\n")

clean <- x %>%
  filter(
    !is.na(PATHWAY_ID),
    !is.na(PATHWAY_NAME),
    !is.na(GENE_SYMBOL),
    PATHWAY_ID!="",
    PATHWAY_NAME!="",
    GENE_SYMBOL!=""
  ) %>%
  mutate(
    GENE_SYMBOL=str_trim(GENE_SYMBOL)
  ) %>%
  distinct(
    PATHWAY_ID,
    PATHWAY_NAME,
    GENE_ID,
    GENE_SYMBOL
  ) %>%
  arrange(
    PATHWAY_ID,
    GENE_SYMBOL,
    GENE_ID
  )

write_csv(
  clean,
  file.path(
    OUT_TAB,
    "20E0_BioPlanet_pathway_gene_universe.csv"
  )
)

# fgsea needs unique gene symbols within each pathway.
gene_sets_tbl <- clean %>%
  distinct(
    PATHWAY_ID,
    PATHWAY_NAME,
    GENE_SYMBOL
  )

# Preserve BioPlanet ID in gene-set name so pathway identity
# remains unambiguous and reversible.
gene_sets_tbl <- gene_sets_tbl %>%
  mutate(
    GENESET_NAME=paste0(
      "BIOPLANET_",
      PATHWAY_ID,
      "__",
      PATHWAY_NAME
    )
  )

gene_sets <- split(
  gene_sets_tbl$GENE_SYMBOL,
  gene_sets_tbl$GENESET_NAME
)

gene_sets <- lapply(
  gene_sets,
  function(z) sort(unique(z))
)

saveRDS(
  gene_sets,
  file.path(
    OUT_GS,
    "20E0_BioPlanet_gene_sets.rds"
  )
)

# ============================================================
# 5. Inventory
# ============================================================

cat("\nSTEP 5/7 — pathway inventory\n")

inventory <- gene_sets_tbl %>%
  group_by(
    PATHWAY_ID,
    PATHWAY_NAME,
    GENESET_NAME
  ) %>%
  summarise(
    n_genes=n_distinct(GENE_SYMBOL),
    .groups="drop"
  ) %>%
  arrange(
    PATHWAY_ID
  )

write_csv(
  inventory,
  file.path(
    OUT_TAB,
    "20E0_BioPlanet_pathway_inventory.csv"
  )
)

cat("Frozen pathways:",nrow(inventory),"\n")
cat(
  "Unique symbols:",
  n_distinct(gene_sets_tbl$GENE_SYMBOL),
  "\n"
)

cat("\nPathway size summary:\n")
print(summary(inventory$n_genes))

cat("\nPathway-size thresholds:\n")

size_qc <- tibble(
  threshold=c(
    "lt10",
    "10_to_500",
    "gt500"
  ),
  n_pathways=c(
    sum(inventory$n_genes<10),
    sum(
      inventory$n_genes>=10 &
      inventory$n_genes<=500
    ),
    sum(inventory$n_genes>500)
  )
)

print(as.data.frame(size_qc),row.names=FALSE)

write_csv(
  size_qc,
  file.path(
    OUT_TAB,
    "20E0_pathway_size_QC.csv"
  )
)

# ============================================================
# 6. Source checksum
# ============================================================

cat("\nSTEP 6/7 — source checksum\n")

sha <- unname(
  tools::md5sum(SRC)
)

# SHA256 generated independently by shasum is preferred.
sha_file <- file.path(
  BASE,"source",
  "BioPlanet_pathway_2019-05-14.sha256"
)

sha256 <- if (file.exists(sha_file)) {
  readLines(sha_file,warn=FALSE)[1]
} else {
  NA_character_
}

source_manifest <- tibble(
  source_file=SRC,
  source_bytes=file.info(SRC)$size,
  source_date_label="2019-05-14",
  source_provider="NCATS BioPlanet",
  source_format="PATHWAY_ID/PATHWAY_NAME/GENE_ID/GENE_SYMBOL long table",
  md5=sha,
  sha256_record=sha256,
  frozen_pathways=nrow(inventory),
  frozen_unique_gene_symbols=
    n_distinct(gene_sets_tbl$GENE_SYMBOL)
)

write_csv(
  source_manifest,
  file.path(
    OUT_TAB,
    "20E0_source_manifest.csv"
  )
)

print(as.data.frame(source_manifest),row.names=FALSE)

# ============================================================
# 7. Formal checks + audit
# ============================================================

cat("\nSTEP 7/7 — formal checks\n")

checks <- tibble(
  check=c(
    "source_present",
    "four_expected_columns",
    "pathways_nonzero",
    "genes_nonzero",
    "gene_sets_equal_inventory",
    "no_empty_gene_sets",
    "unique_genes_within_sets",
    "source_manifest_generated",
    "no_GSEA_performed",
    "no_DE_refit",
    "annotation_unchanged"
  ),
  pass=c(
    file.exists(SRC),
    identical(names(x),expected),
    nrow(inventory)>0,
    n_distinct(gene_sets_tbl$GENE_SYMBOL)>0,
    length(gene_sets)==nrow(inventory),
    all(lengths(gene_sets)>0),
    all(
      vapply(
        gene_sets,
        function(z)
          length(z)==length(unique(z)),
        logical(1)
      )
    ),
    file.exists(
      file.path(
        OUT_TAB,
        "20E0_source_manifest.csv"
      )
    ),
    TRUE,
    TRUE,
    TRUE
  )
)

print(as.data.frame(checks),row.names=FALSE)

stopifnot(all(checks$pass))

audit <- c(
  "PIPELINE_V2 — 20E0 BioPlanet acquisition/freeze",
  "",
  paste0("Freeze date: ",Sys.Date()),
  "",
  "Source:",
  "- NCATS BioPlanet official pathway.csv",
  "- local source label: 2019-05-14",
  "- raw source retained unchanged",
  "- SHA-256 recorded separately at acquisition",
  "",
  "Raw schema:",
  "- PATHWAY_ID",
  "- PATHWAY_NAME",
  "- GENE_ID",
  "- GENE_SYMBOL",
  "",
  paste0(
    "Frozen pathways: ",
    nrow(inventory)
  ),
  paste0(
    "Unique gene symbols: ",
    n_distinct(gene_sets_tbl$GENE_SYMBOL)
  ),
  "",
  "Processing:",
  "- incomplete pathway/gene-symbol records removed",
  "- exact records deduplicated",
  "- unique gene symbols retained within each pathway",
  "- BioPlanet pathway ID preserved in GSEA gene-set name",
  "",
  "Planned downstream GSEA eligibility:",
  "- minSize = 10",
  "- maxSize = 500",
  "",
  "No enrichment performed.",
  "No differential expression refit.",
  "No cell annotation modified."
)

writeLines(
  audit,
  file.path(
    BASE,
    "20E0_audit.txt"
  )
)

cat("\n============================================================\n")
cat("20E0 COMPLETE\n")
cat("============================================================\n")
cat("Official BioPlanet collection frozen\n")
cat("No GSEA performed\n")
cat("PASS\n")
