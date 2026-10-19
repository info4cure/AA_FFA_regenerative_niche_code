# ==============================================================================
# 22F4E4C — HYBRID HFSC–INFLAMMATORY STATE
#
# PRIMARY QUESTION
# Does FFA deviate toward concurrent HFSC-identity and inflammatory-transition
# activation along the frozen Bulge/Isthmus trajectory?
#
# INPUT
# Frozen 22F4E4B external-signature scores.
#
# IMPORTANT
# - no new genes
# - no new scoring
# - no reclustering
# - no pseudotime recalculation
# - pseudobulk = patient x frozen trajectory bin
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(svglite)
})

MODULE <- "22F4E4C"

ROOT <- "PIPELINE_V2/results/disease_biology/22_MULTICOMPARTMENT_DISEASE_AXES/22F_BULGE_STEM_NICHE_INTEGRITY/22F4_BULGE_MONOCLE3_TRAJECTORIES/22F4E_DISEASE_TRAJECTORY/22F4E4_HALLMARK_TRAJECTORY"

INPUT <- file.path(
  ROOT,
  "22F4E4B_EXTERNAL_HFSC_SIGNATURE",
  "objects",
  "22F4E4B_external_HFSC_signature_freeze.rds"
)

OUT <- file.path(
  ROOT,
  "22F4E4C_HYBRID_HFSC_INFLAMMATORY_STATE"
)

TABLE_DIR <- file.path(OUT,"tables")
OBJECT_DIR <- file.path(OUT,"objects")
SVG_DIR <- file.path(OUT,"figures","svg")
PDF_DIR <- file.path(OUT,"figures","pdf")
PNG_DIR <- file.path(OUT,"figures","png")

for (d in c(
  TABLE_DIR,
  OBJECT_DIR,
  SVG_DIR,
  PDF_DIR,
  PNG_DIR
)) {
  dir.create(
    d,
    recursive=TRUE,
    showWarnings=FALSE
  )
}

stopifnot(file.exists(INPUT))

x <- readRDS(INPUT)

cat("\n============================================================\n")
cat("22F4E4C — HYBRID HFSC–INFLAMMATORY STATE\n")
cat("============================================================\n")

# ------------------------------------------------------------------------------
# 1. RECOVER FROZEN SCORES
# ------------------------------------------------------------------------------

scores <- x$signature_scores
meta <- x$metadata

stopifnot(
  ncol(scores)==56,
  nrow(meta)==56,
  identical(
    colnames(scores),
    meta$pseudobulk_id
  ),
  all(
    c(
      "STEM_HFSC",
      "INFLAMMATORY_TRANSITION"
    ) %in% rownames(scores)
  )
)

dat <- tibble(
  pseudobulk_id=colnames(scores),
  HFSC_identity=
    as.numeric(
      scores[
        "STEM_HFSC",
      ]
    ),
  inflammatory_transition=
    as.numeric(
      scores[
        "INFLAMMATORY_TRANSITION",
      ]
    )
) |>
  left_join(
    meta,
    by="pseudobulk_id"
  ) |>
  mutate(
    disease=factor(
      disease,
      levels=c(
        "HC",
        "FFA",
        "AA"
      )
    ),
    trajectory_bin5=
      as.integer(
        trajectory_bin5
      )
  )

stopifnot(
  nrow(dat)==56,
  all(
    is.finite(
      dat$HFSC_identity
    )
  ),
  all(
    is.finite(
      dat$inflammatory_transition
    )
  )
)

# ------------------------------------------------------------------------------
# 2. DISEASE × BIN CENTROIDS
# ------------------------------------------------------------------------------

centroids <- dat |>
  group_by(
    disease,
    trajectory_bin5
  ) |>
  summarise(
    n_pseudobulks=n(),
    n_patients=n_distinct(
      patient_id
    ),
    mean_HFSC=
      mean(
        HFSC_identity
      ),
    mean_inflammatory=
      mean(
        inflammatory_transition
      ),
    median_HFSC=
      median(
        HFSC_identity
      ),
    median_inflammatory=
      median(
        inflammatory_transition
      ),
    .groups="drop"
  ) |>
  arrange(
    disease,
    trajectory_bin5
  )

write.csv(
  centroids,
  file.path(
    TABLE_DIR,
    "22F4E4C_disease_bin_centroids.csv"
  ),
  row.names=FALSE
)

# ------------------------------------------------------------------------------
# 3. WITHIN-DISEASE ASSOCIATION
#
# Descriptive only. Pseudobulks from the same patient are not independent,
# therefore these correlations are NOT used as primary inferential tests.
# ------------------------------------------------------------------------------

correlations <- dat |>
  group_by(disease) |>
  summarise(
    n_pseudobulks=n(),
    n_patients=n_distinct(
      patient_id
    ),
    spearman_rho=
      cor(
        HFSC_identity,
        inflammatory_transition,
        method="spearman"
      ),
    pearson_r=
      cor(
        HFSC_identity,
        inflammatory_transition,
        method="pearson"
      ),
    .groups="drop"
  )

write.csv(
  correlations,
  file.path(
    TABLE_DIR,
    "22F4E4C_within_disease_correlations_descriptive.csv"
  ),
  row.names=FALSE
)

# ------------------------------------------------------------------------------
# 4. PATIENT-LEVEL COACTIVATION SUMMARY
#
# A non-arbitrary continuous measure:
# mean product of the two already-standardized frozen scores.
#
# This is SECONDARY / descriptive. It is not used to define a new cell state.
# ------------------------------------------------------------------------------

dat <- dat |>
  mutate(
    coactivation_product=
      HFSC_identity *
      inflammatory_transition
  )

patient_summary <- dat |>
  group_by(
    disease,
    patient_id
  ) |>
  summarise(
    n_bins=n(),
    mean_HFSC=
      mean(
        HFSC_identity
      ),
    mean_inflammatory=
      mean(
        inflammatory_transition
      ),
    mean_coactivation=
      mean(
        coactivation_product
      ),
    max_coactivation=
      max(
        coactivation_product
      ),
    .groups="drop"
  )

write.csv(
  patient_summary,
  file.path(
    TABLE_DIR,
    "22F4E4C_patient_coactivation_summary.csv"
  ),
  row.names=FALSE
)

# ------------------------------------------------------------------------------
# 5. FIGURE — FUNCTIONAL-STATE PLANE
# ------------------------------------------------------------------------------

disease_cols <- c(
  HC="#8A8A8A",
  FFA="#8FD694",
  AA="#D94B4B"
)

p <- ggplot(
  dat,
  aes(
    x=HFSC_identity,
    y=inflammatory_transition,
    colour=disease
  )
) +

  # Individual patient-bin observations
  geom_point(
    aes(
      size=n_cells
    ),
    alpha=0.32,
    stroke=0,
    show.legend=FALSE
  ) +

  # Disease-specific mean trajectory through B5_1 ... B5_5
  geom_path(
    data=centroids,
    aes(
      x=mean_HFSC,
      y=mean_inflammatory,
      colour=disease,
      group=disease
    ),
    linewidth=1.35,
    lineend="round",
    arrow=grid::arrow(
      length=grid::unit(
        0.11,
        "inches"
      ),
      type="closed"
    ),
    inherit.aes=FALSE
  ) +

  geom_point(
    data=centroids,
    aes(
      x=mean_HFSC,
      y=mean_inflammatory,
      colour=disease
    ),
    size=3.0,
    shape=21,
    fill="white",
    stroke=1.1,
    inherit.aes=FALSE
  ) +

  geom_text(
    data=centroids,
    aes(
      x=mean_HFSC,
      y=mean_inflammatory,
      label=paste0(
        "B",
        trajectory_bin5
      ),
      colour=disease
    ),
    size=3.0,
    fontface="bold",
    nudge_y=0.08,
    show.legend=FALSE,
    inherit.aes=FALSE
  ) +

  geom_vline(
    xintercept=0,
    linewidth=0.35,
    colour="grey75"
  ) +

  geom_hline(
    yintercept=0,
    linewidth=0.35,
    colour="grey75"
  ) +

  scale_colour_manual(
    values=disease_cols,
    breaks=c(
      "HC",
      "FFA",
      "AA"
    ),
    drop=FALSE
  ) +

  scale_size_continuous(
    range=c(
      1.8,
      5.0
    )
  ) +

  coord_equal() +

  labs(
    title=
      "FFA diverges into an inflammation-associated HFSC state",
    subtitle=
      "Frozen external HFSC-identity and inflammatory-transition scores; arrows trace mean B5_1 → B5_5 progression",
    x=
      "HFSC identity\nKRT15 · SOX9 · LEF1",
    y=
      "Inflammatory transition\nSTAT1 · IRF1 · IFITM3",
    colour=NULL
  ) +

  theme_classic(
    base_size=10
  ) +

  theme(
    plot.title=element_text(
      face="bold",
      size=12
    ),
    plot.subtitle=element_text(
      size=9,
      colour="grey30",
      margin=margin(
        b=8
      )
    ),
    axis.title=element_text(
      size=10
    ),
    axis.text=element_text(
      size=9,
      colour="black"
    ),
    legend.position="top",
    legend.justification="left",
    legend.text=element_text(
      size=9
    ),
    plot.margin=margin(
      8,
      10,
      8,
      8
    )
  )

# ------------------------------------------------------------------------------
# 6. SAVE
# ------------------------------------------------------------------------------

stem <-
  "Fig22F4E4C_HFSC_inflammatory_state_plane_preflight"

ggsave(
  file.path(
    SVG_DIR,
    paste0(
      stem,
      ".svg"
    )
  ),
  p,
  width=6.8,
  height=6.2,
  units="in",
  device=svglite::svglite
)

ggsave(
  file.path(
    PDF_DIR,
    paste0(
      stem,
      ".pdf"
    )
  ),
  p,
  width=6.8,
  height=6.2,
  units="in"
)

ggsave(
  file.path(
    PNG_DIR,
    paste0(
      stem,
      ".png"
    )
  ),
  p,
  width=6.8,
  height=6.2,
  units="in",
  dpi=400
)

# ------------------------------------------------------------------------------
# 7. FREEZE
# ------------------------------------------------------------------------------

freeze <- list(
  module=MODULE,
  analysis=
    "Joint HFSC identity × inflammatory-transition functional-state plane",
  source=
    "Frozen 22F4E4B external signature scores",
  new_gene_scoring=FALSE,
  new_pseudotime=FALSE,
  n_pseudobulks=nrow(dat),
  n_patients=n_distinct(
    dat$patient_id
  ),
  data=dat,
  centroids=centroids,
  correlations=correlations,
  patient_summary=patient_summary,
  output_stem=stem
)

saveRDS(
  freeze,
  file.path(
    OBJECT_DIR,
    "22F4E4C_hybrid_HFSC_inflammatory_state_freeze.rds"
  )
)

# ------------------------------------------------------------------------------
# 8. FINAL REPORT
# ------------------------------------------------------------------------------

cat("\n===== DISEASE × BIN CENTROIDS =====\n")
print(
  centroids,
  n=Inf
)

cat("\n===== WITHIN-DISEASE CORRELATIONS — DESCRIPTIVE =====\n")
print(
  correlations,
  n=Inf
)

cat("\n===== PATIENT COACTIVATION SUMMARY =====\n")
print(
  patient_summary,
  n=Inf
)

cat(
  "\nPrimary inference remains the frozen 22F4E4B disease-by-trajectory tests.\n"
)

cat(
  "22F4E4C is a joint-state visualization and secondary coactivation analysis.\n"
)

cat(
  "\nSTATUS=PASS\n"
)
