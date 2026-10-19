# ==============================================================================
# 22F4E4B — EXTERNAL LPP/HFSC SIGNATURE PROJECTION
#
# PURPOSE
#   Test an a priori external HFSC programme derived from Chen et al.
#   Cell & Bioscience 2025 on the frozen Bulge/Isthmus disease trajectory.
#
# PRINCIPLES
#   - NO reclustering
#   - NO pseudotime recalculation
#   - NO new gene-level DE
#   - frozen 22F4E2 pseudobulks/design
#   - same score/model framework as frozen 22F4E4A Hallmark analysis
#   - patient remains biological replicate
#
# EXTERNAL SIGNATURES — frozen before inspection
#   STEM_HFSC:
#       KRT15, SOX9, LEF1
#
#   INFLAMMATORY_TRANSITION:
#       STAT1, IRF1, IFITM3
#
#   REMODELING_ASSOCIATED:
#       VIM, ZEB1
#
# NOTE
#   ZEB1 is absent from the frozen 22F4E2 tested gene universe.
#   It will NOT be substituted. Therefore the third score is VIM-only and
#   must be interpreted as supportive rather than as a full EMT signature.
# ==============================================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(svglite)
})

options(stringsAsFactors = FALSE)

MODULE <- "22F4E4B"

# ------------------------------------------------------------------------------
# 1. PATHS
# ------------------------------------------------------------------------------

ROOT <- "PIPELINE_V2/results/disease_biology/22_MULTICOMPARTMENT_DISEASE_AXES/22F_BULGE_STEM_NICHE_INTEGRITY/22F4_BULGE_MONOCLE3_TRAJECTORIES/22F4E_DISEASE_TRAJECTORY"

UPSTREAM <- file.path(
  ROOT,
  "22F4E2_INTERACTION_MODELLING",
  "objects",
  "22F4E2_disease_trajectory_interaction_freeze.rds"
)

HALLMARK_FREEZE <- file.path(
  ROOT,
  "22F4E4_HALLMARK_TRAJECTORY",
  "22F4E4A_SCORE_INTERACTION",
  "objects",
  "22F4E4A_Hallmark_trajectory_interaction_freeze.rds"
)

OUT <- file.path(
  ROOT,
  "22F4E4_HALLMARK_TRAJECTORY",
  "22F4E4B_EXTERNAL_HFSC_SIGNATURE"
)

TABLE_DIR <- file.path(OUT, "tables")
OBJECT_DIR <- file.path(OUT, "objects")
SVG_DIR <- file.path(OUT, "figures", "svg")
PDF_DIR <- file.path(OUT, "figures", "pdf")
PNG_DIR <- file.path(OUT, "figures", "png")

dir.create(TABLE_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(OBJECT_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(SVG_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(PDF_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(PNG_DIR, recursive=TRUE, showWarnings=FALSE)

stopifnot(file.exists(UPSTREAM))
stopifnot(file.exists(HALLMARK_FREEZE))

x <- readRDS(UPSTREAM)
h <- readRDS(HALLMARK_FREEZE)

cat("\n============================================================\n")
cat("22F4E4B — EXTERNAL LPP/HFSC SIGNATURE PROJECTION\n")
cat("============================================================\n")

# ------------------------------------------------------------------------------
# 2. FROZEN EXTERNAL SIGNATURES
# ------------------------------------------------------------------------------

signatures_requested <- list(
  STEM_HFSC = c(
    "KRT15",
    "SOX9",
    "LEF1"
  ),
  INFLAMMATORY_TRANSITION = c(
    "STAT1",
    "IRF1",
    "IFITM3"
  ),
  REMODELING_ASSOCIATED = c(
    "VIM",
    "ZEB1"
  )
)

display_labels <- c(
  STEM_HFSC = "HFSC identity",
  INFLAMMATORY_TRANSITION = "Inflammatory transition",
  REMODELING_ASSOCIATED = "Remodeling-associated VIM"
)

# ------------------------------------------------------------------------------
# 3. INPUT AUDIT
# ------------------------------------------------------------------------------

dge <- x$dge
metadata <- x$metadata
design <- x$design

stopifnot(
  ncol(dge$counts) == nrow(metadata),
  nrow(design) == nrow(metadata),
  identical(colnames(dge$counts), metadata$pseudobulk_id),
  nrow(metadata) == 56,
  length(unique(metadata$patient_id)) == 16
)

genes_tested <- rownames(dge$counts)

gene_audit <- bind_rows(
  lapply(names(signatures_requested), function(sig) {

    req <- signatures_requested[[sig]]

    tibble(
      signature = sig,
      gene = req,
      requested = TRUE,
      present = req %in% genes_tested
    )
  })
)

write.csv(
  gene_audit,
  file.path(
    TABLE_DIR,
    "22F4E4B_external_signature_gene_audit.csv"
  ),
  row.names=FALSE
)

cat("\n===== EXTERNAL SIGNATURE AUDIT =====\n")
print(gene_audit)

# ------------------------------------------------------------------------------
# 4. RECONSTRUCT 22F4E4A SCORE DEFINITION
# ----------------------------------------------------------------
# ------------------------------------------------------------------------------

logCPM <- edgeR::cpm(
  dge,
  log=TRUE,
  prior.count=2,
  normalized.lib.sizes=TRUE
)

stopifnot(all(is.finite(logCPM)))

# Standardize each gene across the same 56 frozen pseudobulks.
z_logCPM <- t(
  scale(
    t(logCPM),
    center=TRUE,
    scale=TRUE
  )
)

stopifnot(all(is.finite(z_logCPM)))

# ------------------------------------------------------------------------------
# 5. FRAMEWORK RECOVERY AUDIT AGAINST FROZEN 22F4E4A
# ------------------------------------------------------------------------------

recovery_audit <- tibble(
  check = c(
    "same_n_pseudobulks",
    "same_design",
    "same_pseudobulk_order",
    "same_trajectory_grid",
    "finite_standardized_logCPM"
  ),
  passed = c(
    ncol(logCPM) == h$n_pseudobulks,
    isTRUE(all.equal(design, h$design, tolerance=0)),
    identical(colnames(logCPM), colnames(h$Hallmark_scores)),
    isTRUE(all.equal(x$trajectory_grid, h$trajectory_grid, tolerance=0)),
    all(is.finite(z_logCPM))
  )
)

write.csv(
  recovery_audit,
  file.path(
    TABLE_DIR,
    "22F4E4B_framework_recovery_audit.csv"
  ),
  row.names=FALSE
)

cat("\n===== FRAMEWORK RECOVERY AUDIT =====\n")
print(recovery_audit)

if (!all(recovery_audit$passed)) {
  stop("Framework recovery audit failed.")
}

# ------------------------------------------------------------------------------
# 6. CALCULATE EXTERNAL SIGNATURE SCORES
# ------------------------------------------------------------------------------

signatures_present <- lapply(
  signatures_requested,
  function(g) intersect(g, rownames(z_logCPM))
)

signature_scores <- sapply(
  signatures_present,
  function(g) {

    if (length(g) == 0) {
      return(rep(NA_real_, ncol(z_logCPM)))
    }

    colMeans(
      z_logCPM[g, , drop=FALSE],
      na.rm=FALSE
    )
  }
)

signature_scores <- t(signature_scores)

colnames(signature_scores) <- colnames(z_logCPM)

stopifnot(
  nrow(signature_scores) == 3,
  ncol(signature_scores) == 56,
  all(is.finite(signature_scores))
)

score_long <- as.data.frame(t(signature_scores)) |>
  rownames_to_column("pseudobulk_id") |>
  left_join(
    metadata,
    by="pseudobulk_id"
  ) |>
  pivot_longer(
    cols=all_of(names(signatures_requested)),
    names_to="signature",
    values_to="score"
  ) |>
  mutate(
    signature_label = unname(display_labels[signature])
  )

write.csv(
  score_long,
  file.path(
    TABLE_DIR,
    "22F4E4B_external_signature_scores.csv"
  ),
  row.names=FALSE
)

# ------------------------------------------------------------------------------
# 7. FIT SAME LIMMA DESIGN AS 22F4E4A
# ------------------------------------------------------------------------------

fit <- limma::lmFit(
  signature_scores,
  design
)

fit <- limma::eBayes(
  fit,
  robust=TRUE
)

common_coef <- x$common_coef
FFA_coef <- x$FFA_coef
AA_coef <- x$AA_coef
joint_coef <- x$joint_coef

# Helper for multi-coefficient tests.
run_multicoef_test <- function(fit_object, coef_index) {

  tt <- limma::topTable(
    fit_object,
    coef=coef_index,
    number=Inf,
    sort.by="none"
  )

  tibble(
    feature = rownames(tt),
    F = if ("F" %in% colnames(tt)) tt$F else NA_real_,
    PValue = if ("P.Value" %in% colnames(tt)) tt$P.Value else NA_real_
  )
}

global_test <- run_multicoef_test(
  fit,
  joint_coef
)

ffa_test <- run_multicoef_test(
  fit,
  FFA_coef
)

aa_test <- run_multicoef_test(
  fit,
  AA_coef
)

# AA vs FFA trajectory interaction:
# (AA_x_s1 - FFA_x_s1, AA_x_s2 - FFA_x_s2)
contrast_matrix <- matrix(
  0,
  nrow=ncol(design),
  ncol=2,
  dimnames=list(
    colnames(design),
    c(
      "AA_minus_FFA_s1",
      "AA_minus_FFA_s2"
    )
  )
)

contrast_matrix[
  "AA_x_s1",
  "AA_minus_FFA_s1"
] <- 1

contrast_matrix[
  "FFA_x_s1",
  "AA_minus_FFA_s1"
] <- -1

contrast_matrix[
  "AA_x_s2",
  "AA_minus_FFA_s2"
] <- 1

contrast_matrix[
  "FFA_x_s2",
  "AA_minus_FFA_s2"
] <- -1

fit_AA_FFA <- limma::contrasts.fit(
  fit,
  contrast_matrix
)

fit_AA_FFA <- limma::eBayes(
  fit_AA_FFA,
  robust=TRUE
)

aa_ffa_test <- run_multicoef_test(
  fit_AA_FFA,
  1:2
)

interaction_results <- tibble(
  signature = rownames(signature_scores),

  global_F = global_test$F[
    match(
      rownames(signature_scores),
      global_test$feature
    )
  ],
  global_PValue = global_test$PValue[
    match(
      rownames(signature_scores),
      global_test$feature
    )
  ],

  FFA_vs_HC_F = ffa_test$F[
    match(
      rownames(signature_scores),
      ffa_test$feature
    )
  ],
  FFA_vs_HC_PValue = ffa_test$PValue[
    match(
      rownames(signature_scores),
      ffa_test$feature
    )
  ],

  AA_vs_HC_F = aa_test$F[
    match(
      rownames(signature_scores),
      aa_test$feature
    )
  ],
  AA_vs_HC_PValue = aa_test$PValue[
    match(
      rownames(signature_scores),
      aa_test$feature
    )
  ],

  AA_vs_FFA_F = aa_ffa_test$F[
    match(
      rownames(signature_scores),
      aa_ffa_test$feature
    )
  ],
  AA_vs_FFA_PValue = aa_ffa_test$PValue[
    match(
      rownames(signature_scores),
      aa_ffa_test$feature
    )
  ]
) |>
  mutate(
    global_FDR = p.adjust(
      global_PValue,
      method="BH"
    ),
    FFA_vs_HC_FDR = p.adjust(
      FFA_vs_HC_PValue,
      method="BH"
    ),
    AA_vs_HC_FDR = p.adjust(
      AA_vs_HC_PValue,
      method="BH"
    ),
    AA_vs_FFA_FDR = p.adjust(
      AA_vs_FFA_PValue,
      method="BH"
    )
  )

write.csv(
  interaction_results,
  file.path(
    TABLE_DIR,
    "22F4E4B_external_signature_interactions.csv"
  ),
  row.names=FALSE
)

# ------------------------------------------------------------------------------
# 8. RECOVER FROZEN SPLINE BASIS
# ------------------------------------------------------------------------------

trajectory_grid <- h$trajectory_grid
B <- h$recovered_spline_basis

stopifnot(
  nrow(B) == length(common_coef),
  ncol(B) == length(trajectory_grid)
)

coef_mat <- fit$coefficients

HC_curve <- (
  coef_mat[
    ,
    common_coef,
    drop=FALSE
  ]
) %*% B

FFA_curve <- (
  coef_mat[
    ,
    common_coef,
    drop=FALSE
  ] +
  coef_mat[
    ,
    FFA_coef,
    drop=FALSE
  ]
) %*% B

AA_curve <- (
  coef_mat[
    ,
    common_coef,
    drop=FALSE
  ] +
  coef_mat[
    ,
    AA_coef,
    drop=FALSE
  ]
) %*% B

rownames(HC_curve) <- rownames(signature_scores)
rownames(FFA_curve) <- rownames(signature_scores)
rownames(AA_curve) <- rownames(signature_scores)

curve_names <- sprintf(
  "T%03d",
  seq_along(trajectory_grid)
)

colnames(HC_curve) <- curve_names
colnames(FFA_curve) <- curve_names
colnames(AA_curve) <- curve_names

FFA_minus_HC <- FFA_curve - HC_curve
AA_minus_HC <- AA_curve - HC_curve
AA_minus_FFA <- AA_curve - FFA_curve

# ------------------------------------------------------------------------------
# 9. CURVE TABLES
# ------------------------------------------------------------------------------

matrix_to_long <- function(mat, disease_name) {

  as.data.frame(mat) |>
    rownames_to_column("signature") |>
    pivot_longer(
      cols=-signature,
      names_to="grid_id",
      values_to="fitted_score"
    ) |>
    mutate(
      grid_index=as.integer(
        sub(
          "^T",
          "",
          grid_id
        )
      ),
      trajectory_position=
        trajectory_grid[grid_index],
      disease=disease_name,
      signature_label=
        unname(
          display_labels[signature]
        )
    ) |>
    select(
      signature,
      signature_label,
      disease,
      trajectory_position,
      fitted_score
    )
}

curve_long <- bind_rows(
  matrix_to_long(
    HC_curve,
    "HC"
  ),
  matrix_to_long(
    FFA_curve,
    "FFA"
  ),
  matrix_to_long(
    AA_curve,
    "AA"
  )
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
    signature_label=factor(
      signature_label,
      levels=c(
        "HFSC identity",
        "Inflammatory transition",
        "Remodeling-associated VIM"
      )
    )
  )

write.csv(
  curve_long,
  file.path(
    TABLE_DIR,
    "22F4E4B_external_signature_fitted_curves.csv"
  ),
  row.names=FALSE
)

difference_to_long <- function(
  mat,
  contrast_name
) {

  as.data.frame(mat) |>
    rownames_to_column("signature") |>
    pivot_longer(
      cols=-signature,
      names_to="grid_id",
      values_to="difference"
    ) |>
    mutate(
      grid_index=as.integer(
        sub(
          "^T",
          "",
          grid_id
        )
      ),
      trajectory_position=
        trajectory_grid[grid_index],
      contrast=contrast_name
    ) |>
    select(
      signature,
      contrast,
      trajectory_position,
      difference
    )
}

difference_long <- bind_rows(
  difference_to_long(
    FFA_minus_HC,
    "FFA_minus_HC"
  ),
  difference_to_long(
    AA_minus_HC,
    "AA_minus_HC"
  ),
  difference_to_long(
    AA_minus_FFA,
    "AA_minus_FFA"
  )
)

write.csv(
  difference_long,
  file.path(
    TABLE_DIR,
    "22F4E4B_external_signature_difference_curves.csv"
  ),
  row.names=FALSE
)

# ------------------------------------------------------------------------------
# 10. CURVE METRICS
# ------------------------------------------------------------------------------

curve_metrics <- tibble(
  signature=rownames(
    signature_scores
  ),

  max_abs_FFA_minus_HC=
    apply(
      abs(FFA_minus_HC),
      1,
      max
    ),

  max_abs_AA_minus_HC=
    apply(
      abs(AA_minus_HC),
      1,
      max
    ),

  max_abs_AA_minus_FFA=
    apply(
      abs(AA_minus_FFA),
      1,
      max
    ),

  peak_FFA_position=
    trajectory_grid[
      apply(
        abs(FFA_minus_HC),
        1,
        which.max
      )
    ],

  peak_AA_position=
    trajectory_grid[
      apply(
        abs(AA_minus_HC),
        1,
        which.max
      )
    ],

  peak_AA_FFA_position=
    trajectory_grid[
      apply(
        abs(AA_minus_FFA),
        1,
        which.max
      )
    ]
)

write.csv(
  curve_metrics,
  file.path(
    TABLE_DIR,
    "22F4E4B_external_signature_curve_metrics.csv"
  ),
  row.names=FALSE
)

# ------------------------------------------------------------------------------
# 11. FORMAL CHECKS
# ------------------------------------------------------------------------------

stem_present <- gene_audit |>
  filter(
    signature=="STEM_HFSC"
  ) |>
  pull(present)

inflammatory_present <- gene_audit |>
  filter(
    signature==
      "INFLAMMATORY_TRANSITION"
  ) |>
  pull(present)

remodel_present <- gene_audit |>
  filter(
    signature==
      "REMODELING_ASSOCIATED"
  )

formal_checks <- tibble(
  check=c(
    "upstream_exists",
    "Hallmark_freeze_exists",
    "n_pseudobulks_56",
    "n_patients_16",
    "design_identical_to_22F4E4A",
    "pseudobulk_order_identical",
    "trajectory_grid_identical",
    "stem_signature_3_of_3",
    "inflammatory_signature_3_of_3",
    "VIM_present",
    "ZEB1_recorded_absent",
    "all_scores_finite",
    "all_curves_finite",
    "no_new_gene_level_DE"
  ),

  passed=c(
    file.exists(UPSTREAM),
    file.exists(HALLMARK_FREEZE),
    nrow(metadata)==56,
    length(
      unique(
        metadata$patient_id
      )
    )==16,
    isTRUE(
      all.equal(
        design,
        h$design,
        tolerance=0
      )
    ),
    identical(
      colnames(signature_scores),
      colnames(h$Hallmark_scores)
    ),
    isTRUE(
      all.equal(
        trajectory_grid,
        h$trajectory_grid,
        tolerance=0
      )
    ),
    length(stem_present)==3 &&
      all(stem_present),
    length(inflammatory_present)==3 &&
      all(inflammatory_present),
    remodel_present$present[
      remodel_present$gene=="VIM"
    ],
    !remodel_present$present[
      remodel_present$gene=="ZEB1"
    ],
    all(
      is.finite(
        signature_scores
      )
    ),
    all(
      is.finite(
        HC_curve
      )
    ) &&
      all(
        is.finite(
          FFA_curve
        )
      ) &&
      all(
        is.finite(
          AA_curve
        )
      ),
    TRUE
  )
)

write.csv(
  formal_checks,
  file.path(
    TABLE_DIR,
    "22F4E4B_formal_checks.csv"
  ),
  row.names=FALSE
)

cat("\n===== FORMAL CHECKS =====\n")
print(formal_checks)

if (!all(formal_checks$passed)) {
  stop(
    "One or more formal checks failed."
  )
}

# ------------------------------------------------------------------------------
# 12. FIGURE
# ------------------------------------------------------------------------------

# Frozen disease palette used throughout the paper.
disease_cols <- c(
  HC="#8A8A8A",
  FFA="#8FD694",
  AA="#D94B4B"
)

# Frozen boundaries from the existing Bulge trajectory.
cut_T3_T2 <- 0.438
cut_T2_T1 <- 0.844

# Use the same coordinate system in every facet.
curve_plot <- ggplot(
  curve_long,
  aes(
    x=trajectory_position,
    y=fitted_score,
    colour=disease
  )
) +
  annotate(
    "rect",
    xmin=-Inf,
    xmax=cut_T3_T2,
    ymin=-Inf,
    ymax=Inf,
    fill="grey95",
    alpha=0.35
  ) +
  annotate(
    "rect",
    xmin=cut_T3_T2,
    xmax=cut_T2_T1,
    ymin=-Inf,
    ymax=Inf,
    fill="grey90",
    alpha=0.22
  ) +
  annotate(
    "rect",
    xmin=cut_T2_T1,
    xmax=Inf,
    ymin=-Inf,
    ymax=Inf,
    fill="grey95",
    alpha=0.35
  ) +
  geom_vline(
    xintercept=c(
      cut_T3_T2,
      cut_T2_T1
    ),
    linewidth=0.35,
    linetype="dashed",
    colour="grey45"
  ) +
  geom_hline(
    yintercept=0,
    linewidth=0.3,
    colour="grey75"
  ) +
  geom_line(
    linewidth=1.15,
    lineend="round"
  ) +
  facet_wrap(
    ~signature_label,
    ncol=1,
    scales="free_y"
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
  scale_x_continuous(
    limits=c(
      0,
      1
    ),
    breaks=c(
      0,
      0.2,
      0.4,
      0.6,
      0.8,
      1
    ),
    expand=c(
      0,
      0
    )
  ) +
  labs(
    title=
      "External HFSC programmes show disease-dependent kinetics along the shared Bulge trajectory",
    subtitle=
      "A priori LPP-derived HFSC identity and inflammatory-transition signatures projected onto the frozen HC–FFA–AA trajectory",
    x=
      "Frozen Bulge/Isthmus trajectory",
    y=
      "Fitted standardized signature score",
    colour=NULL
  ) +
  theme_classic(
    base_size=10
  ) +
  theme(
    plot.title=element_text(
      face="bold",
      size=12,
      hjust=0
    ),
    plot.subtitle=element_text(
      size=9,
      colour="grey30",
      margin=margin(
        b=8
      )
    ),
    strip.background=element_rect(
      fill="grey96",
      colour="grey75",
      linewidth=0.4
    ),
    strip.text=element_text(
      face="bold",
      size=10
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
    panel.spacing=grid::unit(
      0.65,
      "lines"
    ),
    plot.margin=margin(
      8,
      10,
      8,
      8
    )
  )

# Add T3 / T2 / T1 labels only once above the complete figure.
curve_plot <- curve_plot +
  annotate(
    "text",
    x=cut_T3_T2/2,
    y=Inf,
    label="T3",
    vjust=1.5,
    size=3.2,
    fontface="bold",
    colour="grey30"
  ) +
  annotate(
    "text",
    x=(cut_T3_T2 + cut_T2_T1)/2,
    y=Inf,
    label="T2",
    vjust=1.5,
    size=3.2,
    fontface="bold",
    colour="grey30"
  ) +
  annotate(
    "text",
    x=(cut_T2_T1 + 1)/2,
    y=Inf,
    label="T1",
    vjust=1.5,
    size=3.2,
    fontface="bold",
    colour="grey30"
  )

# ------------------------------------------------------------------------------
# 13. SAVE FIGURE
# ------------------------------------------------------------------------------

output_stem <-
  "Fig22F4E4B_external_HFSC_signature_trajectory_preflight"

ggsave(
  filename=file.path(
    SVG_DIR,
    paste0(
      output_stem,
      ".svg"
    )
  ),
  plot=curve_plot,
  width=7.2,
  height=8.2,
  units="in",
  device=svglite::svglite
)

ggsave(
  filename=file.path(
    PDF_DIR,
    paste0(
      output_stem,
      ".pdf"
    )
  ),
  plot=curve_plot,
  width=7.2,
  height=8.2,
  units="in"
)

ggsave(
  filename=file.path(
    PNG_DIR,
    paste0(
      output_stem,
      ".png"
    )
  ),
  plot=curve_plot,
  width=7.2,
  height=8.2,
  units="in",
  dpi=400
)

# ------------------------------------------------------------------------------
# 14. FREEZE
# ------------------------------------------------------------------------------

freeze <- list(
  module=MODULE,

  analysis=
    "External LPP-derived HFSC programme projection onto frozen Bulge trajectory",

  source=
    "Frozen 22F4E2 pseudobulks/design and frozen 22F4E4A spline basis",

  external_signature_source=
    "Chen et al., Cell & Bioscience 2025",

  new_gene_level_DE=FALSE,

  n_pseudobulks=
    nrow(metadata),

  n_patients=
    length(
      unique(
        metadata$patient_id
      )
    ),

  signatures_requested=
    signatures_requested,

  signatures_present=
    signatures_present,

  gene_audit=
    gene_audit,

  score_definition=
    paste(
      "Mean standardized TMM logCPM across represented",
      "a priori external signature genes"
    ),

  metadata=
    metadata,

  design=
    design,

  common_coef=
    common_coef,

  FFA_coef=
    FFA_coef,

  AA_coef=
    AA_coef,

  joint_coef=
    joint_coef,

  trajectory_grid=
    trajectory_grid,

  recovered_spline_basis=
    B,

  signature_scores=
    signature_scores,

  fit=
    fit,

  interaction_results=
    interaction_results,

  HC_curve=
    HC_curve,

  FFA_curve=
    FFA_curve,

  AA_curve=
    AA_curve,

  FFA_minus_HC=
    FFA_minus_HC,

  AA_minus_HC=
    AA_minus_HC,

  AA_minus_FFA=
    AA_minus_FFA,

  curve_metrics=
    curve_metrics,

  formal_checks=
    formal_checks,

  output_stem=
    output_stem
)

saveRDS(
  freeze,
  file.path(
    OBJECT_DIR,
    "22F4E4B_external_HFSC_signature_freeze.rds"
  )
)

# ------------------------------------------------------------------------------
# 15. FINAL REPORT
# ------------------------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("22F4E4B — FINAL REPORT\n")
cat("============================================================\n")

cat(
  "\nPseudobulks:",
  nrow(metadata),
  "\n"
)

cat(
  "Patients:",
  length(
    unique(
      metadata$patient_id
    )
  ),
  "\n"
)

cat(
  "Disease order: HC -> FFA -> AA\n"
)

cat(
  "Trajectory: frozen 22F4E2 / 22F4E4A\n"
)

cat(
  "Trajectory boundaries: T3/T2 =",
  cut_T3_T2,
  "; T2/T1 =",
  cut_T2_T1,
  "\n"
)

cat(
  "New gene-level DE: NO\n"
)

cat(
  "Patient-aware frozen pseudobulk design: YES\n"
)

cat(
  "\n===== SIGNATURE COVERAGE =====\n"
)

coverage_report <- gene_audit |>
  group_by(signature) |>
  summarise(
    requested=n(),
    represented=sum(present),
    missing=paste(
      gene[!present],
      collapse=", "
    ),
    .groups="drop"
  ) |>
  mutate(
    missing=ifelse(
      missing=="",
      "NONE",
      missing
    )
  )

print(
  coverage_report,
  n=Inf
)

cat(
  "\n===== DISEASE x TRAJECTORY TESTS =====\n"
)

print(
  interaction_results |>
    arrange(global_PValue),
  n=Inf
)

cat(
  "\n===== CURVE METRICS =====\n"
)

print(
  curve_metrics,
  n=Inf
)

cat(
  "\n===== FORMAL CHECKS =====\n"
)

print(
  formal_checks,
  n=Inf
)

cat(
  "\nOutput directory:\n",
  OUT,
  "\n",
  sep=""
)

cat(
  "\nFigure stem:\n",
  output_stem,
  "\n",
  sep=""
)

cat(
  "\nInterpretation guardrail:\n",
  "REMODELING_ASSOCIATED is VIM-only because ZEB1 was not represented ",
  "in the frozen 22F4E2 tested gene universe. It must not be described ",
  "as a validated EMT signature.\n",
  sep=""
)

cat(
  "\nSTATUS=PASS\n"
)
