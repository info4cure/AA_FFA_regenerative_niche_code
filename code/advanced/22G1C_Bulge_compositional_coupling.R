
# ============================================================
# 22G1C — BULGE-CENTERED COMPOSITIONAL COUPLING ATLAS
#
# Patient-level exploratory analysis
# Bulge T1/T2/T3 versus whole-atlas frozen L2 composition
#
# No cell-level pseudoreplication.
# No causal interpretation.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(svglite)
})

cat("\n============================================================\n")
cat("22G1C — BULGE-CENTERED COMPOSITIONAL COUPLING ATLAS\n")
cat("============================================================\n")

# ============================================================
# 1. PATHS
# ============================================================

ROOT <- "PIPELINE_V2"

K3_FILE <- file.path(
  ROOT,
  "results/disease_biology",
  "22_MULTICOMPARTMENT_DISEASE_AXES",
  "22F_BULGE_STEM_NICHE_INTEGRITY",
  "22F3_BULGE_SPATIAL_TERRITORIES",
  "22F3B_K3_PSEUDOBULK",
  "tables",
  "22F3B_patient_K3_composition.csv"
)

L2_FILE <- file.path(
  ROOT,
  "results/disease_biology",
  "20_CELLULAR_COMPOSITION",
  "20C_PATIENT_LEVEL_COMPOSITION",
  "tables",
  "20C_L2_patient_frequencies_within_L1.csv"
)

OUT <- file.path(
  ROOT,
  "results/disease_biology",
  "22_MULTICOMPARTMENT_DISEASE_AXES",
  "22G_BULGE_CENTERED_ECOSYSTEM",
  "22G1_PATIENT_COUPLING_ATLAS",
  "22G1C_COMPOSITIONAL_COUPLING"
)

TAB <- file.path(OUT,"tables")
OBJ <- file.path(OUT,"objects")
AUD <- file.path(OUT,"audit")

SVG <- file.path(OUT,"figures","svg")
PDF <- file.path(OUT,"figures","pdf")
PNG <- file.path(OUT,"figures","png")

for(d in c(TAB,OBJ,AUD,SVG,PDF,PNG)) {
  dir.create(
    d,
    recursive=TRUE,
    showWarnings=FALSE
  )
}

stopifnot(
  file.exists(K3_FILE),
  file.exists(L2_FILE)
)

# ============================================================
# 2. LOAD
# ============================================================

K3 <- read.csv(
  K3_FILE,
  stringsAsFactors=FALSE
)

L2 <- read.csv(
  L2_FILE,
  stringsAsFactors=FALSE
)

cat("\n===== INPUT =====\n")
cat("K3 rows:",nrow(K3),"\n")
cat("L2 rows:",nrow(L2),"\n")
cat("Patients K3:",n_distinct(K3$patient_id),"\n")
cat("Patients L2:",n_distinct(L2$patient_id),"\n")

# ============================================================
# 3. PATIENT MASTER
# ============================================================

patients <- bind_rows(
  K3 |>
    dplyr::select(
      patient_id,
      disease
    ),
  L2 |>
    dplyr::select(
      patient_id,
      disease
    )
) |>
  distinct() |>
  arrange(
    factor(
      disease,
      levels=c("HC","FFA","AA")
    ),
    patient_id
  )

stopifnot(
  nrow(patients) == 16
)

# ============================================================
# 4. COMPLETE BULGE K3 COMPOSITION
# ============================================================

cat("\n===== COMPLETE K3 COMPOSITION =====\n")

K3_complete <- patients |>
  crossing(
    territory=c(
      "T1",
      "T2",
      "T3"
    )
  ) |>
  left_join(
    K3 |>
      dplyr::select(
        patient_id,
        disease,
        territory,
        n_cells,
        total_bulge
      ),
    by=c(
      "patient_id",
      "disease",
      "territory"
    )
  ) |>
  group_by(
    patient_id,
    disease
  ) |>
  mutate(
    total_bulge=max(
      total_bulge,
      na.rm=TRUE
    ),
    n_cells=replace_na(
      n_cells,
      0
    ),
    pct_bulge=
      100*n_cells/total_bulge
  ) |>
  ungroup()

stopifnot(
  nrow(K3_complete) == 48,
  all(
    abs(
      K3_complete |>
        group_by(patient_id) |>
        summarise(
          s=sum(pct_bulge),
          .groups="drop"
        ) |>
        pull(s) - 100
    ) < 1e-6
  )
)

bulge_wide <- K3_complete |>
  dplyr::select(
    patient_id,
    disease,
    territory,
    pct_bulge
  ) |>
  pivot_wider(
    names_from=territory,
    values_from=pct_bulge,
    names_prefix="Bulge_"
  )

print(
  bulge_wide,
  n=Inf
)

# ============================================================
# 5. EXTERNAL L2 UNIVERSE
# ============================================================

cat("\n===== EXTERNAL L2 UNIVERSE =====\n")

l2_clean <- L2 |>
  filter(
    eligible,
    !is.na(frequency),
    L2 != "Bulge_Isthmus"
  ) |>
  mutate(
    feature=paste(
      L1,
      L2,
      sep=" :: "
    )
  )

support <- l2_clean |>
  group_by(
    L1,
    L2,
    feature
  ) |>
  summarise(
    n_patients=n_distinct(patient_id),
    n_HC=n_distinct(
      patient_id[disease=="HC"]
    ),
    n_FFA=n_distinct(
      patient_id[disease=="FFA"]
    ),
    n_AA=n_distinct(
      patient_id[disease=="AA"]
    ),
    total_cells=sum(n),
    .groups="drop"
  ) |>
  mutate(
    eligible_coupling=
      n_patients >= 12 &
      n_HC >= 3 &
      n_FFA >= 3 &
      n_AA >= 3 &
      total_cells >= 40
  )

write.csv(
  support,
  file.path(
    TAB,
    "22G1C_L2_support.csv"
  ),
  row.names=FALSE
)

eligible_features <- support |>
  filter(
    eligible_coupling
  ) |>
  pull(feature)

cat(
  "Eligible external L2 features:",
  length(eligible_features),
  "\n"
)

# ============================================================
# 6. PATIENT × L2 MATRIX
# ============================================================

l2_matrix_long <- l2_clean |>
  filter(
    feature %in%
      eligible_features
  ) |>
  dplyr::select(
    patient_id,
    disease,
    L1,
    L2,
    feature,
    frequency,
    percent,
    n
  )

write.csv(
  l2_matrix_long,
  file.path(
    TAB,
    "22G1C_patient_L2_matrix_long.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 7. GLOBAL BULGE-L2 SPEARMAN COUPLING
# ============================================================

cat("\n===== GLOBAL SPEARMAN COUPLING =====\n")

bulge_axes <- c(
  "Bulge_T3",
  "Bulge_T2",
  "Bulge_T1"
)

global_results <- list()

counter <- 1L

for(feature_i in eligible_features) {

  ext <- l2_matrix_long |>
    filter(
      feature == feature_i
    ) |>
    dplyr::select(
      patient_id,
      disease,
      L1,
      L2,
      feature,
      frequency
    )

  dat <- bulge_wide |>
    inner_join(
      ext,
      by=c(
        "patient_id",
        "disease"
      )
    )

  for(axis_i in bulge_axes) {

    x <- dat[[axis_i]]
    y <- dat$frequency

    ok <- is.finite(x) &
      is.finite(y)

    if(sum(ok) >= 8) {

      ct <- suppressWarnings(
        cor.test(
          x[ok],
          y[ok],
          method="spearman",
          exact=FALSE
        )
      )

      global_results[[counter]] <- tibble(
        feature=feature_i,
        L1=dat$L1[1],
        L2=dat$L2[1],
        bulge_axis=axis_i,
        n=sum(ok),
        rho=unname(ct$estimate),
        p=ct$p.value
      )

      counter <- counter+1L
    }
  }
}

global_results <- bind_rows(
  global_results
) |>
  group_by(
    bulge_axis
  ) |>
  mutate(
    FDR=p.adjust(
      p,
      method="BH"
    )
  ) |>
  ungroup()

write.csv(
  global_results,
  file.path(
    TAB,
    "22G1C_global_Bulge_L2_Spearman.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 8. WITHIN-DISEASE DIRECTIONAL ROBUSTNESS
# ============================================================

cat("\n===== WITHIN-DISEASE ROBUSTNESS =====\n")

within_results <- list()
counter <- 1L

for(feature_i in eligible_features) {

  ext <- l2_matrix_long |>
    filter(
      feature == feature_i
    ) |>
    dplyr::select(
      patient_id,
      disease,
      L1,
      L2,
      feature,
      frequency
    )

  dat <- bulge_wide |>
    inner_join(
      ext,
      by=c(
        "patient_id",
        "disease"
      )
    )

  for(dis_i in c(
    "HC",
    "FFA",
    "AA"
  )) {

    dd <- dat |>
      filter(
        disease == dis_i
      )

    for(axis_i in bulge_axes) {

      x <- dd[[axis_i]]
      y <- dd$frequency

      ok <- is.finite(x) &
        is.finite(y)

      rho <- NA_real_

      if(
        sum(ok) >= 4 &&
        sd(x[ok]) > 0 &&
        sd(y[ok]) > 0
      ) {
        rho <- suppressWarnings(
          cor(
            x[ok],
            y[ok],
            method="spearman"
          )
        )
      }

      within_results[[counter]] <- tibble(
        feature=feature_i,
        L1=dd$L1[1],
        L2=dd$L2[1],
        disease=dis_i,
        bulge_axis=axis_i,
        n=sum(ok),
        rho=rho
      )

      counter <- counter+1L
    }
  }
}

within_results <- bind_rows(
  within_results
)

write.csv(
  within_results,
  file.path(
    TAB,
    "22G1C_within_disease_directional_coupling.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 9. DISEASE EFFECT ON EXTERNAL L2
# ============================================================

cat("\n===== EXTERNAL L2 DISEASE EFFECT =====\n")

disease_results <- list()
counter <- 1L

for(feature_i in eligible_features) {

  dat <- l2_matrix_long |>
    filter(
      feature == feature_i
    )

  # Kruskal-Wallis
  kw <- kruskal.test(
    frequency ~ disease,
    data=dat
  )

  med <- dat |>
    group_by(
      disease
    ) |>
    summarise(
      median=median(
        frequency,
        na.rm=TRUE
      ),
      .groups="drop"
    )

  get_med <- function(d) {
    z <- med$median[
      med$disease==d
    ]
    if(length(z)==0) NA_real_ else z
  }

  # AA vs FFA
  aa <- dat |>
    filter(
      disease=="AA"
    ) |>
    pull(frequency)

  ffa <- dat |>
    filter(
      disease=="FFA"
    ) |>
    pull(frequency)

  wt <- suppressWarnings(
    wilcox.test(
      aa,
      ffa,
      exact=FALSE
    )
  )

  disease_results[[counter]] <- tibble(
    feature=feature_i,
    L1=dat$L1[1],
    L2=dat$L2[1],
    median_HC=get_med("HC"),
    median_FFA=get_med("FFA"),
    median_AA=get_med("AA"),
    delta_AA_FFA=
      get_med("AA")-
      get_med("FFA"),
    KW_p=kw$p.value,
    AA_FFA_p=wt$p.value
  )

  counter <- counter+1L
}

disease_results <- bind_rows(
  disease_results
) |>
  mutate(
    KW_FDR=p.adjust(
      KW_p,
      method="BH"
    ),
    AA_FFA_FDR=p.adjust(
      AA_FFA_p,
      method="BH"
    )
  )

write.csv(
  disease_results,
  file.path(
    TAB,
    "22G1C_external_L2_disease_effects.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 10. INTEGRATED COUPLING TABLE
# ============================================================

cat("\n===== INTEGRATED COUPLING =====\n")

within_wide <- within_results |>
  dplyr::select(
    feature,
    bulge_axis,
    disease,
    rho
  ) |>
  pivot_wider(
    names_from=disease,
    values_from=rho,
    names_prefix="rho_"
  )

integrated <- global_results |>
  left_join(
    within_wide,
    by=c(
      "feature",
      "bulge_axis"
    )
  ) |>
  left_join(
    disease_results,
    by=c(
      "feature",
      "L1",
      "L2"
    )
  ) |>
  mutate(
    abs_rho=abs(rho),

    same_direction_AA_FFA=
      !is.na(rho_AA) &
      !is.na(rho_FFA) &
      sign(rho_AA)==sign(rho_FFA),

    direction_support_n=
      rowSums(
        cbind(
          sign(rho_HC)==sign(rho),
          sign(rho_FFA)==sign(rho),
          sign(rho_AA)==sign(rho)
        ),
        na.rm=TRUE
      ),

    coupling_class=case_when(
      FDR < 0.10 &
        abs_rho >= 0.60 &
        direction_support_n >= 2
        ~ "ROBUST",

      FDR < 0.10 &
        abs_rho >= 0.50
        ~ "SUPPORTED",

      abs_rho >= 0.50
        ~ "EXPLORATORY",

      TRUE
        ~ "WEAK"
    ),

    disease_divergent=
      AA_FFA_FDR < 0.10,

    priority_score=
      abs_rho +
      0.20*(FDR < 0.10) +
      0.15*(direction_support_n >= 2) +
      0.20*(AA_FFA_FDR < 0.10)
  ) |>
  arrange(
    desc(priority_score),
    FDR
  )

write.csv(
  integrated,
  file.path(
    TAB,
    "22G1C_integrated_Bulge_centered_coupling.csv"
  ),
  row.names=FALSE
)

cat("\nTOP COUPLINGS:\n")

print(
  integrated |>
    filter(
      coupling_class != "WEAK"
    ) |>
    dplyr::select(
      L1,
      L2,
      bulge_axis,
      rho,
      FDR,
      rho_HC,
      rho_FFA,
      rho_AA,
      delta_AA_FFA,
      AA_FFA_FDR,
      coupling_class
    ) |>
    slice_head(
      n=40
    ),
  n=Inf
)

# ============================================================
# 11. HEATMAP FEATURE SELECTION
# ============================================================

# One row per external L2.
# Select features showing the strongest coupling to ANY Bulge axis.

display_features <- integrated |>
  group_by(
    feature,
    L1,
    L2
  ) |>
  arrange(
    desc(priority_score),
    desc(abs_rho),
    .by_group=TRUE
  ) |>
  slice_head(
    n=1
  ) |>
  ungroup() |>
  arrange(
    desc(priority_score)
  ) |>
  slice_head(
    n=30
  )

display_order <- display_features |>
  arrange(
    factor(
      L1,
      levels=c(
        "Epithelial",
        "Fibroblast",
        "Immune",
        "Endothelial",
        "Mural",
        "Sebocyte",
        "Melanocyte",
        "Lymphatic EC"
      )
    ),
    desc(priority_score)
  ) |>
  pull(feature)

heat_df <- integrated |>
  filter(
    feature %in% display_order
  ) |>
  mutate(
    feature=factor(
      feature,
      levels=rev(display_order)
    ),

    bulge_axis=factor(
      bulge_axis,
      levels=c(
        "Bulge_T3",
        "Bulge_T2",
        "Bulge_T1"
      ),
      labels=c(
        "T3",
        "T2",
        "T1"
      )
    ),

    sig_label=case_when(
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      FDR < 0.10 ~ "†",
      TRUE ~ ""
    )
  )

write.csv(
  display_features,
  file.path(
    TAB,
    "22G1C_display_feature_selection.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 12. PANEL A — BULGE-CENTERED COUPLING HEATMAP
# ============================================================

cat("\n===== FIGURE A — COUPLING HEATMAP =====\n")

pA <- ggplot(
  heat_df,
  aes(
    x=bulge_axis,
    y=feature,
    fill=rho
  )
) +

  geom_tile(
    colour="white",
    linewidth=0.45
  ) +

  geom_text(
    aes(
      label=sig_label
    ),
    size=3,
    fontface="bold"
  ) +

  scale_fill_gradient2(
    low="#2166AC",
    mid="white",
    high="#B2182B",
    midpoint=0,
    limits=c(-1,1),
    name="Spearman\nrho"
  ) +

  labs(
    title="Bulge-centered compositional coupling across the skin atlas",
    subtitle=paste(
      "Patient-level Spearman correlations;",
      "† FDR<0.10, * FDR<0.05, ** FDR<0.01"
    ),
    x="Bulge territory",
    y=NULL
  ) +

  theme_classic(
    base_size=9
  ) +

  theme(
    plot.title=element_text(
      face="bold",
      size=11
    ),
    plot.subtitle=element_text(
      size=8
    ),
    axis.text.y=element_text(
      size=6.7
    ),
    axis.text.x=element_text(
      face="bold",
      size=8.5
    ),
    axis.ticks.y=element_blank(),
    axis.line.y=element_blank(),
    legend.position="right"
  )

# ============================================================
# 13. PANEL B — DISEASE DIVERGENCE VS BULGE COUPLING
# ============================================================

cat("\n===== FIGURE B — DISEASE DIVERGENCE =====\n")

scatter_df <- integrated |>
  group_by(
    feature,
    L1,
    L2
  ) |>
  slice_max(
    order_by=abs_rho,
    n=1,
    with_ties=FALSE
  ) |>
  ungroup() |>
  mutate(
    highlight=
      abs_rho >= 0.50 |
      AA_FFA_FDR < 0.10,

    label=ifelse(
      highlight,
      L2,
      NA_character_
    )
  )

pB <- ggplot(
  scatter_df,
  aes(
    x=delta_AA_FFA,
    y=abs_rho
  )
) +

  geom_hline(
    yintercept=0.5,
    linetype="dashed",
    linewidth=0.4,
    colour="#888888"
  ) +

  geom_vline(
    xintercept=0,
    linewidth=0.35,
    colour="#888888"
  ) +

  geom_point(
    aes(
      shape=L1,
      alpha=highlight
    ),
    size=2.5
  ) +

  scale_alpha_manual(
    values=c(
      `FALSE`=0.35,
      `TRUE`=0.9
    ),
    guide="none"
  ) +

  labs(
    title="Coupling strength versus AA-FFA compositional divergence",
    subtitle=paste(
      "Positive x: relatively enriched in AA;",
      "negative x: relatively enriched in FFA"
    ),
    x="Median within-L1 frequency difference  (AA - FFA)",
    y="Maximum |rho| with Bulge T1/T2/T3",
    shape="Compartment"
  ) +

  theme_classic(
    base_size=9
  ) +

  theme(
    plot.title=element_text(
      face="bold",
      size=10.5
    ),
    plot.subtitle=element_text(
      size=8
    ),
    legend.position="right"
  )

# ============================================================
# 14. PANEL C — TOP COUPLING LEDGER
# ============================================================

top_ledger <- integrated |>
  filter(
    coupling_class != "WEAK"
  ) |>
  group_by(
    feature
  ) |>
  slice_max(
    order_by=priority_score,
    n=1,
    with_ties=FALSE
  ) |>
  ungroup() |>
  arrange(
    desc(priority_score)
  ) |>
  slice_head(
    n=20
  ) |>
  mutate(
    direction=ifelse(
      rho > 0,
      "Positive",
      "Negative"
    ),
    axis_clean=gsub(
      "Bulge_",
      "",
      bulge_axis
    )
  )

write.csv(
  top_ledger,
  file.path(
    TAB,
    "22G1C_top_candidate_couplings.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 15. SAVE FUNCTION
# ============================================================

save_plot <- function(
  plot,
  stem,
  width,
  height
) {

  ggsave(
    file.path(
      SVG,
      paste0(stem,".svg")
    ),
    plot=plot,
    width=width,
    height=height,
    units="in",
    device=svglite::svglite,
    bg="white"
  )

  ggsave(
    file.path(
      PDF,
      paste0(stem,".pdf")
    ),
    plot=plot,
    width=width,
    height=height,
    units="in",
    bg="white"
  )

  ggsave(
    file.path(
      PNG,
      paste0(stem,".png")
    ),
    plot=plot,
    width=width,
    height=height,
    units="in",
    dpi=400,
    bg="white"
  )
}

# ============================================================
# 16. SAVE FIGURES
# ============================================================

cat("\n===== SAVE FIGURES =====\n")

save_plot(
  pA,
  "Fig22G1C_A_Bulge_centered_coupling_heatmap",
  width=8.5,
  height=9.5
)

save_plot(
  pB,
  "Fig22G1C_B_coupling_vs_AA_FFA_divergence",
  width=8.0,
  height=6.0
)

final_plot <- pA + pB +
  patchwork::plot_layout(
    widths=c(
      1.15,
      1
    )
  ) +
  patchwork::plot_annotation(
    title=
      "Patient-level ecosystem remodeling coupled to Bulge state",
    theme=theme(
      plot.title=element_text(
        face="bold",
        size=14
      )
    )
  )

save_plot(
  final_plot,
  "Fig22G1C_Bulge_centered_ecosystem",
  width=16,
  height=9.5
)

# ============================================================
# 17. FREEZE
# ============================================================

freeze <- list(
  module="22G1C",

  purpose=paste(
    "Patient-level exploratory coupling of frozen Bulge",
    "T1/T2/T3 composition to frozen whole-atlas L2",
    "composition."
  ),

  statistical_unit="patient",

  n_patients=n_distinct(
    patients$patient_id
  ),

  disease_n=
    patients |>
      count(disease),

  Bulge_K3=K3_complete,

  support=support,

  global_coupling=global_results,

  within_disease_coupling=within_results,

  disease_effects=disease_results,

  integrated=integrated,

  top_candidates=top_ledger
)

saveRDS(
  freeze,
  file.path(
    OBJ,
    "22G1C_Bulge_centered_compositional_coupling_freeze.rds"
  )
)

# ============================================================
# 18. FORMAL CHECKS
# ============================================================

cat("\n===== FORMAL CHECKS =====\n")

checks <- tibble(
  check=c(
    "16 patients",
    "48 completed patient-K3 rows",
    "K3 sums to 100 percent",
    "Bulge excluded from external L2",
    "Eligible external features > 0",
    "Three Bulge axes tested",
    "Global coupling table non-empty",
    "Within-disease table non-empty",
    "Disease-effect table non-empty",
    "Integrated table non-empty",
    "SVG heatmap exists",
    "SVG scatter exists",
    "SVG composite exists",
    "Freeze exists"
  ),

  pass=c(
    n_distinct(
      patients$patient_id
    ) == 16,

    nrow(K3_complete) == 48,

    all(
      abs(
        K3_complete |>
          group_by(patient_id) |>
          summarise(
            x=sum(pct_bulge),
            .groups="drop"
          ) |>
          pull(x) - 100
      ) < 1e-6
    ),

    !any(
      l2_matrix_long$L2 ==
        "Bulge_Isthmus"
    ),

    length(
      eligible_features
    ) > 0,

    all(
      bulge_axes %in%
        unique(
          global_results$bulge_axis
        )
    ),

    nrow(global_results) > 0,

    nrow(within_results) > 0,

    nrow(disease_results) > 0,

    nrow(integrated) > 0,

    file.exists(
      file.path(
        SVG,
        "Fig22G1C_A_Bulge_centered_coupling_heatmap.svg"
      )
    ),

    file.exists(
      file.path(
        SVG,
        "Fig22G1C_B_coupling_vs_AA_FFA_divergence.svg"
      )
    ),

    file.exists(
      file.path(
        SVG,
        "Fig22G1C_Bulge_centered_ecosystem.svg"
      )
    ),

    file.exists(
      file.path(
        OBJ,
        "22G1C_Bulge_centered_compositional_coupling_freeze.rds"
      )
    )
  )
)

print(
  checks,
  n=Inf
)

write.csv(
  checks,
  file.path(
    AUD,
    "22G1C_formal_checks.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 19. AUDIT
# ============================================================

audit_lines <- c(
  "============================================================",
  "22G1C — BULGE-CENTERED COMPOSITIONAL COUPLING",
  "============================================================",
  "",
  paste0(
    "Patients: ",
    n_distinct(
      patients$patient_id
    )
  ),
  paste0(
    "Eligible external L2 features: ",
    length(
      eligible_features
    )
  ),
  paste0(
    "Global coupling tests: ",
    nrow(
      global_results
    )
  ),
  paste0(
    "Robust couplings: ",
    sum(
      integrated$coupling_class ==
        "ROBUST",
      na.rm=TRUE
    )
  ),
  paste0(
    "Supported couplings: ",
    sum(
      integrated$coupling_class ==
        "SUPPORTED",
      na.rm=TRUE
    )
  ),
  paste0(
    "Exploratory couplings: ",
    sum(
      integrated$coupling_class ==
        "EXPLORATORY",
      na.rm=TRUE
    )
  ),
  "",
  "Interpretation constraints:",
  "- Patient is the statistical unit.",
  "- Within-disease correlations are directional robustness only.",
  "- Composition is within frozen L1 compartments.",
  "- Bulge_Isthmus is excluded from external features.",
  "- Correlation does not establish directionality or causality.",
  "",
  paste0(
    "FORMAL CHECKS: ",
    sum(checks$pass),
    "/",
    nrow(checks),
    " PASS"
  ),
  "",
  if(all(checks$pass)) {
    "STATUS=PASS"
  } else {
    "STATUS=FAIL"
  }
)

writeLines(
  audit_lines,
  file.path(
    AUD,
    "22G1C_audit.txt"
  )
)

cat(
  paste(
    audit_lines,
    collapse="\n"
  ),
  "\n"
)

if(!all(checks$pass)) {
  stop(
    "22G1C formal checks failed."
  )
}

cat("\n============================================================\n")
cat("22G1C COMPLETE\n")
cat("STATUS=PASS\n")
cat("============================================================\n")

