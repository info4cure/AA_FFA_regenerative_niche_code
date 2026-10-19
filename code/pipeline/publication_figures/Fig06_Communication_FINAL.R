# ============================================================
# FIGURE 6 — DISEASE-ASSOCIATED COMMUNICATION REMODELING
# Publication visualization only
#
# Frozen inputs:
#   21G4C6B — patient communication / cross-modal landscape
#   22E4     — final Bulge growth-survival molecular synthesis
#
# No new statistical inference
# No communication refit
# No imputation
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(svglite)
})

cat("\n============================================================\n")
cat("FIGURE 6 — COMMUNICATION REMODELING\n")
cat("============================================================\n")

# ============================================================
# 1. PATHS
# ============================================================

BASE21 <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "21_CELL_CELL_COMMUNICATION/",
  "21G_MULTICOMPARTMENT_COMMUNICATION/",
  "21G4_CROSSMODAL_INTEGRATION"
)

TAB21 <- file.path(BASE21, "tables")

BASE22 <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "22_MULTICOMPARTMENT_DISEASE_AXES/",
  "22E4_BULGE_GROWTH_SURVIVAL_FINAL_SYNTHESIS"
)

TAB22 <- file.path(BASE22, "tables")

OUT <- file.path(
  "PIPELINE_V2",
  "results",
  "publication_figures",
  "Figure06_Communication"
)

SVG_DIR <- file.path(OUT, "figures", "svg")
PDF_DIR <- file.path(OUT, "figures", "pdf")
PNG_DIR <- file.path(OUT, "figures", "png")
TAB_DIR <- file.path(OUT, "tables")
AUDIT_DIR <- file.path(OUT, "audit")

for (d in c(
  SVG_DIR,
  PDF_DIR,
  PNG_DIR,
  TAB_DIR,
  AUDIT_DIR
)) {
  dir.create(
    d,
    recursive=TRUE,
    showWarnings=FALSE
  )
}

# ============================================================
# 2. INPUT FILES
# ============================================================

LAND_FILE <- file.path(
  TAB21,
  "21G4C6B_signal_effect_landscape.csv"
)

AA_FILE <- file.path(
  TAB21,
  "21G4C6B_AA_vs_HC_patient_signals.csv"
)

FFA_FILE <- file.path(
  TAB21,
  "21G4C6B_FFA_vs_HC_patient_signals.csv"
)

PARENT_FILE <- file.path(
  TAB22,
  "22E4_parent_signal.csv"
)

MOLECULAR_FILE <- file.path(
  TAB22,
  "22E4_publication_molecular_summary.csv"
)

stopifnot(
  file.exists(LAND_FILE),
  file.exists(AA_FILE),
  file.exists(FFA_FILE),
  file.exists(PARENT_FILE),
  file.exists(MOLECULAR_FILE)
)

# ============================================================
# 3. LOAD FROZEN INPUTS
# ============================================================

land <- read_csv(
  LAND_FILE,
  show_col_types=FALSE
)

aa <- read_csv(
  AA_FILE,
  show_col_types=FALSE
)

ffa <- read_csv(
  FFA_FILE,
  show_col_types=FALSE
)

parent <- read_csv(
  PARENT_FILE,
  show_col_types=FALSE
)

molecular <- read_csv(
  MOLECULAR_FILE,
  show_col_types=FALSE
)

cat("Landscape signals :", nrow(land), "\n")
cat("AA patient rows   :", nrow(aa), "\n")
cat("FFA patient rows  :", nrow(ffa), "\n")
cat("Parent rows       :", nrow(parent), "\n")
cat("Molecular rows    :", nrow(molecular), "\n")

# ============================================================
# 4. CONSTANTS / HELPERS
# ============================================================

DISEASE_COLORS <- c(
  HC  = "#8F8D8E",
  FFA = "#4D9E4F",
  AA  = "#A8142A"
)

clean_receiver <- function(x) {

  x <- gsub(
    "_",
    " ",
    x,
    fixed=TRUE
  )

  x <- gsub(
    "IFE Infundibulum ",
    "IFE ",
    x,
    fixed=TRUE
  )

  x
}

theme_pub <- function(base_size=10) {

  theme_classic(
    base_size=base_size
  ) +
    theme(
      axis.text=element_text(
        color="black"
      ),
      axis.title=element_text(
        face="bold"
      ),
      strip.background=element_blank(),
      strip.text=element_text(
        face="bold"
      ),
      plot.title=element_text(
        face="bold",
        hjust=0,
        size=11
      ),
      plot.subtitle=element_text(
        hjust=0,
        size=8.5
      ),
      legend.title=element_text(
        face="bold"
      ),
      plot.margin=margin(
        6, 8, 6, 6
      )
    )
}

save_plot <- function(
  plot,
  stem,
  width,
  height
) {

  ggsave(
    file.path(
      SVG_DIR,
      paste0(stem, ".svg")
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
      PDF_DIR,
      paste0(stem, ".pdf")
    ),
    plot=plot,
    width=width,
    height=height,
    units="in",
    bg="white"
  )

  ggsave(
    file.path(
      PNG_DIR,
      paste0(stem, ".png")
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
# 5. PANEL A
# AA IGF2R-ASSOCIATED FOLLICULAR COMMUNICATION
# ============================================================

aa_igf <- aa |>
  filter(
    communication_program == "IGF2R_CONTEXTUAL",
    receiver %in% c(
      "Bulge_Isthmus",
      "ORS_basal",
      "ORS"
    )
  ) |>
  mutate(
    receiver_display=recode(
      receiver,
      Bulge_Isthmus="Bulge / isthmus",
      ORS_basal="Basal ORS",
      ORS="ORS"
    ),
    receiver_display=factor(
      receiver_display,
      levels=c(
        "Bulge / isthmus",
        "Basal ORS",
        "ORS"
      )
    ),
    disease=factor(
      disease,
      levels=c(
        "HC",
        "AA"
      )
    )
  )

# Frozen communication statistics for Panel A
stats_A <- land |>
  filter(
    contrast=="AA_vs_HC",
    communication_program=="IGF2R_CONTEXTUAL",
    receiver %in% c(
      "Bulge_Isthmus",
      "ORS_basal",
      "ORS"
    )
  ) |>
  mutate(
    receiver_display=recode(
      receiver,
      Bulge_Isthmus="Bulge / isthmus",
      ORS_basal="Basal ORS",
      ORS="ORS"
    ),
    receiver_display=factor(
      receiver_display,
      levels=c(
        "Bulge / isthmus",
        "Basal ORS",
        "ORS"
      )
    ),
    q_label=paste0(
      "FDR=",
      sprintf("%.3f", communication_q_BH)
    ),
    x1=1,
    x2=2,
    y=max(
      aa_igf$mean_incoming_A_axis,
      na.rm=TRUE
    ) + 0.035
  )

pA <- ggplot(
  aa_igf,
  aes(
    x=disease,
    y=mean_incoming_A_axis
  )
) +

  geom_violin(
    aes(fill=disease),
    width=0.82,
    alpha=0.16,
    linewidth=0.40,
    colour="grey35",
    trim=FALSE
  ) +

  geom_boxplot(
    aes(fill=disease),
    width=0.24,
    outlier.shape=NA,
    alpha=0.25,
    linewidth=0.40
  ) +

  geom_point(
    aes(fill=disease),
    position=position_jitter(
      width=0.055,
      height=0
    ),
    shape=21,
    colour="#303030",
    stroke=0.40,
    size=2.7,
    alpha=0.95
  ) +

  facet_wrap(
    ~receiver_display,
    nrow=1
  ) +

  geom_segment(
    data=stats_A,
    aes(
      x=x1,
      xend=x2,
      y=y,
      yend=y
    ),
    inherit.aes=FALSE,
    linewidth=0.40,
    colour="black"
  ) +

  geom_segment(
    data=stats_A,
    aes(
      x=x1,
      xend=x1,
      y=y,
      yend=y-0.008
    ),
    inherit.aes=FALSE,
    linewidth=0.40,
    colour="black"
  ) +

  geom_segment(
    data=stats_A,
    aes(
      x=x2,
      xend=x2,
      y=y,
      yend=y-0.008
    ),
    inherit.aes=FALSE,
    linewidth=0.40,
    colour="black"
  ) +

  geom_text(
    data=stats_A,
    aes(
      x=1.5,
      y=y+0.009,
      label=q_label
    ),
    inherit.aes=FALSE,
    size=3.0,
    vjust=0
  ) +

  scale_fill_manual(
    values=DISEASE_COLORS,
    guide="none",
    drop=FALSE
  ) +

  scale_y_continuous(
    expand=expansion(
      mult=c(0.05,0.18)
    )
  ) +

  labs(
    title=
      "A  Increased IGF2R-associated communication toward AA follicular receivers",
    subtitle=NULL,
    x=NULL,
    y="Mean incoming communication (A-axis)"
  ) +

  theme_pub(10)

# ============================================================
# 6. PANEL B
# AA BULGE GROWTH-SURVIVAL MOLECULAR REWIRING
# ============================================================

parent_B <- parent |>
  transmute(
    level="Parent programme",
    feature="Growth / survival",
    effect=mean_difference,
    FDR=BH_q,
    direction=ifelse(
      mean_difference > 0,
      "Higher in AA",
      "Lower in AA"
    )
  )

molecular_B <- molecular |>
  transmute(
    level=case_when(
      level=="LIGAND" ~ "Ligand",
      level=="EXACT_LR" ~ "Ligand–receptor",
      TRUE ~ level
    ),
    feature=gsub(
      "→",
      " \u2192 ",
      feature,
      fixed=TRUE
    ),
    effect=mean_difference,
    FDR=BH_q,
    direction=ifelse(
      mean_difference > 0,
      "Higher in AA",
      "Lower in AA"
    )
  )

bdat <- bind_rows(
  parent_B,
  molecular_B
) |>
  mutate(
    evidence=case_when(
      FDR < 0.05 ~ "FDR < 0.05",
      FDR < 0.10 ~ "FDR 0.05–0.10",
      TRUE ~ "FDR ≥ 0.10"
    ),
    evidence=factor(
      evidence,
      levels=c(
        "FDR < 0.05",
        "FDR 0.05–0.10",
        "FDR ≥ 0.10"
      )
    ),
    level=factor(
      level,
      levels=c(
        "Parent programme",
        "Ligand",
        "Ligand–receptor"
      )
    ),
    FDR_label=paste0(
      "FDR=",
      sprintf("%.3f", FDR)
    )
  )

feature_order <- c(
  "Growth / survival",
  "HBEGF",
  "EREG",
  "HBEGF → ERBB2",
  "HBEGF → EGFR",
  "EREG → EGFR",
  "EREG → ERBB2",
  "EREG → ERBB3"
)

bdat <- bdat |>
  mutate(
    feature=factor(
      feature,
      levels=rev(feature_order)
    )
  )

pB <- ggplot(
  bdat,
  aes(
    x=effect,
    y=feature
  )
) +

  geom_vline(
    xintercept=0,
    linewidth=0.35,
    colour="grey65",
    linetype="dashed"
  ) +

  geom_segment(
    aes(
      x=0,
      xend=effect,
      y=feature,
      yend=feature
    ),
    linewidth=0.60,
    colour="grey60"
  ) +

  geom_point(
    aes(
      fill=direction,
      alpha=evidence
    ),
    shape=21,
    size=3.8,
    colour="#303030",
    stroke=0.50
  ) +

  geom_text(
    aes(
      label=FDR_label
    ),
    x=Inf,
    hjust=1.05,
    size=2.65,
    colour="grey25"
  ) +

  facet_grid(
    level~.,
    scales="free_y",
    space="free_y"
  ) +

  scale_fill_manual(
    values=c(
      "Higher in AA"=DISEASE_COLORS["AA"],
      "Lower in AA"="#4F6D8A"
    ),
    name="Direction"
  ) +

  scale_alpha_manual(
    values=c(
      "FDR < 0.05"=1.00,
      "FDR 0.05–0.10"=0.45,
      "FDR ≥ 0.10"=0.20
    ),
    breaks=c(
      "FDR < 0.05",
      "FDR 0.05–0.10"
    ),
    name="Evidence"
  ) +

  scale_x_continuous(
    expand=expansion(
      mult=c(0.08,0.30)
    )
  ) +

  labs(
    title=
      "B  AA shows ligand-specific rewiring of Bulge growth-survival communication",
    subtitle=NULL,
    x="AA − HC mean communication difference",
    y=NULL
  ) +

  theme_pub(9.5) +

  theme(
    legend.position="top",
    legend.box="horizontal",
    panel.spacing.y=grid::unit(
      0.45,
      "lines"
    )
  )

# ============================================================
# 7. PANEL C
# FFA DISTRIBUTED GROWTH-SURVIVAL REDUCTION
# ============================================================

fg <- ffa |>
  filter(
    communication_program=="GROWTH_SURVIVAL"
  ) |>
  mutate(
    receiver=clean_receiver(receiver),
    disease=factor(
      disease,
      levels=c(
        "HC",
        "FFA"
      )
    )
  )

ord <- land |>
  filter(
    contrast=="FFA_vs_HC",
    communication_program=="GROWTH_SURVIVAL"
  ) |>
  arrange(
    communication_mean_difference
  ) |>
  transmute(
    receiver=clean_receiver(receiver)
  ) |>
  pull(receiver)

fg$receiver <- factor(
  fg$receiver,
  levels=unique(ord)
)

# Frozen communication statistics for Panel C.
# Display only prespecified communication FDR < 0.10 results.
stats_C <- land |>
  filter(
    contrast=="FFA_vs_HC",
    communication_program=="GROWTH_SURVIVAL",
    communication_q_BH < 0.10
  ) |>
  mutate(
    receiver=clean_receiver(receiver),
    receiver=factor(
      receiver,
      levels=levels(fg$receiver)
    ),
    q_label=paste0(
      "FDR=",
      sprintf("%.3f", communication_q_BH)
    )
  )

# Receiver-specific bracket heights.
# These are plotting coordinates only; no inference is performed.
c_y <- fg |>
  group_by(receiver) |>
  summarise(
    y=max(
      mean_incoming_A_axis,
      na.rm=TRUE
    ) + 0.025,
    .groups="drop"
  )

stats_C <- stats_C |>
  left_join(
    c_y,
    by="receiver"
  ) |>
  mutate(
    x=as.numeric(receiver),
    x1=x-0.22,
    x2=x+0.22
  )

pC <- ggplot(
  fg,
  aes(
    x=receiver,
    y=mean_incoming_A_axis
  )
) +

  geom_boxplot(
    aes(fill=disease),
    position=position_dodge(
      0.66
    ),
    width=0.46,
    alpha=0.20,
    outlier.shape=NA,
    linewidth=0.42
  ) +

  geom_point(
    aes(fill=disease),
    position=position_jitterdodge(
      jitter.width=0.055,
      dodge.width=0.66
    ),
    shape=21,
    size=2.6,
    stroke=0.40,
    colour="#303030"
  ) +

  geom_segment(
    data=stats_C,
    aes(
      x=x1,
      xend=x2,
      y=y,
      yend=y
    ),
    inherit.aes=FALSE,
    linewidth=0.38,
    colour="black"
  ) +

  geom_segment(
    data=stats_C,
    aes(
      x=x1,
      xend=x1,
      y=y,
      yend=y-0.007
    ),
    inherit.aes=FALSE,
    linewidth=0.38,
    colour="black"
  ) +

  geom_segment(
    data=stats_C,
    aes(
      x=x2,
      xend=x2,
      y=y,
      yend=y-0.007
    ),
    inherit.aes=FALSE,
    linewidth=0.38,
    colour="black"
  ) +

  geom_text(
    data=stats_C,
    aes(
      x=x,
      y=y+0.007,
      label=q_label
    ),
    inherit.aes=FALSE,
    size=2.65,
    vjust=0
  ) +

  scale_fill_manual(
    values=DISEASE_COLORS,
    name=NULL,
    drop=FALSE
  ) +

  scale_y_continuous(
    expand=expansion(
      mult=c(0.05,0.18)
    )
  ) +

  labs(
    title=
      "C  Reduced growth-survival communication across FFA epithelial states",
    subtitle=NULL,
    x=NULL,
    y="Mean incoming communication (A-axis)"
  ) +

  theme_pub(9.5) +

  theme(
    axis.text.x=element_text(
      angle=30,
      hjust=1,
      vjust=1
    ),
    legend.position="top"
  )

# ============================================================
# 8. PANEL D
# SUPPORTED COMMUNICATION LANDSCAPE
# ============================================================

d_contrast_labels <- c(
  AA_vs_HC="AA vs HC",
  FFA_vs_HC="FFA vs HC",
  AA_vs_FFA="AA vs FFA"
)

d_program_labels <- c(
  IGF2R_CONTEXTUAL="IGF2R-associated",
  GROWTH_SURVIVAL="Growth / survival",
  CYTOTOXIC_INJURY="Cytotoxic injury",
  IMMUNE_REGULATION="Immune regulation",
  IGF_TROPHIC="IGF trophic"
)

ddat <- land |>
  mutate(
    receiver_display=clean_receiver(receiver),
    program_display=ifelse(
      communication_program %in% names(d_program_labels),
      unname(
        d_program_labels[
          communication_program
        ]
      ),
      gsub(
        "_",
        " ",
        communication_program
      )
    ),
    signal_label=paste0(
      receiver_display,
      " | ",
      program_display
    ),
    contrast_display=recode(
      contrast,
      !!!d_contrast_labels
    ),
    evidence=case_when(
      communication_q_BH < 0.05 ~
        "FDR < 0.05",
      communication_q_BH < 0.10 ~
        "FDR 0.05–0.10",
      TRUE ~
        "FDR ≥ 0.10"
    ),
    evidence=factor(
      evidence,
      levels=c(
        "FDR < 0.05",
        "FDR 0.05–0.10",
        "FDR ≥ 0.10"
      )
    ),
    FDR_label=ifelse(
      communication_q_BH < 0.10,
      sprintf(
        "%.3f",
        communication_q_BH
      ),
      ""
    )
  )

# Order signals by contrast and effect in the frozen landscape.
signal_levels <- ddat |>
  arrange(
    factor(
      contrast,
      levels=c(
        "AA_vs_HC",
        "FFA_vs_HC",
        "AA_vs_FFA"
      )
    ),
    communication_mean_difference
  ) |>
  distinct(
    signal_label
  ) |>
  pull(
    signal_label
  )

ddat <- ddat |>
  mutate(
    signal_label=factor(
      signal_label,
      levels=rev(
        unique(signal_levels)
      )
    ),
    contrast_display=factor(
      contrast_display,
      levels=c(
        "AA vs HC",
        "FFA vs HC",
        "AA vs FFA"
      )
    )
  )

effect_limit <- max(
  abs(
    ddat$communication_mean_difference
  ),
  na.rm=TRUE
)

pD <- ggplot(
  ddat,
  aes(
    x=contrast_display,
    y=signal_label
  )
) +

  geom_point(
    aes(
      fill=communication_mean_difference,
      alpha=evidence
    ),
    shape=21,
    size=5.2,
    colour="#303030",
    stroke=0.55
  ) +

  geom_text(
    aes(
      label=FDR_label
    ),
    size=2.15,
    colour="black",
    fontface="bold"
  ) +

  scale_fill_gradient2(
    low="#3B6FB6",
    mid="white",
    high="#B6383A",
    midpoint=0,
    limits=c(
      -effect_limit,
      effect_limit
    ),
    name="Communication\nchange"
  ) +

  scale_alpha_manual(
    values=c(
      "FDR < 0.05"=1.00,
      "FDR 0.05–0.10"=0.55,
      "FDR ≥ 0.10"=0.16
    ),
    name="Evidence"
  ) +

  labs(
    title=
      "D  Supported communication changes converge on epithelial receiver states",
    subtitle=
      "Numbers inside supported points denote Benjamini–Hochberg FDR",
    x=NULL,
    y=NULL
  ) +

  theme_pub(9.2) +

  theme(
    legend.position="right",
    axis.text.x=element_text(
      face="bold"
    ),
    axis.text.y=element_text(
      size=7.7
    ),
    panel.grid.major.x=element_line(
      colour="grey92",
      linewidth=0.35
    )
  )


# ============================================================
# 8B. PUBLICATION REDESIGN — PANELS B AND D
# Frozen evidence only; no new inference
# ============================================================

cat("\n===== REDESIGN PANELS B AND D =====\n")

E4_DIR <- file.path(
  "PIPELINE_V2",
  "results",
  "disease_biology",
  "22_MULTICOMPARTMENT_DISEASE_AXES",
  "22E4_BULGE_GROWTH_SURVIVAL_FINAL_SYNTHESIS",
  "tables"
)

D4_DIR <- file.path(
  "PIPELINE_V2",
  "results",
  "disease_biology",
  "22_MULTICOMPARTMENT_DISEASE_AXES",
  "22D4_IGF2_IGF2R_FINAL_SYNTHESIS",
  "tables"
)

E4_PARENT <- readr::read_csv(
  file.path(
    E4_DIR,
    "22E4_parent_signal.csv"
  ),
  show_col_types=FALSE
)

E4_LR <- readr::read_csv(
  file.path(
    E4_DIR,
    "22E4_exactLR_synthesis.csv"
  ),
  show_col_types=FALSE
)

D4_IGF2 <- readr::read_csv(
  file.path(
    D4_DIR,
    "22D4_primary_disease_synthesis.csv"
  ),
  show_col_types=FALSE
)

# ------------------------------------------------------------
# PANEL B
# Molecular resolution of AA Bulge growth/survival rewiring
# ------------------------------------------------------------

b_parent <- E4_PARENT |>
  dplyr::filter(
    receiver=="Bulge_Isthmus",
    contrast=="AA_vs_HC",
    mechanistic_class=="GROWTH_SURVIVAL_TROPHIC"
  )

stopifnot(
  nrow(b_parent)==1,
  b_parent$BH_q[1] < 0.05
)

b_lr <- E4_LR |>
  dplyr::filter(
    BH_q_exactLR < 0.10
  ) |>
  dplyr::mutate(
    interaction=paste0(
      ligand,
      " \u2192 ",
      receptor
    ),
    evidence=dplyr::case_when(
      BH_q_exactLR < 0.05 ~ "FDR < 0.05",
      BH_q_exactLR < 0.10 ~ "FDR 0.05\u20130.10",
      TRUE ~ "Other"
    ),
    interaction=factor(
      interaction,
      levels=rev(
        c(
          "EREG \u2192 EGFR",
          "EREG \u2192 ERBB2",
          "EREG \u2192 ERBB3",
          "HBEGF \u2192 EGFR",
          "HBEGF \u2192 ERBB2"
        )
      )
    )
  )

stopifnot(
  nrow(b_lr)==5,
  sum(b_lr$BH_q_exactLR < 0.05)==1
)

pB <- ggplot(
  b_lr,
  aes(
    x=mean_difference,
    y=interaction
  )
) +

  geom_vline(
    xintercept=0,
    linewidth=0.45,
    colour="grey55"
  ) +

  geom_segment(
    aes(
      x=0,
      xend=mean_difference,
      y=interaction,
      yend=interaction
    ),
    linewidth=0.85,
    colour="grey70"
  ) +

  geom_point(
    aes(
      fill=direction,
      shape=evidence
    ),
    size=4.3,
    stroke=1.0,
    colour="#303030"
  ) +

  scale_fill_manual(
    values=c(
      HIGHER_IN_AA="#A8142A",
      LOWER_IN_AA="#3B6FB6"
    ),
    labels=c(
      HIGHER_IN_AA="Higher in AA",
      LOWER_IN_AA="Lower in AA"
    ),
    name=NULL
  ) +

  scale_shape_manual(
    values=c(
      "FDR < 0.05"=21,
      "FDR 0.05\u20130.10"=21
    ),
    name=NULL
  ) +

  geom_text(
    aes(
      label=sprintf(
        "FDR=%.3f",
        BH_q_exactLR
      )
    ),
    hjust=dplyr::if_else(
      b_lr$mean_difference >= 0,
      -0.12,
      1.12
    ),
    size=3.25,
    colour="#303030"
  ) +

  annotate(
    "text",
    x=0,
    y=5.65,
    label=sprintf(
      "Parent growth/survival: \u0394A = %.3f; FDR = %.3f",
      b_parent$mean_difference[1],
      b_parent$BH_q[1]
    ),
    hjust=0.5,
    size=3.45,
    fontface="bold"
  ) +

  coord_cartesian(
    xlim=c(-0.095,0.205),
    clip="off"
  ) +

  labs(
    title="B  Ligand-specific rewiring of AA Bulge growth/survival communication",
    subtitle="Reduced parent communication contains opposing EREG- and HBEGF-associated changes",
    x="AA \u2212 HC mean communication difference",
    y=NULL
  ) +

  theme_classic(
    base_size=10
  ) +

  theme(
    plot.title=element_text(
      face="bold",
      size=11
    ),
    plot.subtitle=element_text(
      size=8.8
    ),
    axis.title.x=element_text(
      face="bold",
      size=9
    ),
    axis.text.y=element_text(
      size=9,
      face="bold"
    ),
    legend.position="top",
    legend.justification="left",
    plot.margin=margin(
      8, 28, 8, 8
    )
  )

# ------------------------------------------------------------
# PANEL D
# Compact mechanistic synthesis
# ------------------------------------------------------------

d_igf <- D4_IGF2 |>
  dplyr::filter(
    contrast=="AA_vs_HC",
    receiver %in% c(
      "Bulge_Isthmus",
      "ORS_basal",
      "ORS"
    ),
    BH_q_primary_family < 0.05
  ) |>
  dplyr::transmute(
    disease="AA",
    receiver=dplyr::recode(
      receiver,
      Bulge_Isthmus="Bulge / isthmus",
      ORS_basal="Basal ORS",
      ORS="ORS"
    ),
    signal="IGF2\u2192IGF2R",
    change="Increased",
    effect=mean_difference,
    FDR=BH_q_primary_family,
    evidence="FDR < 0.05"
  )

d_growth <- b_lr |>
  dplyr::transmute(
    disease="AA",
    receiver="Bulge / isthmus",
    signal=as.character(interaction),
    change=dplyr::if_else(
      direction=="HIGHER_IN_AA",
      "Increased",
      "Decreased"
    ),
    effect=mean_difference,
    FDR=BH_q_exactLR,
    evidence=dplyr::if_else(
      BH_q_exactLR < 0.05,
      "FDR < 0.05",
      "FDR 0.05\u20130.10"
    )
  )

d_ffa <- land |>
  dplyr::filter(
    contrast=="FFA_vs_HC",
    communication_program=="GROWTH_SURVIVAL",
    communication_q_BH < 0.10
  ) |>
  dplyr::transmute(
    disease="FFA",
    receiver=clean_receiver(receiver),
    signal="Growth/survival",
    change="Decreased",
    effect=communication_mean_difference,
    FDR=communication_q_BH,
    evidence="FDR 0.05\u20130.10"
  )

d_table <- dplyr::bind_rows(
  d_igf,
  d_growth,
  d_ffa
) |>
  dplyr::mutate(
    direction_symbol=dplyr::if_else(
      change=="Increased",
      "\u2191",
      "\u2193"
    ),
    effect_label=sprintf(
      "%+.3f",
      effect
    ),
    FDR_label=sprintf(
      "%.3f",
      FDR
    ),
    row_id=dplyr::row_number()
  )

# Explicit publication order
d_table <- d_table |>
  dplyr::arrange(
    factor(
      disease,
      levels=c("AA","FFA")
    ),
    dplyr::case_when(
      signal=="IGF2\u2192IGF2R" ~ 1,
      grepl("EREG",signal) ~ 2,
      grepl("HBEGF",signal) ~ 3,
      TRUE ~ 4
    ),
    FDR
  ) |>
  dplyr::mutate(
    y=rev(
      seq_len(
        dplyr::n()
      )
    )
  )

# Column positions
x_disease  <- 0.02
x_receiver <- 0.14
x_signal   <- 0.42
x_dir      <- 0.70
x_effect   <- 0.79
x_fdr      <- 0.92

pD <- ggplot(
  d_table
) +

  annotate(
    "text",
    x=c(
      x_disease,
      x_receiver,
      x_signal,
      x_dir,
      x_effect,
      x_fdr
    ),
    y=max(d_table$y)+1.35,
    label=c(
      "Disease",
      "Receiver",
      "Molecular / communication signal",
      "Direction",
      "\u0394",
      "FDR"
    ),
    hjust=c(0,0,0,0.5,0.5,0.5),
    fontface="bold",
    size=3.35
  ) +

  geom_hline(
    yintercept=max(d_table$y)+0.80,
    linewidth=0.45,
    colour="grey35"
  ) +

  geom_text(
    aes(
      x=x_disease,
      y=y,
      label=disease
    ),
    hjust=0,
    size=3.25,
    fontface="bold"
  ) +

  geom_text(
    aes(
      x=x_receiver,
      y=y,
      label=receiver
    ),
    hjust=0,
    size=3.05
  ) +

  geom_text(
    aes(
      x=x_signal,
      y=y,
      label=signal
    ),
    hjust=0,
    size=3.05
  ) +

  geom_text(
    aes(
      x=x_dir,
      y=y,
      label=direction_symbol,
      colour=change,
      alpha=evidence
    ),
    size=4.2,
    fontface="bold"
  ) +

  geom_text(
    aes(
      x=x_effect,
      y=y,
      label=effect_label,
      alpha=evidence
    ),
    hjust=0.5,
    size=3.05
  ) +

  geom_text(
    aes(
      x=x_fdr,
      y=y,
      label=FDR_label,
      alpha=evidence
    ),
    hjust=0.5,
    size=3.05,
    fontface="bold"
  ) +

  scale_colour_manual(
    values=c(
      Increased="#A8142A",
      Decreased="#3B6FB6"
    ),
    guide="none"
  ) +

  scale_alpha_manual(
    values=c(
      "FDR < 0.05"=1,
      "FDR 0.05\u20130.10"=0.58
    ),
    guide="none"
  ) +

  scale_x_continuous(
    limits=c(0,1)
  ) +

  scale_y_continuous(
    limits=c(
      0,
      max(d_table$y)+1.8
    ),
    expand=c(0,0)
  ) +

  labs(
    title="D  Molecular resolution of epithelial communication remodeling",
    subtitle="Frozen multiplicity-supported signals; darker entries denote FDR <0.05",
    x=NULL,
    y=NULL
  ) +

  theme_void(
    base_size=10
  ) +

  theme(
    plot.title=element_text(
      face="bold",
      size=11,
      hjust=0
    ),
    plot.subtitle=element_text(
      size=8.8,
      hjust=0
    ),
    plot.margin=margin(
      8,8,8,8
    )
  )

cat(
  "Panel B exact-LR rows:",
  nrow(b_lr),
  "\n"
)

cat(
  "Panel D summary rows:",
  nrow(d_table),
  "\n"
)



# ============================================================
# 9. SAVE INDIVIDUAL PANELS
# ============================================================

cat("\n===== SAVE INDIVIDUAL PANELS =====\n")

save_plot(
  pA,
  "Figure06_A_AA_IGF2R_communication",
  width=10.5,
  height=4.5
)

save_plot(
  pB,
  "Figure06_B_AA_Bulge_molecular_rewiring",
  width=7.8,
  height=5.8
)

save_plot(
  pC,
  "Figure06_C_FFA_growth_survival_communication",
  width=13.0,
  height=5.8
)

save_plot(
  pD,
  "Figure06_D_molecular_summary",
  width=13.0,
  height=7.2
)

# ============================================================
# 10. FINAL FIGURE ASSEMBLY
# ============================================================

cat("\n===== BUILD FINAL FIGURE 6 =====\n")

# Narrative:
#
# A  AA shows increased IGF2-IGF2R-associated communication
#    toward Bulge/Isthmus, basal ORS and ORS receivers.
#
# B  AA Bulge growth/survival communication is reduced overall,
#    but exact-LR decomposition reveals opposing EREG-high and
#    HBEGF-low molecular components.
#
# C  FFA shows distributed reduction of epithelial
#    growth/survival communication.
#
# D  Compact molecular synthesis connects receiver-level
#    communication changes with their supported molecular
#    resolution.
#
# Layout:
#
#   +---------------- A ----------------+--------- B ---------+
#
#   +---------------------- C ------------------------------+
#
#   +---------------------- D ------------------------------+

row_AB <- patchwork::wrap_plots(
  pA,
  pB,
  nrow=1,
  widths=c(
    1.48,
    1.00
  )
)

final_plot <- patchwork::wrap_plots(
  row_AB,
  pC,
  pD,
  ncol=1,
  heights=c(
    1.00,
    0.88,
    1.08
  )
) +
  patchwork::plot_annotation(
    title=
      "AA and FFA exhibit distinct patterns of epithelial trophic communication remodeling",
    subtitle=
      paste(
        "AA shows focal molecular rewiring of follicular communication,",
        "whereas FFA shows distributed attenuation of epithelial growth-survival communication"
      ),
    theme=theme(
      plot.title=element_text(
        face="bold",
        size=15,
        hjust=0
      ),
      plot.subtitle=element_text(
        size=9.5,
        hjust=0,
        margin=margin(
          b=8
        )
      )
    )
  )

# ============================================================
# 11. SAVE FINAL COMPOSITE
# ============================================================

save_plot(
  final_plot,
  "Figure06_FINAL_preflight",
  width=15.5,
  height=16.2
)

# ============================================================
# 12. SAVE PUBLICATION SOURCE TABLES
# ============================================================

cat("\n===== SAVE SOURCE TABLES =====\n")

write.csv(
  aa_igf,
  file.path(
    TAB_DIR,
    "Figure06_A_source.csv"
  ),
  row.names=FALSE
)

# Panel B:
# preserve both parent-level and exact-LR evidence.

write.csv(
  b_parent,
  file.path(
    TAB_DIR,
    "Figure06_B_parent_source.csv"
  ),
  row.names=FALSE
)

write.csv(
  b_lr,
  file.path(
    TAB_DIR,
    "Figure06_B_exactLR_source.csv"
  ),
  row.names=FALSE
)

write.csv(
  fg,
  file.path(
    TAB_DIR,
    "Figure06_C_source.csv"
  ),
  row.names=FALSE
)

write.csv(
  d_table,
  file.path(
    TAB_DIR,
    "Figure06_D_source.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 13. FORMAL CHECKS
# ============================================================

cat("\n===== FORMAL CHECKS =====\n")

checks <- tibble::tibble(
  check=c(
    "AA IGF2R panel nonempty",
    "AA IGF2R has three receivers",
    "AA IGF2R contains HC and AA",

    "Parent growth-survival row present",
    "Parent growth-survival FDR<0.05",

    "Five exact LR interactions FDR<0.10",
    "EREG-EGFR exact LR FDR<0.05",
    "EREG exact-LR effects positive",
    "HBEGF exact-LR effects negative",

    "FFA growth-survival panel nonempty",
    "FFA growth-survival contains HC and FFA",
    "FFA has four communication FDR<0.10 receivers",

    "Panel D summary nonempty",
    "Panel D contains AA and FFA",
    "Panel D contains IGF2-IGF2R",
    "No missing panel-B effects",

    "Panel A SVG exists",
    "Panel B SVG exists",
    "Panel C SVG exists",
    "Panel D SVG exists",

    "Composite SVG exists",
    "Composite PDF exists",
    "Composite PNG exists"
  ),

  passed=c(
    nrow(aa_igf) > 0,

    dplyr::n_distinct(
      aa_igf$receiver_display
    ) == 3,

    all(
      c("HC","AA") %in%
        as.character(
          unique(
            aa_igf$disease
          )
        )
    ),

    nrow(b_parent) == 1,

    all(
      b_parent$BH_q < 0.05
    ),

    nrow(b_lr) == 5 &&
      all(
        b_lr$BH_q_exactLR < 0.10
      ),

    any(
      b_lr$ligand=="EREG" &
        b_lr$receptor=="EGFR" &
        b_lr$BH_q_exactLR < 0.05
    ),

    all(
      b_lr$mean_difference[
        b_lr$ligand=="EREG"
      ] > 0
    ),

    all(
      b_lr$mean_difference[
        b_lr$ligand=="HBEGF"
      ] < 0
    ),

    nrow(fg) > 0,

    all(
      c("HC","FFA") %in%
        as.character(
          unique(
            fg$disease
          )
        )
    ),

    sum(
      land$contrast=="FFA_vs_HC" &
        land$communication_program=="GROWTH_SURVIVAL" &
        land$communication_q_BH < 0.10,
      na.rm=TRUE
    ) == 4,

    nrow(d_table) > 0,

    all(
      c("AA","FFA") %in%
        d_table$disease
    ),

    any(
      d_table$signal=="IGF2\u2192IGF2R"
    ),

    all(
      is.finite(
        b_lr$mean_difference
      )
    ),

    file.exists(
      file.path(
        SVG_DIR,
        "Figure06_A_AA_IGF2R_communication.svg"
      )
    ),

    file.exists(
      file.path(
        SVG_DIR,
        "Figure06_B_AA_Bulge_molecular_rewiring.svg"
      )
    ),

    file.exists(
      file.path(
        SVG_DIR,
        "Figure06_C_FFA_growth_survival_communication.svg"
      )
    ),

    file.exists(
      file.path(
        SVG_DIR,
        "Figure06_D_molecular_summary.svg"
      )
    ),

    file.exists(
      file.path(
        SVG_DIR,
        "Figure06_FINAL_preflight.svg"
      )
    ),

    file.exists(
      file.path(
        PDF_DIR,
        "Figure06_FINAL_preflight.pdf"
      )
    ),

    file.exists(
      file.path(
        PNG_DIR,
        "Figure06_FINAL_preflight.png"
      )
    )
  )
)

print(checks)

write.csv(
  checks,
  file.path(
    AUDIT_DIR,
    "Figure06_formal_checks.csv"
  ),
  row.names=FALSE
)



# ============================================================
# 14. AUDIT / REPORT
# ============================================================

audit_lines <- c(
  "============================================================",
  "FIGURE 6 — COMMUNICATION REMODELING",
  "============================================================",
  "",
  "Publication visualization only.",
  "No new statistical inference.",
  "No communication-network refit.",
  "No imputation.",
  "",
  paste0(
    "AA IGF2R patient rows: ",
    nrow(aa_igf)
  ),
  paste0(
    "FFA growth/survival patient rows: ",
    nrow(fg)
  ),
  paste0(
    "Cross-modal signals: ",
    dplyr::n_distinct(
      ddat$signal_id
    )
  ),
  "",
  "Panel B frozen evidence:",
  paste0(
    "Parent Bulge growth/survival: AA-HC mean difference = ",
    sprintf(
      "%.3f",
      parent$mean_difference[1]
    ),
    "; q=",
    sprintf(
      "%.3f",
      parent$BH_q[1]
    )
  ),
  paste0(
    "EREG ligand: q=",
    sprintf(
      "%.3f",
      molecular$BH_q[
        molecular$level=="LIGAND" &
          molecular$feature=="EREG"
      ]
    )
  ),
  paste0(
    "HBEGF ligand: q=",
    sprintf(
      "%.3f",
      molecular$BH_q[
        molecular$level=="LIGAND" &
          molecular$feature=="HBEGF"
      ]
    )
  ),
  paste0(
    "EREG-EGFR: q=",
    sprintf(
      "%.3f",
      molecular$BH_q[
        molecular$level=="EXACT_LR" &
          molecular$feature=="EREG→EGFR"
      ]
    )
  ),
  "",
  paste0(
    "FORMAL CHECKS: ",
    sum(checks$passed),
    "/",
    nrow(checks),
    " PASS"
  ),
  "",
  if (all(checks$passed)) {
    "STATUS=PASS"
  } else {
    "STATUS=FAIL"
  }
)

writeLines(
  audit_lines,
  file.path(
    AUDIT_DIR,
    "Figure06_audit.txt"
  )
)

cat(
  paste(
    audit_lines,
    collapse="\n"
  ),
  "\n"
)

# ============================================================
# 15. FINAL STATUS
# ============================================================

if (!all(checks$passed)) {
  stop(
    "Figure 6 formal checks failed."
  )
}

cat("\n============================================================\n")
cat("FIGURE 6 COMPLETE\n")
cat("STATUS=PASS\n")
cat("============================================================\n")
