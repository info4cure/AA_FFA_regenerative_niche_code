suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

cat("============================================================\n")
cat("21G4D1 — QUANTITATIVE FIGURE DISCOVERY\n")
cat("============================================================\n")

BASE <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "21_CELL_CELL_COMMUNICATION/",
  "21G_MULTICOMPARTMENT_COMMUNICATION/",
  "21G4_CROSSMODAL_INTEGRATION"
)

TAB <- file.path(BASE, "tables")
FIG <- file.path(BASE, "figures", "21G4D1")
dir.create(FIG, recursive=TRUE, showWarnings=FALSE)

land <- read_csv(
  file.path(TAB, "21G4C6B_signal_effect_landscape.csv"),
  show_col_types=FALSE
)

aa <- read_csv(
  file.path(TAB, "21G4C6B_AA_vs_HC_patient_signals.csv"),
  show_col_types=FALSE
)

ffa <- read_csv(
  file.path(TAB, "21G4C6B_FFA_vs_HC_patient_signals.csv"),
  show_col_types=FALSE
)

ad <- read_csv(
  file.path(TAB, "21G4C6B_patient_AD_geometry.csv"),
  show_col_types=FALSE
)

DISEASE_COLORS <- c(
  HC  = "#8F8D8E",
  FFA = "#4D9E4F",
  AA  = "#A8142A"
)

theme_pub <- function(base_size=11) {
  theme_classic(base_size=base_size) +
    theme(
      axis.title=element_text(face="bold"),
      axis.text=element_text(color="black"),
      strip.background=element_blank(),
      strip.text=element_text(face="bold"),
      legend.title=element_text(face="bold"),
      plot.title=element_text(face="bold", hjust=0.5, size=13),
      plot.subtitle=element_text(hjust=0.5, size=10),
      plot.margin=margin(10,12,10,10)
    )
}

save_plot <- function(p, stem, width, height) {
  ggsave(
    file.path(FIG, paste0(stem, ".pdf")),
    p, width=width, height=height,
    units="in", device="pdf"
  )
  ggsave(
    file.path(FIG, paste0(stem, ".png")),
    p, width=width, height=height,
    units="in", dpi=500, bg="white"
  )
}

clean_receiver <- function(x) {
  x <- gsub("_", " ", x, fixed=TRUE)
  x <- gsub(
    "IFE Infundibulum ",
    "IFE ",
    x,
    fixed=TRUE
  )
  x
}

cat("Frozen datasets loaded\n")
cat("Landscape:", nrow(land), "signals\n")
cat("AA patient rows:", nrow(aa), "\n")
cat("FFA patient rows:", nrow(ffa), "\n")
cat("A-D geometry:", nrow(ad), "rows\n")

# ============================================================
# A — AA IGF2R FOLLICULAR ARCHITECTURE
# ============================================================

aa_igf <- aa %>%
  filter(
    communication_program == "IGF2R_CONTEXTUAL",
    receiver %in% c("Bulge_Isthmus","ORS_basal","ORS")
  ) %>%
  mutate(
    receiver_display=recode(
      receiver,
      Bulge_Isthmus="Bulge / isthmus",
      ORS_basal="Basal ORS",
      ORS="ORS"
    ),
    receiver_display=factor(
      receiver_display,
      levels=c("Bulge / isthmus","Basal ORS","ORS")
    )
  )

pA <- ggplot(
  aa_igf,
  aes(x=disease, y=mean_incoming_A_axis)
) +
  geom_violin(
    aes(fill=disease),
    width=.82, alpha=.16,
    linewidth=.45, color="grey35",
    trim=FALSE
  ) +
  geom_boxplot(
    aes(fill=disease),
    width=.24, outlier.shape=NA,
    alpha=.25, linewidth=.45
  ) +
  geom_point(
    aes(fill=disease),
    position=position_jitter(width=.07),
    shape=21, color="#303030",
    stroke=.45, size=3.2, alpha=.95
  ) +
  facet_wrap(~receiver_display, nrow=1) +
  scale_fill_manual(values=DISEASE_COLORS, guide="none") +
  labs(
    title="AA follicular IGF2R-associated communication",
    subtitle="Patient-level incoming IGF-trophic communication",
    x=NULL,
    y="Mean incoming A-axis"
  ) +
  theme_pub(12)

save_plot(
  pA,
  "Fig21G4D1_A_AA_IGF2R_patient_architecture",
  10.5, 4.5
)

# ============================================================
# B — AA BULGE: IGF2R vs GROWTH/SURVIVAL
# ============================================================

aa_bulge <- aa %>%
  filter(
    receiver == "Bulge_Isthmus",
    communication_program %in%
      c("IGF2R_CONTEXTUAL","GROWTH_SURVIVAL")
  ) %>%
  mutate(
    program_display=recode(
      communication_program,
      IGF2R_CONTEXTUAL="IGF2R-associated\nIGF trophic",
      GROWTH_SURVIVAL="Growth /\nsurvival"
    )
  )

pB <- ggplot(
  aa_bulge,
  aes(
    x=program_display,
    y=mean_incoming_A_axis,
    group=interaction(patient_id,disease)
  )
) +
  geom_line(
    aes(color=disease),
    linewidth=.65, alpha=.38
  ) +
  geom_point(
    aes(fill=disease),
    shape=21, color="#303030",
    stroke=.45, size=3.4
  ) +
  scale_color_manual(values=DISEASE_COLORS, guide="none") +
  scale_fill_manual(values=DISEASE_COLORS, name=NULL) +
  labs(
    title="Divergent trophic architecture in the AA bulge/isthmus niche",
    subtitle="Patient-resolved communication across two biological programmes",
    x=NULL,
    y="Mean incoming A-axis"
  ) +
  theme_pub(12) +
  theme(legend.position="top")

save_plot(
  pB,
  "Fig21G4D1_B_AA_bulge_program_divergence",
  7.5, 5.2
)

# ============================================================
# C — FFA DISTRIBUTED GROWTH/SURVIVAL LOSS
# ============================================================

ffa_growth <- ffa %>%
  filter(communication_program == "GROWTH_SURVIVAL")

ord <- land %>%
  filter(
    contrast=="FFA_vs_HC",
    communication_program=="GROWTH_SURVIVAL"
  ) %>%
  arrange(communication_mean_difference) %>%
  pull(receiver)

ffa_growth <- ffa_growth %>%
  mutate(
    receiver_display=clean_receiver(receiver),
    receiver_display=factor(
      receiver_display,
      levels=clean_receiver(ord)
    )
  )

pC <- ggplot(
  ffa_growth,
  aes(
    x=receiver_display,
    y=mean_incoming_A_axis
  )
) +
  geom_boxplot(
    aes(fill=disease),
    position=position_dodge(width=.68),
    width=.52,
    outlier.shape=NA,
    alpha=.18,
    linewidth=.45
  ) +
  geom_point(
    aes(fill=disease),
    position=position_jitterdodge(
      jitter.width=.10,
      dodge.width=.68
    ),
    shape=21,
    color="#303030",
    stroke=.35,
    size=2.8,
    alpha=.95
  ) +
  scale_fill_manual(values=DISEASE_COLORS, name=NULL) +
  labs(
    title="Distributed epithelial growth-survival communication loss in FFA",
    subtitle="Patient-resolved architecture across follicular and interfollicular epithelial states",
    x=NULL,
    y="Mean incoming A-axis"
  ) +
  theme_pub(11) +
  theme(
    axis.text.x=element_text(angle=38,hjust=1),
    legend.position="top"
  )

save_plot(
  pC,
  "Fig21G4D1_C_FFA_growth_survival_landscape",
  13.2, 6.3
)

# ============================================================
# D — A–D COMMUNICATION GEOMETRY
# ============================================================

ad_main <- ad %>%
  filter(
    signal_id %in% c(
      "21G4PUB_02","21G4PUB_03",
      "21G4PUB_04","21G4PUB_05",
      "21G4PUB_08","21G4PUB_09",
      "21G4PUB_10","21G4PUB_11"
    )
  ) %>%
  mutate(
    panel_label=paste0(
      clean_receiver(receiver),
      "\n",
      gsub("_"," ",communication_program)
    )
  )

pD <- ggplot(
  ad_main,
  aes(
    x=pooled_D_axis,
    y=mean_incoming_A_axis,
    fill=disease
  )
) +
  geom_point(
    shape=21,
    color="#303030",
    stroke=.45,
    size=3.3,
    alpha=.94
  ) +
  facet_grid(
    contrast ~ panel_label,
    scales="free",
    space="free"
  ) +
  scale_fill_manual(values=DISEASE_COLORS, name=NULL) +
  labs(
    title="Patient-level communication geometry",
    subtitle="Incoming communication strength versus active-sender density",
    x="Pooled sender-density axis (D)",
    y="Mean incoming communication axis (A)"
  ) +
  theme_pub(10) +
  theme(
    legend.position="top",
    strip.text=element_text(size=8)
  )

save_plot(
  pD,
  "Fig21G4D1_D_patient_A_D_geometry",
  15, 7.5
)

cat("Figures A-D generated\n")

# ============================================================
# A — AA IGF2R FOLLICULAR ARCHITECTURE
# ============================================================

aa_igf <- aa %>%
  filter(
    communication_program == "IGF2R_CONTEXTUAL",
    receiver %in% c("Bulge_Isthmus","ORS_basal","ORS")
  ) %>%
  mutate(
    receiver_display=recode(
      receiver,
      Bulge_Isthmus="Bulge / isthmus",
      ORS_basal="Basal ORS",
      ORS="ORS"
    ),
    receiver_display=factor(
      receiver_display,
      levels=c("Bulge / isthmus","Basal ORS","ORS")
    )
  )

pA <- ggplot(
  aa_igf,
  aes(x=disease, y=mean_incoming_A_axis)
) +
  geom_violin(
    aes(fill=disease),
    width=.82, alpha=.16,
    linewidth=.45, color="grey35",
    trim=FALSE
  ) +
  geom_boxplot(
    aes(fill=disease),
    width=.24, outlier.shape=NA,
    alpha=.25, linewidth=.45
  ) +
  geom_point(
    aes(fill=disease),
    position=position_jitter(width=.07),
    shape=21, color="#303030",
    stroke=.45, size=3.2, alpha=.95
  ) +
  facet_wrap(~receiver_display, nrow=1) +
  scale_fill_manual(values=DISEASE_COLORS, guide="none") +
  labs(
    title="AA follicular IGF2R-associated communication",
    subtitle="Patient-level incoming IGF-trophic communication",
    x=NULL,
    y="Mean incoming A-axis"
  ) +
  theme_pub(12)

save_plot(
  pA,
  "Fig21G4D1_A_AA_IGF2R_patient_architecture",
  10.5, 4.5
)

# ============================================================
# B — AA BULGE: IGF2R vs GROWTH/SURVIVAL
# ============================================================

aa_bulge <- aa %>%
  filter(
    receiver == "Bulge_Isthmus",
    communication_program %in%
      c("IGF2R_CONTEXTUAL","GROWTH_SURVIVAL")
  ) %>%
  mutate(
    program_display=recode(
      communication_program,
      IGF2R_CONTEXTUAL="IGF2R-associated\nIGF trophic",
      GROWTH_SURVIVAL="Growth /\nsurvival"
    )
  )

pB <- ggplot(
  aa_bulge,
  aes(
    x=program_display,
    y=mean_incoming_A_axis,
    group=interaction(patient_id,disease)
  )
) +
  geom_line(
    aes(color=disease),
    linewidth=.65, alpha=.38
  ) +
  geom_point(
    aes(fill=disease),
    shape=21, color="#303030",
    stroke=.45, size=3.4
  ) +
  scale_color_manual(values=DISEASE_COLORS, guide="none") +
  scale_fill_manual(values=DISEASE_COLORS, name=NULL) +
  labs(
    title="Divergent trophic architecture in the AA bulge/isthmus niche",
    subtitle="Patient-resolved communication across two biological programmes",
    x=NULL,
    y="Mean incoming A-axis"
  ) +
  theme_pub(12) +
  theme(legend.position="top")

save_plot(
  pB,
  "Fig21G4D1_B_AA_bulge_program_divergence",
  7.5, 5.2
)

# ============================================================
# C — FFA DISTRIBUTED GROWTH/SURVIVAL LOSS
# ============================================================

ffa_growth <- ffa %>%
  filter(communication_program == "GROWTH_SURVIVAL")

ord <- land %>%
  filter(
    contrast=="FFA_vs_HC",
    communication_program=="GROWTH_SURVIVAL"
  ) %>%
  arrange(communication_mean_difference) %>%
  pull(receiver)

ffa_growth <- ffa_growth %>%
  mutate(
    receiver_display=clean_receiver(receiver),
    receiver_display=factor(
      receiver_display,
      levels=clean_receiver(ord)
    )
  )

pC <- ggplot(
  ffa_growth,
  aes(
    x=receiver_display,
    y=mean_incoming_A_axis
  )
) +
  geom_boxplot(
    aes(fill=disease),
    position=position_dodge(width=.68),
    width=.52,
    outlier.shape=NA,
    alpha=.18,
    linewidth=.45
  ) +
  geom_point(
    aes(fill=disease),
    position=position_jitterdodge(
      jitter.width=.10,
      dodge.width=.68
    ),
    shape=21,
    color="#303030",
    stroke=.35,
    size=2.8,
    alpha=.95
  ) +
  scale_fill_manual(values=DISEASE_COLORS, name=NULL) +
  labs(
    title="Distributed epithelial growth-survival communication loss in FFA",
    subtitle="Patient-resolved architecture across follicular and interfollicular epithelial states",
    x=NULL,
    y="Mean incoming A-axis"
  ) +
  theme_pub(11) +
  theme(
    axis.text.x=element_text(angle=38,hjust=1),
    legend.position="top"
  )

save_plot(
  pC,
  "Fig21G4D1_C_FFA_growth_survival_landscape",
  13.2, 6.3
)

# ============================================================
# D — A–D COMMUNICATION GEOMETRY
# ============================================================

ad_main <- ad %>%
  filter(
    signal_id %in% c(
      "21G4PUB_02","21G4PUB_03",
      "21G4PUB_04","21G4PUB_05",
      "21G4PUB_08","21G4PUB_09",
      "21G4PUB_10","21G4PUB_11"
    )
  ) %>%
  mutate(
    panel_label=paste0(
      clean_receiver(receiver),
      "\n",
      gsub("_"," ",communication_program)
    )
  )

pD <- ggplot(
  ad_main,
  aes(
    x=pooled_D_axis,
    y=mean_incoming_A_axis,
    fill=disease
  )
) +
  geom_point(
    shape=21,
    color="#303030",
    stroke=.45,
    size=3.3,
    alpha=.94
  ) +
  facet_grid(
    contrast ~ panel_label,
    scales="free",
    space="free"
  ) +
  scale_fill_manual(values=DISEASE_COLORS, name=NULL) +
  labs(
    title="Patient-level communication geometry",
    subtitle="Incoming communication strength versus active-sender density",
    x="Pooled sender-density axis (D)",
    y="Mean incoming communication axis (A)"
  ) +
  theme_pub(10) +
  theme(
    legend.position="top",
    strip.text=element_text(size=8)
  )

save_plot(
  pD,
  "Fig21G4D1_D_patient_A_D_geometry",
  15, 7.5
)

cat("Figures A-D generated\n")

# ============================================================
# E — CROSS-MODAL EFFECT LANDSCAPE
# ============================================================

effect_long <- land %>%
  mutate(
    signal_label=paste0(
      clean_receiver(receiver),
      " | ",
      gsub("_"," ",communication_program)
    ),
    Communication=communication_mean_difference,
    Density=density_mean_difference,
    Abundance=DA_log2_OR,
    Transcription=log10(DE_n_FDR005 + 1)
  ) %>%
  select(
    signal_id, contrast, signal_label,
    Communication, Density,
    Abundance, Transcription
  ) %>%
  pivot_longer(
    cols=c(
      Communication,
      Density,
      Abundance,
      Transcription
    ),
    names_to="evidence_domain",
    values_to="raw_value"
  ) %>%
  group_by(evidence_domain) %>%
  mutate(
    visual_scaled_value=
      as.numeric(scale(raw_value))
  ) %>%
  ungroup()

signal_order <- land %>%
  arrange(
    factor(
      contrast,
      levels=c("AA_vs_HC","FFA_vs_HC","AA_vs_FFA")
    ),
    communication_mean_difference
  ) %>%
  mutate(
    signal_label=paste0(
      clean_receiver(receiver),
      " | ",
      gsub("_"," ",communication_program)
    )
  ) %>%
  pull(signal_label)

effect_long <- effect_long %>%
  mutate(
    signal_label=factor(
      signal_label,
      levels=unique(signal_order)
    ),
    evidence_domain=factor(
      evidence_domain,
      levels=c(
        "Communication",
        "Density",
        "Abundance",
        "Transcription"
      )
    )
  )

pE <- ggplot(
  effect_long,
  aes(
    x=evidence_domain,
    y=signal_label,
    fill=visual_scaled_value
  )
) +
  geom_tile(
    color="white",
    linewidth=.8
  ) +
  geom_point(
    aes(size=abs(visual_scaled_value)),
    shape=21,
    color="#333333",
    stroke=.25
  ) +
  scale_fill_gradient2(
    low="#3B6FB6",
    mid="white",
    high="#B6383A",
    midpoint=0,
    name="Within-domain\nscaled value"
  ) +
  scale_size_continuous(
    range=c(.8,5.8),
    guide="none"
  ) +
  facet_grid(
    contrast ~ .,
    scales="free_y",
    space="free_y"
  ) +
  labs(
    title="Cross-modal architecture of communication-associated disease states",
    subtitle="Communication, sender density, abundance and transcription shown as independent evidence domains",
    x=NULL,
    y=NULL
  ) +
  theme_pub(10.5) +
  theme(
    axis.text.x=element_text(angle=25,hjust=1),
    legend.position="right",
    panel.spacing.y=grid::unit(.7,"lines")
  )

save_plot(
  pE,
  "Fig21G4D1_E_crossmodal_effect_landscape",
  11.8, 9.3
)

# ============================================================
# F — PATIENT COMMUNICATION FINGERPRINT
# ============================================================

finger <- bind_rows(aa,ffa) %>%
  distinct(
    signal_id,
    contrast,
    patient_id,
    disease,
    receiver,
    communication_program,
    mean_incoming_A_axis
  ) %>%
  group_by(signal_id) %>%
  mutate(
    A_visual_z=
      as.numeric(scale(mean_incoming_A_axis))
  ) %>%
  ungroup() %>%
  mutate(
    signal_label=paste0(
      clean_receiver(receiver),
      " | ",
      gsub("_"," ",communication_program)
    )
  )

patient_order <- c(
  sort(unique(finger$patient_id[finger$disease=="HC"])),
  sort(unique(finger$patient_id[finger$disease=="AA"])),
  sort(unique(finger$patient_id[finger$disease=="FFA"]))
)

finger <- finger %>%
  mutate(
    patient_id=factor(
      patient_id,
      levels=patient_order
    )
  )

pF <- ggplot(
  finger,
  aes(
    x=patient_id,
    y=signal_label,
    fill=A_visual_z
  )
) +
  geom_tile(
    color="white",
    linewidth=.65
  ) +
  scale_fill_gradient2(
    low="#3B6FB6",
    mid="white",
    high="#B6383A",
    midpoint=0,
    na.value="grey92",
    name="Within-signal\nA-axis z-score"
  ) +
  facet_grid(
    contrast ~ .,
    scales="free_y",
    space="free_y"
  ) +
  labs(
    title="Patient-resolved communication fingerprints",
    subtitle="Disease-associated communication heterogeneity without clustering or imputation",
    x="Patient",
    y=NULL
  ) +
  theme_pub(10.5) +
  theme(
    axis.text.x=element_text(angle=45,hjust=1),
    legend.position="right",
    panel.spacing.y=grid::unit(.7,"lines")
  )

save_plot(
  pF,
  "Fig21G4D1_F_patient_communication_fingerprint",
  13, 9.4
)

# ============================================================
# CONTACT SHEET — A/B/C only for rapid visual adjudication
# ============================================================

contact <- (
  pA /
  pB /
  pC
) +
  plot_annotation(
    title="21G4D1 quantitative communication discovery"
  )

save_plot(
  contact,
  "Fig21G4D1_contact_A_B_C",
  14, 15
)

generated <- list.files(
  FIG,
  pattern="\\.(pdf|png)$",
  full.names=FALSE
)

cat("\n===== GENERATED FIGURES =====\n")
cat(paste(generated, collapse="\n"), "\n")

stopifnot(
  length(generated) == 14
)

cat("\n============================================================\n")
cat("21G4D1 COMPLETE\n")
cat("============================================================\n")
cat("Six quantitative candidate figures generated\n")
cat("One contact sheet generated\n")
cat("PDF + 500-dpi PNG generated for each\n")
cat("No signal selection changed\n")
cat("No statistical inference\n")
cat("No analytical refit\n")
cat("No manuscript modification\n")
cat("PASS\n")
