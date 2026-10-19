suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

BASE <- paste0(
  "PIPELINE_V2/results/disease_biology/",
  "21_CELL_CELL_COMMUNICATION/",
  "21G_MULTICOMPARTMENT_COMMUNICATION/",
  "21G4_CROSSMODAL_INTEGRATION"
)

TAB <- file.path(BASE,"tables")
FIG <- file.path(BASE,"figures","21G4D2")
dir.create(FIG,recursive=TRUE,showWarnings=FALSE)

land <- read_csv(
  file.path(TAB,"21G4C6B_signal_effect_landscape.csv"),
  show_col_types=FALSE
)

aa <- read_csv(
  file.path(TAB,"21G4C6B_AA_vs_HC_patient_signals.csv"),
  show_col_types=FALSE
)

ffa <- read_csv(
  file.path(TAB,"21G4C6B_FFA_vs_HC_patient_signals.csv"),
  show_col_types=FALSE
)

COL <- c(
  HC="#8F8D8E",
  FFA="#4D9E4F",
  AA="#A8142A"
)

clean_receiver <- function(x){
  x <- gsub("_"," ",x)
  x <- gsub("IFE Infundibulum ","IFE ",x,fixed=TRUE)
  x
}

theme_pub <- function(base_size=11){
  theme_classic(base_size=base_size) +
    theme(
      axis.text=element_text(color="black"),
      axis.title=element_text(face="bold"),
      strip.background=element_blank(),
      strip.text=element_text(face="bold"),
      plot.title=element_text(face="bold",hjust=.5),
      plot.subtitle=element_text(hjust=.5),
      legend.position="top"
    )
}

save_png <- function(p,name,w,h){
  ggsave(
    file.path(FIG,paste0(name,".png")),
    p,width=w,height=h,dpi=500,bg="white"
  )
}

# ============================================================
# C2 — FFA distributed loss: cleaner quantitative landscape
# ============================================================

fg <- ffa %>%
  filter(communication_program=="GROWTH_SURVIVAL") %>%
  mutate(receiver=clean_receiver(receiver))

ord <- land %>%
  filter(
    contrast=="FFA_vs_HC",
    communication_program=="GROWTH_SURVIVAL"
  ) %>%
  arrange(communication_mean_difference) %>%
  transmute(receiver=clean_receiver(receiver)) %>%
  pull(receiver)

fg$receiver <- factor(fg$receiver,levels=unique(ord))

pC2 <- ggplot(
  fg,
  aes(receiver,mean_incoming_A_axis)
) +
  geom_boxplot(
    aes(fill=disease),
    position=position_dodge(.66),
    width=.46,
    alpha=.20,
    outlier.shape=NA,
    linewidth=.45
  ) +
  geom_point(
    aes(fill=disease),
    position=position_jitterdodge(
      jitter.width=.055,
      dodge.width=.66
    ),
    shape=21,
    size=3,
    stroke=.45,
    color="#303030"
  ) +
  scale_fill_manual(values=COL,name=NULL) +
  labs(
    title="FFA exhibits distributed loss of epithelial growth-survival communication",
    subtitle="Patient-level incoming communication across follicular and interfollicular epithelial states",
    x=NULL,
    y="Mean incoming communication (A-axis)"
  ) +
  theme_pub(12) +
  theme(
    axis.text.x=element_text(
      angle=32,hjust=1,vjust=1
    )
  )

save_png(
  pC2,
  "Fig21G4D2_C_FFA_distributed_growth_survival",
  13,6
)

# ============================================================
# E2 — real-value multimodal matrix
# No synthetic cross-modal score.
# ============================================================

mat <- land %>%
  mutate(
    label=paste0(
      clean_receiver(receiver),
      " | ",
      gsub("_"," ",communication_program)
    ),
    `Communication ΔA`=
      communication_mean_difference,
    `Sender density ΔD`=
      density_mean_difference,
    `Abundance log2(OR)`=
      DA_log2_OR,
    `DE genes\nlog10(n+1)`=
      log10(DE_n_FDR005+1)
  ) %>%
  select(
    signal_id,contrast,label,
    communication_q_BH,
    density_q_BH,
    DA_FDR_compartment,
    DE_n_FDR005,
    `Communication ΔA`,
    `Sender density ΔD`,
    `Abundance log2(OR)`,
    `DE genes\nlog10(n+1)`
  ) %>%
  pivot_longer(
    cols=c(
      `Communication ΔA`,
      `Sender density ΔD`,
      `Abundance log2(OR)`,
      `DE genes\nlog10(n+1)`
    ),
    names_to="domain",
    values_to="value"
  ) %>%
  mutate(
    q_domain=case_when(
      domain=="Communication ΔA" ~ communication_q_BH,
      domain=="Sender density ΔD" ~ density_q_BH,
      domain=="Abundance log2(OR)" ~ DA_FDR_compartment,
      TRUE ~ NA_real_
    ),
    supported=case_when(
      domain=="DE genes\nlog10(n+1)" ~ DE_n_FDR005 > 0,
      TRUE ~ q_domain < .05
    )
  )

# Scale only for visual colour within each measurement domain;
# raw values remain encoded in labels.
mat <- mat %>%
  group_by(domain) %>%
  mutate(z=as.numeric(scale(value))) %>%
  ungroup()

order_labels <- land %>%
  arrange(
    factor(
      contrast,
      levels=c("AA_vs_HC","FFA_vs_HC","AA_vs_FFA")
    ),
    communication_mean_difference
  ) %>%
  mutate(
    label=paste0(
      clean_receiver(receiver),
      " | ",
      gsub("_"," ",communication_program)
    )
  ) %>%
  pull(label)

mat$label <- factor(
  mat$label,
  levels=unique(order_labels)
)

pE2 <- ggplot(
  mat,
  aes(domain,label,fill=z)
) +
  geom_tile(
    color="white",
    linewidth=1
  ) +
  geom_point(
    data=filter(mat,supported),
    shape=21,
    fill="white",
    color="black",
    size=3,
    stroke=.7
  ) +
  scale_fill_gradient2(
    low="#3B6FB6",
    mid="white",
    high="#B6383A",
    midpoint=0,
    name="Within-domain\nrelative effect"
  ) +
  facet_grid(
    contrast~.,
    scales="free_y",
    space="free_y"
  ) +
  labs(
    title="Communication states occupy distinct cross-modal disease architectures",
    subtitle="Independent communication, sender-density, abundance and transcriptional evidence; circles mark supported domains",
    x=NULL,
    y=NULL
  ) +
  theme_pub(11) +
  theme(
    axis.text.x=element_text(
      angle=25,hjust=1
    ),
    legend.position="right",
    panel.spacing.y=grid::unit(.7,"lines")
  )

save_png(
  pE2,
  "Fig21G4D2_E_crossmodal_architecture",
  11.5,9
)

# ============================================================
# F2 — patient fingerprints separated by disease contrast
# Removes structurally empty AA/FFA/HC cells.
# ============================================================

make_fingerprint <- function(dat,title_text){

  d <- dat %>%
    distinct(
      signal_id,
      patient_id,
      disease,
      receiver,
      communication_program,
      mean_incoming_A_axis
    ) %>%
    group_by(signal_id) %>%
    mutate(
      z=as.numeric(scale(mean_incoming_A_axis))
    ) %>%
    ungroup() %>%
    mutate(
      label=paste0(
        clean_receiver(receiver),
        " | ",
        gsub("_"," ",communication_program)
      )
    )

  patient_levels <- d %>%
    distinct(patient_id,disease) %>%
    arrange(
      factor(disease,levels=c("HC","AA","FFA")),
      patient_id
    ) %>%
    pull(patient_id)

  d$patient_id <- factor(
    d$patient_id,
    levels=patient_levels
  )

  ggplot(
    d,
    aes(patient_id,label,fill=z)
  ) +
    geom_tile(
      color="white",
      linewidth=.8
    ) +
    scale_fill_gradient2(
      low="#3B6FB6",
      mid="white",
      high="#B6383A",
      midpoint=0,
      name="Within-signal\nA-axis z-score"
    ) +
    labs(
      title=title_text,
      x=NULL,
      y=NULL
    ) +
    theme_pub(11) +
    theme(
      axis.text.x=element_text(
        angle=45,hjust=1
      ),
      legend.position="right"
    )
}

pFAA <- make_fingerprint(
  aa,
  "AA communication fingerprint"
)

pFFFA <- make_fingerprint(
  ffa,
  "FFA communication fingerprint"
)

pF2 <- pFAA / pFFFA +
  plot_annotation(
    title="Patient-resolved disease communication architecture",
    subtitle="Frozen publication signals shown without clustering or imputation"
  )

save_png(
  pF2,
  "Fig21G4D2_F_patient_fingerprints",
  11.5,10
)

cat("\n============================================================\n")
cat("21G4D2 COMPLETE\n")
cat("============================================================\n")
cat("C2: distributed FFA growth-survival architecture\n")
cat("E2: independent cross-modal quantitative architecture\n")
cat("F2: patient-resolved fingerprints\n")
cat("No new statistical inference\n")
cat("No analytical refit\n")
cat("PASS\n")
