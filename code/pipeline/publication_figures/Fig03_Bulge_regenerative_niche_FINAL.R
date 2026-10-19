# ============================================================
# 22F4F6C — INTEGRATED BULGE TRAJECTORY
# Publication visualization only
# No trajectory reinference
# No new differential-expression testing
# ============================================================

suppressPackageStartupMessages({
  library(monocle3)
  library(SingleCellExperiment)
  library(igraph)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(svglite)
})

cat("\n============================================================\n")
cat("22F4F6C — INTEGRATED BULGE TRAJECTORY\n")
cat("============================================================\n")

# ============================================================
# 1. PATHS
# ============================================================

BASE <- file.path(
  "PIPELINE_V2",
  "results",
  "disease_biology",
  "22_MULTICOMPARTMENT_DISEASE_AXES",
  "22F_BULGE_STEM_NICHE_INTEGRITY",
  "22F4_BULGE_MONOCLE3_TRAJECTORIES"
)

A_FILE <- file.path(
  BASE,
  "22F4A_UNROOTED_TOPOLOGY",
  "objects",
  "22F4A_Bulge_Monocle3_unrooted.rds"
)

C_FILE <- file.path(
  BASE,
  "22F4C_TRUNK_PSEUDOTIME",
  "objects",
  "22F4C_trunk_definition_freeze.rds"
)

D2_FILE <- file.path(
  BASE,
  "22F4D_DYNAMIC_EXPRESSION",
  "22F4D2_DYNAMIC_GENE_MODELLING",
  "objects",
  "22F4D2_dynamic_gene_modelling_freeze.rds"
)

D3_FILE <- file.path(
  BASE,
  "22F4D_DYNAMIC_EXPRESSION",
  "22F4D3B_K5_FREEZE",
  "objects",
  "22F4D3B_K5_dynamic_modules_freeze.rds"
)

OUT <- file.path(
  "PIPELINE_V2",
  "results",
  "publication_figures",
  "Figure03_Bulge_regenerative_niche"
)

SVG_DIR <- file.path(OUT, "figures", "svg")
PDF_DIR <- file.path(OUT, "figures", "pdf")
PNG_DIR <- file.path(OUT, "figures", "png")
TAB_DIR <- file.path(OUT, "tables")
OBJ_DIR <- file.path(OUT, "objects")
AUDIT_DIR <- file.path(OUT, "audit")

for (d in c(
  SVG_DIR,
  PDF_DIR,
  PNG_DIR,
  TAB_DIR,
  OBJ_DIR,
  AUDIT_DIR
)) {
  dir.create(
    d,
    recursive=TRUE,
    showWarnings=FALSE
  )
}

stopifnot(
  file.exists(A_FILE),
  file.exists(C_FILE),
  file.exists(D2_FILE),
  file.exists(D3_FILE)
)


# ============================================================
# FIGURE 3 — ADDITIONAL FROZEN INPUTS
# ============================================================

E4A_FILE <- file.path(
  BASE,
  "22F4E_DISEASE_TRAJECTORY",
  "22F4E4_HALLMARK_TRAJECTORY",
  "22F4E4A_SCORE_INTERACTION",
  "objects",
  "22F4E4A_Hallmark_trajectory_interaction_freeze.rds"
)

E4B_FILE <- file.path(
  BASE,
  "22F4E_DISEASE_TRAJECTORY",
  "22F4E4_HALLMARK_TRAJECTORY",
  "22F4E4B_EXTERNAL_HFSC_SIGNATURE",
  "objects",
  "22F4E4B_external_HFSC_signature_freeze.rds"
)

stopifnot(
  file.exists(E4A_FILE),
  file.exists(E4B_FILE)
)

# ============================================================
# 2. LOAD FROZEN INPUTS
# ============================================================

cat("\n===== LOAD FROZEN INPUTS =====\n")

cds <- readRDS(A_FILE)
C <- readRDS(C_FILE)
D2 <- readRDS(D2_FILE)
D3 <- readRDS(D3_FILE)
E4A <- readRDS(E4A_FILE)
E4B <- readRDS(E4B_FILE)

cat("Bulge cells        :", ncol(cds), "\n")
cat("Root node          :", C$root_node, "\n")
cat("End node           :", C$end_node, "\n")
cat("Strict trunk cells :", C$n_strict_trunk, "\n")
cat("Dynamic genes      :", D3$n_dynamic_genes, "\n")
cat("Frozen K           :", D3$K_final, "\n")

stopifnot(
  ncol(cds) == 5420,
  C$root_node == "Y_1",
  C$end_node == "Y_12",
  D3$n_dynamic_genes == 2774,
  D3$K_final == 5
)

# ============================================================
# 3. CELL-LEVEL UMAP DATA
# ============================================================

cat("\n===== BUILD CELL DATA =====\n")

umap <- reducedDims(cds)$UMAP
meta <- as.data.frame(colData(cds))

stopifnot(
  nrow(umap) == nrow(meta),
  identical(
    rownames(umap),
    rownames(meta)
  )
)

cell_df <- tibble(
  cell_id=rownames(umap),
  UMAP_1=umap[,1],
  UMAP_2=umap[,2],
  patient_id=as.character(meta$patient_id),
  disease=as.character(meta$disease),
  territory_K3=as.character(meta$territory_K3)
) |>
  mutate(
    disease=factor(
      disease,
      levels=c("HC","FFA","AA")
    ),
    territory_K3=factor(
      territory_K3,
      levels=c("T3","T2","T1")
    )
  )

cat("Cell rows:", nrow(cell_df), "\n")

# ============================================================
# 4. PRINCIPAL GRAPH COORDINATES
# ============================================================

cat("\n===== RECOVER PRINCIPAL GRAPH =====\n")

g <- principal_graph(cds)[["UMAP"]]

stopifnot(
  igraph::vcount(g) == 40,
  igraph::ecount(g) == 39
)

aux <- principal_graph_aux(cds)[["UMAP"]]

node_coord <- aux$dp_mst

if (is.null(node_coord)) {
  stop("dp_mst not available in principal_graph_aux.")
}

node_coord <- as.matrix(node_coord)

cat(
  "dp_mst dimensions:",
  paste(dim(node_coord), collapse=" x "),
  "\n"
)

# Monocle3 commonly stores dp_mst as 2 x N.
# Standardize to N x 2.

if (
  nrow(node_coord) == 2 &&
  ncol(node_coord) == igraph::vcount(g)
) {

  node_coord <- t(node_coord)

} else if (
  nrow(node_coord) == igraph::vcount(g) &&
  ncol(node_coord) >= 2
) {

  node_coord <- node_coord[,1:2,drop=FALSE]

} else {

  stop(
    paste(
      "Unexpected dp_mst dimensions:",
      paste(dim(node_coord), collapse=" x ")
    )
  )
}

node_names <- igraph::V(g)$name

# Recover names robustly

if (
  !is.null(rownames(node_coord)) &&
  all(node_names %in% rownames(node_coord))
) {

  node_coord <- node_coord[
    node_names,
    ,
    drop=FALSE
  ]

} else {

  rownames(node_coord) <- node_names
}

node_df <- tibble(
  principal_node=node_names,
  UMAP_1=node_coord[,1],
  UMAP_2=node_coord[,2]
)

edges_raw <- igraph::as_data_frame(
  g,
  what="edges"
)

edge_df <- edges_raw |>
  transmute(
    from=as.character(from),
    to=as.character(to)
  ) |>
  left_join(
    node_df |>
      rename(
        from=principal_node,
        x=UMAP_1,
        y=UMAP_2
      ),
    by="from"
  ) |>
  left_join(
    node_df |>
      rename(
        to=principal_node,
        xend=UMAP_1,
        yend=UMAP_2
      ),
    by="to"
  )

stopifnot(
  nrow(node_df) == 40,
  nrow(edge_df) == 39,
  all(complete.cases(
    edge_df[,c("x","y","xend","yend")]
  ))
)

# ============================================================
# 5. TRUNK VERSUS BRANCH EDGES
# ============================================================

cat("\n===== DEFINE FROZEN TRUNK GRAPH =====\n")

trunk_nodes <- as.character(C$trunk_nodes)

edge_df <- edge_df |>
  mutate(
    on_trunk=
      from %in% trunk_nodes &
      to %in% trunk_nodes
  )

node_df <- node_df |>
  mutate(
    on_trunk=principal_node %in% trunk_nodes,
    special_node=case_when(
      principal_node == C$root_node ~ "Root",
      principal_node == C$end_node ~ "Trunk end",
      principal_node == "Y_37" ~ "Y37 branch",
      TRUE ~ NA_character_
    )
  )

cat(
  "Trunk nodes:",
  sum(node_df$on_trunk),
  "\n"
)

cat(
  "Trunk edges:",
  sum(edge_df$on_trunk),
  "\n"
)

# ============================================================
# 6. CELL -> CLOSEST PRINCIPAL NODE
# ============================================================

cat("\n===== MAP CELLS TO PRINCIPAL NODES =====\n")

closest <- aux$pr_graph_cell_proj_closest_vertex

closest <- as.matrix(closest)

if (nrow(closest) != nrow(cell_df)) {
  stop(
    paste(
      "Closest-vertex rows do not match cells:",
      nrow(closest),
      "vs",
      nrow(cell_df)
    )
  )
}

closest_idx <- as.integer(
  closest[,1]
)

# Monocle3 may store either graph indices or node labels.
# Handle both cases.

if (
  all(
    !is.na(closest_idx) &
    closest_idx >= 1 &
    closest_idx <= length(node_names)
  )
) {

  closest_node <- node_names[closest_idx]

} else {

  closest_node <- as.character(
    closest[,1]
  )
}

cell_df$principal_node <- closest_node

# ============================================================
# 7. FROZEN TRAJECTORY POSITION
# ============================================================

cat("\n===== MAP FROZEN TRAJECTORY POSITION =====\n")

node_position <- C$node_profile |>
  select(
    principal_node,
    trajectory_position
  )

cell_df <- cell_df |>
  left_join(
    node_position,
    by="principal_node"
  )

cell_df <- cell_df |>
  mutate(
    on_strict_trunk=
      principal_node %in% trunk_nodes &
      !is.na(trajectory_position)
  )

cat(
  "Cells assigned to frozen trunk position:",
  sum(cell_df$on_strict_trunk),
  "\n"
)

# ============================================================
# 8. PALETTES
# ============================================================

disease_cols <- c(
  HC="#8F8D8E",
  FFA="#4D9E4F",
  AA="#A8142A"
)

state_cols <- c(
  T3="#2F7189",
  T2="#4B4C7A",
  T1="#B36B00"
)

module_cols <- c(
  M1="#4E79A7",
  M2="#59A14F",
  M3="#E15759",
  M4="#F28E2B",
  M5="#B07AA1"
)


# ============================================================
# 9. CANONICAL VISUAL SYSTEM
# ============================================================

K3_COL <- c(
  T3="#2F7189",
  T2="#4B4C7A",
  T1="#B36B00"
)

MODULE_COL <- c(
  M1="#4E79A7",
  M2="#59A14F",
  M3="#E15759",
  M4="#F28E2B",
  M5="#B07AA1"
)

# ============================================================
# 10. SHARED GRAPH OBJECTS
# ============================================================

graph_nontrunk <- edge_df |>
  filter(!on_trunk)

graph_trunk <- edge_df |>
  filter(on_trunk)

special_labels <- node_df |>
  filter(
    principal_node %in% c(
      "Y_1",
      "Y_12",
      "Y_37"
    )
  )

# ============================================================
# 11. GLOBAL K3 LABEL POSITIONS
# ============================================================

# Use medians only as starting coordinates, followed by small
# deterministic offsets for publication placement.

state_label_df <- cell_df |>
  filter(
    !is.na(territory_K3)
  ) |>
  group_by(
    territory_K3
  ) |>
  summarise(
    UMAP_1=median(UMAP_1),
    UMAP_2=median(UMAP_2),
    .groups="drop"
  ) |>
  mutate(
    label=case_when(
      territory_K3=="T3" ~ "T3\nStem-like",
      territory_K3=="T2" ~ "T2\nActivated/remodeling",
      territory_K3=="T1" ~ "T1\nProliferative",
      TRUE ~ as.character(territory_K3)
    ),

    # Minor deterministic offsets to reduce graph overlap.
    UMAP_1=case_when(
      territory_K3=="T3" ~ UMAP_1 - 0.45,
      territory_K3=="T2" ~ UMAP_1 + 0.45,
      territory_K3=="T1" ~ UMAP_1 + 0.35,
      TRUE ~ UMAP_1
    ),

    UMAP_2=case_when(
      territory_K3=="T3" ~ UMAP_2 + 0.40,
      territory_K3=="T2" ~ UMAP_2 + 0.30,
      territory_K3=="T1" ~ UMAP_2 - 0.45,
      TRUE ~ UMAP_2
    )
  )

cat("\n===== K3 LABEL POSITIONS =====\n")
print(state_label_df)

# ============================================================
# 12. PANEL A — GLOBAL TRAJECTORY
# ============================================================

cat("\n===== PANEL A — GLOBAL TRAJECTORY =====\n")

set.seed(220461)

plot_cells_A <- cell_df |>
  slice_sample(
    prop=1
  )

pA <- ggplot() +

  geom_point(
    data=plot_cells_A,
    aes(
      x=UMAP_1,
      y=UMAP_2,
      fill=territory_K3
    ),
    shape=21,
    size=1.25,
    stroke=0.16,
    colour="#3A3A3A",
    alpha=0.95
  ) +

  geom_segment(
    data=graph_nontrunk,
    aes(
      x=x,
      y=y,
      xend=xend,
      yend=yend
    ),
    colour="#4D4D4D",
    linewidth=0.40,
    alpha=0.55,
    lineend="round"
  ) +

  geom_segment(
    data=graph_trunk,
    aes(
      x=x,
      y=y,
      xend=xend,
      yend=yend
    ),
    colour="#111111",
    linewidth=0.72,
    alpha=0.88,
    lineend="round"
  ) +

  geom_label(
    data=state_label_df,
    aes(
      x=UMAP_1,
      y=UMAP_2,
      label=label,
      colour=territory_K3
    ),
    fill=scales::alpha(
      "white",
      0.86
    ),
    label.size=0.18,
    size=3.0,
    fontface="bold",
    lineheight=0.92,
    show.legend=FALSE
  ) +

  geom_text(
    data=special_labels,
    aes(
      x=UMAP_1,
      y=UMAP_2,
      label=principal_node
    ),
    colour="#111111",
    size=2.15,
    fontface="bold",
    nudge_y=0.40,
    check_overlap=TRUE
  ) +

  scale_fill_manual(
    values=K3_COL,
    breaks=c(
      "T3",
      "T2",
      "T1"
    ),
    drop=FALSE
  ) +

  scale_colour_manual(
    values=K3_COL,
    guide="none"
  ) +

  coord_equal() +

  labs(
    title="A  Global Bulge/Isthmus trajectory",
    subtitle="Frozen Monocle3 topology across canonical functional territories",
    x="UMAP 1",
    y="UMAP 2",
    fill="State"
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
      size=8.5
    ),
    axis.title=element_text(
      size=8.5
    ),
    axis.text=element_text(
      size=7.5
    ),
    legend.title=element_text(
      size=8
    ),
    legend.text=element_text(
      size=7
    ),
    legend.position="right"
  )

# ============================================================
# 13. PANEL B — DISEASE SPLIT
# ============================================================

cat("\n===== PANEL B — HC / FFA / AA =====\n")

set.seed(220462)

plot_cells_B <- cell_df |>
  slice_sample(
    prop=1
  )

pB <- ggplot() +

  geom_point(
    data=plot_cells_B,
    aes(
      x=UMAP_1,
      y=UMAP_2,
      fill=territory_K3
    ),
    shape=21,
    size=1.18,
    stroke=0.15,
    colour="#3A3A3A",
    alpha=0.95
  ) +

  geom_segment(
    data=graph_nontrunk,
    aes(
      x=x,
      y=y,
      xend=xend,
      yend=yend
    ),
    colour="#555555",
    linewidth=0.34,
    alpha=0.50
  ) +

  geom_segment(
    data=graph_trunk,
    aes(
      x=x,
      y=y,
      xend=xend,
      yend=yend
    ),
    colour="#111111",
    linewidth=0.55,
    alpha=0.75
  ) +

  facet_wrap(
    ~disease,
    nrow=1,
    drop=FALSE
  ) +

  scale_fill_manual(
    values=K3_COL,
    breaks=c(
      "T3",
      "T2",
      "T1"
    ),
    drop=FALSE
  ) +

  coord_equal(
    xlim=c(-6.5,1.5),
    ylim=c(8.0,15.5),
    expand=FALSE,
    clip="on"
  ) +

  labs(
    title="B  Disease-specific occupancy of the shared Bulge/Isthmus landscape",
    subtitle="Common UMAP window, K3 colours and frozen principal graph across HC, FFA and AA",
    x="UMAP 1",
    y="UMAP 2",
    fill="State"
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
      size=8.2
    ),
    strip.background=element_blank(),
    strip.text=element_text(
      face="bold",
      size=10
    ),
    legend.position="right"
  )


# ============================================================
# 14. PANEL C — REPRESENTATIVE DYNAMIC GENES
# ============================================================

cat("\n===== PANEL C — DYNAMIC GENES =====\n")

# 10 previously frozen representative genes/module.
# This is a display subset only.
# No new gene ranking or module assignment.

rep50 <- D3$representative_genes |>
  arrange(
    module_order,
    order_within_module
  ) |>
  group_by(
    module
  ) |>
  slice_head(
    n=10
  ) |>
  ungroup()

stopifnot(
  nrow(rep50) == 50,
  all(
    rep50$gene %in%
      rownames(D3$row_scaled_matrix)
  )
)

cat("\nDisplayed genes/module:\n")

print(
  rep50 |>
    count(module),
  n=Inf
)

# ============================================================
# 15. HEATMAP MATRIX
# ============================================================

heat_mat <- D3$row_scaled_matrix[
  rep50$gene,
  ,
  drop=FALSE
]

heat_long <- as.data.frame(
  heat_mat
) |>
  rownames_to_column(
    "gene"
  ) |>
  pivot_longer(
    cols=-gene,
    names_to="trajectory_bin",
    values_to="z"
  ) |>
  left_join(
    rep50 |>
      dplyr::select(
        gene,
        module,
        module_order,
        order_within_module
      ),
    by="gene"
  ) |>
  mutate(
    trajectory_index=
      as.integer(
        sub(
          "^T",
          "",
          trajectory_bin
        )
      ),

    trajectory_position=
      (trajectory_index-1)/100
  )

# Preserve M1 -> M5 from TOP to BOTTOM.
# Within each module preserve frozen representative ordering.

gene_order <- rep50 |>
  arrange(
    module_order,
    order_within_module
  ) |>
  pull(gene)

heat_long$gene <- factor(
  heat_long$gene,
  levels=rev(gene_order)
)

# ============================================================
# 16. MODULE GEOMETRY FOR DISPLAY
# ============================================================

gene_display <- rep50 |>
  arrange(
    module_order,
    order_within_module
  ) |>
  mutate(
    row_from_top=row_number(),
    y_numeric=
      51-row_from_top
  )

module_geometry <- gene_display |>
  group_by(
    module,
    module_order
  ) |>
  summarise(
    y_min=min(y_numeric)-0.5,
    y_max=max(y_numeric)+0.5,
    y_mid=mean(y_numeric),
    n_display=n(),
    .groups="drop"
  ) |>
  arrange(
    module_order
  )

cat("\n===== DISPLAY MODULE GEOMETRY =====\n")

print(
  module_geometry,
  n=Inf
)

# Boundaries between M1-M5.
module_boundaries <- module_geometry |>
  arrange(module_order) |>
  slice_head(
    n=4
  ) |>
  pull(y_min)

# ============================================================
# 17. FUNCTIONAL ANNOTATIONS
# Derived from 22F4F6A ORA
# ============================================================

module_annotation <- tibble(
  module=c(
    "M1",
    "M2",
    "M3",
    "M4",
    "M5"
  ),

  functional_title=c(
    "Early / transcriptional adaptation",
    "Interferon / antigen presentation",
    "Tissue remodeling / EMT",
    "Proliferative / cell cycle",
    "Late trajectory program"
  ),

  pathways=c(
    paste(
      "Transcriptional regulation",
      "Hypoxia [FDR 0.061]",
      sep="\n"
    ),

    paste(
      "IFN-gamma response",
      "IFN-alpha response",
      "Antigen processing/presentation",
      sep="\n"
    ),

    paste(
      "Epithelial-mesenchymal transition",
      "Wound healing",
      "Vasculature development",
      sep="\n"
    ),

    paste(
      "E2F targets",
      "G2M checkpoint",
      "Mitotic spindle",
      sep="\n"
    ),

    paste(
      "UV response DN",
      "No robust GO:BP enrichment",
      sep="\n"
    )
  ),

  evidence=c(
    "supportive",
    "strong",
    "strong",
    "very strong",
    "limited"
  )
) |>
  left_join(
    module_geometry |>
      dplyr::select(
        module,
        module_order,
        y_mid
      ),
    by="module"
  ) |>
  arrange(
    module_order
  )

write.csv(
  module_annotation,
  file.path(
    TAB_DIR,
    "22F4F6C_module_functional_annotation.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 18. STATE BAND
# ============================================================

state_ref <- tibble(
  state=factor(
    c(
      "T3",
      "T2",
      "T1"
    ),
    levels=c(
      "T3",
      "T2",
      "T1"
    )
  ),

  xmin=c(
    0,
    0.438,
    0.844
  ),

  xmax=c(
    0.438,
    0.844,
    1
  )
) |>
  mutate(
    midpoint=
      (xmin+xmax)/2
  )

p_state <- ggplot(
  state_ref
) +

  geom_rect(
    aes(
      xmin=xmin,
      xmax=xmax,
      ymin=0,
      ymax=1,
      fill=state
    ),
    colour="white",
    linewidth=0.6
  ) +

  geom_text(
    aes(
      x=midpoint,
      y=0.5,
      label=state
    ),
    fontface="bold",
    size=3.4
  ) +

  scale_fill_manual(
    values=K3_COL,
    guide="none"
  ) +

  scale_x_continuous(
    limits=c(0,1),
    expand=c(0,0)
  ) +

  scale_y_continuous(
    limits=c(0,1),
    expand=c(0,0)
  ) +

  theme_void()

# ============================================================
# 19. DYNAMIC-GENE HEATMAP
# ============================================================

p_heat <- ggplot(
  heat_long,
  aes(
    x=trajectory_position,
    y=gene,
    fill=z
  )
) +

  geom_tile() +

  # Frozen K=5 module bars along the left edge
  geom_segment(
    data=module_geometry,
    aes(
      x=-0.018,
      xend=-0.018,
      y=y_min,
      yend=y_max,
      colour=module
    ),
    inherit.aes=FALSE,
    linewidth=3.2,
    lineend="butt"
  ) +

  geom_text(
    data=module_geometry,
    aes(
      x=-0.030,
      y=y_mid,
      label=module,
      colour=module
    ),
    inherit.aes=FALSE,
    hjust=1,
    fontface="bold",
    size=2.7
  ) +

  geom_vline(
    xintercept=c(
      0.438,
      0.844
    ),
    colour="white",
    linewidth=0.80
  ) +

  geom_hline(
    yintercept=module_boundaries,
    colour="white",
    linewidth=1.0
  ) +

  scale_colour_manual(
    values=MODULE_COL,
    guide="none"
  ) +

  scale_fill_gradient2(
    low="#2166AC",
    mid="white",
    high="#B2182B",
    midpoint=0,
    limits=c(
      -2.5,
      2.5
    ),
    oob=scales::squish,
    name="Trajectory\nz-score"
  ) +

  scale_x_continuous(
    breaks=c(
      0,
      0.438,
      0.844,
      1
    ),
    labels=c(
      "0",
      "0.438",
      "0.844",
      "1"
    ),
    expand=c(0,0)
  ) +

  coord_cartesian(
    xlim=c(-0.055,1),
    expand=FALSE,
    clip="off"
  ) +

  labs(
    title=
      "C  Dynamic genes along the Bulge/Isthmus pseudotrajectory",
    subtitle=
      "50 representative genes from five frozen trajectory-shape modules",
    x="Normalized trajectory position",
    y=NULL
  ) +

  theme_classic(
    base_size=8
  ) +

  theme(
    plot.title=element_text(
      face="bold",
      size=10.5
    ),

    plot.subtitle=element_text(
      size=8.2
    ),

    axis.text.y=element_text(
      size=6.1,
      face="italic"
    ),

    axis.text.x=element_text(
      size=7.5
    ),

    axis.ticks.y=
      element_blank(),

    axis.line.y=
      element_blank(),

    legend.position="right",

    plot.margin=margin(
      4,
      4,
      5,
      4
    )
  )

p_heat_with_state <- p_state / p_heat +
  patchwork::plot_layout(
    heights=c(
      0.055,
      1
    )
  )

# ============================================================
# ============================================================
# 20. MODULE / PATHWAY ANNOTATION PANEL
# ============================================================

# Titles and pathways are deliberately separated vertically:
# title first, supporting pathways immediately below.

annotation_plot_df <- module_annotation |>
  mutate(
    title_label=paste0(
      module,
      "  ",
      functional_title
    ),
    pathway_y=y_mid - 1.15
  )

p_annotation <- ggplot(
  annotation_plot_df
) +

  geom_segment(
    aes(
      x=0,
      xend=0.07,
      y=y_mid,
      yend=y_mid,
      colour=module
    ),
    linewidth=5,
    lineend="butt"
  ) +

  geom_text(
    aes(
      x=0.10,
      y=y_mid,
      label=title_label,
      colour=module
    ),
    hjust=0,
    vjust=0.5,
    fontface="bold",
    size=2.9
  ) +

  geom_text(
    aes(
      x=0.10,
      y=pathway_y,
      label=pathways
    ),
    hjust=0,
    vjust=1,
    size=2.45,
    lineheight=0.98,
    colour="#444444"
  ) +

  scale_colour_manual(
    values=MODULE_COL,
    guide="none"
  ) +

  scale_x_continuous(
    limits=c(0,2.55),
    expand=c(0,0)
  ) +

  scale_y_continuous(
    limits=c(0.5,50.5),
    expand=c(0,0)
  ) +

  coord_cartesian(
    clip="off"
  ) +

  labs(
    title="Functional annotation"
  ) +

  theme_void() +

  theme(
    plot.title=element_text(
      face="bold",
      size=9.5,
      hjust=0
    ),
    plot.margin=margin(
      42,
      8,
      8,
      5
    )
  )

# 21. COMPLETE PANEL C
# ============================================================

pC <- patchwork::wrap_plots(
  list(
    p_heat_with_state,
    p_annotation
  ),
  nrow=1,
  widths=c(
    1.85,
    1
  )
)



# ============================================================
# FIGURE 3 PANEL D — DISEASE-DEPENDENT HALLMARK KINETICS
# Frozen 22F4E4A fitted curves; no new statistical testing.
# ============================================================

cat("\n===== FIGURE 3 PANEL D — HALLMARK KINETICS =====\n")

DISEASE_COL <- c(
  HC="#8A8A8A",
  FFA="#8FD694",
  AA="#D94B4B"
)

hallmark_labels <- c(
  HALLMARK_UV_RESPONSE_DN="UV response DN",
  HALLMARK_WNT_BETA_CATENIN_SIGNALING="WNT beta-catenin signaling",
  HALLMARK_HEDGEHOG_SIGNALING="Hedgehog signaling",
  HALLMARK_ANGIOGENESIS="Angiogenesis"
)

hallmark_order <- c(
  "HALLMARK_UV_RESPONSE_DN",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "HALLMARK_HEDGEHOG_SIGNALING",
  "HALLMARK_ANGIOGENESIS"
)

stopifnot(
  all(hallmark_order %in% rownames(E4A$HC_curve)),
  length(E4A$trajectory_grid)==101
)

hallmark_matrix_long <- function(mat, disease_name) {

  out <- as.data.frame(
    mat[
      hallmark_order,
      ,
      drop=FALSE
    ]
  )

  out$pathway <- rownames(out)

  out |>
    pivot_longer(
      cols=-pathway,
      names_to="grid_id",
      values_to="fitted_score"
    ) |>
    group_by(pathway) |>
    mutate(
      trajectory_position=
        E4A$trajectory_grid[
          seq_len(n())
        ]
    ) |>
    ungroup() |>
    mutate(
      disease=disease_name,
      pathway_label=
        unname(
          hallmark_labels[pathway]
        )
    )
}

hallmark_long <- bind_rows(
  hallmark_matrix_long(
    E4A$HC_curve,
    "HC"
  ),
  hallmark_matrix_long(
    E4A$FFA_curve,
    "FFA"
  ),
  hallmark_matrix_long(
    E4A$AA_curve,
    "AA"
  )
) |>
  mutate(
    disease=factor(
      disease,
      levels=c("HC","FFA","AA")
    ),
    pathway_label=factor(
      pathway_label,
      levels=unname(
        hallmark_labels[
          hallmark_order
        ]
      )
    )
  )

pD <- ggplot(
  hallmark_long,
  aes(
    x=trajectory_position,
    y=fitted_score,
    colour=disease
  )
) +
  geom_hline(
    yintercept=0,
    linewidth=0.25,
    colour="grey78"
  ) +
  geom_vline(
    xintercept=c(
      0.438,
      0.844
    ),
    linetype="dashed",
    linewidth=0.30,
    colour="grey60"
  ) +
  geom_line(
    linewidth=0.90,
    lineend="round"
  ) +
  facet_wrap(
    ~pathway_label,
    ncol=2,
    scales="free_y"
  ) +
  scale_colour_manual(
    values=DISEASE_COL,
    breaks=c("HC","FFA","AA"),
    drop=FALSE
  ) +
  scale_x_continuous(
    limits=c(0,1),
    breaks=c(
      0,
      0.438,
      0.844,
      1
    ),
    labels=c(
      "0",
      "0.438",
      "0.844",
      "1"
    ),
    expand=c(0,0)
  ) +
  labs(
    title="D  Disease-dependent pathway kinetics",
    subtitle="Hallmark programmes along the shared T3 -> T2 -> T1 trajectory",
    x="Normalized trajectory position",
    y="Fitted Hallmark score",
    colour=NULL
  ) +
  theme_classic(
    base_size=8.5
  ) +
  theme(
    plot.title=element_text(
      face="bold",
      size=10.5
    ),
    plot.subtitle=element_text(
      size=7.8
    ),
    strip.background=element_blank(),
    strip.text=element_text(
      face="bold",
      size=7.8
    ),
    axis.title=element_text(
      size=8
    ),
    axis.text=element_text(
      size=7
    ),
    legend.position="top",
    legend.justification="left",
    legend.text=element_text(
      size=7.2
    ),
    panel.spacing=grid::unit(
      0.55,
      "lines"
    )
  )

# ============================================================
# FIGURE 3 PANEL E — EXTERNAL HFSC PROGRAMMES
# Frozen 22F4E4B fitted curves; no new scoring/testing.
# ============================================================

cat("\n===== FIGURE 3 PANEL E — EXTERNAL HFSC PROGRAMMES =====\n")

hfsc_labels <- c(
  STEM_HFSC="HFSC identity",
  INFLAMMATORY_TRANSITION="Inflammatory transition",
  REMODELING_ASSOCIATED="Remodeling-associated VIM"
)

hfsc_order <- c(
  "STEM_HFSC",
  "INFLAMMATORY_TRANSITION",
  "REMODELING_ASSOCIATED"
)

hfsc_matrix_long <- function(mat, disease_name) {

  out <- as.data.frame(
    mat[
      hfsc_order,
      ,
      drop=FALSE
    ]
  )

  out$signature <- rownames(out)

  out |>
    pivot_longer(
      cols=-signature,
      names_to="grid_id",
      values_to="fitted_score"
    ) |>
    group_by(signature) |>
    mutate(
      trajectory_position=
        E4B$trajectory_grid[
          seq_len(n())
        ]
    ) |>
    ungroup() |>
    mutate(
      disease=disease_name,
      signature_label=
        unname(
          hfsc_labels[signature]
        )
    )
}

hfsc_long <- bind_rows(
  hfsc_matrix_long(
    E4B$HC_curve,
    "HC"
  ),
  hfsc_matrix_long(
    E4B$FFA_curve,
    "FFA"
  ),
  hfsc_matrix_long(
    E4B$AA_curve,
    "AA"
  )
) |>
  mutate(
    disease=factor(
      disease,
      levels=c("HC","FFA","AA")
    ),
    signature_label=factor(
      signature_label,
      levels=unname(
        hfsc_labels[
          hfsc_order
        ]
      )
    )
  )

pE <- ggplot(
  hfsc_long,
  aes(
    x=trajectory_position,
    y=fitted_score,
    colour=disease
  )
) +
  geom_hline(
    yintercept=0,
    linewidth=0.25,
    colour="grey78"
  ) +
  geom_vline(
    xintercept=c(
      0.438,
      0.844
    ),
    linetype="dashed",
    linewidth=0.30,
    colour="grey60"
  ) +
  geom_line(
    linewidth=0.95,
    lineend="round"
  ) +
  facet_wrap(
    ~signature_label,
    ncol=1,
    scales="free_y"
  ) +
  scale_colour_manual(
    values=DISEASE_COL,
    breaks=c("HC","FFA","AA"),
    drop=FALSE
  ) +
  scale_x_continuous(
    limits=c(0,1),
    breaks=c(
      0,
      0.438,
      0.844,
      1
    ),
    labels=c(
      "0",
      "0.438",
      "0.844",
      "1"
    ),
    expand=c(0,0)
  ) +
  labs(
    title="E  External HFSC programme projection",
    subtitle="A priori LPP-derived programmes projected onto the frozen trajectory",
    x="Normalized trajectory position",
    y="Fitted standardized score",
    colour=NULL
  ) +
  theme_classic(
    base_size=8.5
  ) +
  theme(
    plot.title=element_text(
      face="bold",
      size=10.5
    ),
    plot.subtitle=element_text(
      size=7.8
    ),
    strip.background=element_rect(
      fill="grey97",
      colour="grey82",
      linewidth=0.3
    ),
    strip.text=element_text(
      face="bold",
      size=7.8
    ),
    axis.title=element_text(
      size=8
    ),
    axis.text=element_text(
      size=7
    ),
    legend.position="top",
    legend.justification="left",
    legend.text=element_text(
      size=7.2
    ),
    panel.spacing=grid::unit(
      0.50,
      "lines"
    )
  )

# ============================================================
# 22. SAVE FUNCTION
# ============================================================

save_plot <- function(
  plot,
  stem,
  width,
  height
) {

  ggsave(
    file.path(
      SVG_DIR,
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
      PDF_DIR,
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
      PNG_DIR,
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
# 23. SAVE INDIVIDUAL PANELS
# ============================================================

cat("\n===== SAVE INDIVIDUAL PANELS =====\n")

save_plot(
  pA,
  "Fig22F4F6C_A_global_trajectory",
  width=7.8,
  height=6.3
)

save_plot(
  pB,
  "Fig22F4F6C_B_disease_split_trajectory",
  width=13.5,
  height=5.2
)

save_plot(
  pC,
  "Fig22F4F6C_C_dynamic_genes_functional_annotation",
  width=13.5,
  height=10.2
)


save_plot(
  pD,
  "Figure03_D_disease_pathway_kinetics",
  width=8.0,
  height=5.5
)

save_plot(
  pE,
  "Figure03_E_external_HFSC_programmes",
  width=6.8,
  height=7.3
)

# ============================================================
# 24. FINAL COMPOSITE
# ============================================================

cat("\n===== BUILD FINAL COMPOSITE =====\n")

# A is deliberately narrower than B/C.
# Spacer preserves a balanced publication layout.


# ============================================================
# FIGURE 3 — FINAL PUBLICATION POLISH
# ============================================================

# Remove explanatory subtitles; detailed interpretation belongs
# in the figure legend rather than inside the artwork.

pA <- pA +
  labs(
    title="A  Shared Bulge/Isthmus trajectory",
    subtitle=NULL
  ) +
  theme(
    plot.title=element_text(
      face="bold",
      size=11
    ),
    plot.subtitle=element_blank(),
    plot.margin=margin(
      4, 4, 4, 4
    )
  )

pB <- pB +
  labs(
    title="B  Disease-specific occupancy",
    subtitle=NULL
  ) +
  theme(
    plot.title=element_text(
      face="bold",
      size=11
    ),
    plot.subtitle=element_blank(),
    plot.margin=margin(
      4, 4, 4, 4
    )
  )

# Panel C retains its original integrated heatmap title.

pD <- pD +
  labs(
    title="D  Disease-dependent pathway kinetics",
    subtitle=NULL
  ) +
  theme(
    plot.title=element_text(
      face="bold",
      size=11
    ),
    plot.subtitle=element_blank(),
    legend.position="none",
    plot.margin=margin(
      4, 5, 4, 4
    )
  )

pE <- pE +
  labs(
    title="E  External HFSC programme projection",
    subtitle=NULL
  ) +
  theme(
    plot.title=element_text(
      face="bold",
      size=11
    ),
    plot.subtitle=element_blank(),
    legend.position="top",
    legend.justification="left",
    legend.box="horizontal",
    legend.margin=margin(
      0, 0, 2, 0
    ),
    plot.margin=margin(
      4, 4, 4, 5
    )
  )

row_AB <- patchwork::wrap_plots(
  list(
    pA,
    pB
  ),
  nrow=1,
  widths=c(
    0.42,
    0.58
  ),
  guides="keep"
)

row_DE <- patchwork::wrap_plots(
  list(
    pD,
    pE
  ),
  nrow=1,
  widths=c(
    0.53,
    0.47
  ),
  guides="collect"
)

final_plot <- patchwork::wrap_plots(
  list(
    row_AB,
    pC,
    row_DE
  ),
  ncol=1,
  heights=c(
    0.74,
    1.25,
    1.12
  ),
  guides="keep"
)

save_plot(
  final_plot,
  "Figure03_FINAL",
  width=15.5,
  height=18.0
)

# ============================================================
# 25. SOURCE TABLES

# ============================================================

cat("\n===== SAVE SOURCE TABLES =====\n")

write.csv(
  state_label_df,
  file.path(
    TAB_DIR,
    "22F4F6C_global_K3_label_positions.csv"
  ),
  row.names=FALSE
)

write.csv(
  cell_df,
  file.path(
    TAB_DIR,
    "22F4F6C_cell_UMAP_K3_source.csv"
  ),
  row.names=FALSE
)

write.csv(
  edge_df,
  file.path(
    TAB_DIR,
    "22F4F6C_principal_graph_edges.csv"
  ),
  row.names=FALSE
)

write.csv(
  rep50,
  file.path(
    TAB_DIR,
    "22F4F6C_50_representative_dynamic_genes.csv"
  ),
  row.names=FALSE
)

write.csv(
  heat_long,
  file.path(
    TAB_DIR,
    "22F4F6C_dynamic_gene_heatmap_source.csv"
  ),
  row.names=FALSE
)

write.csv(
  module_geometry,
  file.path(
    TAB_DIR,
    "22F4F6C_module_display_geometry.csv"
  ),
  row.names=FALSE
)

write.csv(
  module_annotation,
  file.path(
    TAB_DIR,
    "22F4F6C_module_functional_annotation_final.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 26. MODULE LEDGER
# ============================================================

module_ledger <- D3$module_summary |>
  dplyr::select(
    module,
    n_genes,
    centroid_peak_position,
    centroid_trough_position,
    centroid_amplitude
  ) |>
  left_join(
    module_annotation |>
      dplyr::select(
        module,
        functional_title,
        pathways,
        evidence
      ),
    by="module"
  ) |>
  arrange(
    match(
      module,
      paste0("M",1:5)
    )
  )

write.csv(
  module_ledger,
  file.path(
    TAB_DIR,
    "22F4F6C_module_ledger.csv"
  ),
  row.names=FALSE
)

cat("\n===== MODULE LEDGER =====\n")

print(
  module_ledger,
  n=Inf
)

# ============================================================
# 27. FREEZE
# ============================================================

freeze <- list(
  module="22F4F6C",

  purpose=paste(
    "Integrated publication visualization of frozen",
    "Bulge/Isthmus trajectory topology, disease-specific",
    "K3 occupancy, dynamic transcriptional modules and",
    "functional annotations."
  ),

  analysis_status=
    "PUBLICATION POLISH ONLY — no reinference",

  n_cells=nrow(cell_df),

  n_patients=
    dplyr::n_distinct(
      cell_df$patient_id
    ),

  disease_counts=
    cell_df |>
      count(
        disease,
        name="n_cells"
      ),

  K3_counts=
    cell_df |>
      count(
        territory_K3,
        name="n_cells"
      ),

  graph=list(
    root_node=C$root_node,
    end_node=C$end_node,
    n_nodes=nrow(node_df),
    n_edges=nrow(edge_df),
    trunk_nodes=C$trunk_nodes,
    n_strict_trunk=C$n_strict_trunk
  ),

  state_reference=state_ref,

  dynamic_modules=list(
    K=D3$K_final,
    n_dynamic_genes=D3$n_dynamic_genes,
    displayed_genes=rep50,
    module_geometry=module_geometry,
    module_summary=D3$module_summary,
    functional_annotation=module_annotation
  ),

  canonical_colours=list(
    K3=K3_COL,
    modules=MODULE_COL
  )
)

saveRDS(
  freeze,
  file.path(
    OBJ_DIR,
    "22F4F6C_integrated_trajectory_figure_freeze.rds"
  )
)

# ============================================================
# 28. FORMAL CHECKS
# ============================================================

cat("\n===== FORMAL CHECKS =====\n")

per_module_display <- rep50 |>
  count(
    module,
    name="n"
  )

checks <- tibble(
  check=c(
    "Bulge cells = 5420",
    "Patients = 16",
    "Principal graph nodes = 40",
    "Principal graph edges = 39",
    "Strict trunk cells = 4295",
    "Root = Y_1",
    "End = Y_12",
    "Dynamic genes = 2774",
    "K = 5",
    "Displayed genes = 50",
    "10 displayed genes/module",
    "All M1-M5 annotated",
    "T3/T2 boundary = 0.438",
    "T2/T1 boundary = 0.844",
    "Canonical T3 colour",
    "Canonical T2 colour",
    "Canonical T1 colour",
    "Panel A SVG exists",
    "Panel B SVG exists",
    "Panel C SVG exists",
    "Figure03 composite SVG exists",
    "Figure03 composite PDF exists",
    "Figure03 composite PNG exists"
  ),

  pass=c(
    nrow(cell_df) == 5420,

    dplyr::n_distinct(
      cell_df$patient_id
    ) == 16,

    nrow(node_df) == 40,

    nrow(edge_df) == 39,

    C$n_strict_trunk == 4295,

    identical(
      C$root_node,
      "Y_1"
    ),

    identical(
      C$end_node,
      "Y_12"
    ),

    D3$n_dynamic_genes == 2774,

    D3$K_final == 5,

    nrow(rep50) == 50,

    nrow(per_module_display) == 5 &&
      all(per_module_display$n == 10),

    all(
      paste0("M",1:5) %in%
        module_annotation$module
    ),

    isTRUE(
      all.equal(
        state_ref$xmax[
          state_ref$state=="T3"
        ],
        0.438
      )
    ),

    isTRUE(
      all.equal(
        state_ref$xmax[
          state_ref$state=="T2"
        ],
        0.844
      )
    ),

    identical(
      unname(K3_COL["T3"]),
      "#2F7189"
    ),

    identical(
      unname(K3_COL["T2"]),
      "#4B4C7A"
    ),

    identical(
      unname(K3_COL["T1"]),
      "#B36B00"
    ),

    file.exists(
      file.path(
        SVG_DIR,
        "Fig22F4F6C_A_global_trajectory.svg"
      )
    ),

    file.exists(
      file.path(
        SVG_DIR,
        "Fig22F4F6C_B_disease_split_trajectory.svg"
      )
    ),

    file.exists(
      file.path(
        SVG_DIR,
        "Fig22F4F6C_C_dynamic_genes_functional_annotation.svg"
      )
    ),

    file.exists(
      file.path(
        SVG_DIR,
        "Figure03_FINAL.svg"
      )
    ),

    file.exists(
      file.path(
        PDF_DIR,
        "Figure03_FINAL.pdf"
      )
    ),

    file.exists(
      file.path(
        PNG_DIR,
        "Figure03_FINAL.png"
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
    AUDIT_DIR,
    "22F4F6C_formal_checks.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 29. AUDIT
# ============================================================

audit_lines <- c(
  "============================================================",
  "22F4F6C — INTEGRATED BULGE TRAJECTORY FIGURE",
  "============================================================",
  "",
  paste0(
    "Cells: ",
    nrow(cell_df)
  ),
  paste0(
    "Patients: ",
    dplyr::n_distinct(
      cell_df$patient_id
    )
  ),
  paste0(
    "Principal graph: ",
    nrow(node_df),
    " nodes / ",
    nrow(edge_df),
    " edges"
  ),
  paste0(
    "Strict trunk cells: ",
    C$n_strict_trunk
  ),
  paste0(
    "Dynamic genes: ",
    D3$n_dynamic_genes
  ),
  paste0(
    "Displayed genes: ",
    nrow(rep50),
    " (10/module)"
  ),
  "",
  "Frozen trajectory:",
  "T3: 0.000-0.438",
  "T2: 0.438-0.844",
  "T1: 0.844-1.000",
  "",
  "Functional module display:",
  "M1 — Early / transcriptional adaptation; Hypoxia FDR=0.061",
  "M2 — Interferon / antigen presentation",
  "M3 — Tissue remodeling / EMT / wound healing",
  "M4 — Proliferative / cell cycle",
  "M5 — Late trajectory program; UV response DN",
  "",
  paste(
    "No trajectory topology, pseudotime, differential-expression",
    "model or K=5 module assignment was re-estimated."
  ),
  "",
  paste0(
    "FORMAL CHECKS: ",
    sum(checks$pass),
    "/",
    nrow(checks),
    " PASS"
  ),
  "",
  if (all(checks$pass)) {
    "STATUS=PASS"
  } else {
    "STATUS=FAIL"
  }
)

writeLines(
  audit_lines,
  file.path(
    AUDIT_DIR,
    "22F4F6C_audit.txt"
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
# 30. FINAL STATUS
# ============================================================

if (!all(checks$pass)) {
  stop(
    "22F4F6C formal checks failed."
  )
}

cat("\n============================================================\n")
cat("22F4F6C COMPLETE\n")
cat("STATUS=PASS\n")
cat("============================================================\n")

