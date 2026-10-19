# ============================================================
# 22F4F5 — BULGE TRAJECTORY + DYNAMIC GENES
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
cat("22F4F5 — BULGE TRAJECTORY + DYNAMIC GENES\n")
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
  BASE,
  "22F4F_PUBLICATION_FIGURES",
  "22F4F5_TRAJECTORY_DYNAMIC_GENES"
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
# 2. LOAD FROZEN INPUTS
# ============================================================

cat("\n===== LOAD FROZEN INPUTS =====\n")

cds <- readRDS(A_FILE)
C <- readRDS(C_FILE)
D2 <- readRDS(D2_FILE)
D3 <- readRDS(D3_FILE)

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
# 9. PANEL A
# FULL UMAP + PRINCIPAL GRAPH + K3 STATES
# ============================================================

cat("\n===== PANEL A =====\n")

set.seed(22045)

plot_cells_A <- cell_df |>
  slice_sample(
    prop=1
  )

special_labels_A <- node_df |>
  filter(
    principal_node %in% c(
      "Y_1",
      "Y_12",
      "Y_37"
    )
  )

pA <- ggplot() +

  # ----------------------------------------------------------
  # Original K3 UMAP appearance
  # ----------------------------------------------------------

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

  # ----------------------------------------------------------
  # Full Monocle graph: subtle context
  # ----------------------------------------------------------

  geom_segment(
    data=edge_df |>
      filter(!on_trunk),
    aes(
      x=x,
      y=y,
      xend=xend,
      yend=yend
    ),
    colour="#4D4D4D",
    linewidth=0.42,
    alpha=0.60,
    lineend="round"
  ) +

  # ----------------------------------------------------------
  # Frozen strict trunk
  # ----------------------------------------------------------

  geom_segment(
    data=edge_df |>
      filter(on_trunk),
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

  # ----------------------------------------------------------
  # Only biologically relevant node labels.
  # No white graph-node circles.
  # ----------------------------------------------------------

  geom_text(
    data=special_labels_A,
    aes(
      x=UMAP_1,
      y=UMAP_2,
      label=principal_node
    ),
    colour="#111111",
    size=2.5,
    fontface="bold",
    nudge_y=0.40,
    check_overlap=TRUE
  ) +

  scale_fill_manual(
    values=state_cols,
    breaks=c(
      "T3",
      "T2",
      "T1"
    ),
    drop=FALSE
  ) +

  coord_equal() +

  labs(
    title="Bulge/Isthmus principal trajectory",
    subtitle="Frozen Monocle3 topology over the canonical K3 state landscape",
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
# 10. PANEL B
# DISEASE OCCUPANCY ON SAME GRAPH
# ============================================================

cat("\n===== PANEL B =====\n")

set.seed(22046)

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
      fill=disease
    ),
    shape=21,
    size=0.85,
    stroke=0.12,
    colour="#3A3A3A",
    alpha=0.88
  ) +

  geom_segment(
    data=edge_df,
    aes(
      x=x,
      y=y,
      xend=xend,
      yend=yend
    ),
    colour="black",
    linewidth=0.62,
    alpha=0.78,
    lineend="round"
  ) +

  facet_wrap(
    ~disease,
    nrow=1
  ) +

  scale_fill_manual(
    values=disease_cols,
    drop=FALSE
  ) +

  coord_equal() +

  labs(
    title="Disease occupancy of the shared trajectory landscape",
    x="UMAP 1",
    y="UMAP 2",
    fill="Disease"
  ) +

  theme_classic(
    base_size=9
  ) +

  theme(
    plot.title=element_text(
      face="bold",
      size=10.5
    ),
    strip.background=element_blank(),
    strip.text=element_text(
      face="bold",
      size=10
    ),
    legend.position="none"
  )

# ============================================================
# 11. PANEL C
# 75 REPRESENTATIVE DYNAMIC GENES
# ============================================================

cat("\n===== PANEL C =====\n")

rep_genes <- D3$representative_genes |>
  arrange(
    module_order,
    order_within_module
  )

stopifnot(
  nrow(rep_genes) == 75,
  all(rep_genes$gene %in% rownames(D3$row_scaled_matrix))
)

heat_mat <- D3$row_scaled_matrix[
  rep_genes$gene,
  ,
  drop=FALSE
]

heat_long <- as.data.frame(heat_mat) |>
  rownames_to_column("gene") |>
  pivot_longer(
    cols=-gene,
    names_to="trajectory_bin",
    values_to="z"
  ) |>
  left_join(
    rep_genes |>
      select(
        gene,
        module,
        module_order,
        display_order,
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
      (trajectory_index - 1) / 100
  )

gene_levels <- rev(
  rep_genes$gene
)

heat_long$gene <- factor(
  heat_long$gene,
  levels=gene_levels
)

state_reference <- D3$state_reference

state_labels <- state_reference |>
  mutate(
    midpoint=(xmin+xmax)/2
  )

module_bounds_display <- rep_genes |>
  group_by(
    module,
    module_order
  ) |>
  summarise(
    first=min(display_order),
    last=max(display_order),
    .groups="drop"
  ) |>
  arrange(module_order)

pC <- ggplot(
  heat_long,
  aes(
    x=trajectory_position,
    y=gene,
    fill=z
  )
) +

  geom_tile() +

  geom_vline(
    xintercept=c(
      0.438,
      0.844
    ),
    colour="white",
    linewidth=0.65,
    linetype="solid"
  ) +

  scale_fill_gradient2(
    low="#2166AC",
    mid="white",
    high="#B2182B",
    midpoint=0,
    limits=c(-2.5,2.5),
    oob=scales::squish,
    name="Trajectory\nz-score"
  ) +

  scale_x_continuous(
    breaks=c(
      0,
      0.219,
      0.438,
      0.641,
      0.844,
      1
    ),
    labels=c(
      "0",
      "",
      "0.438",
      "",
      "0.844",
      "1"
    ),
    expand=c(0,0)
  ) +

  coord_cartesian(
    xlim=c(0,1),
    expand=FALSE
  ) +

  labs(
    title="Dynamic genes along the Bulge/Isthmus trajectory",
    subtitle="75 representative genes from five frozen trajectory-shape modules",
    x="Trajectory position  (T3 -> T2 -> T1)",
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
      size=8.5
    ),
    axis.text.y=element_text(
      size=5.7,
      face="italic"
    ),
    axis.text.x=element_text(
      size=7.5
    ),
    axis.ticks.y=element_blank(),
    axis.line.y=element_blank(),
    legend.position="right",
    legend.title=element_text(
      size=8
    ),
    legend.text=element_text(
      size=7
    ),
    plot.margin=margin(
      5.5,
      5.5,
      5.5,
      5.5
    )
  )

# ============================================================
# 11B. PANEL C TOP STATE STRIP
# ============================================================

state_strip <- ggplot(
  state_reference,
  aes(
    xmin=xmin,
    xmax=xmax,
    ymin=0,
    ymax=1,
    fill=state
  )
) +

  geom_rect(
    colour="white",
    linewidth=0.5
  ) +

  geom_text(
    data=state_labels,
    aes(
      x=midpoint,
      y=0.5,
      label=state
    ),
    inherit.aes=FALSE,
    fontface="bold",
    size=3.2
  ) +

  scale_fill_manual(
    values=state_cols,
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

  theme_void() +

  theme(
    plot.margin=margin(
      0,
      5.5,
      0,
      5.5
    )
  )

pC_full <- state_strip / pC +
  patchwork::plot_layout(
    heights=c(0.055,1)
  )

# ============================================================
# 12. PANEL D
# FROZEN MODULE CENTROID CURVES
# ============================================================

cat("\n===== PANEL D =====\n")

centroid_df <- D3$centroid_curves |>
  mutate(
    module=factor(
      module,
      levels=paste0(
        "M",
        1:5
      )
    )
  )

pD <- ggplot(
  centroid_df,
  aes(
    x=trajectory_position,
    y=z,
    colour=module
  )
) +

  annotate(
    "rect",
    xmin=0,
    xmax=0.438,
    ymin=-Inf,
    ymax=Inf,
    fill="#2166AC",
    alpha=0.035
  ) +

  annotate(
    "rect",
    xmin=0.438,
    xmax=0.844,
    ymin=-Inf,
    ymax=Inf,
    fill="#FDAE61",
    alpha=0.035
  ) +

  annotate(
    "rect",
    xmin=0.844,
    xmax=1,
    ymin=-Inf,
    ymax=Inf,
    fill="#B2182B",
    alpha=0.035
  ) +

  geom_vline(
    xintercept=c(
      0.438,
      0.844
    ),
    linetype="dashed",
    linewidth=0.4,
    colour="#777777"
  ) +

  geom_hline(
    yintercept=0,
    linewidth=0.3,
    colour="#BDBDBD"
  ) +

  geom_line(
    linewidth=1
  ) +

  scale_colour_manual(
    values=module_cols,
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
    expand=c(0,0)
  ) +

  labs(
    title="Five frozen trajectory-shape modules",
    subtitle="Centroid expression dynamics across T3 -> T2 -> T1",
    x="Trajectory position  (T3 -> T2 -> T1)",
    y="Mean trajectory z-score",
    colour="Module"
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
      size=8.5
    ),
    legend.position="right"
  )

# ============================================================
# 13. SAVE FUNCTION
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
# 14. SAVE INDIVIDUAL PANELS
# ============================================================

cat("\n===== SAVE INDIVIDUAL PANELS =====\n")

save_plot(
  pA,
  "Fig22F4F5_A_UMAP_principal_graph",
  7.2,
  6.3
)

save_plot(
  pB,
  "Fig22F4F5_B_disease_trajectory_occupancy",
  12.5,
  4.8
)

save_plot(
  pC_full,
  "Fig22F4F5_C_dynamic_gene_heatmap",
  9.2,
  10.5
)

save_plot(
  pD,
  "Fig22F4F5_D_dynamic_module_centroids",
  7.2,
  5.2
)

# ============================================================
# 15. PUBLICATION COMPOSITE
# ============================================================

cat("\n===== BUILD PUBLICATION COMPOSITE =====\n")

top_row <- patchwork::wrap_plots(
  list(
    pA,
    pB
  ),
  nrow=1,
  widths=c(
    0.82,
    1.55
  ),
  guides="collect"
)

bottom_row <- patchwork::wrap_plots(
  list(
    pC_full,
    pD
  ),
  nrow=1,
  widths=c(
    1.35,
    1
  ),
  guides="keep"
)

final_plot <- patchwork::wrap_plots(
  list(
    top_row,
    bottom_row
  ),
  ncol=1,
  heights=c(
    0.82,
    1.35
  )
) +
  patchwork::plot_annotation(
    title="Bulge/Isthmus trajectory and dynamic transcriptional remodeling",
    tag_levels="A",
    theme=theme(
      plot.title=element_text(
        face="bold",
        size=14,
        hjust=0
      ),
      plot.tag=element_text(
        face="bold",
        size=13
      )
    )
  )

save_plot(
  final_plot,
  "Fig22F4F5_Bulge_trajectory_dynamic_genes",
  16,
  12.5
)

# ============================================================
# 16. SOURCE TABLES
# ============================================================

cat("\n===== SAVE SOURCE TABLES =====\n")

write.csv(
  cell_df,
  file.path(
    TAB_DIR,
    "22F4F5_cell_UMAP_trajectory_source.csv"
  ),
  row.names=FALSE
)

write.csv(
  edge_df,
  file.path(
    TAB_DIR,
    "22F4F5_principal_graph_edges.csv"
  ),
  row.names=FALSE
)

write.csv(
  node_df,
  file.path(
    TAB_DIR,
    "22F4F5_principal_graph_nodes.csv"
  ),
  row.names=FALSE
)

write.csv(
  rep_genes,
  file.path(
    TAB_DIR,
    "22F4F5_representative_dynamic_genes.csv"
  ),
  row.names=FALSE
)

write.csv(
  heat_long,
  file.path(
    TAB_DIR,
    "22F4F5_representative_gene_trajectory_matrix_long.csv"
  ),
  row.names=FALSE
)

write.csv(
  centroid_df,
  file.path(
    TAB_DIR,
    "22F4F5_dynamic_module_centroid_curves.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 17. FREEZE
# ============================================================

freeze <- list(
  module="22F4F5",
  purpose=paste(
    "Publication visualization of the frozen Bulge/Isthmus",
    "Monocle3 topology and K=5 dynamic transcriptional modules.",
    "No trajectory reinference or new differential-expression testing."
  ),
  n_cells=nrow(cell_df),
  n_graph_nodes=nrow(node_df),
  n_graph_edges=nrow(edge_df),
  root_node=C$root_node,
  end_node=C$end_node,
  trunk_nodes=C$trunk_nodes,
  state_reference=D3$state_reference,
  n_dynamic_genes=D3$n_dynamic_genes,
  K_final=D3$K_final,
  representative_genes=rep_genes,
  module_summary=D3$module_summary
)

saveRDS(
  freeze,
  file.path(
    OBJ_DIR,
    "22F4F5_trajectory_dynamic_genes_freeze.rds"
  )
)

# ============================================================
# 18. FORMAL CHECKS
# ============================================================

cat("\n===== FORMAL CHECKS =====\n")

checks <- tibble(
  check=c(
    "Bulge cell count = 5420",
    "principal graph nodes = 40",
    "principal graph edges = 39",
    "root node = Y_1",
    "end node = Y_12",
    "dynamic genes = 2774",
    "K dynamic modules = 5",
    "representative genes = 75",
    "T3/T2 boundary = 0.438",
    "T2/T1 boundary = 0.844",
    "SVG composite exists",
    "PDF composite exists",
    "PNG composite exists"
  ),
  pass=c(
    nrow(cell_df) == 5420,
    nrow(node_df) == 40,
    nrow(edge_df) == 39,
    identical(C$root_node,"Y_1"),
    identical(C$end_node,"Y_12"),
    D3$n_dynamic_genes == 2774,
    D3$K_final == 5,
    nrow(rep_genes) == 75,
    isTRUE(
      all.equal(
        D3$state_reference$xmax[
          D3$state_reference$state=="T3"
        ],
        0.438
      )
    ),
    isTRUE(
      all.equal(
        D3$state_reference$xmax[
          D3$state_reference$state=="T2"
        ],
        0.844
      )
    ),
    file.exists(
      file.path(
        SVG_DIR,
        "Fig22F4F5_Bulge_trajectory_dynamic_genes.svg"
      )
    ),
    file.exists(
      file.path(
        PDF_DIR,
        "Fig22F4F5_Bulge_trajectory_dynamic_genes.pdf"
      )
    ),
    file.exists(
      file.path(
        PNG_DIR,
        "Fig22F4F5_Bulge_trajectory_dynamic_genes.png"
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
    "22F4F5_formal_checks.csv"
  ),
  row.names=FALSE
)

# ============================================================
# 19. AUDIT
# ============================================================

audit_lines <- c(
  "============================================================",
  "22F4F5 — BULGE TRAJECTORY + DYNAMIC GENES",
  "============================================================",
  "",
  paste0(
    "Cells: ",
    nrow(cell_df)
  ),
  paste0(
    "Principal graph nodes: ",
    nrow(node_df)
  ),
  paste0(
    "Principal graph edges: ",
    nrow(edge_df)
  ),
  paste0(
    "Frozen root: ",
    C$root_node
  ),
  paste0(
    "Frozen end: ",
    C$end_node
  ),
  paste0(
    "Strict trunk cells: ",
    C$n_strict_trunk
  ),
  paste0(
    "Strict trunk fraction: ",
    round(
      100*C$strict_trunk_fraction,
      1
    ),
    "%"
  ),
  paste0(
    "Dynamic genes: ",
    D3$n_dynamic_genes
  ),
  paste0(
    "Dynamic modules: K=",
    D3$K_final
  ),
  paste0(
    "Representative genes displayed: ",
    nrow(rep_genes)
  ),
  "",
  "Frozen state sequence:",
  "T3: 0.000–0.438",
  "T2: 0.438–0.844",
  "T1: 0.844–1.000",
  "",
  paste(
    "No topology, pseudotime model, DE model,",
    "module assignment or state boundary was",
    "re-estimated in 22F4F5."
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
    "22F4F5_audit.txt"
  )
)

cat(
  paste(
    audit_lines,
    collapse="\n"
  ),
  "\n"
)

if (!all(checks$pass)) {
  stop(
    "22F4F5 formal checks failed."
  )
}

cat("\n============================================================\n")
cat("22F4F5 COMPLETE\n")
cat("STATUS=PASS\n")
cat("============================================================\n")

