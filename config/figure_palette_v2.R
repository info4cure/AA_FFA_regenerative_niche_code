# ============================================================
# PIPELINE_V2 — CANONICAL FIGURE PALETTE
# ============================================================
#
# Version: V2_20260930
# Status : FROZEN — disease + L1
#
# PURPOSE
# -------
# Canonical colour specification for all downstream figures.
#
# RULES
# -----
# 1. Disease colours are invariant across the project.
# 2. L1 compartment colours are invariant across the project.
# 3. L2 identities inherit the chromatic family of their
#    parent L1 compartment.
# 4. L2 colours will be individually frozen after visual
#    discrimination QC.
# 5. Figure scripts must source this file rather than define
#    local disease/L1 palettes.
#
# IMPORTANT
# ---------
# This file is a graphical specification only.
# It does NOT alter biological annotation or frozen identities.
# ============================================================


# ------------------------------------------------------------
# 1. Disease palette
# ------------------------------------------------------------

DISEASE_COLORS <- c(
  HC  = "#8F8D8E",
  FFA = "#4D9E4F",
  AA  = "#A8142A"
)


# ------------------------------------------------------------
# 2. L1 compartment palette
# ------------------------------------------------------------

L1_COLORS <- c(
  Epithelial     = "#8A034D",
  Fibroblast     = "#4F8A6B",
  Immune         = "#75559A",
  Endothelial    = "#258B8B",
  Mural          = "#B7643E",
  Sweat_gland    = "#3975B7",
  Melanocyte     = "#85562F",
  Sebocyte       = "#C5942D",
  Neural_Schwann = "#A87900"
)


# ------------------------------------------------------------
# 3. Canonical ordering
# ------------------------------------------------------------

DISEASE_ORDER <- c(
  "HC",
  "FFA",
  "AA"
)

L1_ORDER <- c(
  "Epithelial",
  "Fibroblast",
  "Immune",
  "Endothelial",
  "Mural",
  "Sweat_gland",
  "Melanocyte",
  "Sebocyte",
  "Neural_Schwann"
)


# ------------------------------------------------------------
# 4. Validation
# ------------------------------------------------------------

stopifnot(
  identical(
    names(DISEASE_COLORS),
    DISEASE_ORDER
  )
)

stopifnot(
  identical(
    names(L1_COLORS),
    L1_ORDER
  )
)

stopifnot(
  !anyDuplicated(DISEASE_COLORS)
)

stopifnot(
  !anyDuplicated(L1_COLORS)
)

stopifnot(
  all(
    grepl(
      "^#[0-9A-Fa-f]{6}$",
      c(
        DISEASE_COLORS,
        L1_COLORS
      )
    )
  )
)


# ------------------------------------------------------------
# 5. L2 palette
# ------------------------------------------------------------
#
# Intentionally not defined yet.
#
# L2 colours will be derived from the corresponding L1
# chromatic family and frozen only after visual QC.
#
# Example principle:
#
# Epithelial L2 -> burgundy / magenta family
# Fibroblast L2 -> green family
# Immune L2 -> violet family
# Endothelial L2 -> teal family
# Mural L2 -> terracotta family
# Sweat-gland L2 -> blue family
# Melanocyte L2 -> brown family
# Sebocyte L2 -> ochre/gold family
# Neural/Schwann -> dark yellow/gold family
#
# ============================================================
# END CANONICAL PALETTE
# ============================================================
