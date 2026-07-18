# Numeric-value type refinement, ported from
# generate_forcefield_inputs_correct_type_bond_angle_dihedrals_DFT_05_13_2024.py.
#
# The Chen-Manz atom-type combination alone can lump together bonds/angles/
# dihedrals whose actual geometric values differ (e.g. two P-O bonds in a
# phosphonate MOF with the same atom-type combination but a genuinely
# different bond length because one P-O is part of a P=O-like linkage and
# the other bridges to the metal). This splits each existing Chen-Manz type
# into finer sub-types so two instances only share a type if their
# equilibrium values also agree:
#
#  - Bonds: same Chen-Manz bond type AND bond length within 1% of a common
#    anchor length (greedy clustering from the first unassigned instance).
#  - Angles: same Chen-Manz angle type AND angle value equal after rounding
#    to 2 decimal places in radians (~0.01 rad bins, matching the original)
#    AND the same (sorted) pair of flanking bond sub-types.
#  - Dihedrals: same Chen-Manz dihedral type AND the same (sorted) pair of
#    flanking angle sub-types, AND (for non-linear instances) the same
#    |dihedral value| after rounding to 2 decimals in radians. Instances
#    already flagged "linear" (a flanking angle within ~1.72 deg of 180)
#    are grouped into their own sub-type per Chen-Manz type rather than
#    further split by value, since a linear torsion doesn't have a
#    well-defined equilibrium dihedral angle.
#
# This must run AFTER add_type_columns() (needs the base Chen-Manz type
# strings) and BEFORE classify_dihedral_rotatability() / dihedral pruning,
# matching the original pipeline's actual step order: typing ->
# DFT-value refinement -> ring/rotatable classification -> dihedral
# pruning. It overwrites bond_type/angle_type/dihedral_type in place with
# the refined label, so every downstream step (type-index assignment,
# summary tables, pruning) picks up the finer grouping automatically.
#
# Two deliberate deviations from the original, consistent with earlier
# fixes in this app:
#  - The original's greedy bond-length clustering has an off-by-one bug
#    that silently drops the LAST unassigned instance in each Chen-Manz
#    bucket from the output entirely (`if j < len-1` skips creating a
#    group for it). This implementation assigns every instance a sub-type.
#  - The original's linear-dihedral "rescue" branch indexes a loop
#    variable from the wrong loop (`bonded_pos_with_types_new[j]` where j
#    is the outer dihedral-instance index, not the inner bond-scan index
#    k) -- almost certainly a bug. Rather than replicate it, linear
#    dihedral instances are simply not split further by value here.

split_bonds_by_length <- function(bonds, tol_frac = 0.01) {
  if (is.null(bonds) || nrow(bonds) == 0) return(bonds)
  new_type <- character(nrow(bonds))
  for (ty in unique(bonds$bond_type)) {
    remaining <- which(bonds$bond_type == ty)
    sub <- 0L
    while (length(remaining) > 0) {
      anchor_len <- bonds$distance[remaining[1]]
      match_mask <- abs(bonds$distance[remaining] - anchor_len) < tol_frac * anchor_len
      sub <- sub + 1L
      new_type[remaining[match_mask]] <- paste0(ty, " #", sub)
      remaining <- remaining[!match_mask]
    }
  }
  bonds$bond_type <- new_type
  bonds
}

split_angles_by_value <- function(angles, bonds, round_digits = 2) {
  if (is.null(angles) || nrow(angles) == 0) return(angles)

  bond_lookup <- new.env(hash = TRUE, parent = emptyenv())
  if (!is.null(bonds) && nrow(bonds) > 0) {
    keys <- paste(pmin(bonds$atom1, bonds$atom2), pmax(bonds$atom1, bonds$atom2))
    for (i in seq_len(nrow(bonds))) {
      if (!exists(keys[i], envir = bond_lookup, inherits = FALSE)) {
        assign(keys[i], bonds$bond_type[i], envir = bond_lookup)
      }
    }
  }
  get_bond_type <- function(a, b) {
    k <- paste(min(a, b), max(a, b))
    if (exists(k, envir = bond_lookup, inherits = FALSE)) get(k, envir = bond_lookup, inherits = FALSE) else NA_character_
  }

  n <- nrow(angles)
  flank1 <- vapply(seq_len(n), function(r) get_bond_type(angles$atom1[r], angles$atom2[r]), character(1))
  flank2 <- vapply(seq_len(n), function(r) get_bond_type(angles$atom2[r], angles$atom3[r]), character(1))
  flank_key <- mapply(function(a, b) paste(sort(c(a, b)), collapse = " :: "), flank1, flank2)
  value_key <- round(angles$angle_deg * pi / 180, round_digits)
  composite <- paste(angles$angle_type, value_key, flank_key, sep = " || ")

  new_type <- character(n)
  for (ty in unique(angles$angle_type)) {
    idx <- which(angles$angle_type == ty)
    uniq_composite <- unique(composite[idx])
    sub_map <- setNames(seq_along(uniq_composite), uniq_composite)
    new_type[idx] <- paste0(ty, " #", sub_map[composite[idx]])
  }
  angles$angle_type <- new_type
  angles
}

split_dihedrals_by_value <- function(dihedrals, angles, round_digits = 2) {
  if (is.null(dihedrals) || nrow(dihedrals) == 0) return(dihedrals)

  linear_tol_deg <- 0.03 * 180 / pi # ~1.72 deg, matches the 0.03 rad tolerance used elsewhere

  angle_lookup <- new.env(hash = TRUE, parent = emptyenv())
  if (!is.null(angles) && nrow(angles) > 0) {
    keys <- paste(angles$atom2, pmin(angles$atom1, angles$atom3), pmax(angles$atom1, angles$atom3))
    for (i in seq_len(nrow(angles))) {
      if (!exists(keys[i], envir = angle_lookup, inherits = FALSE)) {
        assign(keys[i], angles$angle_type[i], envir = angle_lookup)
      }
    }
  }
  get_angle_type <- function(center, a, b) {
    k <- paste(center, min(a, b), max(a, b))
    if (exists(k, envir = angle_lookup, inherits = FALSE)) get(k, envir = angle_lookup, inherits = FALSE) else NA_character_
  }

  n <- nrow(dihedrals)
  flank1 <- vapply(seq_len(n), function(r) get_angle_type(dihedrals$atom2[r], dihedrals$atom1[r], dihedrals$atom3[r]), character(1))
  flank2 <- vapply(seq_len(n), function(r) get_angle_type(dihedrals$atom3[r], dihedrals$atom2[r], dihedrals$atom4[r]), character(1))
  flank_key <- mapply(function(a, b) paste(sort(c(a, b)), collapse = " :: "), flank1, flank2)

  is_linear <- (abs(abs(dihedrals$flank_angle1_deg) - 180) < linear_tol_deg) |
               (abs(abs(dihedrals$flank_angle2_deg) - 180) < linear_tol_deg)
  value_key <- ifelse(is_linear, "linear",
                       as.character(round(abs(dihedrals$dihedral_deg) * pi / 180, round_digits)))

  composite <- paste(dihedrals$dihedral_type, value_key, flank_key, sep = " || ")

  new_type <- character(n)
  for (ty in unique(dihedrals$dihedral_type)) {
    idx <- which(dihedrals$dihedral_type == ty)
    uniq_composite <- unique(composite[idx])
    sub_map <- setNames(seq_along(uniq_composite), uniq_composite)
    new_type[idx] <- paste0(ty, " #", sub_map[composite[idx]])
  }
  dihedrals$dihedral_type <- new_type
  dihedrals
}
