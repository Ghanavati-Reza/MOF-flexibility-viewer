# Dihedral-type pruning, ported from dihedral_pruning_materials_1_30_2024.py.
#
# Idea: two different dihedral TYPES can end up describing torsion around
# the exact same set of middle bonds (e.g. atom j's neighbors include two
# different atom types, atom k's neighbors include two different atom
# types -> up to 4 distinct 4-atom-type combinations, i.e. up to 4 distinct
# "types", all sitting on the same physical j-k bond). Since a forcefield
# only needs ONE torsion term per physical bond, types that share an
# identical middle-bond set are redundant -- keep just one representative
# per group and drop the rest (both the type and its instances).
#
# Middle-bond identity is by atom NUMBER only (not periodic image), matching
# the original: two instances whose center/corner atoms are the same pair
# of atom numbers count as "the same middle bond" even if reached through
# different images.
#
# Representative choice within a redundant group: prefer the type whose
# flanking bond angles are furthest from 180 deg (least likely to be a
# numerically ill-conditioned/near-linear torsion), penalized by instance
# count. This ports the original's `(pi - max_angle) / count` metric.
#
# Two deliberate deviations from the original, both worth calling out:
#  - The original accumulates `dihedrals_angle_1`/`dihedrals_angle_2` across
#    ALL types processed so far without resetting them per type, so its
#    "average flanking angle" is actually a running average over every type
#    seen up to that point, not that type's own instances. That looks like
#    a bug (a missing `= []` reset each loop iteration), so this
#    implementation computes the average over each type's own instances only.
#  - The original breaks ties between equally-good candidates with a seeded
#    random pick. This picks the first (lowest type index) instead, so the
#    result is deterministic given the same input and tolerance.

middle_bond_set_key <- function(dihedrals, type_idx) {
  sub <- dihedrals[dihedrals$dihedral_type_index == type_idx, , drop = FALSE]
  pairs <- mapply(function(a, b) paste(sort(c(a, b)), collapse = "_"), sub$atom2, sub$atom3)
  paste(sort(unique(pairs)), collapse = ";")
}

dihedral_prune_metric <- function(dihedrals, type_idx) {
  sub <- dihedrals[dihedrals$dihedral_type_index == type_idx, , drop = FALSE]
  m1 <- mean(sub$flank_angle1_deg)
  m2 <- mean(sub$flank_angle2_deg)
  (180 - max(m1, m2)) / nrow(sub)
}

# Returns the pruning result: which type indices survive, plus (for
# transparency) how many redundant types were collapsed into each survivor.
prune_dihedral_types <- function(dihedrals, dihedral_types_summary) {
  if (is.null(dihedrals) || nrow(dihedrals) == 0 ||
      is.null(dihedral_types_summary) || nrow(dihedral_types_summary) == 0) {
    return(list(kept_indices = integer(0), merged_count = integer(0)))
  }

  idx <- dihedral_types_summary$Index
  keys <- vapply(idx, function(i) middle_bond_set_key(dihedrals, i), character(1))
  metrics <- vapply(idx, function(i) dihedral_prune_metric(dihedrals, i), numeric(1))

  kept <- integer(0)
  merged_count <- integer(0)
  processed <- rep(FALSE, length(idx))

  for (i in seq_along(idx)) {
    if (processed[i]) next
    same_group <- which(keys == keys[i])
    processed[same_group] <- TRUE
    best_local <- same_group[which.max(metrics[same_group])]
    kept <- c(kept, idx[best_local])
    merged_count <- c(merged_count, length(same_group) - 1L)
  }

  ord <- order(kept)
  list(kept_indices = kept[ord], merged_count = merged_count[ord])
}

# Applies the pruning result to produce the two extra "after pruning" lists.
apply_dihedral_pruning <- function(dihedrals, dihedral_types_summary) {
  pruning <- prune_dihedral_types(dihedrals, dihedral_types_summary)

  if (length(pruning$kept_indices) == 0) {
    return(list(dihedrals_pruned = data.frame(), dihedral_types_summary_pruned = data.frame()))
  }

  dihedrals_pruned <- dihedrals[dihedrals$dihedral_type_index %in% pruning$kept_indices, , drop = FALSE]

  dihedral_types_summary_pruned <- dihedral_types_summary[
    dihedral_types_summary$Index %in% pruning$kept_indices, , drop = FALSE
  ]
  dihedral_types_summary_pruned$Merged <- pruning$merged_count[
    match(dihedral_types_summary_pruned$Index, pruning$kept_indices)
  ]

  list(dihedrals_pruned = dihedrals_pruned, dihedral_types_summary_pruned = dihedral_types_summary_pruned)
}
