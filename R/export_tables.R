# Builds the same friendly-column-named data frame shown in each DT tab, so
# both the on-screen tables and the "Download All (Excel)" workbook use one
# shared definition instead of two copies that could drift apart.

.empty_msg <- function(msg) data.frame(Message = msg)

export_bonds <- function(res) {
  if (is.null(res$bonds) || nrow(res$bonds) == 0) return(.empty_msg("No bonds found"))
  df <- res$bonds[, c("atom1", "element1", "image1", "atom2", "element2", "image2",
                       "bond", "distance", "bond_type_index", "bond_type")]
  colnames(df) <- c("Atom 1 #", "Elem 1", "Image 1", "Atom 2 #", "Elem 2", "Image 2",
                     "Bond", "Distance (Å)", "Type #", "Bond Type")
  df
}

export_angles <- function(res) {
  if (is.null(res$angles) || nrow(res$angles) == 0) return(.empty_msg("No angles found"))
  df <- res$angles[, c("atom1", "element1", "image1", "atom2", "element2", "image2",
                        "atom3", "element3", "image3",
                        "angle_label", "angle_deg", "ring_label", "angle_type_index", "angle_type")]
  colnames(df) <- c("Atom 1 #", "Elem 1", "Image 1", "Atom 2 (center) #", "Elem 2", "Image 2",
                     "Atom 3 #", "Elem 3", "Image 3",
                     "Angle", "Angle (deg)", "Ring", "Type #", "Angle Type")
  df
}

export_dihedrals <- function(res) {
  if (is.null(res$dihedrals) || nrow(res$dihedrals) == 0) return(.empty_msg("No dihedrals found"))
  df <- res$dihedrals[, c("atom1", "element1", "image1", "atom2", "element2", "image2",
                           "atom3", "element3", "image3", "atom4", "element4", "image4",
                           "dihedral_label", "dihedral_deg", "ring_excluded",
                           "dihedral_type_index", "dihedral_type", "rotatability")]
  df$ring_excluded <- ifelse(df$ring_excluded, "Yes", "No")
  colnames(df) <- c("Atom 1 #", "Elem 1", "Image 1", "Atom 2 #", "Elem 2", "Image 2",
                     "Atom 3 #", "Elem 3", "Image 3", "Atom 4 #", "Elem 4", "Image 4",
                     "Dihedral", "Dihedral (deg)", "Ring-Excluded", "Type #", "Dihedral Type", "Rotatability")
  df
}

export_atom_types <- function(res) {
  if (is.null(res$elements) || length(res$elements) == 0) return(.empty_msg("No atoms found"))
  data.frame(Atom = seq_along(res$elements), Element = res$elements, `Chen-Manz atom type` = res$atom_types,
             check.names = FALSE)
}

export_bond_types <- function(res) {
  if (is.null(res$bond_types_summary) || nrow(res$bond_types_summary) == 0) return(.empty_msg("No bonds found"))
  df <- res$bond_types_summary[, c("Index", "bond", "Type", "Count", "Mean", "Min", "Max")]
  colnames(df) <- c("Type #", "Elements", "Bond Type", "Count", "Mean Distance (Å)", "Min (Å)", "Max (Å)")
  df
}

export_angle_types <- function(res) {
  if (is.null(res$angle_types_summary) || nrow(res$angle_types_summary) == 0) return(.empty_msg("No angles found"))
  df <- res$angle_types_summary[, c("Index", "angle_label", "Type", "Count", "Mean", "Min", "Max")]
  colnames(df) <- c("Type #", "Elements", "Angle Type", "Count", "Mean Angle (deg)", "Min (deg)", "Max (deg)")
  df
}

export_dihedral_types <- function(res) {
  if (is.null(res$dihedral_types_summary) || nrow(res$dihedral_types_summary) == 0) return(.empty_msg("No dihedrals found"))
  df <- res$dihedral_types_summary[, c("Index", "dihedral_label", "Type", "Count", "Mean", "Min", "Max", "rotatability")]
  colnames(df) <- c("Type #", "Elements", "Dihedral Type", "Count", "Mean |Dihedral| (deg)", "Min (deg)", "Max (deg)", "Rotatability")
  df
}

export_ub_bonds <- function(res) {
  if (is.null(res$ub_bonds) || nrow(res$ub_bonds) == 0) return(.empty_msg("No 4-membered rings found -- no Urey-Bradley diagonals to add"))
  df <- res$ub_bonds[, c("atom1", "element1", "image1", "atom2", "element2", "image2",
                          "bond", "distance", "ub_type")]
  colnames(df) <- c("Atom 1 #", "Elem 1", "Image 1", "Atom 2 #", "Elem 2", "Image 2",
                     "Diagonal", "Distance (Å)", "UB Type")
  df
}

export_bonds_dedup <- function(res) {
  if (is.null(res$bonds_deduplicated) || nrow(res$bonds_deduplicated) == 0) return(.empty_msg("No bonds found"))
  df <- res$bonds_deduplicated[, c("atom1", "element1", "image1", "atom2", "element2", "image2",
                                    "bond", "distance", "bond_type_index", "bond_type", "duplicate_count")]
  colnames(df) <- c("Atom 1 #", "Elem 1", "Image 1", "Atom 2 #", "Elem 2", "Image 2",
                     "Bond", "Distance (Å)", "Type #", "Bond Type", "Duplicate Count")
  df
}

export_dihedrals_pruned <- function(res) {
  if (is.null(res$dihedrals_pruned) || nrow(res$dihedrals_pruned) == 0) return(.empty_msg("No dihedrals found"))
  df <- res$dihedrals_pruned[, c("atom1", "element1", "image1", "atom2", "element2", "image2",
                                  "atom3", "element3", "image3", "atom4", "element4", "image4",
                                  "dihedral_label", "dihedral_deg", "duplicate_count",
                                  "dihedral_type_index", "dihedral_type", "rotatability")]
  colnames(df) <- c("Atom 1 #", "Elem 1", "Image 1", "Atom 2 #", "Elem 2", "Image 2",
                     "Atom 3 #", "Elem 3", "Image 3", "Atom 4 #", "Elem 4", "Image 4",
                     "Dihedral", "Dihedral (deg)", "Duplicate Count", "Type #", "Dihedral Type", "Rotatability")
  df
}

export_dihedral_types_pruned <- function(res) {
  if (is.null(res$dihedral_types_summary_pruned) || nrow(res$dihedral_types_summary_pruned) == 0) return(.empty_msg("No dihedrals found"))
  df <- res$dihedral_types_summary_pruned[, c("Index", "dihedral_label", "Type", "Count", "Mean", "Min", "Max",
                                                "rotatability", "Merged")]
  colnames(df) <- c("Type #", "Elements", "Dihedral Type", "Count", "Mean |Dihedral| (deg)", "Min (deg)", "Max (deg)",
                     "Rotatability", "Redundant Types Merged")
  df
}

# Sheet name -> export function, in the order sheets should appear in the workbook.
EXPORT_SHEET_BUILDERS <- list(
  "Bonds" = export_bonds,
  "Angles" = export_angles,
  "Dihedrals" = export_dihedrals,
  "AtomTypes" = export_atom_types,
  "BondTypes" = export_bond_types,
  "AngleTypes" = export_angle_types,
  "DihedralTypes" = export_dihedral_types,
  "UreyBradley" = export_ub_bonds,
  "BondsDedup" = export_bonds_dedup,
  "DihedralsPruned" = export_dihedrals_pruned,
  "DihedralTypesPruned" = export_dihedral_types_pruned
)

build_export_workbook_sheets <- function(res) {
  lapply(EXPORT_SHEET_BUILDERS, function(f) f(res))
}
