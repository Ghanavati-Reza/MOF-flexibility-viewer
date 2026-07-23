# Double-counted-image removal, ported from
# code_for_counting_double_counted_images_stretch_dihedral_1_30_2024.py.
#
# Two instances within the same type are "duplicates" if they involve the
# same atom NUMBERS in the same or reverse order. For dihedrals this
# matches the original exactly: atom identity only, periodic image is
# ignored, since a torsion between the same four named atoms is the same
# physical interaction for parameter-fitting purposes regardless of which
# translated copy produced this particular instance. Only one
# representative instance per duplicate group is kept.
#
# For bonds, working through the original's image-difference check shows
# it only ever matches truly identical (atom1, atom2, image) rows (a safety
# net for accidental exact duplicates from earlier steps) -- two bonds
# between the same atom pair via genuinely different periodic images are
# NOT duplicates of each other (e.g. each atom in a periodic chain/sheet is
# legitimately bonded to more than one image of its neighbor), so bond
# dedup here only removes literal exact-duplicate rows.

dedupe_bonds <- function(bonds) {
  if (is.null(bonds) || nrow(bonds) == 0) {
    return(list(deduped = bonds, duplicate_count = integer(0)))
  }
  key <- paste(bonds$atom1, bonds$atom2, bonds$img_a, bonds$img_b, bonds$img_c)
  first_idx <- match(unique(key), key)
  dup_count <- as.integer(table(key)[key[first_idx]])
  deduped <- bonds[first_idx, , drop = FALSE]
  deduped$duplicate_count <- dup_count
  list(deduped = deduped, duplicate_count = dup_count)
}

dedupe_dihedrals <- function(dihedrals) {
  if (is.null(dihedrals) || nrow(dihedrals) == 0) {
    return(list(deduped = dihedrals, duplicate_count = integer(0)))
  }
  fwd <- paste(dihedrals$atom1, dihedrals$atom2, dihedrals$atom3, dihedrals$atom4)
  rev <- paste(dihedrals$atom4, dihedrals$atom3, dihedrals$atom2, dihedrals$atom1)
  # canonical key: within a type, group by atom-number sequence regardless
  # of forward/reverse direction, so a dihedral and its reverse-listed
  # duplicate collapse together
  canon_key <- ifelse(fwd <= rev, fwd, rev)
  key <- paste(dihedrals$dihedral_type, canon_key)

  first_idx <- match(unique(key), key)
  dup_count <- as.integer(table(key)[key[first_idx]])
  deduped <- dihedrals[first_idx, , drop = FALSE]
  deduped$duplicate_count <- dup_count
  list(deduped = deduped, duplicate_count = dup_count)
}
