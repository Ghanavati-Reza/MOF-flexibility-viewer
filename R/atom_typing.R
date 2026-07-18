# Chen-Manz forcefield-precursor atom typing, ported from
# create_atom_type_dir_POSCAR_1_30_2024.py.
# Ref: T. Chen and T. A. Manz, "A collection of forcefield precursors for
# metal-organic frameworks", RSC Adv., 2019, 9, 36492-36507.
#
# Each atom's type string encodes its atomic number, and for every bonded
# neighbor: that neighbor's atomic number plus the (sorted) atomic numbers
# of the neighbor's *other* neighbors (i.e. the 2nd coordination shell,
# with the bond back to the original atom removed once).

chen_manz_atom_types <- function(elements, neighbor_list) {
  n <- length(elements)
  atomic_num <- unname(ATOMIC_NUM[elements])

  neighbor_atomic_num <- vector("list", n)
  for (p in seq_len(n)) {
    nb <- neighbor_list[[p]]
    if (length(nb$idx) == 0) {
      neighbor_atomic_num[[p]] <- integer(0)
    } else {
      an <- atomic_num[nb$idx]
      neighbor_atomic_num[[p]] <- sort(an)
    }
  }

  types <- character(n)
  for (p in seq_len(n)) {
    nb <- neighbor_list[[p]]
    if (length(nb$idx) == 0) {
      types[p] <- paste0(atomic_num[p], "[]")
      next
    }
    an_p <- atomic_num[nb$idx]
    nb_idx_sorted <- nb$idx[order(an_p)]

    sng <- character(length(nb_idx_sorted))
    for (qi in seq_along(nb_idx_sorted)) {
      q <- nb_idx_sorted[qi]
      q_an <- atomic_num[q]
      q_neighbors_an <- neighbor_atomic_num[[q]]
      rm_pos <- match(atomic_num[p], q_neighbors_an)
      remaining <- if (!is.na(rm_pos)) q_neighbors_an[-rm_pos] else q_neighbors_an
      part <- paste0("(", paste(remaining, collapse = ","), ")")
      sng[qi] <- paste0(sprintf("%02d", q_an), "-", part)
    }
    sng <- sort(sng)
    sng <- sub("^0", "", sng)
    types[p] <- paste0(atomic_num[p], "[", paste(sng, collapse = ","), "]")
  }
  types
}
