# Two features ported from the SAVESTEPS pipeline, reimplemented to work
# directly off this app's in-memory tables instead of the original's
# intermediate text-file format:
#
# 1. bond_type / angle_type / dihedral_type: canonical labels built from the
#    Chen-Manz atom types (see atom_typing.R), the same "type" concept used
#    by bonded_atoms_type.py / angle_atoms_type.py / finding_dihedral_materials.py.
#    Since this app works from a single already-optimized geometry (not an
#    ensemble of DFT/AIMD snapshots), each instance's own measured
#    distance/angle/dihedral value IS its equilibrium value -- there's no
#    redundant-calculation problem to solve, so unlike
#    generate_forcefield_inputs_correct_type_bond_angle_dihedrals_DFT.py this
#    does not further split a type by numeric-value tolerance. The per-type
#    summary tables just report the count and value range (mean/min/max)
#    actually observed among that type's instances.
#
# 2. Rotatable / Nonrotatable / Linear dihedral-type classification, ported
#    from Finding_cycles_in_materials_05_13_2024.py: a dihedral's middle
#    bond is "in a ring" if there's some OTHER path (through the periodic
#    bond graph) connecting its two atoms besides the direct bond itself.
#    The original does this with a bidirectional BFS bounded by
#    `2 * number_of_atoms` hops; here it's a single-sided BFS bounded by a
#    fixed, generous hop count instead. A depth bound scaled by atom count
#    breaks down for the tiny hand-built unit cells this app ships as
#    examples (e.g. diamond's 2-atom primitive cell), so a fixed cap that
#    comfortably covers any chemically plausible ring is used instead, with
#    a large node-visit safety cap to bound worst-case runtime. Exceeding
#    the cap without finding a cycle is treated as "not a ring" -- the same
#    fallback behavior as the original's depth-bounded search.

canonical_pair <- function(a, b) {
  s <- sort(c(a, b))
  paste(s, collapse = " :: ")
}

add_type_columns <- function(bonds, angles, dihedrals, atom_types) {
  if (!is.null(bonds) && nrow(bonds) > 0) {
    bonds$bond_type <- mapply(canonical_pair, atom_types[bonds$atom1], atom_types[bonds$atom2])
  }

  if (!is.null(angles) && nrow(angles) > 0) {
    corners <- mapply(canonical_pair, atom_types[angles$atom1], atom_types[angles$atom3])
    angles$angle_type <- paste0(atom_types[angles$atom2], " | ", corners)
  }

  if (!is.null(dihedrals) && nrow(dihedrals) > 0) {
    t1 <- atom_types[dihedrals$atom1]; t2 <- atom_types[dihedrals$atom2]
    t3 <- atom_types[dihedrals$atom3]; t4 <- atom_types[dihedrals$atom4]
    fwd <- paste(t1, t2, t3, t4, sep = " - ")
    rev <- paste(t4, t3, t2, t1, sep = " - ")
    dihedrals$dihedral_type <- ifelse(fwd <= rev, fwd, rev)
  }

  list(bonds = bonds, angles = angles, dihedrals = dihedrals)
}

# BFS: is (atom_k, img_jk) reachable from (atom_j, [0,0,0]) via any path
# OTHER than the direct bond itself? TRUE => that bond lies on a ring/cycle.
is_bond_in_ring <- function(atom_j, atom_k, img_jk, neighbor_list,
                             max_hops = 50, max_visited = 20000) {
  key_of <- function(atom, img) paste(atom, img[1], img[2], img[3], sep = "_")
  target_key <- key_of(atom_k, img_jk)

  visited <- new.env(hash = TRUE, parent = emptyenv())
  assign(key_of(atom_j, c(0, 0, 0)), TRUE, envir = visited)

  queue_atom <- atom_j; queue_img <- matrix(c(0, 0, 0), nrow = 1); queue_hop <- 0L
  qi <- 1L; n_visited <- 1L
  is_origin <- TRUE

  while (qi <= length(queue_atom)) {
    cur_atom <- queue_atom[qi]; cur_img <- queue_img[qi, ]; cur_hop <- queue_hop[qi]
    qi <- qi + 1L

    nb <- neighbor_list[[cur_atom]]
    if (length(nb$idx) > 0 && cur_hop < max_hops) {
      for (m in seq_along(nb$idx)) {
        nb_atom <- nb$idx[m]
        nb_img <- nb$image[m, ] + cur_img

        # never traverse the exact direct edge being tested, on the first hop
        if (is_origin && nb_atom == atom_k && all(nb_img == img_jk)) next

        k <- key_of(nb_atom, nb_img)
        if (!exists(k, envir = visited, inherits = FALSE)) {
          if (k == target_key) return(TRUE)
          assign(k, TRUE, envir = visited)
          n_visited <- n_visited + 1L
          if (n_visited > max_visited) return(FALSE)
          queue_atom <- c(queue_atom, nb_atom)
          queue_img <- rbind(queue_img, nb_img)
          queue_hop <- c(queue_hop, cur_hop + 1L)
        }
      }
    }
    is_origin <- FALSE
  }
  FALSE
}

# Adds in_ring / is_linear / rotatability columns to the dihedrals table.
# Classification is done per dihedral_type (so add_type_columns() must be
# called first), following Finding_cycles_in_materials' rule: a type is
# "Nonrotatable" if ANY of its instances sits on a ring, "Rotatable" only if
# ALL its non-linear instances are bridges, "Linear" if a flanking angle is
# always within ~1.72 deg (0.03 rad) of 180 deg.
classify_dihedral_rotatability <- function(dihedrals, neighbor_list) {
  n <- nrow(dihedrals)
  if (is.null(dihedrals) || n == 0) return(dihedrals)

  linear_tol_deg <- 0.03 * 180 / pi # ~1.72 deg, matches the original's 0.03 rad tolerance
  is_linear <- (abs(abs(dihedrals$flank_angle1_deg) - 180) < linear_tol_deg) |
               (abs(abs(dihedrals$flank_angle2_deg) - 180) < linear_tol_deg)

  # memoize ring lookups per unique middle bond -- many dihedral instances
  # (e.g. all H-C-C-H torsions in ethane) share the same central bond
  ring_cache <- new.env(hash = TRUE, parent = emptyenv())
  in_ring <- logical(n)
  for (r in seq_len(n)) {
    if (is_linear[r]) { in_ring[r] <- NA; next }
    j <- dihedrals$atom2[r]; k <- dihedrals$atom3[r]
    img_jk <- c(dihedrals$img3_a[r], dihedrals$img3_b[r], dihedrals$img3_c[r])
    ck <- paste(j, k, img_jk[1], img_jk[2], img_jk[3], sep = "_")
    if (!exists(ck, envir = ring_cache, inherits = FALSE)) {
      assign(ck, is_bond_in_ring(j, k, img_jk, neighbor_list), envir = ring_cache)
    }
    in_ring[r] <- get(ck, envir = ring_cache, inherits = FALSE)
  }
  dihedrals$in_ring <- in_ring
  dihedrals$is_linear <- is_linear

  rotatability <- character(n)
  for (ty in unique(dihedrals$dihedral_type)) {
    idx <- which(dihedrals$dihedral_type == ty)
    c1 <- 0; c2 <- 0
    for (r in idx) {
      if (isTRUE(is_linear[r])) {
        c1 <- -1; c2 <- -1
      } else {
        if (isFALSE(in_ring[r])) c1 <- 1
        if (isTRUE(in_ring[r])) c2 <- 1
      }
    }
    label <- if (c1 == 1 && c2 == 1) "Nonrotatable"
      else if (c1 == 0 && c2 == 1) "Nonrotatable"
      else if (c1 == 1 && c2 == 0) "Rotatable"
      else if (c1 == -1 && c2 == -1) "Linear"
      else "Unknown"
    rotatability[idx] <- label
  }
  dihedrals$rotatability <- rotatability

  dihedrals
}

# Generic "instances -> unique types" summary: count + value range per type.
# The returned "Index" column is the row number in this summary (sorted by
# descending count) -- e.g. Index 1 is whichever type has the most
# instances. Instance tables carry a matching *_type_index column (see
# attach_type_index()) so a row in the Bonds/Angles/Dihedrals table can be
# cross-referenced directly against its Bond/Angle/Dihedral Type row.
summarize_types <- function(df, type_col, value_col, extra_cols = character(0)) {
  if (is.null(df) || nrow(df) == 0) return(data.frame())
  types <- unique(df[[type_col]])
  rows <- lapply(types, function(ty) {
    sub <- df[df[[type_col]] == ty, , drop = FALSE]
    vals <- sub[[value_col]]
    base <- data.frame(
      Type = ty,
      Count = nrow(sub),
      Mean = round(mean(vals), 3),
      Min = round(min(vals), 3),
      Max = round(max(vals), 3)
    )
    for (ec in extra_cols) base[[ec]] <- sub[[ec]][1]
    base
  })
  out <- do.call(rbind, rows)
  out <- out[order(-out$Count), ]
  out$Index <- seq_len(nrow(out))
  rownames(out) <- NULL
  out
}

# Maps each instance's type string to the row number ("Index") of that type
# in its corresponding summary table (from summarize_types()).
attach_type_index <- function(df, type_col, summary_df, index_col_name) {
  if (is.null(df) || nrow(df) == 0) return(df)
  df[[index_col_name]] <- summary_df$Index[match(df[[type_col]], summary_df$Type)]
  df
}
