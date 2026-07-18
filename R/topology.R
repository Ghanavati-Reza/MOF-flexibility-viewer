# Bond / angle / dihedral perception for a (possibly periodic) structure.
#
# Bonding rule: two atoms are bonded if the distance between atom i and any
# of the 27 nearest periodic images of atom j (offsets -1..1 along each
# lattice vector) is <= (radius_i + radius_j) * bond_tolerance. Every
# qualifying image is kept as its own bond (not just the closest one) --
# this matters for symmetric/small-cell motifs (common in MOFs) where an
# atom is genuinely bonded to two different periodic images of the same
# neighbor at equal distance; keeping only the single nearest image would
# silently drop one of the two real bonds. For non-periodic input
# (lattice = NULL) this reduces to ordinary distance-based bonding.
#
# Angles: every pair of bonds sharing a central atom.
# Dihedrals: for every bond j-k, every neighbor i of j (other than k) paired
# with every neighbor l of k (other than j) -> dihedral i-j-k-l. Each bond
# is visited once (atom1 < atom2), so each physical dihedral is produced
# exactly once.

find_bonds <- function(elements, coords, lattice, tolerance = 1.0) {
  n <- length(elements)
  if (n < 2) stop("Structure has fewer than 2 atoms; cannot compute bonds.")

  radii <- ATOM_TYPING_RADII[elements]
  missing_el <- unique(elements[is.na(radii)])
  if (length(missing_el) > 0) {
    stop("Unrecognized element symbol(s): ", paste(missing_el, collapse = ", "),
         " -- cannot look up bonding radii.")
  }

  if (!is.null(lattice)) {
    grid <- expand.grid(na = -1:1, nb = -1:1, nc = -1:1)
    images <- as.matrix(grid[, c("na", "nb", "nc")])
    image_shift <- images %*% lattice
  } else {
    images <- matrix(c(0, 0, 0), nrow = 1)
    image_shift <- matrix(c(0, 0, 0), nrow = 1)
  }
  n_img <- nrow(images)

  bond_atom1 <- integer(0); bond_atom2 <- integer(0)
  bond_image <- matrix(nrow = 0, ncol = 3)
  bond_dist <- numeric(0)

  for (i in seq_len(n - 1)) {
    js <- (i + 1):n
    base_j <- coords[js, , drop = FALSE]
    cutoff <- (radii[i] + radii[js]) * tolerance
    for (im in seq_len(n_img)) {
      shifted <- sweep(base_j, 2, image_shift[im, ], "+")
      d <- sqrt(rowSums(sweep(shifted, 2, coords[i, ], "-")^2))
      hit <- which(d <= cutoff & d > 1e-6)
      if (length(hit) > 0) {
        bond_atom1 <- c(bond_atom1, rep(i, length(hit)))
        bond_atom2 <- c(bond_atom2, js[hit])
        bond_image <- rbind(bond_image, matrix(images[im, ], nrow = length(hit), ncol = 3, byrow = TRUE))
        bond_dist <- c(bond_dist, d[hit])
      }
    }
  }

  data.frame(
    atom1 = bond_atom1, atom2 = bond_atom2,
    img_a = bond_image[, 1], img_b = bond_image[, 2], img_c = bond_image[, 3],
    distance = round(bond_dist, 4)
  )
}

build_neighbor_list <- function(bonds, n) {
  neighbor_list <- vector("list", n)
  for (i in seq_len(n)) {
    neighbor_list[[i]] <- list(idx = integer(0), image = matrix(nrow = 0, ncol = 3), dist = numeric(0))
  }
  if (nrow(bonds) == 0) return(neighbor_list)
  for (k in seq_len(nrow(bonds))) {
    i <- bonds$atom1[k]; j <- bonds$atom2[k]
    img <- c(bonds$img_a[k], bonds$img_b[k], bonds$img_c[k])
    neighbor_list[[i]]$idx <- c(neighbor_list[[i]]$idx, j)
    neighbor_list[[i]]$image <- rbind(neighbor_list[[i]]$image, img)
    neighbor_list[[i]]$dist <- c(neighbor_list[[i]]$dist, bonds$distance[k])
    neighbor_list[[j]]$idx <- c(neighbor_list[[j]]$idx, i)
    neighbor_list[[j]]$image <- rbind(neighbor_list[[j]]$image, -img)
    neighbor_list[[j]]$dist <- c(neighbor_list[[j]]$dist, bonds$distance[k])
  }
  neighbor_list
}

image_str <- function(a, b, c) sprintf("(%d,%d,%d)", a, b, c)

find_angles <- function(coords, lattice, neighbor_list) {
  rows <- list()
  n <- length(neighbor_list)
  for (j in seq_len(n)) {
    nb <- neighbor_list[[j]]
    m <- length(nb$idx)
    if (m < 2) next
    pos_j <- coords[j, ]
    for (a in seq_len(m - 1)) {
      for (b in (a + 1):m) {
        k_idx <- nb$idx[a]; k_img <- nb$image[a, ]
        l_idx <- nb$idx[b]; l_img <- nb$image[b, ]
        pos_k <- coords[k_idx, ] + (if (!is.null(lattice)) as.numeric(k_img %*% lattice) else c(0, 0, 0))
        pos_l <- coords[l_idx, ] + (if (!is.null(lattice)) as.numeric(l_img %*% lattice) else c(0, 0, 0))
        if (sqrt(sum((pos_k - pos_l)^2)) < 1e-6) next
        ang <- vec_angle_deg(pos_k, pos_j, pos_l)
        rows[[length(rows) + 1]] <- data.frame(
          atom1 = k_idx, img1_a = k_img[1], img1_b = k_img[2], img1_c = k_img[3],
          atom2 = j,
          atom3 = l_idx, img3_a = l_img[1], img3_b = l_img[2], img3_c = l_img[3],
          angle_deg = round(ang, 2)
        )
      }
    }
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

find_dihedrals <- function(coords, lattice, bonds, neighbor_list) {
  rows <- list()
  zero3 <- c(0, 0, 0)
  if (nrow(bonds) == 0) return(data.frame())

  for (bi in seq_len(nrow(bonds))) {
    j <- bonds$atom1[bi]; k <- bonds$atom2[bi]
    img_jk <- c(bonds$img_a[bi], bonds$img_b[bi], bonds$img_c[bi])
    nb_j <- neighbor_list[[j]]; nb_k <- neighbor_list[[k]]

    # exclude the neighbor-list entry that IS this central bond itself
    # (bonds df has at most one row per atom pair, so at most one match each side)
    j_back <- if (length(nb_j$idx) == 0) logical(0) else
      (nb_j$idx == k & apply(nb_j$image, 1, function(r) all(r == img_jk)))
    k_back <- if (length(nb_k$idx) == 0) logical(0) else
      (nb_k$idx == j & apply(nb_k$image, 1, function(r) all(r == -img_jk)))
    j_keep <- !j_back
    k_keep <- !k_back
    if (!any(j_keep) || !any(k_keep)) next

    i_idx_list <- nb_j$idx[j_keep]; i_img_list <- nb_j$image[j_keep, , drop = FALSE]
    l_idx_list <- nb_k$idx[k_keep]; l_img_list <- nb_k$image[k_keep, , drop = FALSE]

    pos_j <- coords[j, ]
    pos_k <- coords[k, ] + (if (!is.null(lattice)) as.numeric(img_jk %*% lattice) else zero3)

    for (a in seq_along(i_idx_list)) {
      i_idx <- i_idx_list[a]; img_i <- i_img_list[a, ]
      pos_i <- coords[i_idx, ] + (if (!is.null(lattice)) as.numeric(img_i %*% lattice) else zero3)

      for (b in seq_along(l_idx_list)) {
        l_idx <- l_idx_list[b]; img_l_rel_k <- l_img_list[b, ]
        img_l <- img_jk + img_l_rel_k
        pos_l <- coords[l_idx, ] + (if (!is.null(lattice)) as.numeric(img_l %*% lattice) else zero3)

        pts <- rbind(pos_i, pos_j, pos_k, pos_l)
        dmat <- as.matrix(dist(pts))
        if (any(dmat[upper.tri(dmat)] < 1e-6)) next

        dh <- tryCatch(dihedral_deg(pos_i, pos_j, pos_k, pos_l), error = function(e) NA_real_)
        if (is.na(dh)) next

        # flanking bond angles (i-j-k and j-k-l): needed to flag "linear"
        # dihedrals where a flanking angle is near 180 deg and the torsion
        # is geometrically ill-defined.
        flank1 <- vec_angle_deg(pos_i, pos_j, pos_k)
        flank2 <- vec_angle_deg(pos_j, pos_k, pos_l)

        rows[[length(rows) + 1]] <- data.frame(
          atom1 = i_idx, img1_a = img_i[1], img1_b = img_i[2], img1_c = img_i[3],
          atom2 = j,     img2_a = 0L, img2_b = 0L, img2_c = 0L,
          atom3 = k,     img3_a = img_jk[1], img3_b = img_jk[2], img3_c = img_jk[3],
          atom4 = l_idx, img4_a = img_l[1], img4_b = img_l[2], img4_c = img_l[3],
          dihedral_deg = round(dh, 2),
          flank_angle1_deg = round(flank1, 2),
          flank_angle2_deg = round(flank2, 2)
        )
      }
    }
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

# Runs the full pipeline and adds human-readable element / Chen-Manz
# atom-type label columns to each table.
analyze_structure <- function(struct, tolerance = 1.0) {
  elements <- struct$elements
  coords <- struct$coords
  lattice <- struct$lattice

  bonds <- find_bonds(elements, coords, lattice, tolerance)
  neighbor_list <- build_neighbor_list(bonds, length(elements))
  angles <- find_angles(coords, lattice, neighbor_list)
  dihedrals <- find_dihedrals(coords, lattice, bonds, neighbor_list)
  atom_types <- chen_manz_atom_types(elements, neighbor_list)

  if (nrow(bonds) > 0) {
    bonds$element1 <- elements[bonds$atom1]
    bonds$element2 <- elements[bonds$atom2]
    bonds$bond <- paste0(bonds$element1, "-", bonds$element2)
    bonds$type1 <- atom_types[bonds$atom1]
    bonds$type2 <- atom_types[bonds$atom2]
    # atom1 is always the reference atom (home image [0,0,0]); atom2's image
    # is relative to it. Both are shown so every atom's image sits right
    # after its element, even when it's trivially [0,0,0].
    bonds$image1 <- image_str(0L, 0L, 0L)
    bonds$image2 <- image_str(bonds$img_a, bonds$img_b, bonds$img_c)
  }

  if (nrow(angles) > 0) {
    angles$element1 <- elements[angles$atom1]
    angles$element2 <- elements[angles$atom2]
    angles$element3 <- elements[angles$atom3]
    angles$angle_label <- paste0(angles$element1, "-", angles$element2, "-", angles$element3)
    angles$image1 <- image_str(angles$img1_a, angles$img1_b, angles$img1_c)
    angles$image2 <- image_str(0L, 0L, 0L) # center atom is always the home reference
    angles$image3 <- image_str(angles$img3_a, angles$img3_b, angles$img3_c)
  }

  if (nrow(dihedrals) > 0) {
    dihedrals$element1 <- elements[dihedrals$atom1]
    dihedrals$element2 <- elements[dihedrals$atom2]
    dihedrals$element3 <- elements[dihedrals$atom3]
    dihedrals$element4 <- elements[dihedrals$atom4]
    dihedrals$dihedral_label <- paste0(dihedrals$element1, "-", dihedrals$element2, "-",
                                        dihedrals$element3, "-", dihedrals$element4)
    dihedrals$image1 <- image_str(dihedrals$img1_a, dihedrals$img1_b, dihedrals$img1_c)
    dihedrals$image2 <- image_str(0L, 0L, 0L) # atom2 is always the home reference
    dihedrals$image3 <- image_str(dihedrals$img3_a, dihedrals$img3_b, dihedrals$img3_c)
    dihedrals$image4 <- image_str(dihedrals$img4_a, dihedrals$img4_b, dihedrals$img4_c)
  }

  typed <- add_type_columns(bonds, angles, dihedrals, atom_types)
  bonds <- typed$bonds; angles <- typed$angles; dihedrals <- typed$dihedrals

  # Numeric-value refinement: split each Chen-Manz type into finer
  # sub-types by bond length / angle value / dihedral value (see
  # numeric_subtyping.R). Must run in this order (bonds -> angles ->
  # dihedrals) since each level looks up its flanking lower-level
  # sub-types, and must run before rotatability classification and
  # pruning so those operate on the refined types -- this is the actual
  # step order in the original pipeline.
  bonds <- split_bonds_by_length(bonds)
  angles <- split_angles_by_value(angles, bonds)
  dihedrals <- split_dihedrals_by_value(dihedrals, angles)

  dihedrals <- classify_dihedral_rotatability(dihedrals, neighbor_list)

  bond_types_summary <- summarize_types(bonds, "bond_type", "distance", extra_cols = "bond")
  angle_types_summary <- summarize_types(angles, "angle_type", "angle_deg", extra_cols = "angle_label")
  dihedral_types_summary <- summarize_types(dihedrals, "dihedral_type", "dihedral_deg",
                                             extra_cols = c("dihedral_label", "rotatability"))

  bonds <- attach_type_index(bonds, "bond_type", bond_types_summary, "bond_type_index")
  angles <- attach_type_index(angles, "angle_type", angle_types_summary, "angle_type_index")
  dihedrals <- attach_type_index(dihedrals, "dihedral_type", dihedral_types_summary, "dihedral_type_index")

  # group instances by type (Type # ascending) so all bond_type 1 rows come
  # first, then all bond_type 2 rows, etc. -- same for angles/dihedrals
  if (!is.null(bonds) && nrow(bonds) > 0) bonds <- bonds[order(bonds$bond_type_index), ]
  if (!is.null(angles) && nrow(angles) > 0) angles <- angles[order(angles$angle_type_index), ]
  if (!is.null(dihedrals) && nrow(dihedrals) > 0) dihedrals <- dihedrals[order(dihedrals$dihedral_type_index), ]

  pruned <- apply_dihedral_pruning(dihedrals, dihedral_types_summary)

  list(
    elements = elements, coords = coords, lattice = lattice, name = struct$name,
    atom_types = atom_types, bonds = bonds, angles = angles, dihedrals = dihedrals,
    neighbor_list = neighbor_list,
    bond_types_summary = bond_types_summary,
    angle_types_summary = angle_types_summary,
    dihedral_types_summary = dihedral_types_summary,
    dihedrals_pruned = pruned$dihedrals_pruned,
    dihedral_types_summary_pruned = pruned$dihedral_types_summary_pruned
  )
}
