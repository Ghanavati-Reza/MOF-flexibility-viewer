# Rotatable -> Hindered reclassification, ported (without geometry-file
# generation) from
# code_for_generating_geoms_for_rotatable_dihedrals_materials_05_13_2024.py.
#
# For each dihedral TYPE currently classified "Rotatable": take one
# representative instance, find the smaller of the two branches hanging off
# the middle bond (via BFS, excluding the bond itself), and rotate that
# branch around the middle-bond axis in 10-degree increments from -180 to
# +170 (36 positions, matching the original). After each rotation, re-type
# every atom (same Chen-Manz scheme) and compare to the original typing. If
# ANY rotation changes ANY atom's type -- i.e. the rotation would break or
# remake a bond -- the type is reclassified "Hindered" and we move on to
# the next rotatable type (no need to check the remaining angles once one
# has failed). This app does not generate/save the displaced geometries
# themselves, only performs the connectivity check.
#
# A branch that turns out to be part of the (possibly infinite) periodic
# framework rather than a finite pendant group would make the BFS run
# forever; a node-count cap guards against that, and such a dihedral type
# is left as "Rotatable" (inconclusive) rather than crashing or hanging.
#
# Performance: only the (usually small) rotated fragment's coordinates
# change between the optimized and displaced geometries, so bonds between
# two UNCHANGED atoms can't possibly differ -- only bonds touching a
# fragment atom need to be recomputed. find_bonds_incremental() does that
# instead of a full O(n^2) re-search, which matters here since this check
# runs up to 36 times per rotatable dihedral type.

find_bonds_incremental <- function(elements, coords, lattice, tolerance, moved_atoms, original_bonds) {
  n <- length(elements)
  unchanged <- original_bonds[!(original_bonds$atom1 %in% moved_atoms | original_bonds$atom2 %in% moved_atoms), , drop = FALSE]

  radii <- ATOM_TYPING_RADII[elements]
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
  bond_image <- matrix(nrow = 0, ncol = 3); bond_dist <- numeric(0)

  for (i in moved_atoms) {
    others <- setdiff(seq_len(n), i)
    if (length(others) == 0) next
    base_j <- coords[others, , drop = FALSE]
    cutoff <- (radii[i] + radii[others]) * tolerance
    for (im in seq_len(n_img)) {
      shifted <- sweep(base_j, 2, image_shift[im, ], "+")
      d <- sqrt(rowSums(sweep(shifted, 2, coords[i, ], "-")^2))
      hit <- which(d <= cutoff & d > 1e-6)
      if (length(hit) > 0) {
        a1 <- rep(i, length(hit)); a2 <- others[hit]
        img <- matrix(images[im, ], nrow = length(hit), ncol = 3, byrow = TRUE)
        swap <- a1 > a2
        tmp1 <- ifelse(swap, a2, a1); tmp2 <- ifelse(swap, a1, a2)
        img[swap, ] <- -img[swap, , drop = FALSE]
        bond_atom1 <- c(bond_atom1, tmp1); bond_atom2 <- c(bond_atom2, tmp2)
        bond_image <- rbind(bond_image, img); bond_dist <- c(bond_dist, d[hit])
      }
    }
  }
  new_bonds <- data.frame(
    atom1 = bond_atom1, atom2 = bond_atom2,
    img_a = bond_image[, 1], img_b = bond_image[, 2], img_c = bond_image[, 3],
    distance = round(bond_dist, 4)
  )
  new_bonds <- unique(new_bonds) # a pair with both endpoints moved gets computed from each side
  rbind(unchanged, new_bonds)
}

bfs_branch <- function(start_atom, start_image, exclude_atom, exclude_image, neighbor_list, max_nodes = 500) {
  key_of <- function(atom, img) paste(atom, img[1], img[2], img[3], sep = "_")
  visited <- new.env(hash = TRUE, parent = emptyenv())
  assign(key_of(start_atom, start_image), TRUE, envir = visited)

  atoms <- start_atom
  images <- list(start_image)
  queue_atom <- start_atom; queue_img <- matrix(start_image, nrow = 1)
  qi <- 1L; is_origin <- TRUE; overflow <- FALSE

  while (qi <= length(queue_atom)) {
    cur_atom <- queue_atom[qi]; cur_img <- queue_img[qi, ]
    qi <- qi + 1L
    nb <- neighbor_list[[cur_atom]]
    if (length(nb$idx) > 0) {
      for (m in seq_along(nb$idx)) {
        nb_atom <- nb$idx[m]
        nb_img <- nb$image[m, ] + cur_img
        if (is_origin && nb_atom == exclude_atom && all(nb_img == exclude_image)) next
        k <- key_of(nb_atom, nb_img)
        if (!exists(k, envir = visited, inherits = FALSE)) {
          assign(k, TRUE, envir = visited)
          atoms <- c(atoms, nb_atom)
          images[[length(images) + 1]] <- nb_img
          if (length(atoms) > max_nodes) { overflow <- TRUE; break }
          queue_atom <- c(queue_atom, nb_atom)
          queue_img <- rbind(queue_img, nb_img)
        }
      }
    }
    is_origin <- FALSE
    if (overflow) break
  }
  list(atoms = atoms, images = images, overflow = overflow)
}

# Rotates the given fragment atoms (atom numbers + their images, expressed
# relative to the structure's home cell) around `axis_unit` through
# `pivot_pos`, by `rotation_amount` radians (Rodrigues' rotation formula).
# Returns a full coordinate matrix (only the fragment rows differ from the
# input).
rotate_fragment <- function(coords, lattice, frag_atoms, frag_images, pivot_pos, axis_unit, rotation_amount) {
  displaced <- coords
  cos_t <- cos(rotation_amount); sin_t <- sin(rotation_amount)
  zero3 <- c(0, 0, 0)
  for (i in seq_along(frag_atoms)) {
    a <- frag_atoms[i]; img <- frag_images[[i]]
    img_shift <- if (!is.null(lattice)) as.numeric(img %*% lattice) else zero3
    abs_pos <- coords[a, ] + img_shift
    v <- abs_pos - pivot_pos
    k_dot_v <- sum(axis_unit * v)
    k_cross_v <- cross3(axis_unit, v)
    rotated <- v * cos_t + k_cross_v * sin_t + axis_unit * k_dot_v * (1 - cos_t)
    new_abs_pos <- rotated + pivot_pos
    displaced[a, ] <- new_abs_pos - img_shift # back to a home-cell-equivalent coordinate
  }
  displaced
}

# Returns the Type strings (from dihedral_types_summary) that should be
# reclassified "Hindered".
check_hindered_dihedral_types <- function(dihedrals_active, dihedral_types_summary, elements, coords, lattice,
                                           neighbor_list, atom_types, tolerance, bonds) {
  hindered_types <- character(0)
  if (is.null(dihedral_types_summary) || nrow(dihedral_types_summary) == 0) return(hindered_types)
  rotatable_rows <- dihedral_types_summary[dihedral_types_summary$rotatability == "Rotatable", , drop = FALSE]
  if (nrow(rotatable_rows) == 0) return(hindered_types)

  zero3 <- c(0, 0, 0)

  for (ty in rotatable_rows$Type) {
    inst <- dihedrals_active[dihedrals_active$dihedral_type == ty, , drop = FALSE][1, ]
    atom_j <- inst$atom2; atom_k <- inst$atom3
    img_jk <- c(inst$img3_a, inst$img3_b, inst$img3_c)

    left <- bfs_branch(atom_j, zero3, atom_k, img_jk, neighbor_list)
    right <- bfs_branch(atom_k, img_jk, atom_j, zero3, neighbor_list)
    # An overflowing branch just confirms it's the LARGER side (definitely
    # bigger than max_nodes, hence bigger than any finite branch) -- that's
    # still a conclusive comparison, so only skip when BOTH sides overflow
    # (genuinely ambiguous / both reach into the infinite framework).
    if (left$overflow && right$overflow) next

    use_left <- if (left$overflow) FALSE else if (right$overflow) TRUE else (length(left$atoms) <= length(right$atoms))
    if (use_left) {
      pivot_atom <- atom_j; pivot_image <- zero3
      other_atom <- atom_k; other_image <- img_jk
      frag_atoms <- left$atoms[-1]; frag_images <- left$images[-1]
    } else {
      pivot_atom <- atom_k; pivot_image <- img_jk
      other_atom <- atom_j; other_image <- zero3
      frag_atoms <- right$atoms[-1]; frag_images <- right$images[-1]
    }
    if (length(frag_atoms) == 0) next # nothing to rotate

    pivot_shift <- if (!is.null(lattice)) as.numeric(pivot_image %*% lattice) else zero3
    other_shift <- if (!is.null(lattice)) as.numeric(other_image %*% lattice) else zero3
    pivot_pos <- coords[pivot_atom, ] + pivot_shift
    other_pos <- coords[other_atom, ] + other_shift
    axis_vec <- pivot_pos - other_pos
    axis_len <- sqrt(sum(axis_vec^2))
    if (axis_len < 1e-8) next
    axis_unit <- axis_vec / axis_len

    old_dihedral_rad <- inst$dihedral_deg * pi / 180

    is_hindered <- FALSE
    for (new_angle_deg in seq(-180, 170, by = 10)) {
      rotation_amount <- (new_angle_deg * pi / 180) - old_dihedral_rad
      displaced_coords <- rotate_fragment(coords, lattice, frag_atoms, frag_images, pivot_pos, axis_unit, rotation_amount)

      displaced_bonds <- find_bonds_incremental(elements, displaced_coords, lattice, tolerance, frag_atoms, bonds)
      displaced_neighbor_list <- build_neighbor_list(displaced_bonds, length(elements))
      displaced_types <- chen_manz_atom_types(elements, displaced_neighbor_list)

      if (!identical(displaced_types, atom_types)) {
        is_hindered <- TRUE
        break # connectivity changed -- no need to check the remaining angles for this type
      }
    }
    if (is_hindered) hindered_types <- c(hindered_types, ty)
  }
  hindered_types
}
