# 3- and 4-membered ring detection, ported from
# finding_3_member_rings_1_30_2024.py / finding_4_member_rings_1_30_2024.py,
# and their downstream use in generate_bonds_materials_1_30_2024.py /
# generate_angles_materials_1_30_2024.py, matching the paper's Step 4:
#   (1) Urey-Bradley stretches are added for the two diagonals of each
#       4-membered ring.
#   (2) Angles in 3- and 4-membered rings are flagged and excluded from the
#       "active" angle set (angle-bending potential in the paper; here,
#       excluded from the Angle Types summary / type-index).
#   (3) Dihedrals containing an angle from a 3- or 4-membered ring are
#       removed from the active internal coordinate list (here: excluded
#       from the Dihedral Types summary, dihedral-type numeric refinement,
#       and pruning -- but still shown, flagged, in the raw Dihedrals table
#       for transparency).
#
# A 3-membered ring: an angle instance (corner1-center-corner2) where
# corner1 and corner2 are ALSO directly bonded to each other.
#
# A 4-membered ring: two DIFFERENT angle instances that share the same pair
# of corner atoms but have different center atoms -- i.e. corner_lo -
# center_A - corner_hi - center_B - corner_lo. Detected by comparing, for
# each pair of angles sharing a corner-atom pair, whether the periodic
# image offset between the corners (measured through each angle's own
# center) agrees -- this confirms it's a genuine closed ring, not a
# coincidental atom-number match through unrelated periodic copies.

find_3_member_rings <- function(angles, neighbor_list) {
  n <- nrow(angles)
  if (n == 0) return(logical(0))
  is_3ring <- logical(n)
  for (r in seq_len(n)) {
    a1 <- angles$atom1[r]; a3 <- angles$atom3[r]
    needed_img <- c(angles$img3_a[r] - angles$img1_a[r],
                     angles$img3_b[r] - angles$img1_b[r],
                     angles$img3_c[r] - angles$img1_c[r])
    nb <- neighbor_list[[a1]]
    if (length(nb$idx) > 0) {
      is_3ring[r] <- any(nb$idx == a3 & apply(nb$image, 1, function(x) all(x == needed_img)))
    }
  }
  is_3ring
}

find_4_member_rings <- function(angles) {
  n <- nrow(angles)
  empty_rings <- data.frame(
    corner_lo = integer(0), corner_hi = integer(0),
    corner_off_a = integer(0), corner_off_b = integer(0), corner_off_c = integer(0),
    center1 = integer(0), center1_row = integer(0),
    center2 = integer(0), center2_row = integer(0),
    center_off_a = integer(0), center_off_b = integer(0), center_off_c = integer(0)
  )
  if (n < 2) return(list(is_4ring = logical(n), rings = empty_rings))

  lo <- pmin(angles$atom1, angles$atom3)
  hi <- pmax(angles$atom1, angles$atom3)
  hi_is_atom1 <- angles$atom1 == hi
  hi_img <- cbind(ifelse(hi_is_atom1, angles$img1_a, angles$img3_a),
                   ifelse(hi_is_atom1, angles$img1_b, angles$img3_b),
                   ifelse(hi_is_atom1, angles$img1_c, angles$img3_c))
  lo_img <- cbind(ifelse(hi_is_atom1, angles$img3_a, angles$img1_a),
                   ifelse(hi_is_atom1, angles$img3_b, angles$img1_b),
                   ifelse(hi_is_atom1, angles$img3_c, angles$img1_c))
  corner_off <- hi_img - lo_img # hi's image relative to lo, via this angle's own center

  corner_key <- paste(lo, hi)
  offset_key <- paste(corner_off[, 1], corner_off[, 2], corner_off[, 3])

  is_4ring <- logical(n)
  ring_rows <- list()

  for (ck in unique(corner_key)) {
    idx <- which(corner_key == ck)
    if (length(idx) < 2) next
    for (a in seq_along(idx)) {
      for (b in seq_along(idx)) {
        if (a >= b) next
        ri <- idx[a]; rj <- idx[b]
        if (angles$atom2[ri] == angles$atom2[rj]) next # same center: not a 4-ring pair
        if (offset_key[ri] != offset_key[rj]) next
        is_4ring[ri] <- TRUE
        is_4ring[rj] <- TRUE
        center_off <- hi_img[ri, ] - hi_img[rj, ] # center_rj relative to center_ri, via corner_hi bridge
        ring_rows[[length(ring_rows) + 1]] <- data.frame(
          corner_lo = lo[ri], corner_hi = hi[ri],
          corner_off_a = corner_off[ri, 1], corner_off_b = corner_off[ri, 2], corner_off_c = corner_off[ri, 3],
          center1 = angles$atom2[ri], center1_row = ri,
          center2 = angles$atom2[rj], center2_row = rj,
          center_off_a = center_off[1], center_off_b = center_off[2], center_off_c = center_off[3]
        )
      }
    }
  }
  rings <- if (length(ring_rows) > 0) do.call(rbind, ring_rows) else empty_rings
  list(is_4ring = is_4ring, rings = rings)
}

# For each dihedral instance, checks whether EITHER flanking angle
# (i-j-k or j-k-l) is part of a 3- or 4-membered ring.
compute_dihedral_ring_exclusion <- function(dihedrals, angles) {
  n <- nrow(dihedrals)
  if (n == 0) return(logical(0))

  lookup <- new.env(hash = TRUE, parent = emptyenv())
  if (nrow(angles) > 0) {
    keys <- paste(angles$atom2, pmin(angles$atom1, angles$atom3), pmax(angles$atom1, angles$atom3))
    for (i in seq_len(nrow(angles))) {
      if (!exists(keys[i], envir = lookup, inherits = FALSE)) {
        assign(keys[i], angles$ring_flag[i], envir = lookup)
      }
    }
  }
  get_ring_flag <- function(center, a, b) {
    k <- paste(center, min(a, b), max(a, b))
    if (exists(k, envir = lookup, inherits = FALSE)) get(k, envir = lookup, inherits = FALSE) else 0L
  }

  flag1 <- vapply(seq_len(n), function(r) get_ring_flag(dihedrals$atom2[r], dihedrals$atom1[r], dihedrals$atom3[r]), integer(1))
  flag2 <- vapply(seq_len(n), function(r) get_ring_flag(dihedrals$atom3[r], dihedrals$atom2[r], dihedrals$atom4[r]), integer(1))
  (flag1 != 0) | (flag2 != 0)
}

# Urey-Bradley diagonal "stretches" for each detected 4-membered ring: the
# two diagonals corner_lo--corner_hi and center1--center2. Diagonals that
# happen to already be real bonds are dropped (matching the original: "adds
# the diagonal bonds ... if those are not listed in the stretch list").
# These are NOT real chemical bonds -- not drawn in the 3D viewer, not
# counted in the Bonds tab.
build_ub_bonds <- function(rings, bonds, coords, lattice, elements, atom_types) {
  if (is.null(rings) || nrow(rings) == 0) return(data.frame())

  d1 <- data.frame(atom1 = rings$corner_lo, atom2 = rings$corner_hi,
                    img_a = rings$corner_off_a, img_b = rings$corner_off_b, img_c = rings$corner_off_c)
  d2 <- data.frame(atom1 = rings$center1, atom2 = rings$center2,
                    img_a = rings$center_off_a, img_b = rings$center_off_b, img_c = rings$center_off_c)
  cand <- rbind(d1, d2)

  canon <- function(df) {
    swap <- df$atom1 > df$atom2
    data.frame(
      atom1 = ifelse(swap, df$atom2, df$atom1),
      atom2 = ifelse(swap, df$atom1, df$atom2),
      img_a = ifelse(swap, -df$img_a, df$img_a),
      img_b = ifelse(swap, -df$img_b, df$img_b),
      img_c = ifelse(swap, -df$img_c, df$img_c)
    )
  }
  cand <- unique(canon(cand))

  if (nrow(bonds) > 0) {
    real <- canon(bonds[, c("atom1", "atom2", "img_a", "img_b", "img_c")])
    real_key <- paste(real$atom1, real$atom2, real$img_a, real$img_b, real$img_c)
    cand_key <- paste(cand$atom1, cand$atom2, cand$img_a, cand$img_b, cand$img_c)
    cand <- cand[!(cand_key %in% real_key), , drop = FALSE]
  }
  if (nrow(cand) == 0) return(data.frame())

  dist <- numeric(nrow(cand))
  for (i in seq_len(nrow(cand))) {
    p1 <- coords[cand$atom1[i], ]
    img <- c(cand$img_a[i], cand$img_b[i], cand$img_c[i])
    p2 <- coords[cand$atom2[i], ] + (if (!is.null(lattice)) as.numeric(img %*% lattice) else c(0, 0, 0))
    dist[i] <- sqrt(sum((p2 - p1)^2))
  }
  cand$distance <- round(dist, 4)
  cand$element1 <- elements[cand$atom1]
  cand$element2 <- elements[cand$atom2]
  cand$bond <- paste0(cand$element1, "-", cand$element2)
  cand$type1 <- atom_types[cand$atom1]
  cand$type2 <- atom_types[cand$atom2]
  cand$ub_type <- mapply(canonical_pair, cand$type1, cand$type2)
  cand$image1 <- image_str(0L, 0L, 0L)
  cand$image2 <- image_str(cand$img_a, cand$img_b, cand$img_c)
  cand
}
