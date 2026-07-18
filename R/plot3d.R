# Interactive 3D structure viewer built with rgl (ball-and-stick style):
# atoms as CPK-colored spheres, bonds as solid cylinders (drawn to the true
# bonded periodic image, even if that puts the far end outside the drawn
# cell), and an optional unit-cell wireframe. Returns an htmlwidget
# (rgl::rglwidget()) suitable for rgl::renderRglwidget() in Shiny.

build_structure_widget <- function(analysis, show_cell = TRUE, show_labels = FALSE,
                                    bond_radius = 0.12, atom_scale = 0.28) {
  elements <- analysis$elements
  coords <- analysis$coords
  lattice <- analysis$lattice
  bonds <- analysis$bonds

  radii <- unname(ATOM_TYPING_RADII[elements])
  radii[is.na(radii)] <- 1.5
  sphere_radii <- radii * atom_scale
  cols <- vapply(elements, element_color, character(1))

  rgl::open3d(useNULL = TRUE)

  result <- tryCatch({
    rgl::bg3d(color = "white")

    rgl::spheres3d(coords, radius = sphere_radii, color = cols, specular = "black")

    if (!is.null(bonds) && nrow(bonds) > 0) {
      for (k in seq_len(nrow(bonds))) {
        i <- bonds$atom1[k]; j <- bonds$atom2[k]
        img <- c(bonds$img_a[k], bonds$img_b[k], bonds$img_c[k])
        p1 <- coords[i, ]
        p2 <- coords[j, ] + (if (!is.null(lattice)) as.numeric(img %*% lattice) else c(0, 0, 0))
        if (sqrt(sum((p2 - p1)^2)) < 1e-6) next
        cyl <- rgl::cylinder3d(rbind(p1, p2), radius = bond_radius, sides = 12)
        rgl::shade3d(cyl, color = "#4d4d4d", specular = "black")
      }
    }

    if (show_cell && !is.null(lattice)) {
      origin <- c(0, 0, 0)
      a_v <- lattice[1, ]; b_v <- lattice[2, ]; c_v <- lattice[3, ]
      verts <- list(origin, a_v, b_v, c_v, a_v + b_v, a_v + c_v, b_v + c_v, a_v + b_v + c_v)
      edges <- list(
        c(1, 2), c(1, 3), c(1, 4), c(2, 5), c(2, 6), c(3, 5),
        c(3, 7), c(4, 6), c(4, 7), c(5, 8), c(6, 8), c(7, 8)
      )
      seg <- do.call(rbind, lapply(edges, function(e) rbind(verts[[e[1]]], verts[[e[2]]])))
      rgl::segments3d(seg, color = "#8888cc", lwd = 1.5)
    }

    if (show_labels) {
      rgl::text3d(coords[, 1], coords[, 2], coords[, 3] + sphere_radii + 0.2,
                  texts = paste0(elements, seq_along(elements)), color = "black", cex = 0.8)
    }

    rgl::aspect3d("iso")
    rgl::view3d(theta = 30, phi = 20, zoom = 0.8)

    rgl::rglwidget()
  }, finally = {
    rgl::close3d()
  })

  result
}
