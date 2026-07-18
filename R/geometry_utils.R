# Basic vector geometry helpers used by topology.R

cell_params_to_vectors <- function(a, b, c, alpha, beta, gamma) {
  # alpha, beta, gamma in degrees. Standard crystallographic convention:
  # a-vector along x, b-vector in the xy-plane, c-vector completes the set.
  al <- alpha * pi / 180
  be <- beta * pi / 180
  ga <- gamma * pi / 180

  ax <- a; ay <- 0; az <- 0
  bx <- b * cos(ga); by <- b * sin(ga); bz <- 0
  cx <- c * cos(be)
  cy <- c * (cos(al) - cos(be) * cos(ga)) / sin(ga)
  cz_sq <- c^2 - cx^2 - cy^2
  cz <- if (cz_sq > 0) sqrt(cz_sq) else 0

  matrix(c(ax, ay, az, bx, by, bz, cx, cy, cz), nrow = 3, byrow = TRUE)
}

frac_to_cart <- function(frac_xyz, lattice) {
  # frac_xyz: N x 3 matrix of fractional coords, lattice: 3x3 matrix with
  # rows a_vector, b_vector, c_vector. Returns N x 3 Cartesian coords.
  frac_xyz %*% lattice
}

vec_angle_deg <- function(p1, p2, p3) {
  # angle at p2 formed by p1-p2-p3, in degrees
  v1 <- p1 - p2
  v2 <- p3 - p2
  cos_theta <- sum(v1 * v2) / (sqrt(sum(v1^2)) * sqrt(sum(v2^2)))
  cos_theta <- max(-1, min(1, cos_theta))
  acos(cos_theta) * 180 / pi
}

dihedral_deg <- function(p1, p2, p3, p4) {
  # dihedral angle i-j-k-l (p1-p2-p3-p4), in degrees, standard sign convention
  b1 <- p2 - p1
  b2 <- p3 - p2
  b3 <- p4 - p3
  n1 <- cross3(b1, b2)
  n2 <- cross3(b2, b3)
  m1 <- cross3(n1, b2 / sqrt(sum(b2^2)))
  x <- sum(n1 * n2)
  y <- sum(m1 * n2)
  atan2(y, x) * 180 / pi
}

cross3 <- function(u, v) {
  c(u[2]*v[3] - u[3]*v[2], u[3]*v[1] - u[1]*v[3], u[1]*v[2] - u[2]*v[1])
}
