# Parsers for the structure file formats the app accepts. Every parser
# returns a common list:
#   list(elements = character vector (length N),
#        coords   = N x 3 Cartesian coordinate matrix (Angstrom),
#        lattice  = 3x3 matrix (rows a,b,c) or NULL if non-periodic,
#        name     = structure label)

strip_esd <- function(x) sub("\\(.*\\)", "", x)

# Tokenize a CIF data line, respecting single/double quoted strings that may
# contain spaces.
tokenize_cif_line <- function(line) {
  m <- gregexpr("'[^']*'|\"[^\"]*\"|\\S+", line)
  toks <- regmatches(line, m)[[1]]
  gsub("^['\"]|['\"]$", "", toks)
}

parse_cif <- function(text_lines, name = "structure") {
  lines <- text_lines
  # drop comments-only lines but keep inline data
  get_tag_value <- function(tag) {
    hit <- grep(paste0("^", tag, "(\\s|$)"), trimws(lines))
    if (length(hit) == 0) return(NA_character_)
    toks <- tokenize_cif_line(trimws(lines[hit[1]]))
    if (length(toks) >= 2) strip_esd(toks[2]) else NA_character_
  }

  a <- as.numeric(get_tag_value("_cell_length_a"))
  b <- as.numeric(get_tag_value("_cell_length_b"))
  c <- as.numeric(get_tag_value("_cell_length_c"))
  alpha <- as.numeric(get_tag_value("_cell_angle_alpha"))
  beta  <- as.numeric(get_tag_value("_cell_angle_beta"))
  gamma <- as.numeric(get_tag_value("_cell_angle_gamma"))

  if (any(is.na(c(a, b, c, alpha, beta, gamma)))) {
    stop("CIF file is missing one or more _cell_length_*/_cell_angle_* tags.")
  }
  lattice <- cell_params_to_vectors(a, b, c, alpha, beta, gamma)

  # locate the atom_site loop block
  loop_starts <- grep("^loop_", trimws(lines))
  atom_site_loop <- NULL
  for (ls in loop_starts) {
    j <- ls + 1
    tags <- character(0)
    while (j <= length(lines) && grepl("^_", trimws(lines[j]))) {
      tags <- c(tags, trimws(lines[j]))
      j <- j + 1
    }
    if (any(grepl("^_atom_site_", tags))) {
      atom_site_loop <- list(tags = tags, data_start = j)
      break
    }
  }
  if (is.null(atom_site_loop)) stop("No _atom_site_ loop found in CIF file.")

  tags <- atom_site_loop$tags
  col_index <- function(tag_name) {
    idx <- which(tags == tag_name)
    if (length(idx) == 0) NA_integer_ else idx[1]
  }
  i_symbol <- col_index("_atom_site_type_symbol")
  i_label  <- col_index("_atom_site_label")
  i_fx <- col_index("_atom_site_fract_x")
  i_fy <- col_index("_atom_site_fract_y")
  i_fz <- col_index("_atom_site_fract_z")
  i_cx <- col_index("_atom_site_Cartn_x")
  i_cy <- col_index("_atom_site_Cartn_y")
  i_cz <- col_index("_atom_site_Cartn_z")

  use_fractional <- !any(is.na(c(i_fx, i_fy, i_fz)))
  if (!use_fractional && any(is.na(c(i_cx, i_cy, i_cz)))) {
    stop("CIF atom_site loop has neither fractional nor Cartesian coordinate columns.")
  }

  j <- atom_site_loop$data_start
  rows <- list()
  while (j <= length(lines)) {
    raw <- trimws(lines[j])
    if (raw == "" || grepl("^loop_", raw) || grepl("^_", raw) || grepl("^#", raw)) break
    rows[[length(rows) + 1]] <- tokenize_cif_line(raw)
    j <- j + 1
  }
  if (length(rows) == 0) stop("atom_site loop in CIF file has no data rows.")

  ncols <- length(tags)
  elements <- character(length(rows))
  coord_mat <- matrix(NA_real_, nrow = length(rows), ncol = 3)
  for (k in seq_along(rows)) {
    r <- rows[[k]]
    if (length(r) < ncols) r <- c(r, rep(NA_character_, ncols - length(r)))
    sym <- if (!is.na(i_symbol)) r[i_symbol] else r[i_label]
    sym <- gsub("[0-9+\\-].*$", "", sym) # strip trailing charge/index e.g. "Zn1" -> "Zn", "O2-" -> "O"
    sym <- paste0(toupper(substr(sym, 1, 1)), tolower(substr(sym, 2, nchar(sym))))
    elements[k] <- sym
    if (use_fractional) {
      coord_mat[k, ] <- c(as.numeric(strip_esd(r[i_fx])),
                           as.numeric(strip_esd(r[i_fy])),
                           as.numeric(strip_esd(r[i_fz])))
    } else {
      coord_mat[k, ] <- c(as.numeric(strip_esd(r[i_cx])),
                           as.numeric(strip_esd(r[i_cy])),
                           as.numeric(strip_esd(r[i_cz])))
    }
  }

  coords <- if (use_fractional) frac_to_cart(coord_mat, lattice) else coord_mat
  list(elements = elements, coords = coords, lattice = lattice, name = name)
}

parse_poscar <- function(text_lines, name = "structure") {
  lines <- text_lines
  scale <- as.numeric(trimws(lines[2]))
  a_vec <- as.numeric(strsplit(trimws(lines[3]), "\\s+")[[1]]) * scale
  b_vec <- as.numeric(strsplit(trimws(lines[4]), "\\s+")[[1]]) * scale
  c_vec <- as.numeric(strsplit(trimws(lines[5]), "\\s+")[[1]]) * scale
  lattice <- rbind(a_vec, b_vec, c_vec)
  rownames(lattice) <- NULL

  line6 <- strsplit(trimws(lines[6]), "\\s+")[[1]]
  is_vasp5 <- suppressWarnings(any(is.na(as.numeric(line6))))
  if (is_vasp5) {
    element_syms <- line6
    counts <- as.integer(strsplit(trimws(lines[7]), "\\s+")[[1]])
    next_line <- 8
  } else {
    stop("POSCAR/CONTCAR without an element-symbol line (VASP4 format) is not supported; add element symbols on line 6 (VASP5 format).")
  }

  sel_dyn <- grepl("^\\s*[Ss]", lines[next_line])
  if (sel_dyn) next_line <- next_line + 1

  coord_mode_line <- trimws(lines[next_line])
  is_direct <- grepl("^[DdFf]", coord_mode_line) # Direct / Fractional
  next_line <- next_line + 1

  natoms <- sum(counts)
  elements <- rep(element_syms, counts)
  coord_mat <- matrix(NA_real_, nrow = natoms, ncol = 3)
  for (i in seq_len(natoms)) {
    toks <- strsplit(trimws(lines[next_line + i - 1]), "\\s+")[[1]]
    coord_mat[i, ] <- as.numeric(toks[1:3])
  }

  coords <- if (is_direct) frac_to_cart(coord_mat, lattice) else coord_mat
  list(elements = elements, coords = coords, lattice = lattice, name = name)
}

parse_xyz_with_lattice <- function(text_lines, name = "structure") {
  lines <- text_lines
  natoms <- as.integer(trimws(strsplit(trimws(lines[1]), "\\s+")[[1]][1]))
  comment <- lines[2]

  m <- regmatches(comment, regexpr("\\[.*\\]", comment))
  if (length(m) == 0 || m == "") {
    stop("XYZ comment line does not contain a 'unitcell [ {..},{..},{..} ]' lattice block.")
  }
  nums <- as.numeric(regmatches(m, gregexpr("-?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?", m))[[1]])
  if (length(nums) < 9) stop("Could not parse 3 lattice vectors from the XYZ comment line.")
  lattice <- matrix(nums[1:9], nrow = 3, byrow = TRUE)

  elements <- character(natoms)
  coord_mat <- matrix(NA_real_, nrow = natoms, ncol = 3)
  for (i in seq_len(natoms)) {
    toks <- strsplit(trimws(lines[2 + i]), "\\s+")[[1]]
    elements[i] <- toks[1]
    coord_mat[i, ] <- as.numeric(toks[2:4])
  }
  list(elements = elements, coords = coord_mat, lattice = lattice, name = name)
}

parse_xyz_plain <- function(text_lines, name = "structure") {
  lines <- text_lines
  natoms <- as.integer(trimws(strsplit(trimws(lines[1]), "\\s+")[[1]][1]))
  elements <- character(natoms)
  coord_mat <- matrix(NA_real_, nrow = natoms, ncol = 3)
  for (i in seq_len(natoms)) {
    toks <- strsplit(trimws(lines[2 + i]), "\\s+")[[1]]
    elements[i] <- toks[1]
    coord_mat[i, ] <- as.numeric(toks[2:4])
  }
  list(elements = elements, coords = coord_mat, lattice = NULL, name = name)
}

# Dispatcher. `format_choice` is one of:
# "cif", "poscar", "xyz_lattice", "xyz_plain", "auto"
parse_structure_file <- function(filepath, original_filename, format_choice = "auto") {
  text_lines <- readLines(filepath, warn = FALSE)
  name <- tools::file_path_sans_ext(basename(original_filename))
  ext <- tolower(tools::file_ext(original_filename))
  base_upper <- toupper(basename(original_filename))

  fmt <- format_choice
  if (fmt == "auto") {
    if (ext == "cif") {
      fmt <- "cif"
    } else if (grepl("POSCAR|CONTCAR", base_upper)) {
      fmt <- "poscar"
    } else if (ext %in% c("xyz", "txt") && length(text_lines) >= 2 &&
               grepl("unitcell|jmolscript", text_lines[2], ignore.case = TRUE)) {
      fmt <- "xyz_lattice"
    } else if (ext %in% c("xyz", "txt")) {
      fmt <- "xyz_plain"
    } else {
      stop("Could not auto-detect the file format from its name/extension. Please pick a format explicitly.")
    }
  }

  switch(fmt,
    cif = parse_cif(text_lines, name),
    poscar = parse_poscar(text_lines, name),
    xyz_lattice = parse_xyz_with_lattice(text_lines, name),
    xyz_plain = parse_xyz_plain(text_lines, name),
    stop("Unknown structure format: ", fmt)
  )
}
