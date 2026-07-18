# Element table and bonding radii, ported from the SAVESTEPS python pipeline
# (bonded_atoms_type_1_30_2024.py / create_atom_type_dir_POSCAR_1_30_2024.py).
# Two atoms are considered bonded when their minimum-image distance is <=
# (radius[i] + radius[j]) * bond_tolerance.

ELEMENTS <- c(
  "H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca",
  "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
  "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
  "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
  "Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm",
  "Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds"
)

ATOMIC_NUM <- setNames(seq_along(ELEMENTS), ELEMENTS)

ATOM_TYPING_RADII <- setNames(c(
  0.38,1.11,0.86,0.53,1.01,0.88,0.86,0.89,0.82,1.3,1.15,1.28,1.53,1.38,1.28,1.2,1.17,1.75,1.44,1.17,
  1.62,1.65,1.51,1.53,1.53,1.43,1.31,1.33,1.31,1.41,1.4,1.35,1.39,1.4,1.39,2.09,1.65,1.3,1.84,1.73,1.66,
  1.57,1.53,1.58,1.63,1.68,1.56,1.56,1.53,1.64,1.64,1.65,1.58,2.16,1.85,1.52,1.91,1.98,1.75,1.92,1.98,
  1.89,1.83,1.79,1.82,1.79,1.63,1.8,1.84,1.8,1.86,1.73,1.61,1.33,1.29,1.55,1.5,1.66,1.68,1.88,1.73,1.72,
  1.72,1.86,1.88,2.58,2.18,2.08,2.06,1.97,1.79,1.76,1.73,1.71,1.69,1.68,1.68,1.68,1.68,1.68,1.68,1.68,
  1.68,1.78,1.78,1.78,1.78,1.78,1.78,1.78
), ELEMENTS)
# Note: the source (SAVESTEPS python pipeline) radii table has only 109
# entries for 110 elements -- it silently omits Ds (element 110), a
# synthetic, radioactive element with a half-life of seconds that no real
# MOF will ever contain. The line above appends a 110th entry (1.78,
# matching the placeholder value already used for its neighbors
# Hs/Mt/Rf/Db/Sg/Bh) purely so ELEMENTS and ATOM_TYPING_RADII stay the same
# length; it has no effect on any real structure.

# CPK-style colors for common elements; anything not listed falls back to a
# deterministic color generated from the element's atomic number.
CPK_COLORS <- c(
  H = "#FFFFFF", C = "#909090", N = "#3050F8", O = "#FF0D0D", F = "#90E050",
  Cl = "#1FF01F", Br = "#A62929", I = "#940094", S = "#FFFF30", P = "#FF8000",
  B = "#FFB5B5", Si = "#F0C8A0",
  Zn = "#7D80B0", Cu = "#C88033", Fe = "#E06633", Ni = "#50D050", Co = "#F090A0",
  Mn = "#9C7AC7", Cr = "#8A99C7", Mg = "#8AFF00", Ca = "#3DFF00", Na = "#AB5CF2",
  K = "#8F40D4", Al = "#BFA6A6", Ti = "#BFC2C7", Zr = "#94E0E0", V = "#A6A6AB",
  Cd = "#FFD98F", Ag = "#C0C0C0", Au = "#FFD123", Pt = "#D0D0E0", Pd = "#006985",
  La = "#70D4FF", Ce = "#FFFFC7"
)

element_color <- function(element_symbol) {
  known <- CPK_COLORS[element_symbol]
  if (!is.na(known)) return(unname(known))
  z <- ATOMIC_NUM[[element_symbol]]
  if (is.null(z) || is.na(z)) return("#CCCCCC")
  grDevices::hsv(h = (z * 47 %% 360) / 360, s = 0.55, v = 0.85)
}
