# "About" tab content: citation + an honest map of which SAVESTEPS protocol
# steps (as numbered in the paper below) this app implements vs. does not.
#
# Ghanavati, R.; Escobosa, A.C.; Manz, T.A. An automated protocol to
# construct flexibility parameters for classical forcefields: applications
# to metal-organic frameworks. RSC Advances 2024, 14, 22714-22762,
# doi:10.1039/D4RA01859A.

step_row <- function(step, text, done) {
  tags$tr(
    tags$td(step, style = "white-space:nowrap; vertical-align:top; font-weight:bold;"),
    tags$td(if (done) "✅" else "—", style = "text-align:center; vertical-align:top;"),
    tags$td(text)
  )
}

about_tab_ui <- fluidRow(
  column(
    width = 10, offset = 1,
    tags$h4("About this app"),
    tags$p(
      "This app implements a portion of the SAVESTEPS protocol for perceiving and typing a MOF's ",
      "(or molecule's) bonds, angles, and dihedrals from an already-optimized geometry, including ",
      "numeric-value type refinement (splitting a type further when instances' bond lengths / angle ",
      "values / dihedral values actually differ), 3-/4-membered ring detection (with Urey-Bradley ",
      "diagonal stretches and exclusion of ring-adjacent angles/dihedrals from the active set), ",
      "dihedral rotatable/nonrotatable/linear/hindered classification (including a rotation-based ",
      "connectivity check for hindered dihedrals), double-counted-image removal, and dihedral-type ",
      "pruning, plus an interactive 3D structure viewer."
    ),
    tags$h4("Please cite"),
    tags$blockquote(
      "Ghanavati, R.; Escobosa, A.C.; Manz, T.A. An automated protocol to construct flexibility ",
      "parameters for classical forcefields: applications to metal-organic frameworks. ",
      tags$i("RSC Advances"), " 2024, ", tags$b("14"), ", 22714-22762, ",
      tags$a(href = "https://doi.org/10.1039/D4RA01859A", target = "_blank", "doi:10.1039/D4RA01859A"), "."
    ),
    tags$blockquote(
      "Ghanavati, R.; Escobosa, A.C.; Manz, T.A. Correction to “An automated protocol to construct ",
      "flexibility parameters for classical forcefields: applications to metal-organic frameworks”. ",
      tags$i("RSC Advances"), " 2025, ",
      tags$a(href = "https://doi.org/10.1039/d5ra90084k", target = "_blank", "doi:10.1039/D5RA90084K"), "."
    ),
    tags$h4("Which steps of the protocol this app covers"),
    tags$p(
      "The paper's Section 3 (\"Overview of the SAVESTEPS approach\") lays out the full pipeline as ",
      "13 numbered steps. This app covers part of Steps 3-6 and part of Step 8; everything else -- ",
      "structure verification, the actual quantum chemistry, and the force-constant fitting -- is not ",
      "implemented here and would need the full SAVESTEPS Python package or equivalent tooling."
    ),
    tags$table(
      class = "table table-striped table-condensed",
      tags$thead(tags$tr(tags$th("Step"), tags$th("Here?"), tags$th("What it does"))),
      tags$tbody(
        step_row("Step 1", "Check the starting structure for misbonded atoms and other chemical errors; reject bad structures.", FALSE),
        step_row("Step 2", "Run a quantum chemistry geometry optimization to get the ground-state structure.", FALSE),
        step_row("Step 3", "Type atoms (Chen-Manz scheme); type bonds/angles/dihedrals requiring (i) same atom-type combination, (ii) equilibrium value matching within a tolerance, and (iii) same combination of the lower-level (bond/angle) types.", TRUE),
        step_row("Step 4a", "Classify a dihedral as “linear” when a flanking bond angle is within ~0.03 rad of π.", TRUE),
        step_row("Step 4b", "Add Urey-Bradley stretches for 4-membered ring diagonals; flag angles in 3-/4-membered rings and exclude them from the active angle set; exclude dihedrals containing those flagged angles from the active set.", TRUE),
        step_row("Step 5", "Classify a dihedral instance as rotatable/nonrotatable by whether its middle bond sits on a ring; a type is nonrotatable if any instance is.", TRUE),
        step_row("Step 6", "Dihedral-type pruning: when two+ types sit on an identical set of middle bonds, keep only one representative type. Also removes double-counted instances (same atom numbers, different periodic image) down to one representative before pruning.", TRUE),
        step_row("Step 7", "Run quantum chemistry: finite-displacement (“Hessian”) geometries/forces plus AIMD training and validation geometries/forces.", FALSE),
        step_row("Step 8", "For each rotatable dihedral type, rotate the smaller branch off the middle bond in 10° steps through a full circle and re-type the structure at each step; reclassify as “hindered” if any step would break or form a bond. (Connectivity check only -- this app does not generate or save the displaced geometries themselves, since that's for feeding a quantum chemistry calculation this app doesn't run.)", TRUE),
        step_row("Step 9", "Single-point quantum chemistry energies for the Step 8 dihedral-displaced geometries.", FALSE),
        step_row("Step 10", "Project each rotatable dihedral's energy curve onto orthonormal torsion modes; smart-select the significant ones.", FALSE),
        step_row("Step 11", "Assemble the actual potential model (angle-bending, CADT/ADDT dihedral torsion, stretches, optional cross terms).", FALSE),
        step_row("Step 12", "Fit all force constants via regularized (LASSO) linear least-squares on the training set.", FALSE),
        step_row("Step 13", "Validate the fitted model (R-squared, RMSE) against a held-out AIMD validation set.", FALSE)
      )
    ),
    tags$p(
      style = "color:#888; font-size:0.9em;",
      "Note: the Hindered check can only detect hindrance that shows up as a spurious bond forming or ",
      "an existing bond breaking during the rotation scan (the same criterion the original protocol ",
      "uses) -- it won't catch purely steric strain that never crosses the bonding-distance cutoff. A ",
      "branch that reaches into the (possibly infinite) periodic framework rather than terminating in a ",
      "finite pendant group is left as Rotatable (inconclusive) rather than checked."
    )
  )
)
