# MOF / crystal structure analyzer: reads an optimized geometry (CIF, VASP
# POSCAR/CONTCAR, or XYZ) and reports bonds, angles, dihedrals, and an
# interactive 3D picture of the structure.
#
# Run with:  shiny::runApp("MOF_Shiny_App")
#
# Copyright (C) 2026 Reza Ghanavati
#
# Several algorithms in this app (bond/angle/dihedral perception, Chen-Manz
# atom typing, ring-based dihedral rotatability classification, and
# dihedral-type pruning) are reimplementations of methods from the
# SAVESTEPS Python package, Copyright (C) 2024 Reza Ghanavati,
# Alma C. Escobosa, Thomas A. Manz, licensed under GPL-3.0. See:
# Ghanavati, R.; Escobosa, A.C.; Manz, T.A. An automated protocol to
# construct flexibility parameters for classical forcefields: applications
# to metal-organic frameworks. RSC Advances 2024, 14, 22714-22762,
# doi:10.1039/D4RA01859A (correction: RSC Advances 2025, doi:10.1039/D5RA90084K).
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version. This program is distributed WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License, included in
# this repository as LICENSE, for more details.

library(shiny)
library(rgl)
library(DT)

options(rgl.useNULL = TRUE) # never try to open a real OS window; render to the browser only

# Explicitly source helpers (Shiny auto-sources R/ too, but this makes the
# app self-contained when run with plain Rscript as well).
invisible(lapply(list.files(file.path(getwd(), "R"), pattern = "\\.R$", full.names = TRUE), source))

ui <- fluidPage(
  theme = shinythemes::shinytheme("flatly"),
  titlePanel("MOF / Crystal Structure: Bonds, Angles & Dihedrals"),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      fileInput("struct_file", "Upload optimized geometry",
                accept = NULL, placeholder = "CIF, POSCAR/CONTCAR, or XYZ"),
      selectInput("format_choice", "File format", choices = c(
        "Auto-detect from filename" = "auto",
        "CIF" = "cif",
        "VASP POSCAR/CONTCAR (VASP5, Direct or Cartesian)" = "poscar",
        "XYZ with lattice (unitcell comment line)" = "xyz_lattice",
        "XYZ molecule (no lattice / non-periodic)" = "xyz_plain"
      )),
      selectInput("example_choice", "...or load a bundled example", choices = c(
        "None" = "none",
        "Ethane (XYZ molecule, non-periodic)" = "ethane",
        "Benzene (XYZ molecule, non-periodic)" = "benzene",
        "Acetic acid (XYZ molecule, non-periodic)" = "acetic_acid",
        "Carbon chain (CIF, 1D periodic)" = "chain",
        "Diamond (CIF, 3D periodic)" = "diamond",
        "Graphene sheet (XYZ+lattice, 2D periodic)" = "graphene",
        "Zn-N2 pcu MOF-like net (POSCAR, 3D periodic)" = "zn_mof",
        "JIVFUQ: real Zn-phosphonate MOF (XYZ+lattice, 3D periodic)" = "jivfuq"
      )),
      sliderInput("tolerance", "Bond distance tolerance (× sum of atomic radii)",
                  min = 0.7, max = 1.5, value = 1.0, step = 0.05),
      checkboxInput("show_cell", "Show unit cell wireframe", value = TRUE),
      checkboxInput("show_labels", "Show atom labels", value = FALSE),
      actionButton("analyze_btn", "Analyze structure", class = "btn-primary"),
      hr(),
      uiOutput("summary_box"),
      hr(),
      downloadButton("dl_bonds", "Bonds CSV"),
      downloadButton("dl_angles", "Angles CSV"),
      downloadButton("dl_dihedrals", "Dihedrals CSV"),
      downloadButton("dl_dihedrals_pruned", "Dihedrals (Pruned) CSV"),
      hr(),
      helpText(tags$b("Note:"), "the Bond/Angle/Dihedral Type tabs report the",
               "equilibrium value (mean/min/max) directly from this single input",
               "structure -- they assume it is already geometry-optimized, not an",
               "ensemble of DFT/AIMD snapshots."),
      helpText(tags$b("Pruning:"), "the \"(Pruned)\" tabs drop dihedral types that",
               "sit on the exact same set of middle bonds as another type -- only",
               "one representative type is needed per physical bond, chosen for",
               "having the least near-linear (most well-defined) flanking angles."),
      hr(),
      helpText(tags$b("Please cite:"), tags$br(),
               "Ghanavati, R.; Escobosa, A.C.; Manz, T.A. An automated protocol to",
               "construct flexibility parameters for classical forcefields:",
               "applications to metal-organic frameworks.", tags$i("RSC Advances"),
               "2024,", tags$b("14"), ", 22714-22762,",
               tags$a(href = "https://doi.org/10.1039/D4RA01859A", target = "_blank",
                      "doi:10.1039/D4RA01859A"), ".",
               tags$br(), tags$br(),
               "Ghanavati, R.; Escobosa, A.C.; Manz, T.A. Correction to “An",
               "automated protocol to construct flexibility parameters for classical",
               "forcefields: applications to metal-organic frameworks”.",
               tags$i("RSC Advances"), "2025,",
               tags$a(href = "https://doi.org/10.1039/d5ra90084k", target = "_blank",
                      "doi:10.1039/D5RA90084K"), ".",
               tags$br(), tags$br(),
               "See the \"About\" tab for which steps of that protocol this app covers.")
    ),
    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel("3D Structure", rgl::rglwidgetOutput("structure_plot", height = "650px")),
        tabPanel("Bonds", DTOutput("bonds_table")),
        tabPanel("Angles", DTOutput("angles_table")),
        tabPanel("Dihedrals", DTOutput("dihedrals_table")),
        tabPanel("Atom Types", DTOutput("atom_types_table")),
        tabPanel("Bond Types", DTOutput("bond_types_table")),
        tabPanel("Angle Types", DTOutput("angle_types_table")),
        tabPanel("Dihedral Types", DTOutput("dihedral_types_table")),
        tabPanel("Dihedrals (Pruned)", DTOutput("dihedrals_pruned_table")),
        tabPanel("Dihedral Types (Pruned)", DTOutput("dihedral_types_pruned_table")),
        tabPanel("About", about_tab_ui)
      )
    )
  )
)

server <- function(input, output, session) {

  analysis <- eventReactive(input$analyze_btn, {
    req(input$example_choice != "none" || !is.null(input$struct_file))

    if (input$example_choice != "none") {
      example_map <- list(
        ethane      = list(file = "ethane_molecule.xyz",      fmt = "xyz_plain"),
        benzene     = list(file = "benzene.xyz",               fmt = "xyz_plain"),
        acetic_acid = list(file = "acetic_acid.xyz",            fmt = "xyz_plain"),
        chain       = list(file = "periodic_carbon_chain.cif", fmt = "cif"),
        diamond     = list(file = "diamond_carbon.cif",         fmt = "cif"),
        graphene    = list(file = "graphene_sheet.xyz",         fmt = "xyz_lattice"),
        zn_mof      = list(file = "POSCAR_Zn_N2_pcu_mof",       fmt = "poscar"),
        jivfuq      = list(file = "JIVFUQ.xyz",                 fmt = "xyz_lattice")
      )
      ex <- example_map[[input$example_choice]]
      path <- file.path(getwd(), "example_data", ex$file)
      fname <- ex$file
      fmt <- ex$fmt
    } else {
      path <- input$struct_file$datapath
      fname <- input$struct_file$name
      fmt <- input$format_choice
    }

    withProgress(message = "Parsing structure...", value = 0.2, {
      struct <- tryCatch(
        parse_structure_file(path, fname, fmt),
        error = function(e) {
          showNotification(paste("Parse error:", conditionMessage(e)), type = "error", duration = NULL)
          NULL
        }
      )
      req(struct)

      incProgress(0.3, message = "Finding bonds, angles & dihedrals...")
      res <- tryCatch(
        analyze_structure(struct, tolerance = input$tolerance),
        error = function(e) {
          showNotification(paste("Analysis error:", conditionMessage(e)), type = "error", duration = NULL)
          NULL
        }
      )
      req(res)
      incProgress(0.5, message = "Done")
      res
    })
  })

  output$summary_box <- renderUI({
    res <- tryCatch(analysis(), error = function(e) NULL)
    if (is.null(res)) return(helpText('Upload a structure (or pick an example) and click "Analyze structure".'))
    n_bonds <- if (is.null(res$bonds)) 0 else nrow(res$bonds)
    n_angles <- if (is.null(res$angles)) 0 else nrow(res$angles)
    n_dihedrals <- if (is.null(res$dihedrals)) 0 else nrow(res$dihedrals)
    tagList(
      tags$b(res$name), tags$br(),
      sprintf("%d atoms, %d element(s)", length(res$elements), length(unique(res$elements))), tags$br(),
      sprintf("%d bonds, %d angles, %d dihedrals", n_bonds, n_angles, n_dihedrals), tags$br(),
      if (is.null(res$lattice)) "Non-periodic (molecule)" else "Periodic structure (minimum-image bonding)"
    )
  })

  output$structure_plot <- rgl::renderRglwidget({
    res <- analysis()
    req(res)
    build_structure_widget(res, show_cell = input$show_cell, show_labels = input$show_labels)
  })

  output$bonds_table <- renderDT({
    res <- analysis(); req(res)
    if (is.null(res$bonds) || nrow(res$bonds) == 0) return(datatable(data.frame(Message = "No bonds found")))
    df <- res$bonds[, c("atom1", "element1", "image1", "atom2", "element2", "image2",
                         "bond", "distance", "bond_type_index", "bond_type")]
    colnames(df) <- c("Atom 1 #", "Elem 1", "Image 1", "Atom 2 #", "Elem 2", "Image 2",
                       "Bond", "Distance (Å)", "Type #", "Bond Type")
    datatable(df, filter = "top", options = list(pageLength = 15))
  })

  output$angles_table <- renderDT({
    res <- analysis(); req(res)
    if (is.null(res$angles) || nrow(res$angles) == 0) return(datatable(data.frame(Message = "No angles found")))
    df <- res$angles[, c("atom1", "element1", "image1", "atom2", "element2", "image2",
                          "atom3", "element3", "image3",
                          "angle_label", "angle_deg", "angle_type_index", "angle_type")]
    colnames(df) <- c("Atom 1 #", "Elem 1", "Image 1", "Atom 2 (center) #", "Elem 2", "Image 2",
                       "Atom 3 #", "Elem 3", "Image 3",
                       "Angle", "Angle (deg)", "Type #", "Angle Type")
    datatable(df, filter = "top", options = list(pageLength = 15))
  })

  output$dihedrals_table <- renderDT({
    res <- analysis(); req(res)
    if (is.null(res$dihedrals) || nrow(res$dihedrals) == 0) return(datatable(data.frame(Message = "No dihedrals found")))
    df <- res$dihedrals[, c("atom1", "element1", "image1", "atom2", "element2", "image2",
                             "atom3", "element3", "image3", "atom4", "element4", "image4",
                             "dihedral_label", "dihedral_deg", "dihedral_type_index", "dihedral_type", "rotatability")]
    colnames(df) <- c("Atom 1 #", "Elem 1", "Image 1", "Atom 2 #", "Elem 2", "Image 2",
                       "Atom 3 #", "Elem 3", "Image 3", "Atom 4 #", "Elem 4", "Image 4",
                       "Dihedral", "Dihedral (deg)", "Type #", "Dihedral Type", "Rotatability")
    datatable(df, filter = "top", options = list(pageLength = 15))
  })

  output$atom_types_table <- renderDT({
    res <- analysis(); req(res)
    df <- data.frame(Atom = seq_along(res$elements), Element = res$elements, `Chen-Manz atom type` = res$atom_types,
                      check.names = FALSE)
    datatable(df, filter = "top", options = list(pageLength = 20))
  })

  output$bond_types_table <- renderDT({
    res <- analysis(); req(res)
    if (is.null(res$bond_types_summary) || nrow(res$bond_types_summary) == 0)
      return(datatable(data.frame(Message = "No bonds found")))
    df <- res$bond_types_summary[, c("Index", "bond", "Type", "Count", "Mean", "Min", "Max")]
    colnames(df) <- c("Type #", "Elements", "Bond Type", "Count", "Mean Distance (Å)", "Min (Å)", "Max (Å)")
    datatable(df, filter = "top", options = list(pageLength = 15))
  })

  output$angle_types_table <- renderDT({
    res <- analysis(); req(res)
    if (is.null(res$angle_types_summary) || nrow(res$angle_types_summary) == 0)
      return(datatable(data.frame(Message = "No angles found")))
    df <- res$angle_types_summary[, c("Index", "angle_label", "Type", "Count", "Mean", "Min", "Max")]
    colnames(df) <- c("Type #", "Elements", "Angle Type", "Count", "Mean Angle (deg)", "Min (deg)", "Max (deg)")
    datatable(df, filter = "top", options = list(pageLength = 15))
  })

  output$dihedral_types_table <- renderDT({
    res <- analysis(); req(res)
    if (is.null(res$dihedral_types_summary) || nrow(res$dihedral_types_summary) == 0)
      return(datatable(data.frame(Message = "No dihedrals found")))
    df <- res$dihedral_types_summary[, c("Index", "dihedral_label", "Type", "Count", "Mean", "Min", "Max", "rotatability")]
    colnames(df) <- c("Type #", "Elements", "Dihedral Type", "Count", "Mean Dihedral (deg)", "Min (deg)", "Max (deg)", "Rotatability")
    datatable(df, filter = "top", options = list(pageLength = 15))
  })

  output$dihedrals_pruned_table <- renderDT({
    res <- analysis(); req(res)
    if (is.null(res$dihedrals_pruned) || nrow(res$dihedrals_pruned) == 0)
      return(datatable(data.frame(Message = "No dihedrals found")))
    df <- res$dihedrals_pruned[, c("atom1", "element1", "image1", "atom2", "element2", "image2",
                                    "atom3", "element3", "image3", "atom4", "element4", "image4",
                                    "dihedral_label", "dihedral_deg", "dihedral_type_index", "dihedral_type", "rotatability")]
    colnames(df) <- c("Atom 1 #", "Elem 1", "Image 1", "Atom 2 #", "Elem 2", "Image 2",
                       "Atom 3 #", "Elem 3", "Image 3", "Atom 4 #", "Elem 4", "Image 4",
                       "Dihedral", "Dihedral (deg)", "Type #", "Dihedral Type", "Rotatability")
    datatable(df, filter = "top", options = list(pageLength = 15))
  })

  output$dihedral_types_pruned_table <- renderDT({
    res <- analysis(); req(res)
    if (is.null(res$dihedral_types_summary_pruned) || nrow(res$dihedral_types_summary_pruned) == 0)
      return(datatable(data.frame(Message = "No dihedrals found")))
    df <- res$dihedral_types_summary_pruned[, c("Index", "dihedral_label", "Type", "Count", "Mean", "Min", "Max",
                                                  "rotatability", "Merged")]
    colnames(df) <- c("Type #", "Elements", "Dihedral Type", "Count", "Mean Dihedral (deg)", "Min (deg)", "Max (deg)",
                       "Rotatability", "Redundant Types Merged")
    datatable(df, filter = "top", options = list(pageLength = 15))
  })

  output$dl_bonds <- downloadHandler(
    filename = function() paste0(analysis()$name, "_bonds.csv"),
    content = function(file) write.csv(analysis()$bonds, file, row.names = FALSE)
  )
  output$dl_angles <- downloadHandler(
    filename = function() paste0(analysis()$name, "_angles.csv"),
    content = function(file) write.csv(analysis()$angles, file, row.names = FALSE)
  )
  output$dl_dihedrals <- downloadHandler(
    filename = function() paste0(analysis()$name, "_dihedrals.csv"),
    content = function(file) write.csv(analysis()$dihedrals, file, row.names = FALSE)
  )
  output$dl_dihedrals_pruned <- downloadHandler(
    filename = function() paste0(analysis()$name, "_dihedrals_pruned.csv"),
    content = function(file) write.csv(analysis()$dihedrals_pruned, file, row.names = FALSE)
  )
}

shinyApp(ui, server)
