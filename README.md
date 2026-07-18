# MOF Flexibility Viewer

An R Shiny app that reads an optimized geometry of a molecule or a periodic
structure (e.g. a MOF) and reports its bonds, angles, and dihedrals, an
interactive 3D picture, and dihedral rotatable/nonrotatable classification
with type-based pruning.

**Live demo:** https://rghanavati.shinyapps.io/mof-flexibility-viewer/

## Features

- **Input formats:** CIF, VASP POSCAR/CONTCAR (VASP5, Direct or Cartesian),
  XYZ with an embedded lattice (`unitcell [...]` comment line), and plain
  (non-periodic) XYZ molecules.
- **Periodic bond perception:** bonds are found via a minimum-image search
  over the 27 nearest periodic cell images, using a covalent-radii cutoff
  with an adjustable tolerance. Every qualifying image is kept (not just
  the nearest), which matters for symmetric/small-cell motifs common in
  MOFs.
- **Atom typing:** the Chen & Manz second-neighbor-shell atom typing scheme
  (see citation below).
- **Bond / angle / dihedral instance lists**, each atom shown as
  `# / element / periodic image`, grouped and indexed by type.
- **Bond / Angle / Dihedral Type summaries:** unique types with instance
  counts and equilibrium value ranges (mean/min/max), read directly from
  the input geometry.
- **Dihedral rotatability classification:** each dihedral instance is
  classified Rotatable / Nonrotatable (its middle bond is/isn't part of a
  ring, found via periodic-graph BFS) / Linear (a flanking bond angle is
  within ~1.7° of 180°, making the torsion ill-defined).
- **Dihedral-type pruning:** when two or more dihedral types sit on an
  identical set of middle bonds, only one representative type is kept.
- **Interactive 3D structure viewer** (ball-and-stick, via `rgl`), with an
  optional unit-cell wireframe and atom labels.
- **CSV export** for bonds, angles, dihedrals, and pruned dihedrals.
- Seven bundled example structures spanning all four formats, from simple
  molecules to a real MOF (see `example_data/`).

## Running locally

Requires R with the `shiny`, `rgl`, `DT`, and `shinythemes` packages
installed:

```r
install.packages(c("shiny", "rgl", "DT", "shinythemes"))
shiny::runApp("path/to/mof-flexibility-viewer")
```

Pick a bundled example from the sidebar, or upload your own CIF/POSCAR/XYZ,
then click **Analyze structure**.

## Relationship to SAVESTEPS

This app implements part of the structure-typing stage of the **SAVESTEPS**
protocol for constructing classical-forcefield flexibility parameters for
MOFs. See the in-app **About** tab for a step-by-step breakdown of exactly
which of the protocol's 13 steps this app covers (atom/bond/angle/dihedral
typing, dihedral rotatability classification, and dihedral-type pruning)
versus which require the full SAVESTEPS pipeline (structure verification,
quantum chemistry calculations, and force-constant fitting).

### Please cite

> Ghanavati, R.; Escobosa, A.C.; Manz, T.A. An automated protocol to
> construct flexibility parameters for classical forcefields: applications
> to metal-organic frameworks. *RSC Advances* 2024, *14*, 22714-22762,
> [doi:10.1039/D4RA01859A](https://doi.org/10.1039/D4RA01859A).
>
> Ghanavati, R.; Escobosa, A.C.; Manz, T.A. Correction to "An automated
> protocol to construct flexibility parameters for classical forcefields:
> applications to metal-organic frameworks". *RSC Advances* 2025,
> [doi:10.1039/D5RA90084K](https://doi.org/10.1039/d5ra90084k).

## License

GPL-3.0. See [LICENSE](LICENSE). Several algorithms here are
reimplementations of methods from the SAVESTEPS Python package (Copyright
2024 Reza Ghanavati, Alma C. Escobosa, Thomas A. Manz, GPL-3.0).
