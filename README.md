# Homomeric Polypeptide Steered Molecular Dynamics

Steered MD (SMD) pipeline that pulls decapeptides (10-residue homopolymers of
glycine, valine, serine, lysine, and glutamate) along their end-to-end axis
in explicit solvent, then extracts force–extension curves, PMF estimates
(via Jarzynski's equality), and Kuhn lengths for each polymer.

## Pipeline stages

1. **Structure generation** — AlphaFold-style folding of each
   `poly<residue>_10` sequence produces `fold_out/fold_poly*_10_model_0.{cif,pdb}`.
   `convert_cif_pdb.py` batch-converts any `.cif` in a directory to `.pdb`
   (via `Bio.PDB`), giving the working structures in `fold_pdbs/`.

2. **Axis alignment** — the chain must run along +z before force is applied.
   - `rotntrans.tcl` (VMD): computes principal axes of the protein selection,
     rotates the long axis onto z and the second axis onto y (via the VMD
     `Orient` package), then translates so the N atom of residue 1 sits at
     the origin. Writes `<name>_cleaned.pdb`.
   - `rotntrans.py` (PyMOL): alternative/complementary alignment — picks the
     N atoms of residues 1 and 10 (chain A if present), pair-fits them onto
     a synthetic z-axis reference (`Z_LENGTH = 100 Å`), reports the angle to
     +z before/after, and writes `aligned_pdbs/aligned_<name>.pdb`.

3. **Solvation & system setup (CHARMM-GUI)** — each aligned structure is
   submitted to CHARMM-GUI's Solution Builder to generate a solvated,
   ionized, periodic system: PSF/PDB topology, CHARMM36/CGenFF parameter
   files (`toppar/`), and NAMD-ready equilibration scripts
   (`step3_pbcsetup.*`, `step4_equilibration.*`). This step's raw output
   (`charmm_out/`) is **not** tracked in this repo (~440 MB of standard
   force-field parameter files and ion-placement logs, regenerable from
   CHARMM-GUI) — only the pipeline scripts that consume it are.
   - Before submission: re-check z-alignment of the CHARMM-GUI output.
   - After equilibration: re-check alignment again (drift can occur) and
     re-run the rotate/translate step if needed.

4. **Steered MD run** (`smd.namd` + `smd.tcl`, NAMD + Tcl forces):
   - `smd.namd` reads the CHARMM-GUI PSF/PDB and the equilibration
     restart/xsc files, loads the full CHARMM36 parameter set, sets PME
     electrostatics, Langevin thermostat at 300 K, 2 fs timestep, and
     enables `tclForces` pointing at `smd.tcl`. Runs 1,750,000 steps.
   - `smd.tcl` implements the SMD harmonic-constraint protocol: the N atom
     of residue 1 is held by a stationary harmonic spring (`k = 7.2
     kcal/mol/Å²`), while the N atom of residue 10 is pulled by a spring
     whose z-anchor moves at constant velocity `v = 0.00002 Å/timestep`.
     Every `Tclfreq` (50) steps it appends `time, z-position, z-force` of
     the pulled atom to `da_smd_tcl_<runid>.out` — one such run is repeated
     per replicate (10 replicates per residue type, per `analysis.py`).

5. **Analysis** (`analysis.py`):
   - `parse_constraints` reads `c1z`/`c2z` back out of each protein's
     `smd.tcl` to reconstruct the constraint-point trajectory.
   - `load_trajectories` loads the 10 replicate `.out` files per residue
     type into `(n_timesteps, n_traj)` arrays of time/position/force.
   - `compute_extension_and_work` reconstructs the moving constraint
     separation `c(t) = (c2z - c1z) + v·t` and the cumulative
     work `W = v·dt·Σ(-F)`.
   - `compute_pmf` produces three PMF estimates from Jarzynski's equality:
     the exponential average (numerically stabilized), the pure work
     average, and the second-order cumulant expansion.
   - `compute_kuhn_length` fits the low-extension linear regime of
     force vs. extension to a freely-jointed-chain model and back out the
     Kuhn length `a_Kuhn = sqrt(3·k_BT·slope / N_residues)`.
   - `compute_all_results` runs the above for every residue type, then
     `plot_individual`, `plot_force_extension_by_group`,
     `plot_pmf_by_group`, `plot_bar_comparisons`, `print_summary_table`,
     and `plot_correlations` render the comparative figures/tables used in
     the report.

6. **Write-up** — `proj1.ipynb` is the exploratory notebook that develops
   the analysis; `report.tex`/`report.pdf` is the final write-up.

## Layout

```
fold_out/        AlphaFold structure predictions (.cif + converted .pdb)
fold_pdbs/       Working .pdb per homopolymer, pre-alignment
aligned_pdbs/    Post-alignment .pdb (z-axis aligned), CHARMM-GUI input
convert_cif_pdb.py   .cif → .pdb batch converter
rotntrans.tcl        VMD principal-axis alignment
rotntrans.py         PyMOL N1→N10 alignment onto z-axis
smd.namd             NAMD SMD run configuration
smd.tcl              Tcl SMD force protocol (dual-spring pulling)
analysis.py          Force/PMF/Kuhn-length analysis + figures
proj1.ipynb          Exploratory analysis notebook
report.tex/.pdf      Final report
project.md           Working checklist per residue type
```

`charmm_out/` (CHARMM-GUI raw output per residue type) and the background
dissertation reference PDF are excluded from version control — see
`.gitignore`.

## Dependencies

- VMD (with the `Orient` package) or PyMOL, for alignment
- CHARMM-GUI Solution Builder access, for solvated system setup
- NAMD, for running the SMD simulation
- Python: `numpy`, `scipy`, `matplotlib`, `biopython`

## Reproducing a residue type end-to-end

```
1. Fold sequence  → fold_out/fold_poly<res>_10_model_0.cif
2. python convert_cif_pdb.py            # → .pdb
3. pymol -cq rotntrans.py                # or vmd -e rotntrans.tcl
   → aligned_pdbs/aligned_deca<res>.pdb
4. Submit aligned .pdb to CHARMM-GUI Solution Builder
   → charmm_out/deca<res>/charmm-gui-*/...
5. Run equilibration scripts from the CHARMM-GUI bundle
6. Re-check z-alignment; re-run step 3 on equilibrated structure if needed
7. namd2 smd.namd +setup runid=<1..10> seed=<seed>   # per replicate
8. python analysis.py                    # aggregate + plot
```
