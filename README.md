# Europa tidal heating

Numerical code, input files and figure data accompanying the manuscript:

M. Kervazo, M. Běhounková, G. Tobie, G. Choblet and C. Dumoulin, *Impact of melt accumulation on tidal heat production in Europa's rocky mantle*, Planetary Research.

The Fortran code computes the viscoelastic tidal response and tidal dissipation of a radially stratified model of Europa, including the effects of partial melting in its rocky mantle.

## Repository contents

- `src/`: Fortran source files, Makefile and script to run the simulations.
- `input/`: simulation input files (`input.nml`) and radial viscosity and melt-rheology parameter files (`data/`). The simulations are grouped by eccentricity and model configuration.
- `figure_data/`: numerical data associated with Figures 1–5 of the manuscript.

## Compilation and execution

A Fortran compiler (`gfortran`) and `make` are required. From the repository root, run:

```sh
sh src/run_all.sh
```

The script compiles the code, copies the input files to `work/` and runs the simulations in their respective directories. The original input files are not modified. Use a fresh working directory for a new set of calculations. To select the compiler or run several cases in parallel, use, for example:

```sh
FC=gfortran JOBS=4 sh src/run_all.sh
```

Simulation outputs are written to `Melt/` within each calculation directory. For cases with `map = .true.`, heating maps are also written to `Htide_Output/Tot/`.

## Main output files

| File | Columns |
| --- | --- |
| `melt_ptide_pmant_pice_phaselag.dat` | Melt fraction, total/mantle/ice tidal power (GW), Re/Im ratios of surface radial functions y5, y1 and y3, ice thickness (m), ice viscosity (Pa s), molten layer thickness parameter (m). |
| `melt_y5_y1_y3_k2_h2_l2.dat` | Melt fraction, real and imaginary parts of surface y5, y1 and y3, real k2, h2 and l2. |
| `Melt_ecc_Ptot_Pmant_bice_etaice_bast.dat` | Melt fraction, eccentricity, total and mantle tidal power (GW), ice thickness (m), ice viscosity (Pa s), molten layer thickness parameter (m), melt volume (km³). |
| `R_rho_eta_mu_y1_y5_htide_phi_21Gyr_phi.dat` | Radius (m), density (kg/m³), viscosity (Pa s), shear modulus (Pa), complex y1 and y5, volumetric tidal heating (W/m³), melt fraction. |
| `Htide_colatt_lont_radius_tot.dat` | Radius (m), colatitude (degrees), longitude (degrees), volumetric heating (W/m³), including the prescribed radiogenic contribution when applicable. |

Phase lags can be derived from the complex surface radial functions. Use the melt fractions recorded in the output files.

## Code and data provenance

The tidal-response code is based on the implementation of Gabriel Tobie (Takeuchi & Saito, 1972; Tobie et al., 2005), subsequently adapted to partially molten interiors. The radial rheological profiles are based on Běhounková et al. (2021). The melt-rheology parameterization is available at https://github.com/marie-behounkova/rheomelt-parameterization.

This is a research version of the code accompanying the manuscript. The plotting notebooks are not included in this version. Additional documentation may be provided in future releases.

## Citation

Please cite Tobie et al., 2005, the associated manuscript and the archived version of the repository when using the code or data.
