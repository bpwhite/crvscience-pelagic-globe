# Preprints for PelagicGlobe (`crvscience-pelagic-globe`)

**Author:** Bryan P. White (ORCID: [0000-0003-4191-3511](https://orcid.org/0000-0003-4191-3511))  
**Hub:** https://crvscience.com/apps

This directory contains the self-contained LaTeX source (`.tex`) and single-column PDF (`.pdf`) editions for both:

1. **Companion Application Preprint (`APP-PREPRINT-22`)**
   - **Title:** PelagicGlobe: 3D Ocean Currents Globe, Conformal Spilhaus World Ocean Projection, and NOAA ENSO Workstation in Rust
   - **Files:** `pelagic-globe-3d-ocean-spilhaus-enso.tex`, `pelagic-globe-3d-ocean-spilhaus-enso.pdf`
   - **Canonical Zenodo DOI:** https://doi.org/10.5281/zenodo.24442089
   - **Figshare Mirror:** https://figshare.com/articles/preprint/PelagicGlobe_3D_Ocean_Currents_Globe_Conformal_Spilhaus_World_Ocean_Projection_and_NOAA_ENSO_Workstation_in_Rust/35084812
   - **Web Preprint:** https://crvscience.com/apps/preprints/pelagic-globe-3d-ocean-spilhaus-enso

2. **CRV Science Public Apps Network Overview (`APP-PREPRINT-01`)**
   - **Title:** The CRV Science App Network: Twenty-Two Open-Source Rust Programs for Science, Data, and Simulation
   - **Files:** `crv-science-public-apps-network.tex`, `crv-science-public-apps-network.pdf`
   - **Canonical Zenodo DOI:** https://doi.org/10.5281/zenodo.24222398
   - **Figshare Mirror:** https://figshare.com/articles/preprint/The_CRV_Science_App_Network_Seventeen_Open-Source_Rust_Programs_for_Science_Data_and_Simulation/35013736
   - **Web Preprint:** https://crvscience.com/apps/preprints/crv-science-public-apps-network

## Compiling Standalone LaTeX

Both `.tex` files inline their bibliography so each compiles standalone with standard `pdflatex`:

```bash
pdflatex pelagic-globe-3d-ocean-spilhaus-enso.tex
pdflatex crv-science-public-apps-network.tex
```
