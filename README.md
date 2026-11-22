# PelagicGlobe (`pelagic-globe` v1.1.0)

Native `#![forbid(unsafe_code)]` Rust 3D ocean circulation, conformal **Spilhaus World Ocean Map in a Square** (`ESRI:54099` / `PROJ +proj=spilhaus`), and **NOAA ENSO (El Niño / La Niña)** desktop workstation built on `eframe`/`egui` `0.31`, `wgpu` `24.0` WGSL shaders, `glam` `0.29`, `rusqlite` `0.32` (WAL mode), and `clap` `4.5`.

## Core Capabilities

- **High-Fidelity 3D Rotating Ocean Globe & Seamless Spilhaus Morph (`wgpu` `24.0` WGSL)**:
  - Parameterizes a `160 × 160` ocean surface mesh directly in conformal Spilhaus square space `(u, v) ∈ [-K(1/2), +K(1/2)]²` mapped onto the 3D sphere via the exact Bulirsch AGM Jacobi elliptic inverse (`sn, cn, dn`), enabling tear-free continuous 3D morphing (`0%` globe ↔ `100%` Spilhaus World Ocean Square) along the continental East Asia–Bering Strait–South America seam.
  - Renders bathymetric shelf-to-abyss shading, Fresnel sun-glint, Rayleigh limb glow, and an on-globe Spilhaus conformal isoline & cut-seam layer.
- **Exact Conformal Spilhaus Projection Engine (`src/spilhaus.rs`)**:
  - Implements the oblique Adams World in a Square II formulation (`lon_0 = 66.94970198°E`, `lat_0 = 49.56371678°S`, `azimuth = 40.17823482°`, `K(1/2) = 1.8540746773013719`) documented by [Šavrič, Burrows & Kennedy (ArcGIS StoryMap)](https://storymaps.arcgis.com/stories/756bcae18d304a1eac140f19f4d5cb3d), [PROJ 9.9 `+proj=spilhaus`](https://proj.org/en/stable/operations/projections/spilhaus.html), and [Torben Jansen (`guyou-peirce-adams-spilhaus`)](https://observablehq.com/@toja/guyou-peirce-adams-spilhaus).
  - Supports on-globe conformal grid overlays, synchronized 2D docked Spilhaus square and `3×3` periodic world-ocean mosaic views, and standalone SVG export (`--export-spilhaus-svg`).
- **Complete 1950–2026 NOAA CPC ERSSTv6 El Niño / La Niña Dataset (`src/enso.rs`)**:
  - Embeds the seasonal 3-month running mean Oceanic Niño Index (ONI) across Niño 3.4 (`5°N–5°S, 170°W–120°W`) from `1950 DJF` through `2026 ASO (+2.5°C)` from [NOAA Climate Prediction Center ONI v6](https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/enso/oni/v6/), alongside Niño 1+2, Niño 3, Niño 3.4, and Niño 4 region overlays ([UCAR Climate Data Guide](https://climatedataguide.ucar.edu/climate-data/nino-sst-indices-nino-12-3-34-4-oni-and-tni)).
  - Couples ONI forcing to equatorial Pacific `20°C` isotherm thermocline depth (`120°E–80°W`, flattening toward `~150 m` in the eastern Pacific during strong El Niño events per [NOAA PMEL](https://www.pmel.noaa.gov/elnino/el-nino-and-sea-level/)) and coastally trapped Kelvin wave warming along the North American West Coast ([NOAA Fisheries](https://www.fisheries.noaa.gov/feature-story/kelvin-waves-explained)).
- **Calm Lagrangian Ocean Current Ribbons (`src/currents.rs`)**:
  - Simulates 4,200 smooth, non-strobing RK2 Lagrangian drifter streamlines across 14 major global currents grounded in [NASA/ESR OSCAR `0.25°`](https://catalog.data.gov/dataset/ocean-surface-current-analyses-real-time-oscar-surface-currents-final-0-25-degree-version-) and [NOAA Ocean Service](https://oceanservice.noaa.gov/education/tutorial_currents/04currents3.html) observations, including the [Antarctic Circumpolar Current (`135 Sv`)](https://en.wikipedia.org/wiki/Antarctic_Circumpolar_Current), [Gulf Stream (`30–150 Sv`, up to `2.5 m/s`)](https://en.wikipedia.org/wiki/Sverdrup), [Agulhas Current (`70–100 Sv`, peak `2.45 m/s`)](https://en.wikipedia.org/wiki/Agulhas_Current), [Kuroshio Current (`65 Sv`)](https://en.wikipedia.org/wiki/Kuroshio_Current), [Cromwell Equatorial Undercurrent (`30 Sv`, `1.5 m/s`)](https://en.wikipedia.org/wiki/Cromwell_Current), [North Equatorial Countercurrent (`20–25 Sv`)](https://en.wikipedia.org/wiki/Equatorial_Counter_Current), and [Indonesian Throughflow (`15 Sv`)](https://en.wikipedia.org/wiki/Indonesian_Throughflow).

## Quickstart

```bash
# Run mathematical & simulation unit tests
cargo test

# Launch the interactive 3D wgpu + egui desktop globe workstation
cargo run --release

# Run headless verification & inspect Spilhaus / ENSO diagnostics
cargo run --release -- --self-test
cargo run --release -- --episode super-2026 --project-coord=-125.0,44.6
cargo run --release -- --dump-oni-csv oni_1950_2026.csv
cargo run --release -- --episode super-2015 --export-spilhaus-svg spilhaus_enso_2015.svg
```

## Satellite Raster Option (v1.1.0)

Switch seamlessly between the `0 ms` analytical coupled model and **Gridded Satellite Rasters** (`src/raster.rs`):

- **Live NOAA CoastWatch ERDDAP Fetch (`ureq` HTTPS + `rusqlite` WAL Tile Cache)**:
  - **SST & Daily Anomaly (`0.25°` AVHRR OISST v2.1)**: Queries [NOAA NCEI `ncdcOisst21Agg` / `ncdcOisst21NrtAgg`](https://coastwatch.pfeg.noaa.gov/erddap/griddap/ncdcOisst21Agg.html) (`sst`, `anom`, `ice` on `0.125°–359.875°E`, automatically wrapped to `-180°..+180°`).
  - **Geostrophic Satellite Surface Currents (`0.25°` Multi-Mission Altimetry)**: Queries [NOAA/NESDIS `nesdisSSH1day`](https://coastwatch.pfeg.noaa.gov/erddap/griddap/nesdisSSH1day.graph) (`sla`, `ugos`, `vgos` from Jason-2/3, Sentinel-3A, CryoSat-2, SARAL) and historical [ESR/NASA `jplOscar`](https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplOscar.graph) (`u`, `v`).
- **Local Gridded Raster File Import**:
  - Load offline ERDDAP `.json` tables, gridded `.csv` (`lat,lon,sst,anom,u,v`), or equirectangular `.png`/`.jpg` satellite maps onto both the 3D `wgpu` globe and conformal Spilhaus World Ocean Square.

```bash
# Fetch live NOAA OISST v2.1 + NESDIS geostrophic currents (2.0 deg stride = 8) & cache in SQLite
cargo run --release -- --fetch-erddap 2023-12-15 --erddap-stride 8

# Import a local ERDDAP JSON, CSV, or equirectangular satellite PNG/JPEG raster
cargo run --release -- --import-raster oisst_subset.json

# Render offscreen 1600x1000 PNG screenshots (1=3D Globe 1997, 2=2D Map 2015, 3=Spilhaus 1973)
cargo run --release -- --headless-screenshot validation/headless_screenshot.png --screenshot-view 1
cargo run --release -- --headless-screenshot validation/headless_screenshot_2.png --screenshot-view 2
cargo run --release -- --headless-screenshot validation/headless_screenshot_3.png --screenshot-view 3
```

## Validation, Screenshots & Preprints

- **Automated Test Suite**: `cargo test --all-targets` (`12/12` checks passing across Spilhaus round-trip accuracy, historical ENSO flavors, 20°C thermocline scaling, ERDDAP parser, and CLI arguments).
- **Validation Artifacts**: See [validation/pelagic_globe_validation_summary.json](validation/pelagic_globe_validation_summary.json) and the three `1600×1000` headless GPU renders under `validation/`.
- **Companion Preprints**: See [`Preprints/`](Preprints/README.md) (`APP-PREPRINT-22`, Zenodo Concept DOI [`10.5281/zenodo.24442089`](https://doi.org/10.5281/zenodo.24442089)) and citation metadata in [`CITATION.cff`](CITATION.cff).

## License & Attribution

Licensed under the **GNU General Public License v3.0 or later (`GPL-3.0-or-later`)**. See [LICENSE](LICENSE) for details.

- **Repository**: <https://github.com/bpwhite/crvscience-pelagic-globe>
- **Zenodo Concept DOI**: <https://doi.org/10.5281/zenodo.24442089>
- **CRV Science Open-Source Apps**: <https://crvscience.com/apps>
- **Copyright**: Copyright CRV Science Network (<https://crvscience.com>)
