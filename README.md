# PelagicGlobe (`pelagic-globe` v1.1.0)

Native `#![forbid(unsafe_code)]` Rust 3D ocean circulation, conformal **Spilhaus World Ocean Map in a Square** (`ESRI:54099` / `PROJ +proj=spilhaus`), and **NOAA ENSO (El Niño / La Niña)** desktop workstation built on `eframe`/`egui` `0.31`, `wgpu` `24.0` WGSL shaders, `glam` `0.29`, `rusqlite` `0.32` (WAL mode), and `clap` `4.5`.

## Documentation

- [User Guide](docs/USER_GUIDE.md) — Workspace walkthrough, interactive controls, live ERDDAP ingestion, bookmarking, and CSV/SVG exports.
- [Scientific Methods & Equations](docs/SCIENTIFIC_METHODS.md) — Mathematical formulation of the conformal Spilhaus projection, Jacobian velocity rotation, spherical RK2 Lagrangian advection, historical ENSO flavors, and tropical Pacific `20°C` thermocline masking.

## Scientific Summary

- **High-Fidelity 3D Rotating Ocean Globe & Seamless Spilhaus Morph**:
  - Parameterizes a conformal Spilhaus square mesh `(u, v) ∈ [-1, 1]²` mapped onto the 3D sphere via the Bulirsch AGM Jacobi elliptic inverse (`sn, cn, dn`), enabling continuous `0% ↔ 100%` morphing along the East Asia–Bering Strait–South America cut seam.
  - Renders solar directional shading, ocean Blinn-Phong specular sun-glint, Rayleigh atmospheric Fresnel limb glow, bathymetric/thermal color ramps, and an optional synchronized bottom-right Spilhaus inset (`Docked Square`, `3×3 Periodic Mosaic`, or `Hidden`).
- **Exact Conformal Spilhaus Engine (`src/spilhaus.rs`)**:
  - Implements the oblique Adams World in a Square II formulation (`lon_0 = 66.94970198°E`, `lat_0 = 49.56371678°S`, `azimuth = 40.17823482°`, `K(1/2) = 1.854074677301372`) with sub-arcsecond (`< 10⁻⁷ deg`) round-trip inversion, unclamped `3×3` periodic mosaic continuation, and local forward Jacobian `∂(u, v)/∂(East, North)` rotation for exported SVG current vectors.
- **Complete 1950–2026 NOAA CPC ERSSTv6 El Niño / La Niña Dataset (`src/enso.rs`)**:
  - Embeds the 3-month running mean Oceanic Niño Index (ONI) in Niño 3.4 from `1950 DJF` through `2026 ASO (+2.5°C)` across 11 landmark historical episodes (`super-2026` through `neutral-1990` at `1990 DJF +0.2°C`).
  - Differentiates Canonical Eastern Pacific (`1997–98`), Central Pacific Modoki (`2009–10`), Mixed Super + NE Pacific Warm Blob (`2015–16`), and Global Heatwave (`2023–24`) flavors, and evaluates tropical Pacific `20°C` isotherm thermocline depth (`120°E–80°W`, `|lat| ≤ 20°`, masked cleanly outside the tropical Pacific).
- **Spherical Surface Lagrangian Drifters & 14 Major Currents (`src/currents.rs`)**:
  - Advects 4,200 RK2 drifters with `1 / cos(lat)` spherical zonal scaling across 14 verified ocean currents (`acc`, `gulf_stream`, `kuroshio`, `agulhas`, `necc_pacific`, `cromwell_euc`, `sec_pacific`, `nec_pacific`, `humboldt`, `california`, `itf`, `eac`, `brazil`, `benguela`).
  - Excludes `ThermalRegime::SubsurfaceEquatorial` (`cromwell_euc`, `91–200 m`) from surface velocity sampling, eliminates shared-vertex polyline double-counting by selecting the closest segment per current, normalizes tangent directions with `cos(lat)` metric weighting, and splits dateline-crossing (`|Δlon| > 180°`) tracks cleanly across `±180°`.
- **Live NOAA CoastWatch ERDDAP & Local Satellite Raster Ingestion (`src/raster.rs`)**:
  - Fetches live HTTPS NOAA OISST v2.1 (`ncdcOisst21Agg` / `ncdcOisst21NrtAgg`, `0.25°`) and NESDIS altimetry currents (`nesdisSSH1day`, `ugos`/`vgos`, 2017–present) with selectable grid stride (`2`–`16`) and SQLite WAL tile caching, or imports local ERDDAP `.json`, gridded `.csv` (`lat,lon,sst,anom,u,v`), and equirectangular `.png`/`.jpg` maps.

## Controls & Hotkeys

All keyboard shortcuts are guarded by `!ctx.wants_keyboard_input()` so typing in bookmark or export text fields never triggers viewport actions.

| Input | Action | Notes |
| --- | --- | --- |
| `Left Drag` | Orbit 3D globe / Pan 2D map | Disables auto-rotate on globe orbit |
| `Scroll Wheel` | Zoom in / out (`0.55×`–`2.8×`) | Smoothly scales globe, 2D map, or Spilhaus square |
| `Double Click` | Reset camera & zoom | Restores default Pacific camera framing |
| `Hover Cursor` | Live telemetry readout | Displays `lon/lat`, surface class, SST, `ΔSST`, `u/v`, `\|V\|`, and `20°C` depth |
| `Space` | Play / Pause timeline | Steps `1950 DJF → 2026 ASO` chronologically |
| `Left` / `Right` | Step 1 season back / forward | Advances or rewinds 3-month ONI window |
| `1` / `2` / `3` | Switch view mode | `1` = 3D Globe, `2` = 2D World Map, `3` = Spilhaus Square |
| `R` | Reset camera | Resets yaw (`2.55 rad`), pitch (`0.14 rad`), pan, and zoom |
| `F1` / `Esc` | Toggle / Close About dialog | Displays GPL-3.0-or-later and repository metadata |

## CLI Flags

| Flag | Default | Description |
| --- | --- | --- |
| `--episode <slug>` | — | Load a historical ENSO preset (`super-2026`, `strong-2023`, `lanina-2020`, `super-2015`, `modoki-2009`, `lanina-1999`, `super-1997`, `lanina-1988`, `super-1982`, `lanina-1973`, `neutral-1990`) |
| `--oni <celsius>` | — | Override active ONI value (`°C`, hyphen values supported) |
| `--self-test` | `false` | Run headless Spilhaus, mesh, and thermocline diagnostic checks and exit |
| `--project-coord <lon,lat>` | — | Project `lon,lat` into normalized Spilhaus `(u, v)` and round-trip inverse |
| `--dump-oni-csv <path>` | — | Write embedded 1950–2026 NOAA CPC ERSSTv6 ONI table as CSV and exit |
| `--export-spilhaus-svg <path>` | — | Export Jacobian-rotated Spilhaus World Ocean SVG (`ESRI:54099`) and exit |
| `--fetch-erddap <date\|last>` | — | Fetch live NOAA OISST v2.1 + NESDIS currents for `YYYY-MM-DD` or `last` |
| `--erddap-stride <2..16>` | `8` | Subsampling stride on `0.25°` ERDDAP grids (`4` = `1.0°`, `8` = `2.0°`) |
| `--import-raster <path>` | — | Import local ERDDAP `.json`, gridded `.csv` (`lat,lon,sst,anom,u,v`), or `.png`/`.jpg` |
| `--headless-screenshot <path>` | — | Render offscreen `1600×1000` RGBA PNG without opening a GUI window |
| `--screenshot-view <1\|2\|3>` | `1` | Headless capture preset (`1` = 3D Globe 1997, `2` = 2D Map 2015, `3` = Spilhaus 1973) |

## Quickstart

```bash
# Run full unit & CLI test suite
cargo test --all-targets

# Launch the interactive 3D wgpu + egui desktop globe workstation
cargo run --release

# Run headless diagnostics, CSV/SVG export, and offscreen GPU captures
cargo run --release -- --self-test
cargo run --release -- --episode super-2026 --project-coord=-125.0,44.6
cargo run --release -- --dump-oni-csv oni_1950_2026.csv
cargo run --release -- --episode super-2015 --export-spilhaus-svg spilhaus_enso_2015.svg
cargo run --release -- --fetch-erddap 2023-12-15 --erddap-stride 8
cargo run --release -- --headless-screenshot validation/headless_screenshot.png --screenshot-view 1
```

## Validation & Quality Metrics

- **Automated Tests (`cargo test --all-targets`)**: `18/18` passing (`17` library unit tests + `1` CLI parser test; `0` failed, `0` ignored), covering Spilhaus origin/seam/round-trip accuracy, Jacobian velocity vector rotation & unclamped `3×3` mosaic continuity, historical ENSO flavor spatial separation, tropical Pacific `20°C` thermocline masking & `neutral-1990` preset consistency, Cromwell EUC surface exclusion (`u < -0.35 m/s` at `(-150°, -3°)`), ERDDAP `0–360°` & CSV `u,v` raster ingestion (`has_velocity == true`), WGSL shader layer dispatch, and SQLite WAL bookmark CRUD.
- **Validation Artifacts**: See [validation/pelagic_globe_validation_summary.json](validation/pelagic_globe_validation_summary.json), `1600×1000` headless renders in `validation/`, companion preprint [`Preprints/`](Preprints/README.md) (`APP-PREPRINT-22`, Zenodo Concept DOI [`10.5281/zenodo.24442089`](https://doi.org/10.5281/zenodo.24442089)), and [`CITATION.cff`](CITATION.cff).

## License & Attribution

Licensed under the **GNU General Public License v3.0 or later (`GPL-3.0-or-later`)**. See [LICENSE](LICENSE) for full terms.

- **Repository**: <https://github.com/bpwhite/crvscience-pelagic-globe>
- **Zenodo Concept DOI**: <https://doi.org/10.5281/zenodo.24442089>
- **CRV Science Open-Source Apps**: <https://crvscience.com/apps>
- **Copyright**: Copyright CRV Science Network (<https://crvscience.com>)
