# PelagicGlobe User Guide

*Version 1.1.0 — Licensed under GNU General Public License v3.0 or later (`GPL-3.0-or-later`).*

## 1. Workspace Layout

PelagicGlobe organizes its interactive oceanography workstation into four docked panels:

1. **Top Menu Bar**: Displays the active `Year Season | ONI (°C) (Phase)`, quick view buttons (`3D Globe [1]`, `World Map (2D) [2]`, `Spilhaus Square [3]`, `Reset View [R]`), one-click `Save SQLite Bookmark`, and `About [F1]`.
2. **Left Controls Panel (`controls_panel`)**: Houses the `3D <-> Spilhaus Morph` slider, overlay toggles, `Spilhaus Inset Overlay` mode selector (`Docked Square`, `3x3 Mosaic`, `Hidden`), `Field Source` switcher, `ERDDAP Stride` selector, `Data Layer` selector, 11 NOAA CPC historical presets, SQLite bookmarks manager, CSV/SVG export controls, and the 14-current catalog + metadata inspector card.
3. **Central Interactive Viewport**: Renders the 3D sphere, continuous Spilhaus morph, or 2D equirectangular map with Lagrangian drifter ribbons, Niño monitoring boxes, atmospheric Fresnel rim glow, top-left colorbar HUD badge, live cursor telemetry card, and bottom-right Spilhaus inset.
4. **Bottom Chronological & Thermocline Footer**: Combines the interactive 1950–2026 ONI sparkline scrubber (`Play/Pause [Space]`, `< Prev`, `Next >`, Year/Season sliders) on the left with a real-time Equatorial Pacific `20°C` thermocline cross-section (`120°E → 80°W`) on the right.

## 2. View Modes & Camera Navigation

- **3D Globe (`1`)**: Drag to rotate yaw/pitch; scroll wheel zooms (`0.55×`–`2.8×`). Dragging automatically disables auto-rotate.
- **World Map (2D) (`2`)**: Displays an equirectangular `-180°..+180°` projection with drag panning and scroll zoom.
- **Spilhaus Square (`3`)**: Morphs cleanly to `1.0` in conformal `ESRI:54099` coordinates (`[-1, 1]²`) where all world oceans connect without cutting across the Pacific, Atlantic, or Indian basins.
- **Spilhaus Inset Overlay**: In the left panel, choose between:
  - `Docked Square`: A `184×184 px` synchronized single Spilhaus square in the bottom-right corner.
  - `3x3 Mosaic`: A `222×222 px` periodic `3×3` Adams/Spilhaus world-ocean mosaic highlighting the central primary square in gold.
  - `Hidden`: Maximizes unobstructed canvas area.
- **Reset View (`R` or Double-Click)**: Restores default Pacific-centered camera orientation and `1.0×` zoom.

## 3. Cursor Hover Telemetry & Current Metadata Inspector

- **Live Cursor Telemetry**: Hover anywhere over the 3D globe, 2D map, or Spilhaus square to display a crosshair and floating telemetry readout (`Cursor: lat, lon [Ocean/Land/Polar Ice]`, `SST (°C)` & `ΔSST`, `|V|` & `u/v` in `m/s`, and `20°C Thermocline` depth in meters). Outside the tropical Pacific (`|lat| > 20°` or Atlantic/Indian basins), thermocline depth explicitly reports `Masked (outside Trop. Pac.)`.
- **Major Ocean Currents Inspector**: Click any of the 14 currents in the left panel to center the camera on that jet and open its metadata inspector card displaying **Basin**, **Thermal Regime**, **Volume Transport (`Sv`)**, **Peak Speed (`m/s`)**, **ENSO Response**, and a clickable primary literature/agency link. Note that the **Cromwell Current (`cromwell_euc`)** is tagged as `SubsurfaceEquatorial (91–200 m)` and is excluded from surface drifter velocities.

## 4. Analytical Model vs. Live NOAA ERDDAP & Local Rasters

1. **Analytical Model (Default)**: Zero-latency 1950–2026 spatial SST anomaly, seasonal climatology, and coupled circulation field.
2. **Live NOAA ERDDAP Fetch**:
   - Choose an `ERDDAP Stride` (`2 = 0.50°`, `4 = 1.00°`, `8 = 2.00°`, `12 = 3.00°`, `16 = 4.00°`).
   - Click `Fetch NOAA Satellite Raster (YYYY-MM)` to download NOAA OISST v2.1 (`ncdcOisst21Agg` / `ncdcOisst21NrtAgg`, available September 1981–present) and NESDIS altimetric geostrophic currents (`nesdisSSH1day`, 2017–present) on a background thread.
   - Downloaded grids are cached locally in `pelagic_globe_bookmarks.sqlite3` keyed by date and stride.
3. **CLI Local Raster Import**:
   - Load offline ERDDAP `.json`, CSV (`lat,lon,sst,anom,u,v`), or equirectangular `.png`/`.jpg` maps via `--import-raster <path>`.

## 5. Bookmarks & Scientific Exports

- **SQLite WAL Bookmarks**: Enter an optional custom label in the `Saved SQLite Bookmarks` collapsible section and click `Save` (or click `Save SQLite Bookmark` in the top bar). Click any saved `#id` button to restore its exact state, or click `Delete` to remove it from the database.
- **GUI & CLI Exports**: Under `Export CSV & Spilhaus SVG`, customize the filename prefix and click:
  - `Export ONI CSV` → writes `<prefix>_oni_1950_2026.csv` containing all 1950–2026 seasonal ONI entries and phase labels.
  - `Export Spilhaus SVG` → writes `<prefix>_spilhaus_<year>_<season>.svg` containing coastlines and local-Jacobian-rotated surface velocity arrows.
