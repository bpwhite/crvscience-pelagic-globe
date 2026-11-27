// Copyright (C) 2026 CRV Science Network (https://crvscience.com)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://gnu.org>.

#![forbid(unsafe_code)]

//! Native `eframe`/`egui` 0.31 desktop GUI (strictly using `CornerRadius::same(u8)` / `CornerRadius::ZERO`).

use eframe::egui::{self, Color32, CornerRadius, Pos2, Rect, RichText, Sense, Stroke, Vec2};

use crate::coastlines::{WorldLandMask, CONTINENT_POLYLINES};
use crate::currents::{sample_velocity_ms, MAJOR_CURRENTS};
use crate::enso::{
    build_enso_comparative_digest, classify_historical_flavor, climatological_sst_c,
    equatorial_20c_isotherm_depth_historical_m, latest_valid_season_for_year, lookup_oni,
    lookup_oni_clamped, spatial_sst_anomaly_historical, tropical_pacific_20c_isotherm_depth_m,
    EnsoPhase, END_YEAR, ENSO_PRESETS, NINO_BOXES, SEASONS, START_YEAR,
};
use crate::raster::{
    erddap_oisst_supported, fetch_live_noaa_erddap, FieldSource, SatelliteRasterGrid,
};
use crate::spilhaus::{
    lonlat_to_spilhaus_norm, sample_boundary_seam, spilhaus_norm_to_lonlat,
    spilhaus_norm_to_lonlat_unclamped,
};
use crate::storage::BookmarkStore;
use crate::{export_spilhaus_svg, DataLayer, SimulationSnapshot, SpilhausInsetMode};

const SPILHAUS_GRID_N: usize = 96;
const EQUIRECT_COLS: usize = 120;
const EQUIRECT_ROWS: usize = 60;
const INSET_GRID_N: usize = 48;
const MOSAIC_GRID_N: usize = 72;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MapViewMode {
    Globe3D,
    WorldMap2D,
    SpilhausSquare,
}

#[derive(Clone, Copy)]
struct MapGridVertex {
    lon_deg: f64,
    lat_deg: f64,
    cos_lat: f32,
    sin_lat: f32,
    lon_rad: f32,
    su: f32,
    sv: f32,
    is_land: bool,
    is_ice: bool,
}

#[derive(Clone, Copy)]
struct ProjectedCoastPt {
    lon_deg: f64,
    cos_lat: f32,
    sin_lat: f32,
    lon_rad: f32,
    su: f32,
    sv: f32,
}

#[derive(Clone, Copy)]
struct ProjectedScreenVert {
    pos: Pos2,
    z3: f32,
    lit: f32,
    spec: f32,
    fresnel: f32,
    valid: bool,
}

#[derive(Clone, Copy)]
struct HoverTelemetry {
    lon_deg: f64,
    lat_deg: f64,
    sst_c: f64,
    anom_c: f64,
    u_ms: f64,
    v_ms: f64,
    thermo_m: Option<f64>,
    is_land: bool,
    is_ice: bool,
}

fn is_polar_ice_cell(lon_deg: f64, lat_deg: f64) -> bool {
    lat_deg < -65.0 || (lat_deg > 60.0 && (-60.0..=-14.0).contains(&lon_deg))
}

fn make_grid_vertex(lon_deg: f64, lat_deg: f64, su: f32, sv: f32, land_mask: &WorldLandMask) -> MapGridVertex {
    let is_land = land_mask.is_land(lon_deg, lat_deg);
    let is_ice = is_land && is_polar_ice_cell(lon_deg, lat_deg);
    let lon_rad = (lon_deg as f32).to_radians();
    let lat_rad = (lat_deg as f32).to_radians();
    MapGridVertex {
        lon_deg,
        lat_deg,
        cos_lat: lat_rad.cos(),
        sin_lat: lat_rad.sin(),
        lon_rad,
        su,
        sv,
        is_land,
        is_ice,
    }
}

fn build_spilhaus_grid(n: usize, land_mask: &WorldLandMask) -> Vec<MapGridVertex> {
    let mut verts = Vec::with_capacity((n + 1) * (n + 1));
    for r in 0..=n {
        let v = -0.996 + 1.992 * (r as f64 / n as f64);
        for c in 0..=n {
            let u = -0.996 + 1.992 * (c as f64 / n as f64);
            let (lon_deg, lat_deg) = spilhaus_norm_to_lonlat(u, v);
            verts.push(make_grid_vertex(lon_deg, lat_deg, u as f32, v as f32, land_mask));
        }
    }
    verts
}

fn build_spilhaus_mosaic_grid(n: usize, land_mask: &WorldLandMask) -> Vec<MapGridVertex> {
    let mut verts = Vec::with_capacity((n + 1) * (n + 1));
    for r in 0..=n {
        let v = -2.994 + 5.988 * (r as f64 / n as f64);
        for c in 0..=n {
            let u = -2.994 + 5.988 * (c as f64 / n as f64);
            let (lon_deg, lat_deg) = spilhaus_norm_to_lonlat_unclamped(u, v);
            verts.push(make_grid_vertex(lon_deg, lat_deg, u as f32, v as f32, land_mask));
        }
    }
    verts
}

fn build_equirect_grid(cols: usize, rows: usize, land_mask: &WorldLandMask) -> Vec<MapGridVertex> {
    let mut verts = Vec::with_capacity((cols + 1) * (rows + 1));
    for r in 0..=rows {
        let lat_deg = 85.0 - 170.0 * (r as f64 / rows as f64);
        for c in 0..=cols {
            let lon_deg = -180.0 + 360.0 * (c as f64 / cols as f64);
            verts.push(make_grid_vertex(lon_deg, lat_deg, 0.0, 0.0, land_mask));
        }
    }
    verts
}

fn build_ocean_seeds(land_mask: &WorldLandMask) -> Vec<(f64, f64)> {
    let mut seeds = Vec::with_capacity(6000);
    for i in 0..16000_usize {
        let u = (((i * 73 + 17) % 9973) as f64 / 9973.0) * 1.96 - 0.98;
        let v = (((i * 199 + 41) % 9967) as f64 / 9967.0) * 1.96 - 0.98;
        let (lon, lat) = spilhaus_norm_to_lonlat(u, v);
        if (-82.0..=82.0).contains(&lat) && !land_mask.is_land(lon, lat) {
            seeds.push((lon, lat));
        }
    }
    if seeds.is_empty() {
        seeds.push((-145.0, 0.0));
    }
    seeds
}

fn wrap_lon_180_preserving_endpoints(lon_deg: f64) -> f64 {
    if (-180.0..=180.0).contains(&lon_deg) {
        lon_deg
    } else {
        (lon_deg + 180.0).rem_euclid(360.0) - 180.0
    }
}

fn lerp_rgb(a: [f32; 3], b: [f32; 3], t: f32) -> [f32; 3] {
    let u = t.clamp(0.0, 1.0);
    [
        a[0] + (b[0] - a[0]) * u,
        a[1] + (b[1] - a[1]) * u,
        a[2] + (b[2] - a[2]) * u,
    ]
}

fn to_color32(rgb: [f32; 3], lit: f32) -> Color32 {
    Color32::from_rgb(
        (rgb[0] * lit).round().clamp(0.0, 255.0) as u8,
        (rgb[1] * lit).round().clamp(0.0, 255.0) as u8,
        (rgb[2] * lit).round().clamp(0.0, 255.0) as u8,
    )
}

fn land_color(_lon_deg: f64, lat_deg: f64, is_ice: bool) -> [f32; 3] {
    if is_ice {
        [196.0, 212.0, 224.0]
    } else if (15.0..=34.0).contains(&lat_deg.abs()) {
        [114.0, 104.0, 80.0]
    } else {
        [58.0, 86.0, 66.0]
    }
}

fn diverged_anom_rgb(anom: f32) -> [f32; 3] {
    if anom <= -2.0 {
        let t = ((anom + 3.5) / 1.5).clamp(0.0, 1.0);
        lerp_rgb([16.0, 48.0, 148.0], [28.0, 148.0, 232.0], t)
    } else if anom <= -0.8 {
        let t = ((anom + 2.0) / 1.2).clamp(0.0, 1.0);
        lerp_rgb([28.0, 148.0, 232.0], [64.0, 212.0, 224.0], t)
    } else if anom < 0.0 {
        let t = ((anom + 0.8) / 0.8).clamp(0.0, 1.0);
        lerp_rgb([64.0, 212.0, 224.0], [11.0, 36.0, 70.0], t)
    } else if anom <= 0.8 {
        let t = (anom / 0.8).clamp(0.0, 1.0);
        lerp_rgb([11.0, 36.0, 70.0], [246.0, 198.0, 68.0], t)
    } else if anom <= 2.0 {
        let t = ((anom - 0.8) / 1.2).clamp(0.0, 1.0);
        lerp_rgb([246.0, 198.0, 68.0], [248.0, 102.0, 40.0], t)
    } else {
        let t = ((anom - 2.0) / 1.8).clamp(0.0, 1.0);
        lerp_rgb([248.0, 102.0, 40.0], [178.0, 24.0, 34.0], t)
    }
}

pub struct PelagicGlobeApp {
    pub state: SimulationSnapshot,
    store: Option<BookmarkStore>,
    saved_bookmarks: Vec<(i64, String, SimulationSnapshot)>,
    land_mask: WorldLandMask,
    view_mode: MapViewMode,
    show_world_map: bool,
    show_graticule: bool,
    pan_2d: Vec2,
    spilhaus_grid: Vec<MapGridVertex>,
    equirect_grid: Vec<MapGridVertex>,
    inset_grid: Vec<MapGridVertex>,
    mosaic_grid: Vec<MapGridVertex>,
    coastline_rings: Vec<Vec<ProjectedCoastPt>>,
    seam_cache: Vec<(f64, f64)>,
    ocean_seeds: Vec<(f64, f64)>,
    drifters: Vec<(f64, f64)>,
    active_raster: Option<SatelliteRasterGrid>,
    erddap_rx: Option<std::sync::mpsc::Receiver<(String, Result<SatelliteRasterGrid, String>)>>,
    erddap_in_flight: Option<String>,
    timeline_playing: bool,
    timeline_accum_sec: f32,
    bookmark_name_draft: String,
    export_prefix: String,
    status_note: String,
    pub show_about_window: bool,
    pub command_palette_open: bool,
    command_palette_just_opened: bool,
    command_palette_query: String,
}

impl PelagicGlobeApp {
    pub fn new(initial: SimulationSnapshot) -> Self {
        Self::new_with_raster(initial, None)
    }

    pub fn new_headless(
        initial: SimulationSnapshot,
        raster: Option<SatelliteRasterGrid>,
    ) -> Self {
        Self::new_internal(initial, raster, false)
    }

    pub fn new_with_raster(
        initial: SimulationSnapshot,
        raster: Option<SatelliteRasterGrid>,
    ) -> Self {
        Self::new_internal(initial, raster, true)
    }

    fn new_internal(
        initial: SimulationSnapshot,
        raster: Option<SatelliteRasterGrid>,
        open_bookmark_store: bool,
    ) -> Self {
        let store = if open_bookmark_store {
            BookmarkStore::open("pelagic_globe_bookmarks.sqlite3").ok()
        } else {
            None
        };
        let saved_bookmarks = store
            .as_ref()
            .and_then(|s| s.list_bookmarks().ok())
            .unwrap_or_default();
        let land_mask = WorldLandMask::build();
        let spilhaus_grid = build_spilhaus_grid(SPILHAUS_GRID_N, &land_mask);
        let equirect_grid = build_equirect_grid(EQUIRECT_COLS, EQUIRECT_ROWS, &land_mask);
        let inset_grid = build_spilhaus_grid(INSET_GRID_N, &land_mask);
        let mosaic_grid = build_spilhaus_mosaic_grid(MOSAIC_GRID_N, &land_mask);

        let coastline_rings = CONTINENT_POLYLINES
            .iter()
            .map(|ring| {
                ring.iter()
                    .map(|&(lat_deg, lon_deg)| {
                        let (su, sv) = lonlat_to_spilhaus_norm(lon_deg, lat_deg);
                        let lat_rad = (lat_deg as f32).to_radians();
                        ProjectedCoastPt {
                            lon_deg,
                            cos_lat: lat_rad.cos(),
                            sin_lat: lat_rad.sin(),
                            lon_rad: (lon_deg as f32).to_radians(),
                            su: su as f32,
                            sv: sv as f32,
                        }
                    })
                    .collect()
            })
            .collect();

        let seam_cache = sample_boundary_seam(64);
        let ocean_seeds = build_ocean_seeds(&land_mask);
        let mut drifters = Vec::with_capacity(initial.particle_count);
        for i in 0..initial.particle_count {
            let idx = (i * 17) % ocean_seeds.len();
            drifters.push(ocean_seeds[idx]);
        }

        let view_mode = if initial.spilhaus_morph >= 0.99 {
            MapViewMode::SpilhausSquare
        } else {
            MapViewMode::Globe3D
        };

        let status_note = match &raster {
            Some(grid) => format!(
                "Loaded satellite raster ({}, {} ocean cells).",
                grid.source_label,
                grid.valid_cell_count()
            ),
            None => format!(
                "{} {} | ONI {:+.1}°C ({}) — {}",
                initial.year,
                SEASONS[initial.season_index.min(11)],
                initial.oni_celsius,
                EnsoPhase::from_oni(initial.oni_celsius).label(),
                classify_historical_flavor(initial.year, initial.oni_celsius).label()
            ),
        };

        Self {
            state: initial,
            store,
            saved_bookmarks,
            land_mask,
            view_mode,
            show_world_map: true,
            show_graticule: true,
            pan_2d: Vec2::ZERO,
            spilhaus_grid,
            equirect_grid,
            inset_grid,
            mosaic_grid,
            coastline_rings,
            seam_cache,
            ocean_seeds,
            drifters,
            active_raster: raster,
            erddap_rx: None,
            erddap_in_flight: None,
            timeline_playing: false,
            timeline_accum_sec: 0.0,
            bookmark_name_draft: String::new(),
            export_prefix: "pelagic_globe".to_string(),
            status_note,
            show_about_window: false,
            command_palette_open: false,
            command_palette_just_opened: false,
            command_palette_query: String::new(),
        }
    }

    fn apply_historical_date(
        &mut self,
        year: i32,
        season_index: usize,
        focus_pacific: bool,
        custom_note: Option<String>,
    ) {
        let y = year.clamp(START_YEAR, END_YEAR);
        let (sidx, oni) = lookup_oni_clamped(y, season_index);
        self.state.year = y;
        self.state.season_index = sidx;
        self.state.oni_celsius = oni;
        self.state.field_source = FieldSource::AnalyticalModel;
        if focus_pacific {
            self.state.yaw_rad = 2.40;
            self.state.pitch_rad = 0.08;
            self.state.auto_rotate = false;
        }
        let flavor = classify_historical_flavor(y, oni);
        let phase = EnsoPhase::from_oni(oni);
        self.status_note = custom_note.unwrap_or_else(|| {
            format!(
                "{} {} | ONI {:+.1}°C ({}) — {}",
                y,
                SEASONS[sidx],
                oni,
                phase.label(),
                flavor.label()
            )
        });
    }

    fn step_season(&mut self, delta: isize) {
        self.timeline_playing = false;
        let (ny, ns) = if delta < 0 {
            if self.state.season_index > 0 {
                (self.state.year, self.state.season_index - 1)
            } else if self.state.year > START_YEAR {
                (
                    self.state.year - 1,
                    latest_valid_season_for_year(self.state.year - 1),
                )
            } else {
                (END_YEAR, latest_valid_season_for_year(END_YEAR))
            }
        } else {
            let max_s = latest_valid_season_for_year(self.state.year);
            if self.state.season_index < max_s {
                (self.state.year, self.state.season_index + 1)
            } else if self.state.year < END_YEAR {
                (self.state.year + 1, 0)
            } else {
                (START_YEAR, 0)
            }
        };
        self.apply_historical_date(ny, ns, false, None);
    }

    fn reset_camera(&mut self) {
        self.pan_2d = Vec2::ZERO;
        self.state.zoom = 1.0;
        self.state.yaw_rad = 2.55;
        self.state.pitch_rad = 0.14;
        self.state.auto_rotate = false;
    }

    fn trigger_erddap_fetch(&mut self) {
        if self.erddap_in_flight.is_some() {
            return;
        }
        let year = self.state.year;
        let month = (self.state.season_index + 1).clamp(1, 12);
        if !erddap_oisst_supported(year, month) {
            self.state.field_source = FieldSource::AnalyticalModel;
            self.status_note = format!(
                "NOAA OISST v2.1 satellite archive begins Sep 1981 — showing NOAA CPC ERSSTv6 historical model for {} {}.",
                year, SEASONS[self.state.season_index.min(11)]
            );
            return;
        }

        let iso = if year >= 2026 && month >= 10 {
            "last".to_string()
        } else {
            format!("{year:04}-{month:02}-15")
        };
        let stride = self.state.erddap_stride.clamp(2, 16);
        let cache_key = format!("{iso}_s{stride}");

        if let Some(cached) = self.store.as_ref().and_then(|s| {
            s.get_cached_raster(&cache_key)
                .ok()
                .flatten()
                .or_else(|| {
                    if stride == 8 {
                        s.get_cached_raster(&iso).ok().flatten()
                    } else {
                        None
                    }
                })
        }) {
            let cells = cached.valid_cell_count();
            let has_vel = cached.has_velocity;
            self.active_raster = Some(cached);
            self.state.field_source = FieldSource::SatelliteRaster;
            self.status_note = if has_vel {
                format!("Loaded cached NOAA ERDDAP raster ({iso}, stride {stride}, {cells} cells, OISST + NESDIS currents).")
            } else {
                format!("Loaded cached NOAA OISST v2.1 raster ({iso}, stride {stride}, {cells} cells; analytical currents pre-2017).")
            };
            return;
        }

        let iso_thread = iso.clone();
        let (tx, rx) = std::sync::mpsc::channel();
        self.erddap_in_flight = Some(cache_key);
        self.erddap_rx = Some(rx);
        self.status_note = format!("Fetching NOAA ERDDAP {iso} (stride {stride}) in background...");
        std::thread::spawn(move || {
            let res = fetch_live_noaa_erddap(&iso_thread, stride);
            let _ = tx.send((iso_thread, res));
        });
    }

    fn sample_effective_velocity(&self, lon_deg: f64, lat_deg: f64) -> (f64, f64) {
        if self.state.field_source == FieldSource::SatelliteRaster {
            if let Some(grid) = &self.active_raster {
                if grid.has_velocity {
                    if let Some((_sst, _anom, u, v, _sla)) = grid.sample_bilinear(lon_deg, lat_deg) {
                        return (u, v);
                    }
                }
            }
        }
        sample_velocity_ms(lon_deg, lat_deg, self.state.oni_celsius)
    }

    fn sample_effective_sst_anom(&self, lon_deg: f64, lat_deg: f64) -> (f64, f64) {
        if self.state.field_source == FieldSource::SatelliteRaster {
            if let Some(grid) = &self.active_raster {
                if let Some((sst, anom, _u, _v, _sla)) = grid.sample_bilinear(lon_deg, lat_deg) {
                    return (sst, anom);
                }
            }
        }
        let anom = spatial_sst_anomaly_historical(
            lon_deg,
            lat_deg,
            self.state.year,
            self.state.season_index,
            self.state.oni_celsius,
        );
        let base = climatological_sst_c(lon_deg, lat_deg, self.state.season_index);
        ((base + anom).clamp(-2.0, 33.5), anom)
    }

    fn sample_ocean_telemetry(&self, lon_deg: f64, lat_deg: f64) -> HoverTelemetry {
        let is_land = self.land_mask.is_land(lon_deg, lat_deg);
        let is_ice = is_land && is_polar_ice_cell(lon_deg, lat_deg);
        let (sst_c, anom_c) = self.sample_effective_sst_anom(lon_deg, lat_deg);
        let (u_ms, v_ms) = self.sample_effective_velocity(lon_deg, lat_deg);
        let thermo_m = if is_land {
            None
        } else {
            tropical_pacific_20c_isotherm_depth_m(
                lon_deg,
                lat_deg,
                self.state.year,
                self.state.oni_celsius,
            )
        };
        HoverTelemetry {
            lon_deg,
            lat_deg,
            sst_c,
            anom_c,
            u_ms,
            v_ms,
            thermo_m,
            is_land,
            is_ice,
        }
    }

    fn screen_pos_to_lonlat(&self, pos: Pos2, rect: Rect) -> Option<(f64, f64)> {
        if !rect.contains(pos) {
            return None;
        }
        if self.view_mode == MapViewMode::WorldMap2D {
            let map_rect = self.world_map_rect(rect);
            if !map_rect.contains(pos) {
                return None;
            }
            let tx = f64::from((pos.x - map_rect.left()) / map_rect.width().max(1.0));
            let ty = f64::from((pos.y - map_rect.top()) / map_rect.height().max(1.0));
            let lon = (-180.0 + tx * 360.0).clamp(-180.0, 180.0);
            let lat = (85.0 - ty * 170.0).clamp(-85.0, 85.0);
            return Some((lon, lat));
        }

        let morph = self.effective_morph();
        let center = rect.center();
        let radius = rect.width().min(rect.height()) * 0.42 * self.state.zoom;
        if radius < 1.0 {
            return None;
        }
        let bx = (pos.x - center.x) / radius;
        let by = (pos.y - center.y) / radius;

        if morph >= 0.99 {
            let su = f64::from(bx / 0.96);
            let sv = f64::from(-by / 0.96);
            if su.abs() <= 1.0 && sv.abs() <= 1.0 {
                return Some(spilhaus_norm_to_lonlat(su, sv));
            }
            return None;
        }

        if morph < 0.04 {
            let r2 = bx * bx + by * by;
            if r2 > 1.0 {
                return None;
            }
            let x3 = bx;
            let y3 = -by;
            let z3 = (1.0 - r2).max(0.0).sqrt();
            let cos_p = self.state.pitch_rad.cos();
            let sin_p = self.state.pitch_rad.sin();
            let sin_lat = (y3 * cos_p + z3 * sin_p).clamp(-1.0, 1.0);
            let cos_lat_cos_l = z3 * cos_p - y3 * sin_p;
            let cos_lat_sin_l = x3;
            let lat_deg = f64::from(sin_lat.asin().to_degrees());
            let l_rad = f64::from(cos_lat_sin_l.atan2(cos_lat_cos_l));
            let lon_deg = (l_rad - f64::from(self.state.yaw_rad)).to_degrees();
            let lon_wrapped = (lon_deg + 180.0).rem_euclid(360.0) - 180.0;
            return Some((lon_wrapped, lat_deg));
        }

        // Intermediate morph: find nearest visible projected surface grid vertex
        let mut best: Option<(f32, f64, f64)> = None;
        let max_d2 = (radius * 0.075).powi(2);
        for v in &self.spilhaus_grid {
            let pv = self.project_components_3d_spilhaus(
                v.cos_lat,
                v.sin_lat,
                v.lon_rad,
                v.su,
                v.sv,
                morph,
                center,
                radius,
            );
            if !pv.valid {
                continue;
            }
            let d2 = pv.pos.distance_sq(pos);
            if d2 <= max_d2 && best.map(|(bd2, _, _)| d2 < bd2).unwrap_or(true) {
                best = Some((d2, v.lon_deg, v.lat_deg));
            }
        }
        best.map(|(_, lon, lat)| (lon, lat))
    }

    fn effective_morph(&self) -> f32 {
        match self.view_mode {
            MapViewMode::SpilhausSquare => 1.0,
            MapViewMode::WorldMap2D => 0.0,
            MapViewMode::Globe3D => self.state.spilhaus_morph.clamp(0.0, 1.0),
        }
    }

    fn ocean_color(&self, lon_deg: f64, lat_deg: f64) -> [f32; 3] {
        let (sst_val, anom_val) = self.sample_effective_sst_anom(lon_deg, lat_deg);
        match self.state.data_layer {
            DataLayer::EnsoSstAnomaly => diverged_anom_rgb(anom_val as f32),
            DataLayer::CurrentVelocity => {
                let (vx, vy) = self.sample_effective_velocity(lon_deg, lat_deg);
                let t = ((vx * vx + vy * vy).sqrt() as f32 / 1.4).clamp(0.0, 1.0);
                if t < 0.5 {
                    lerp_rgb([8.0, 24.0, 54.0], [24.0, 154.0, 188.0], t * 2.0)
                } else {
                    lerp_rgb([24.0, 154.0, 188.0], [250.0, 204.0, 21.0], (t - 0.5) * 2.0)
                }
            }
            DataLayer::AbsoluteSst => {
                let t = ((sst_val as f32 + 2.0) / 34.0).clamp(0.0, 1.0);
                if t < 0.5 {
                    lerp_rgb([16.0, 32.0, 98.0], [20.0, 148.0, 162.0], t * 2.0)
                } else {
                    lerp_rgb([20.0, 148.0, 162.0], [246.0, 108.0, 62.0], (t - 0.5) * 2.0)
                }
            }
            DataLayer::ThermoclineDepth => {
                match tropical_pacific_20c_isotherm_depth_m(
                    lon_deg,
                    lat_deg,
                    self.state.year,
                    self.state.oni_celsius,
                ) {
                    Some(depth_f64) => {
                        let depth = depth_f64 as f32;
                        let t = ((depth - 20.0) / 170.0).clamp(0.0, 1.0);
                        if t < 0.5 {
                            lerp_rgb([28.0, 34.0, 108.0], [18.0, 158.0, 126.0], t * 2.0)
                        } else {
                            lerp_rgb([18.0, 158.0, 126.0], [244.0, 164.0, 44.0], (t - 0.5) * 2.0)
                        }
                    }
                    None => [14.0, 26.0, 44.0],
                }
            }
        }
    }

    fn world_map_rect(&self, rect: Rect) -> Rect {
        let pad_w = (rect.width() - 28.0).max(160.0);
        let pad_h = (rect.height() - 24.0).max(80.0);
        let (base_w, base_h) = if pad_w >= 2.0 * pad_h {
            (pad_h * 2.0, pad_h)
        } else {
            (pad_w, pad_w * 0.5)
        };
        let z = self.state.zoom.clamp(0.55, 2.8);
        let w = base_w * z;
        let h = base_h * z;
        let center = rect.center() + self.pan_2d;
        Rect::from_center_size(center, Vec2::new(w, h))
    }

    fn project_components_3d_spilhaus(
        &self,
        cos_lat: f32,
        sin_lat: f32,
        lon_rad: f32,
        su: f32,
        sv: f32,
        morph: f32,
        center: Pos2,
        radius: f32,
    ) -> ProjectedScreenVert {
        let cos_p = self.state.pitch_rad.cos();
        let sin_p = self.state.pitch_rad.sin();
        let l = lon_rad + self.state.yaw_rad;
        let (sin_l, cos_l) = l.sin_cos();

        let x3 = cos_lat * sin_l;
        let y3 = sin_lat * cos_p - cos_lat * cos_l * sin_p;
        let z3 = sin_lat * sin_p + cos_lat * cos_l * cos_p;

        let valid = !(morph < 0.04 && z3 < -0.04);
        let sx = su * 0.96;
        let sy = -sv * 0.96;
        let bx = x3 * (1.0 - morph) + sx * morph;
        let by = -y3 * (1.0 - morph) + sy * morph;

        // Solar directional diffuse + Blinn-Phong specular glint + Rayleigh Fresnel limb glow
        let lx = 0.54_f32;
        let ly = 0.46_f32;
        let lz = 0.704_f32;
        let raw_ndotl = (x3 * lx + y3 * ly + z3 * lz).clamp(0.0, 1.0);
        let lit = (0.30 + 0.70 * raw_ndotl) * (1.0 - morph) + morph;

        // Half-vector of normalized light + view (0, 0, 1)
        let hx = 0.292_f32;
        let hy = 0.249_f32;
        let hz = 0.923_f32;
        let ndoth = (x3 * hx + y3 * hy + z3 * hz).max(0.0);
        let spec = ndoth.powi(36) * 0.34 * (1.0 - morph);
        let fresnel = (1.0 - z3.max(0.0)).powf(3.2) * 0.32 * (1.0 - morph);

        ProjectedScreenVert {
            pos: Pos2::new(center.x + bx * radius, center.y + by * radius),
            z3,
            lit,
            spec,
            fresnel,
            valid,
        }
    }

    fn project_lonlat(&self, lon_deg: f64, lat_deg: f64, rect: Rect) -> Option<(Pos2, f32)> {
        if self.view_mode == MapViewMode::WorldMap2D {
            let map_rect = self.world_map_rect(rect);
            let lon = wrap_lon_180_preserving_endpoints(lon_deg) as f32;
            let lat = lat_deg.clamp(-85.0, 85.0) as f32;
            let x = map_rect.left() + ((lon + 180.0) / 360.0) * map_rect.width();
            let y = map_rect.top() + ((85.0 - lat) / 170.0) * map_rect.height();
            Some((Pos2::new(x, y), 1.0))
        } else {
            let morph = self.effective_morph();
            let center = rect.center();
            let radius = rect.width().min(rect.height()) * 0.42 * self.state.zoom;
            let (su, sv) = if morph > 0.001 {
                let (u, v) = lonlat_to_spilhaus_norm(lon_deg, lat_deg);
                (u as f32, v as f32)
            } else {
                (0.0, 0.0)
            };
            let lat_r = (lat_deg as f32).to_radians();
            let p = self.project_components_3d_spilhaus(
                lat_r.cos(),
                lat_r.sin(),
                (lon_deg as f32).to_radians(),
                su,
                sv,
                morph,
                center,
                radius,
            );
            if p.valid {
                Some((p.pos, p.z3))
            } else {
                None
            }
        }
    }

    fn step_drifters(&mut self, dt: f64) {
        if !self.state.show_currents {
            return;
        }
        let scale = f64::from(self.state.flow_speed) * dt * 9.0;
        let tick = (self.state.yaw_rad.abs() * 100.0) as usize;
        let seed_len = self.ocean_seeds.len();
        for i in 0..self.drifters.len() {
            let p = self.drifters[i];
            let (u1, v1) = self.sample_effective_velocity(p.0, p.1);
            let inv_cos1 = 1.0 / p.1.to_radians().cos().max(0.15);
            let mid_lon = p.0 + 0.5 * scale * u1 * inv_cos1;
            let mid_lat = (p.1 + 0.5 * scale * v1).clamp(-82.0, 82.0);
            let (u2, v2) = self.sample_effective_velocity(mid_lon, mid_lat);
            let inv_cos2 = 1.0 / mid_lat.to_radians().cos().max(0.15);
            let next_lon = (p.0 + scale * u2 * inv_cos2 + 180.0).rem_euclid(360.0) - 180.0;
            let next_lat = (p.1 + scale * v2).clamp(-82.0, 82.0);
            if self.land_mask.is_land(next_lon, next_lat) || ((i + tick) % 311) == 0 {
                let idx = (i * 37 + tick * 13) % seed_len;
                self.drifters[i] = self.ocean_seeds[idx];
            } else {
                self.drifters[i] = (next_lon, next_lat);
            }
        }
    }

    fn render_surface_mesh_3d_or_spilhaus(&self, painter: &egui::Painter, rect: Rect) {
        let n = SPILHAUS_GRID_N;
        let cols = n + 1;
        let morph = self.effective_morph();
        let center = rect.center();
        let radius = rect.width().min(rect.height()) * 0.42 * self.state.zoom;
        let jump_limit = radius * 0.35;

        let mut projected = Vec::with_capacity(self.spilhaus_grid.len());
        let mut colors = Vec::with_capacity(self.spilhaus_grid.len());
        for v in &self.spilhaus_grid {
            let pv = self.project_components_3d_spilhaus(
                v.cos_lat,
                v.sin_lat,
                v.lon_rad,
                v.su,
                v.sv,
                morph,
                center,
                radius,
            );
            let is_land_vis = v.is_land && self.show_world_map;
            let base_rgb = if is_land_vis {
                land_color(v.lon_deg, v.lat_deg, v.is_ice)
            } else {
                self.ocean_color(v.lon_deg, v.lat_deg)
            };
            let glint = if is_land_vis { 0.0 } else { pv.spec * 96.0 };
            let limb = if is_land_vis {
                pv.fresnel * 0.35
            } else {
                pv.fresnel
            };
            let shaded = [
                base_rgb[0] * pv.lit + glint + 64.0 * limb,
                base_rgb[1] * pv.lit + glint + 168.0 * limb,
                base_rgb[2] * pv.lit + glint * 1.05 + 245.0 * limb,
            ];
            colors.push(to_color32(shaded, 1.0));
            projected.push(pv);
        }

        let mut quads: Vec<(f32, u32, u32, u32, u32)> = Vec::with_capacity(n * n);
        for row in 0..n {
            for col in 0..n {
                let i00 = row * cols + col;
                let i10 = i00 + 1;
                let i01 = (row + 1) * cols + col;
                let i11 = i01 + 1;
                let p00 = projected[i00];
                let p10 = projected[i10];
                let p11 = projected[i11];
                let p01 = projected[i01];
                if !p00.valid && !p10.valid && !p11.valid && !p01.valid {
                    continue;
                }
                let min_z = p00.z3.min(p10.z3).min(p11.z3).min(p01.z3);
                if morph < 0.04 && min_z < -0.04 {
                    continue;
                }
                if p00.pos.distance(p10.pos) > jump_limit
                    || p10.pos.distance(p11.pos) > jump_limit
                    || p11.pos.distance(p01.pos) > jump_limit
                    || p01.pos.distance(p00.pos) > jump_limit
                {
                    continue;
                }
                let avg_z = 0.25 * (p00.z3 + p10.z3 + p11.z3 + p01.z3);
                quads.push((avg_z, i00 as u32, i10 as u32, i11 as u32, i01 as u32));
            }
        }

        if morph < 0.99 {
            quads.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        }

        let mut mesh = egui::epaint::Mesh::default();
        for (i, pv) in projected.iter().enumerate() {
            mesh.colored_vertex(pv.pos, colors[i]);
        }
        for &(_, i00, i10, i11, i01) in &quads {
            mesh.add_triangle(i00, i10, i11);
            mesh.add_triangle(i00, i11, i01);
        }
        painter.add(egui::Shape::mesh(mesh));
    }

    fn render_surface_mesh_2d(&self, painter: &egui::Painter, rect: Rect) {
        let map_rect = self.world_map_rect(rect);
        painter.rect_filled(map_rect, CornerRadius::ZERO, Color32::from_rgb(8, 22, 42));

        let cols = EQUIRECT_COLS;
        let rows = EQUIRECT_ROWS;
        let stride = cols + 1;
        let mut mesh = egui::epaint::Mesh::default();
        for r in 0..=rows {
            let ty = r as f32 / rows as f32;
            let y = map_rect.top() + ty * map_rect.height();
            for c in 0..=cols {
                let tx = c as f32 / cols as f32;
                let x = map_rect.left() + tx * map_rect.width();
                let v = &self.equirect_grid[r * stride + c];
                let rgb = if v.is_land && self.show_world_map {
                    land_color(v.lon_deg, v.lat_deg, v.is_ice)
                } else {
                    self.ocean_color(v.lon_deg, v.lat_deg)
                };
                mesh.colored_vertex(Pos2::new(x, y), to_color32(rgb, 1.0));
            }
        }
        for r in 0..rows {
            for c in 0..cols {
                let i00 = (r * stride + c) as u32;
                let i10 = i00 + 1;
                let i01 = ((r + 1) * stride + c) as u32;
                let i11 = i01 + 1;
                mesh.add_triangle(i00, i10, i11);
                mesh.add_triangle(i00, i11, i01);
            }
        }
        painter.add(egui::Shape::mesh(mesh));
        painter.rect_stroke(
            map_rect,
            CornerRadius::ZERO,
            Stroke::new(1.2_f32, Color32::from_rgb(71, 128, 168)),
            egui::StrokeKind::Outside,
        );
    }

    fn render_graticule(&self, painter: &egui::Painter, rect: Rect, jump_limit: f32) {
        if !self.show_graticule {
            return;
        }
        let grat_col = Color32::from_rgba_unmultiplied(148, 210, 240, 46);
        let eq_col = Color32::from_rgba_unmultiplied(125, 211, 252, 85);

        for lat_deg in (-60..=60).step_by(30) {
            let stroke = Stroke::new(if lat_deg == 0 { 1.1_f32 } else { 0.8_f32 }, if lat_deg == 0 { eq_col } else { grat_col });
            let mut prev: Option<Pos2> = None;
            for step in 0..=72 {
                let lon = -180.0 + 360.0 * (step as f64 / 72.0);
                if let Some((p, _)) = self.project_lonlat(lon, f64::from(lat_deg), rect) {
                    if let Some(pa) = prev {
                        if pa.distance(p) <= jump_limit {
                            painter.line_segment([pa, p], stroke);
                        }
                    }
                    prev = Some(p);
                } else {
                    prev = None;
                }
            }
        }

        for lon_deg in (-150..=180).step_by(30) {
            let stroke = Stroke::new(0.8_f32, grat_col);
            let mut prev: Option<Pos2> = None;
            for step in 0..=34 {
                let lat = -85.0 + 170.0 * (step as f64 / 34.0);
                if let Some((p, _)) = self.project_lonlat(f64::from(lon_deg), lat, rect) {
                    if let Some(pa) = prev {
                        if pa.distance(p) <= jump_limit {
                            painter.line_segment([pa, p], stroke);
                        }
                    }
                    prev = Some(p);
                } else {
                    prev = None;
                }
            }
        }
    }

    fn render_coastlines(&self, painter: &egui::Painter, rect: Rect, jump_limit: f32) {
        if !self.show_world_map {
            return;
        }
        let casing_stroke = Stroke::new(2.0_f32, Color32::from_rgba_unmultiplied(10, 18, 28, 185));
        let shore_stroke = Stroke::new(1.1_f32, Color32::from_rgb(218, 226, 206));

        let morph = self.effective_morph();
        let center = rect.center();
        let radius = rect.width().min(rect.height()) * 0.42 * self.state.zoom;
        let map_rect = self.world_map_rect(rect);

        let project_pt = |pt: &ProjectedCoastPt| -> Option<Pos2> {
            if self.view_mode == MapViewMode::WorldMap2D {
                let lon = wrap_lon_180_preserving_endpoints(pt.lon_deg) as f32;
                let lat_deg = pt.sin_lat.asin().to_degrees().clamp(-85.0, 85.0);
                let x = map_rect.left() + ((lon + 180.0) / 360.0) * map_rect.width();
                let y = map_rect.top() + ((85.0 - lat_deg) / 170.0) * map_rect.height();
                Some(Pos2::new(x, y))
            } else {
                let pv = self.project_components_3d_spilhaus(
                    pt.cos_lat,
                    pt.sin_lat,
                    pt.lon_rad,
                    pt.su,
                    pt.sv,
                    morph,
                    center,
                    radius,
                );
                if pv.valid {
                    Some(pv.pos)
                } else {
                    None
                }
            }
        };

        for ring in &self.coastline_rings {
            for seg in ring.windows(2) {
                if (seg[0].lon_deg - seg[1].lon_deg).abs() > 180.0 {
                    continue;
                }
                if let (Some(a), Some(b)) = (project_pt(&seg[0]), project_pt(&seg[1])) {
                    if a.distance(b) <= jump_limit {
                        painter.line_segment([a, b], casing_stroke);
                        painter.line_segment([a, b], shore_stroke);
                    }
                }
            }
        }
    }

    fn render_selected_current_track(&self, painter: &egui::Painter, rect: Rect, jump_limit: f32) {
        let Some(cur_id) = self.state.selected_current_id.as_deref() else {
            return;
        };
        let Some(cur) = MAJOR_CURRENTS.iter().find(|c| c.id == cur_id) else {
            return;
        };
        let track_stroke = Stroke::new(2.4_f32, Color32::from_rgb(250, 204, 21));
        for pair in cur.waypoints.windows(2) {
            let (lon0, lat0) = pair[0];
            let (lon1, lat1) = pair[1];
            let mut lon1_u = lon1;
            if (lon1_u - lon0).abs() > 180.0 {
                if lon1_u < lon0 {
                    lon1_u += 360.0;
                } else {
                    lon1_u -= 360.0;
                }
            }
            let sub_steps = 8_usize;
            for s in 0..sub_steps {
                let t0 = s as f64 / sub_steps as f64;
                let t1 = (s + 1) as f64 / sub_steps as f64;
                let a_lon = lon0 + (lon1_u - lon0) * t0;
                let a_lat = lat0 + (lat1 - lat0) * t0;
                let b_lon = lon0 + (lon1_u - lon0) * t1;
                let b_lat = lat0 + (lat1 - lat0) * t1;
                let wa = wrap_lon_180_preserving_endpoints(a_lon);
                let wb = wrap_lon_180_preserving_endpoints(b_lon);

                if (wb - wa).abs() > 180.0 {
                    let seam_u = if b_lon > a_lon { 180.0 } else { -180.0 };
                    let frac = ((seam_u - a_lon) / (b_lon - a_lon)).clamp(0.0, 1.0);
                    let seam_lat = a_lat + (b_lat - a_lat) * frac;
                    let seam_a = if b_lon > a_lon { 180.0 } else { -180.0 };
                    let seam_b = -seam_a;
                    if let (Some((pa, _)), Some((ps, _))) = (
                        self.project_lonlat(wa, a_lat, rect),
                        self.project_lonlat(seam_a, seam_lat, rect),
                    ) {
                        if pa.distance(ps) <= jump_limit {
                            painter.line_segment([pa, ps], track_stroke);
                        }
                    }
                    if let (Some((ps, _)), Some((pb, _))) = (
                        self.project_lonlat(seam_b, seam_lat, rect),
                        self.project_lonlat(wb, b_lat, rect),
                    ) {
                        if ps.distance(pb) <= jump_limit {
                            painter.line_segment([ps, pb], track_stroke);
                        }
                    }
                } else if let (Some((pa, _)), Some((pb, _))) = (
                    self.project_lonlat(wa, a_lat, rect),
                    self.project_lonlat(wb, b_lat, rect),
                ) {
                    if pa.distance(pb) <= jump_limit {
                        painter.line_segment([pa, pb], track_stroke);
                    }
                }
            }
        }
        if let Some((mid, _)) = cur
            .waypoints
            .get(cur.waypoints.len() / 2)
            .and_then(|&(lon, lat)| self.project_lonlat(lon, lat, rect))
        {
            painter.text(
                mid + Vec2::new(0.0, -10.0),
                egui::Align2::CENTER_CENTER,
                cur.name,
                egui::FontId::proportional(11.0),
                Color32::from_rgb(254, 240, 138),
            );
        }
    }

    fn render_nino_boxes(&self, painter: &egui::Painter, rect: Rect, jump_limit: f32) {
        if !self.state.show_nino_boxes {
            return;
        }
        for nbox in NINO_BOXES {
            let c = Color32::from_rgb(nbox.color_rgb[0], nbox.color_rgb[1], nbox.color_rgb[2]);
            let stroke = Stroke::new(1.5_f32, c);
            let lon_e = if nbox.lon_east < nbox.lon_west {
                nbox.lon_east + 360.0
            } else {
                nbox.lon_east
            };

            let draw_edge = |lon_a: f64, lat_a: f64, lon_b: f64, lat_b: f64, steps: usize| {
                let mut prev: Option<(f64, Pos2)> = None;
                for s in 0..=steps {
                    let t = s as f64 / steps as f64;
                    let raw_lon = lon_a + (lon_b - lon_a) * t;
                    let wrapped_lon = (raw_lon + 180.0).rem_euclid(360.0) - 180.0;
                    let lat = lat_a + (lat_b - lat_a) * t;
                    if let Some((p, _)) = self.project_lonlat(wrapped_lon, lat, rect) {
                        if let Some((plon, pa)) = prev {
                            if (wrapped_lon - plon).abs() <= 180.0 && pa.distance(p) <= jump_limit {
                                painter.line_segment([pa, p], stroke);
                            }
                        }
                        prev = Some((wrapped_lon, p));
                    } else {
                        prev = None;
                    }
                }
            };

            draw_edge(nbox.lon_west, nbox.lat_max, lon_e, nbox.lat_max, 14);
            draw_edge(lon_e, nbox.lat_max, lon_e, nbox.lat_min, 4);
            draw_edge(lon_e, nbox.lat_min, nbox.lon_west, nbox.lat_min, 14);
            draw_edge(nbox.lon_west, nbox.lat_min, nbox.lon_west, nbox.lat_max, 4);

            let mid_lon = ((0.5 * (nbox.lon_west + lon_e)) + 180.0).rem_euclid(360.0) - 180.0;
            let mid_lat = 0.5 * (nbox.lat_min + nbox.lat_max);
            if let Some((p, _)) = self.project_lonlat(mid_lon, mid_lat, rect) {
                let short_lbl = match nbox.id {
                    "nino12" => "Nino 1+2",
                    "nino3" => "Nino 3",
                    "nino34" => "Nino 3.4",
                    "nino4" => "Nino 4",
                    _ => nbox.id,
                };
                painter.text(
                    p,
                    egui::Align2::CENTER_CENTER,
                    short_lbl,
                    egui::FontId::proportional(10.0),
                    c,
                );
            }
        }
    }

    fn render_geographic_labels(&self, painter: &egui::Painter, rect: Rect) {
        if !self.show_world_map {
            return;
        }
        const GEO_LABELS: &[(&str, f64, f64, bool)] = &[
            ("PACIFIC OCEAN", -148.0, 6.0, false),
            ("ATLANTIC OCEAN", -34.0, 12.0, false),
            ("INDIAN OCEAN", 78.0, -18.0, false),
            ("SOUTHERN OCEAN", 10.0, -60.0, false),
            ("N. AMERICA", -102.0, 44.0, true),
            ("S. AMERICA", -60.0, -14.0, true),
            ("AFRICA", 20.0, 6.0, true),
            ("EURASIA", 62.0, 50.0, true),
            ("AUSTRALIA", 134.0, -25.0, true),
            ("ANTARCTICA", 0.0, -78.0, true),
        ];

        for &(label, lon, lat, is_land) in GEO_LABELS {
            if let Some((pos, z3)) = self.project_lonlat(lon, lat, rect) {
                if self.view_mode == MapViewMode::Globe3D && self.effective_morph() < 0.5 && z3 < 0.22 {
                    continue;
                }
                let col = if is_land {
                    Color32::from_rgba_unmultiplied(226, 232, 214, 175)
                } else {
                    Color32::from_rgba_unmultiplied(125, 211, 252, 130)
                };
                painter.text(
                    pos,
                    egui::Align2::CENTER_CENTER,
                    label,
                    egui::FontId::proportional(if is_land { 10.0 } else { 10.5 }),
                    col,
                );
            }
        }
    }

    fn render_docked_spilhaus_inset(&self, painter: &egui::Painter, rect: Rect) {
        let is_mosaic = match self.state.spilhaus_inset_mode {
            SpilhausInsetMode::Hidden => return,
            SpilhausInsetMode::DockedSquare => false,
            SpilhausInsetMode::Mosaic3x3 => true,
        };
        let side = if is_mosaic { 222.0_f32 } else { 184.0_f32 };
        if rect.width() < side + 40.0 || rect.height() < side + 40.0 {
            return;
        }
        let inset = Rect::from_min_size(
            Pos2::new(rect.right() - side - 14.0, rect.bottom() - side - 14.0),
            Vec2::splat(side),
        );
        painter.rect_filled(
            inset,
            CornerRadius::same(6),
            Color32::from_rgba_unmultiplied(8, 19, 35, 238),
        );

        let pad = 6.0_f32;
        let inner = inset.shrink(pad);
        let (grid, n, coord_scale): (&[MapGridVertex], usize, f32) = if is_mosaic {
            (&self.mosaic_grid, MOSAIC_GRID_N, 3.0)
        } else {
            (&self.inset_grid, INSET_GRID_N, 1.0)
        };
        let stride = n + 1;
        let mut mesh = egui::epaint::Mesh::default();
        for v in grid {
            let px = inner.center().x + (v.su / coord_scale) * (inner.width() * 0.5);
            let py = inner.center().y - (v.sv / coord_scale) * (inner.height() * 0.5);
            let rgb = if v.is_land && self.show_world_map {
                land_color(v.lon_deg, v.lat_deg, v.is_ice)
            } else {
                self.ocean_color(v.lon_deg, v.lat_deg)
            };
            mesh.colored_vertex(Pos2::new(px, py), to_color32(rgb, 1.0));
        }
        for r in 0..n {
            for c in 0..n {
                let i00 = (r * stride + c) as u32;
                let i10 = i00 + 1;
                let i01 = ((r + 1) * stride + c) as u32;
                let i11 = i01 + 1;
                mesh.add_triangle(i00, i10, i11);
                mesh.add_triangle(i00, i11, i01);
            }
        }
        painter.add(egui::Shape::mesh(mesh));

        if self.show_world_map {
            let inset_shore = Stroke::new(
                if is_mosaic { 0.65_f32 } else { 0.8_f32 },
                Color32::from_rgba_unmultiplied(226, 232, 240, 215),
            );
            let tile_w = inner.width() / coord_scale;
            let inset_jump = tile_w * 0.25;
            for ring in &self.coastline_rings {
                for seg in ring.windows(2) {
                    if (seg[0].lon_deg - seg[1].lon_deg).abs() > 180.0 {
                        continue;
                    }
                    let ax = inner.center().x + (seg[0].su / coord_scale) * (inner.width() * 0.5);
                    let ay = inner.center().y - (seg[0].sv / coord_scale) * (inner.height() * 0.5);
                    let bx = inner.center().x + (seg[1].su / coord_scale) * (inner.width() * 0.5);
                    let by = inner.center().y - (seg[1].sv / coord_scale) * (inner.height() * 0.5);
                    let a = Pos2::new(ax, ay);
                    let b = Pos2::new(bx, by);
                    if a.distance(b) <= inset_jump {
                        painter.line_segment([a, b], inset_shore);
                    }
                }
            }
        }

        if self.state.show_currents {
            for &(lon, lat) in self.drifters.iter().step_by(4) {
                let (u, v) = lonlat_to_spilhaus_norm(lon, lat);
                let px = inner.center().x + (u as f32 / coord_scale) * (inner.width() * 0.5);
                let py = inner.center().y - (v as f32 / coord_scale) * (inner.height() * 0.5);
                painter.circle_filled(
                    Pos2::new(px, py),
                    if is_mosaic { 0.85 } else { 1.0 },
                    Color32::from_rgba_unmultiplied(125, 211, 252, 210),
                );
            }
        }

        if is_mosaic {
            let tile_stroke = Stroke::new(0.7_f32, Color32::from_rgba_unmultiplied(148, 210, 240, 95));
            for k in 1..3 {
                let f = k as f32 / 3.0;
                let gx = inner.left() + f * inner.width();
                let gy = inner.top() + f * inner.height();
                painter.line_segment(
                    [Pos2::new(gx, inner.top()), Pos2::new(gx, inner.bottom())],
                    tile_stroke,
                );
                painter.line_segment(
                    [Pos2::new(inner.left(), gy), Pos2::new(inner.right(), gy)],
                    tile_stroke,
                );
            }
            let primary_rect = Rect::from_center_size(inner.center(), inner.size() / 3.0);
            painter.rect_stroke(
                primary_rect,
                CornerRadius::ZERO,
                Stroke::new(1.1_f32, Color32::from_rgb(250, 204, 21)),
                egui::StrokeKind::Inside,
            );
            painter.rect_stroke(
                inner,
                CornerRadius::ZERO,
                Stroke::new(1.0_f32, Color32::from_rgba_unmultiplied(56, 189, 248, 175)),
                egui::StrokeKind::Inside,
            );
        } else {
            painter.rect_stroke(
                inner,
                CornerRadius::ZERO,
                Stroke::new(1.1_f32, Color32::from_rgb(250, 204, 21)),
                egui::StrokeKind::Inside,
            );
        }
        painter.rect_stroke(
            inset,
            CornerRadius::same(6),
            Stroke::new(1.2_f32, Color32::from_rgb(56, 189, 248)),
            egui::StrokeKind::Outside,
        );
        painter.text(
            Pos2::new(inset.center().x, inset.top() + 9.0),
            egui::Align2::CENTER_CENTER,
            if is_mosaic {
                "Spilhaus 3x3 Mosaic (ESRI:54099)"
            } else {
                "Spilhaus Square (ESRI:54099)"
            },
            egui::FontId::proportional(9.5),
            Color32::from_rgb(226, 232, 240),
        );
    }

    fn render_hud_badge_and_colorbar(&self, painter: &egui::Painter, rect: Rect) {
        if rect.width() < 320.0 || rect.height() < 180.0 {
            return;
        }
        let card = Rect::from_min_size(
            Pos2::new(rect.left() + 12.0, rect.top() + 12.0),
            Vec2::new(294.0, 74.0),
        );
        painter.rect_filled(
            card,
            CornerRadius::same(6),
            Color32::from_rgba_unmultiplied(8, 19, 35, 214),
        );
        painter.rect_stroke(
            card,
            CornerRadius::same(6),
            Stroke::new(1.0_f32, Color32::from_rgba_unmultiplied(56, 189, 248, 145)),
            egui::StrokeKind::Inside,
        );

        let phase = EnsoPhase::from_oni(self.state.oni_celsius);
        let head_col = if self.state.oni_celsius >= 0.5 {
            Color32::from_rgb(251, 146, 60)
        } else if self.state.oni_celsius <= -0.5 {
            Color32::from_rgb(56, 189, 248)
        } else {
            Color32::from_rgb(203, 213, 225)
        };
        painter.text(
            Pos2::new(card.left() + 10.0, card.top() + 14.0),
            egui::Align2::LEFT_CENTER,
            format!(
                "{} {} | ONI {:+.1}°C ({})",
                self.state.year,
                SEASONS[self.state.season_index.min(11)],
                self.state.oni_celsius,
                phase.label()
            ),
            egui::FontId::proportional(11.0),
            head_col,
        );

        let sub_text = if self.state.field_source == FieldSource::SatelliteRaster {
            match &self.active_raster {
                Some(grid) => {
                    let ts = if grid.timestamp_utc.is_empty() {
                        "loaded"
                    } else {
                        grid.timestamp_utc.as_str()
                    };
                    format!("NOAA OISST v2.1 ({ts})")
                }
                None => "Satellite Raster (not yet loaded)".to_string(),
            }
        } else {
            format!(
                "Model: {}",
                classify_historical_flavor(self.state.year, self.state.oni_celsius).label()
            )
        };
        painter.text(
            Pos2::new(card.left() + 10.0, card.top() + 30.0),
            egui::Align2::LEFT_CENTER,
            sub_text,
            egui::FontId::proportional(10.0),
            Color32::from_rgb(226, 232, 240),
        );

        let bar = Rect::from_min_size(
            Pos2::new(card.left() + 10.0, card.top() + 42.0),
            Vec2::new(card.width() - 20.0, 10.0),
        );
        let steps = 24_usize;
        let seg_w = bar.width() / steps as f32;
        for i in 0..steps {
            let t = (i as f32 + 0.5) / steps as f32;
            let rgb = match self.state.data_layer {
                DataLayer::EnsoSstAnomaly => diverged_anom_rgb(-3.5 + t * 7.3),
                DataLayer::CurrentVelocity => {
                    if t < 0.5 {
                        lerp_rgb([8.0, 24.0, 54.0], [24.0, 154.0, 188.0], t * 2.0)
                    } else {
                        lerp_rgb([24.0, 154.0, 188.0], [250.0, 204.0, 21.0], (t - 0.5) * 2.0)
                    }
                }
                DataLayer::AbsoluteSst => {
                    if t < 0.5 {
                        lerp_rgb([16.0, 32.0, 98.0], [20.0, 148.0, 162.0], t * 2.0)
                    } else {
                        lerp_rgb([20.0, 148.0, 162.0], [246.0, 108.0, 62.0], (t - 0.5) * 2.0)
                    }
                }
                DataLayer::ThermoclineDepth => {
                    if t < 0.5 {
                        lerp_rgb([28.0, 34.0, 108.0], [18.0, 158.0, 126.0], t * 2.0)
                    } else {
                        lerp_rgb([18.0, 158.0, 126.0], [244.0, 164.0, 44.0], (t - 0.5) * 2.0)
                    }
                }
            };
            let x0 = bar.left() + (i as f32) * seg_w;
            painter.rect_filled(
                Rect::from_min_size(Pos2::new(x0, bar.top()), Vec2::new(seg_w + 0.5, bar.height())),
                CornerRadius::ZERO,
                to_color32(rgb, 1.0),
            );
        }
        painter.rect_stroke(
            bar,
            CornerRadius::ZERO,
            Stroke::new(0.8_f32, Color32::from_rgb(100, 116, 139)),
            egui::StrokeKind::Inside,
        );

        let (lbl_min, lbl_mid, lbl_max) = match self.state.data_layer {
            DataLayer::EnsoSstAnomaly => ("-3.5°C", "0°C", "+3.8°C"),
            DataLayer::CurrentVelocity => ("0 m/s", "0.7 m/s", "1.4+ m/s"),
            DataLayer::AbsoluteSst => ("-2°C", "15°C", "32°C"),
            DataLayer::ThermoclineDepth => ("20 m", "105 m", "190 m"),
        };
        let tick_y = bar.bottom() + 9.0;
        let tick_col = Color32::from_rgb(186, 200, 218);
        painter.text(
            Pos2::new(bar.left(), tick_y),
            egui::Align2::LEFT_CENTER,
            lbl_min,
            egui::FontId::proportional(9.5),
            tick_col,
        );
        painter.text(
            Pos2::new(bar.center().x, tick_y),
            egui::Align2::CENTER_CENTER,
            lbl_mid,
            egui::FontId::proportional(9.5),
            tick_col,
        );
        painter.text(
            Pos2::new(bar.right(), tick_y),
            egui::Align2::RIGHT_CENTER,
            lbl_max,
            egui::FontId::proportional(9.5),
            tick_col,
        );
    }

    pub fn ui_step(&mut self, ctx: &egui::Context, interactive: bool) {
        if interactive {
            if let Some(rx) = &self.erddap_rx {
                if let Ok((cache_key, res)) = rx.try_recv() {
                    self.erddap_rx = None;
                    let key_to_store = self
                        .erddap_in_flight
                        .take()
                        .unwrap_or_else(|| cache_key.clone());
                    match res {
                        Ok(grid) => {
                            let cells = grid.valid_cell_count();
                            let has_vel = grid.has_velocity;
                            if let Some(store) = &self.store {
                                let _ = store.put_cached_raster(&key_to_store, &grid);
                            }
                            self.active_raster = Some(grid);
                            self.state.field_source = FieldSource::SatelliteRaster;
                            self.status_note = if has_vel {
                                format!("Loaded NOAA ERDDAP raster ({cache_key}, {cells} cells, OISST + NESDIS currents).")
                            } else {
                                format!("Loaded NOAA OISST v2.1 raster ({cache_key}, {cells} cells; analytical currents pre-2017).")
                            };
                        }
                        Err(err) => {
                            self.status_note = err;
                        }
                    }
                }
            }

            let wants_kb = ctx.wants_keyboard_input();
            let mut key_play_pause = false;
            let mut key_step: isize = 0;
            let mut key_view: Option<MapViewMode> = None;
            let mut key_reset = false;
            let mut toggle_palette = false;
            ctx.input_mut(|i| {
                let cmd_k = egui::KeyboardShortcut::new(egui::Modifiers::COMMAND, egui::Key::K);
                if i.consume_shortcut(&cmd_k) {
                    toggle_palette = true;
                }
                if i.key_pressed(egui::Key::F1) {
                    self.show_about_window = !self.show_about_window;
                }
                if i.key_pressed(egui::Key::Escape) {
                    if self.command_palette_open {
                        self.command_palette_open = false;
                        self.command_palette_just_opened = false;
                    } else if self.show_about_window {
                        self.show_about_window = false;
                    }
                }
                if !wants_kb && !toggle_palette {
                    if i.key_pressed(egui::Key::Space) {
                        key_play_pause = true;
                    }
                    if i.key_pressed(egui::Key::ArrowLeft) {
                        key_step = -1;
                    } else if i.key_pressed(egui::Key::ArrowRight) {
                        key_step = 1;
                    }
                    if i.key_pressed(egui::Key::Num1) {
                        key_view = Some(MapViewMode::Globe3D);
                    } else if i.key_pressed(egui::Key::Num2) {
                        key_view = Some(MapViewMode::WorldMap2D);
                    } else if i.key_pressed(egui::Key::Num3) {
                        key_view = Some(MapViewMode::SpilhausSquare);
                    }
                    if i.key_pressed(egui::Key::R) {
                        key_reset = true;
                    }
                }
            });
            if toggle_palette {
                self.command_palette_open = !self.command_palette_open;
                self.command_palette_just_opened = self.command_palette_open;
            }
            if key_play_pause {
                self.timeline_playing = !self.timeline_playing;
                self.timeline_accum_sec = 0.0;
            }
            if key_step != 0 {
                self.step_season(key_step);
            }
            if let Some(vm) = key_view {
                self.view_mode = vm;
                match vm {
                    MapViewMode::Globe3D => self.state.spilhaus_morph = 0.0,
                    MapViewMode::WorldMap2D => self.pan_2d = Vec2::ZERO,
                    MapViewMode::SpilhausSquare => self.state.spilhaus_morph = 1.0,
                }
            }
            if key_reset {
                self.reset_camera();
            }

            let dt = ctx.input(|i| i.stable_dt).min(0.05);
            if self.timeline_playing {
                self.timeline_accum_sec += dt;
                if self.timeline_accum_sec >= 0.25 {
                    self.timeline_accum_sec -= 0.25;
                    let max_s = latest_valid_season_for_year(self.state.year);
                    let (ny, ns) = if self.state.season_index < max_s {
                        (self.state.year, self.state.season_index + 1)
                    } else if self.state.year < END_YEAR {
                        (self.state.year + 1, 0)
                    } else {
                        (START_YEAR, 0)
                    };
                    self.apply_historical_date(ny, ns, false, None);
                }
            }
            if self.view_mode == MapViewMode::Globe3D
                && self.state.auto_rotate
                && self.state.spilhaus_morph < 0.5
            {
                self.state.yaw_rad += dt * 0.16;
            }
            self.step_drifters(f64::from(dt));
            ctx.request_repaint();
        }

        egui::TopBottomPanel::top("menu_bar").show(ctx, |ui| {
            ui.horizontal_wrapped(|ui| {
                ui.label(
                    RichText::new("PELAGIC GLOBE v1.1.0")
                        .strong()
                        .color(Color32::from_rgb(125, 211, 252)),
                );
                ui.separator();
                let phase = EnsoPhase::from_oni(self.state.oni_celsius);
                ui.label(format!(
                    "{} {} | ONI {:+.1}°C ({})",
                    self.state.year,
                    SEASONS[self.state.season_index.min(11)],
                    self.state.oni_celsius,
                    phase.label()
                ));
                ui.separator();

                if ui
                    .selectable_label(self.view_mode == MapViewMode::Globe3D && self.state.spilhaus_morph < 0.99, "3D Globe [1]")
                    .on_hover_text("Interactive rotating 3D ocean sphere (Hotkey: 1)")
                    .clicked()
                {
                    self.view_mode = MapViewMode::Globe3D;
                    self.state.spilhaus_morph = 0.0;
                }
                if ui
                    .selectable_label(self.view_mode == MapViewMode::WorldMap2D, "World Map (2D) [2]")
                    .on_hover_text("Equirectangular 2D global map (Hotkey: 2)")
                    .clicked()
                {
                    self.view_mode = MapViewMode::WorldMap2D;
                    self.pan_2d = Vec2::ZERO;
                }
                if ui
                    .selectable_label(self.view_mode == MapViewMode::SpilhausSquare, "Spilhaus Square [3]")
                    .on_hover_text("Conformal Spilhaus World Ocean Square ESRI:54099 (Hotkey: 3)")
                    .clicked()
                {
                    self.view_mode = MapViewMode::SpilhausSquare;
                    self.state.spilhaus_morph = 1.0;
                }
                if ui
                    .button("Reset View [R]")
                    .on_hover_text("Reset yaw, pitch, 2D pan, and zoom (Hotkey: R)")
                    .clicked()
                {
                    self.reset_camera();
                }

                ui.separator();
                if ui
                    .button("Save SQLite Bookmark")
                    .on_hover_text("Persist current camera, episode, and layer snapshot into WAL SQLite")
                    .clicked()
                {
                    if let Some(store) = &self.store {
                        let custom = self.bookmark_name_draft.trim();
                        let label = if custom.is_empty() {
                            format!(
                                "{} {} ({:+.1}C)",
                                self.state.year,
                                SEASONS[self.state.season_index.min(11)],
                                self.state.oni_celsius
                            )
                        } else {
                            custom.to_string()
                        };
                        if store.save_bookmark(&label, &self.state).is_ok() {
                            self.saved_bookmarks = store.list_bookmarks().unwrap_or_default();
                            self.bookmark_name_draft.clear();
                            self.status_note =
                                format!("Saved bookmark '{label}' to WAL SQLite store.");
                        }
                    }
                }
                ui.separator();
                if ui.button("About [F1]").clicked() {
                    self.show_about_window = !self.show_about_window;
                }
                if ui
                    .selectable_label(self.command_palette_open, "Palette [Ctrl+K]")
                    .on_hover_text("Fuzzy quick-jump palette for Globe/Spilhaus morph, ENSO flavors, thermocline layers, and exports")
                    .clicked()
                {
                    self.command_palette_open = !self.command_palette_open;
                    self.command_palette_just_opened = self.command_palette_open;
                }
            });
        });

        egui::SidePanel::left("controls_panel")
            .default_width(316.0)
            .show(ctx, |ui| {
                egui::Frame::NONE
                    .fill(Color32::from_rgb(12, 22, 36))
                    .corner_radius(CornerRadius::same(8))
                    .inner_margin(10.0)
                    .show(ui, |ui| {
                        egui::ScrollArea::vertical().show(ui, |ui| {
                            ui.heading("Spilhaus & Globe Layer");
                            if ui
                                .add(
                                    egui::Slider::new(&mut self.state.spilhaus_morph, 0.0..=1.0)
                                        .text("3D <-> Spilhaus Morph"),
                                )
                                .on_hover_text("Continuously morph between the 3D sphere (0.0) and conformal Spilhaus World Ocean Square (1.0)")
                                .changed()
                            {
                                self.view_mode = if self.state.spilhaus_morph >= 0.99 {
                                    MapViewMode::SpilhausSquare
                                } else {
                                    MapViewMode::Globe3D
                                };
                            }
                            ui.checkbox(&mut self.show_world_map, "World Map & Coastlines")
                                .on_hover_text("Toggle Natural-Earth continental coastlines, ice shelves, and ocean basin labels");
                            ui.checkbox(&mut self.show_graticule, "Lat/Lon Graticule")
                                .on_hover_text("Show 30-degree latitude/longitude graticule and highlighted Equator");
                            ui.checkbox(
                                &mut self.state.show_spilhaus_overlay,
                                "Spilhaus Seam & Conformal Grid",
                            )
                            .on_hover_text("Highlight the East Asia - Bering Strait - South America Spilhaus boundary cut seam");
                            ui.checkbox(&mut self.state.auto_rotate, "Auto-Rotate 3D Globe")
                                .on_hover_text("Slowly spin the 3D globe in azimuth (disabled automatically when dragging)");
                            ui.checkbox(&mut self.state.show_currents, "Lagrangian Current Ribbons")
                                .on_hover_text("Animate 4,200 RK2 spherical surface drifters across global ocean currents");
                            ui.checkbox(
                                &mut self.state.show_nino_boxes,
                                "NOAA Niño 1+2 / 3 / 3.4 / 4 Boxes",
                            )
                            .on_hover_text("Overlay official NOAA equatorial Pacific monitoring regions");
                            ui.add(
                                egui::Slider::new(&mut self.state.flow_speed, 0.25..=2.5)
                                    .text("Advection Speed"),
                            )
                            .on_hover_text("Multiplier for Lagrangian surface drifter RK2 integration step");

                            ui.add_space(4.0);
                            ui.label(RichText::new("Spilhaus Inset Overlay").strong());
                            ui.horizontal_wrapped(|ui| {
                                ui.selectable_value(
                                    &mut self.state.spilhaus_inset_mode,
                                    SpilhausInsetMode::DockedSquare,
                                    "Docked Square",
                                )
                                .on_hover_text("Render a single synchronized ESRI:54099 Spilhaus inset in the bottom-right corner");
                                ui.selectable_value(
                                    &mut self.state.spilhaus_inset_mode,
                                    SpilhausInsetMode::Mosaic3x3,
                                    "3x3 Mosaic",
                                )
                                .on_hover_text("Render a periodic 3x3 conformal Adams/Spilhaus world-ocean mosaic inset");
                                ui.selectable_value(
                                    &mut self.state.spilhaus_inset_mode,
                                    SpilhausInsetMode::Hidden,
                                    "Hidden",
                                )
                                .on_hover_text("Hide the bottom-right Spilhaus inset panel");
                            });

                            ui.separator();
                            ui.heading("Field Source (Model vs. Satellite Raster)");
                            ui.selectable_value(
                                &mut self.state.field_source,
                                FieldSource::AnalyticalModel,
                                FieldSource::AnalyticalModel.label(),
                            )
                            .on_hover_text("Use zero-latency 1950-2026 NOAA CPC ERSSTv6 coupled ENSO + gyre model");
                            if ui
                                .selectable_label(
                                    self.state.field_source == FieldSource::SatelliteRaster,
                                    FieldSource::SatelliteRaster.label(),
                                )
                                .on_hover_text("Display fetched NOAA CoastWatch ERDDAP OISST v2.1 + NESDIS geostrophic raster")
                                .clicked()
                            {
                                self.state.field_source = FieldSource::SatelliteRaster;
                                if self.active_raster.is_none() {
                                    self.status_note = "No satellite raster loaded yet — click 'Fetch NOAA Satellite Raster' below.".to_string();
                                }
                            }
                            ui.horizontal(|ui| {
                                ui.label("ERDDAP Stride:")
                                    .on_hover_text("Subsampling stride on 0.25 deg NOAA grids (2 = 0.5 deg, 4 = 1.0 deg, 8 = 2.0 deg)");
                                egui::ComboBox::from_id_salt("erddap_stride_combo")
                                    .selected_text(format!(
                                        "{} ({:.2}°)",
                                        self.state.erddap_stride,
                                        self.state.erddap_stride as f32 * 0.25
                                    ))
                                    .show_ui(ui, |ui| {
                                        for (s, desc) in [
                                            (2_usize, "2 (0.50° — High Detail)"),
                                            (4, "4 (1.00° — Standard)"),
                                            (8, "8 (2.00° — Fast)"),
                                            (12, "12 (3.00° — Coarse)"),
                                            (16, "16 (4.00° — Minimal)"),
                                        ] {
                                            ui.selectable_value(&mut self.state.erddap_stride, s, desc);
                                        }
                                    });
                            });
                            let fetch_month = (self.state.season_index + 1).clamp(1, 12);
                            let fetch_label = match &self.erddap_in_flight {
                                Some(iso) => format!("Fetching NOAA ERDDAP ({iso})..."),
                                None => format!(
                                    "Fetch NOAA Satellite Raster ({:04}-{:02})",
                                    self.state.year, fetch_month
                                ),
                            };
                            if ui
                                .add_enabled(
                                    self.erddap_in_flight.is_none(),
                                    egui::Button::new(fetch_label),
                                )
                                .on_hover_text("Query live NOAA OISST v2.1 (Sep 1981+) and NESDIS currents (2017+) over HTTPS with SQLite caching")
                                .clicked()
                            {
                                self.trigger_erddap_fetch();
                            }

                            ui.separator();
                            ui.heading("Data Layer");
                            for layer in [
                                DataLayer::EnsoSstAnomaly,
                                DataLayer::CurrentVelocity,
                                DataLayer::AbsoluteSst,
                                DataLayer::ThermoclineDepth,
                            ] {
                                let tip = match layer {
                                    DataLayer::EnsoSstAnomaly => "Diverging sea-surface temperature anomaly relative to climatology",
                                    DataLayer::CurrentVelocity => "Mixed-layer surface current speed magnitude (m/s), excluding subsurface EUC",
                                    DataLayer::AbsoluteSst => "Seasonal climatology + historical ENSO/marine-heatwave SST (deg C)",
                                    DataLayer::ThermoclineDepth => "20 deg C isotherm depth across tropical Pacific (|lat| <= 20 deg, 120E-80W); masked elsewhere",
                                };
                                ui.selectable_value(
                                    &mut self.state.data_layer,
                                    layer,
                                    layer.label(),
                                )
                                .on_hover_text(tip);
                            }

                            ui.separator();
                            ui.heading("NOAA CPC El Niño / La Niña Presets");
                            for preset in ENSO_PRESETS {
                                let is_active = self.state.field_source
                                    == FieldSource::AnalyticalModel
                                    && self.state.year == preset.year
                                    && self.state.season_index == preset.season_index
                                    && (self.state.oni_celsius - preset.oni).abs() < 0.05;
                                let badge_col = if preset.oni >= 0.5 {
                                    Color32::from_rgb(251, 146, 60)
                                } else if preset.oni <= -0.5 {
                                    Color32::from_rgb(56, 189, 248)
                                } else {
                                    Color32::from_rgb(148, 163, 184)
                                };
                                ui.horizontal(|ui| {
                                    ui.colored_label(badge_col, format!("{:+.1}", preset.oni));
                                    if ui
                                        .selectable_label(is_active, preset.title)
                                        .on_hover_text(preset.summary)
                                        .clicked()
                                    {
                                        self.timeline_playing = false;
                                        self.apply_historical_date(
                                            preset.year,
                                            preset.season_index,
                                            true,
                                            Some(preset.summary.to_string()),
                                        );
                                    }
                                });
                            }

                            ui.separator();
                            ui.collapsing(
                                format!("Saved SQLite Bookmarks ({})", self.saved_bookmarks.len()),
                                |ui| {
                                    ui.horizontal(|ui| {
                                        ui.add(
                                            egui::TextEdit::singleline(&mut self.bookmark_name_draft)
                                                .desired_width(170.0)
                                                .hint_text("Optional bookmark label"),
                                        )
                                        .on_hover_text("Enter an optional name before saving a bookmark");
                                        if ui
                                            .button("Save")
                                            .on_hover_text("Save current simulation view to local WAL SQLite store")
                                            .clicked()
                                        {
                                            if let Some(store) = &self.store {
                                                let custom = self.bookmark_name_draft.trim();
                                                let label = if custom.is_empty() {
                                                    format!(
                                                        "{} {} ({:+.1}C)",
                                                        self.state.year,
                                                        SEASONS[self.state.season_index.min(11)],
                                                        self.state.oni_celsius
                                                    )
                                                } else {
                                                    custom.to_string()
                                                };
                                                if store.save_bookmark(&label, &self.state).is_ok() {
                                                    self.saved_bookmarks =
                                                        store.list_bookmarks().unwrap_or_default();
                                                    self.bookmark_name_draft.clear();
                                                    self.status_note =
                                                        format!("Saved bookmark '{label}' to WAL SQLite store.");
                                                }
                                            }
                                        }
                                    });
                                    if self.saved_bookmarks.is_empty() {
                                        ui.small("No saved bookmarks yet — click Save above.");
                                    } else {
                                        let mut restored: Option<(String, SimulationSnapshot)> = None;
                                        let mut deleted_id: Option<(i64, String)> = None;
                                        for (id, name, snap) in &self.saved_bookmarks {
                                            ui.horizontal(|ui| {
                                                if ui
                                                    .button(format!("#{id} {name}"))
                                                    .on_hover_text("Restore saved camera, layer, and historical ENSO date")
                                                    .clicked()
                                                {
                                                    restored = Some((name.clone(), snap.clone()));
                                                }
                                                if ui
                                                    .small_button("Delete")
                                                    .on_hover_text("Remove this bookmark from SQLite")
                                                    .clicked()
                                                {
                                                    deleted_id = Some((*id, name.clone()));
                                                }
                                            });
                                        }
                                        if let Some((del_id, del_name)) = deleted_id {
                                            if let Some(store) = &self.store {
                                                if store.delete_bookmark(del_id).is_ok() {
                                                    self.saved_bookmarks =
                                                        store.list_bookmarks().unwrap_or_default();
                                                    self.status_note =
                                                        format!("Deleted bookmark #{del_id} '{del_name}'.");
                                                }
                                            }
                                        } else if let Some((name, snap)) = restored {
                                            self.state = snap;
                                            self.state.field_source = FieldSource::AnalyticalModel;
                                            self.timeline_playing = false;
                                            self.view_mode = if self.state.spilhaus_morph >= 0.99 {
                                                MapViewMode::SpilhausSquare
                                            } else {
                                                MapViewMode::Globe3D
                                            };
                                            self.status_note = format!("Restored bookmark '{name}'.");
                                        }
                                    }
                                },
                            );

                            ui.separator();
                            ui.heading("Export CSV & Spilhaus SVG");
                            ui.horizontal(|ui| {
                                ui.label("Prefix:");
                                ui.add(
                                    egui::TextEdit::singleline(&mut self.export_prefix)
                                        .desired_width(150.0),
                                )
                                .on_hover_text("Filename prefix written to working directory");
                            });
                            ui.horizontal(|ui| {
                                if ui
                                    .button("Export ONI CSV")
                                    .on_hover_text("Export complete 1950-2026 NOAA CPC ERSSTv6 ONI table to CSV")
                                    .clicked()
                                {
                                    let prefix = if self.export_prefix.trim().is_empty() {
                                        "pelagic_globe"
                                    } else {
                                        self.export_prefix.trim()
                                    };
                                    let path = format!("{prefix}_oni_1950_2026.csv");
                                    match std::fs::write(&path, crate::enso::export_oni_csv()) {
                                        Ok(()) => {
                                            self.status_note = format!("Exported NOAA ONI CSV -> {path}");
                                        }
                                        Err(e) => {
                                            self.status_note = format!("CSV export failed: {e}");
                                        }
                                    }
                                }
                                if ui
                                    .button("Export Spilhaus SVG")
                                    .on_hover_text("Export Jacobian-rotated Spilhaus World Ocean SVG for active episode")
                                    .clicked()
                                {
                                    let prefix = if self.export_prefix.trim().is_empty() {
                                        "pelagic_globe"
                                    } else {
                                        self.export_prefix.trim()
                                    };
                                    let season = SEASONS[self.state.season_index.min(11)];
                                    let path = format!(
                                        "{prefix}_spilhaus_{}_{}.svg",
                                        self.state.year, season
                                    );
                                    let svg = export_spilhaus_svg(
                                        self.state.year,
                                        self.state.season_index,
                                        self.state.oni_celsius,
                                    );
                                    match std::fs::write(&path, svg) {
                                        Ok(()) => {
                                            self.status_note = format!("Exported Spilhaus SVG -> {path}");
                                        }
                                        Err(e) => {
                                            self.status_note = format!("SVG export failed: {e}");
                                        }
                                    }
                                }
                            });
                            ui.horizontal(|ui| {
                                if ui
                                    .button("Export Digest CSV")
                                    .on_hover_text("Export Niño-3.4 & 20°C thermocline comparative digest across all historical ENSO flavors")
                                    .clicked()
                                {
                                    self.export_enso_digest_file(true);
                                }
                                if ui
                                    .button("Export Digest JSON")
                                    .on_hover_text("Export Niño-3.4 & 20°C thermocline comparative digest as structured JSON")
                                    .clicked()
                                {
                                    self.export_enso_digest_file(false);
                                }
                            });
                            ui.collapsing("ENSO Flavor & Thermocline Digest", |ui| {
                                let digest = build_enso_comparative_digest();
                                egui::Grid::new("pelagic_enso_digest_grid")
                                    .striped(true)
                                    .spacing([6.0, 2.0])
                                    .show(ui, |ui| {
                                        ui.small(RichText::new("Episode").strong());
                                        ui.small(RichText::new("ONI").strong());
                                        ui.small(RichText::new("N3.4").strong());
                                        ui.small(RichText::new("E.20°C").strong());
                                        ui.end_row();
                                        for row in digest.rows.iter().take(6) {
                                            ui.small(&row.slug);
                                            ui.small(format!("{:+.1}°", row.oni_celsius));
                                            ui.small(format!("{:+.2}°", row.nino34_sst_anom_c));
                                            ui.small(format!("{:.0}m", row.east_pac_thermocline_m));
                                            ui.end_row();
                                        }
                                    });
                            });

                            ui.separator();
                            ui.heading("Major Ocean Currents (14)");
                            for cur in MAJOR_CURRENTS {
                                let selected =
                                    self.state.selected_current_id.as_deref() == Some(cur.id);
                                if ui
                                    .selectable_label(
                                        selected,
                                        format!("{} ({})", cur.name, cur.basin),
                                    )
                                    .on_hover_text(format!("{} | {} — click to inspect & focus camera", cur.transport_sv, cur.enso_response))
                                    .clicked()
                                {
                                    self.state.selected_current_id = Some(cur.id.to_string());
                                    self.state.yaw_rad = -(cur.focus_lon_deg as f32).to_radians();
                                    self.state.pitch_rad =
                                        (cur.focus_lat_deg as f32).to_radians() * 0.55;
                                }
                            }

                            if let Some(cur_id) = self.state.selected_current_id.as_deref() {
                                if let Some(cur) = MAJOR_CURRENTS.iter().find(|c| c.id == cur_id) {
                                    ui.add_space(6.0);
                                    egui::Frame::NONE
                                        .fill(Color32::from_rgb(15, 32, 54))
                                        .stroke(Stroke::new(1.0_f32, Color32::from_rgb(56, 189, 248)))
                                        .corner_radius(CornerRadius::same(6))
                                        .inner_margin(8.0)
                                        .show(ui, |ui| {
                                            ui.label(
                                                RichText::new(cur.name)
                                                    .strong()
                                                    .color(Color32::from_rgb(254, 240, 138)),
                                            );
                                            ui.small(format!("Basin: {} | Regime: {:?}", cur.basin, cur.regime));
                                            ui.small(format!("Transport: {}", cur.transport_sv));
                                            ui.small(format!("Peak Speed: {:.2} m/s", cur.peak_speed_ms));
                                            if cur.regime == crate::currents::ThermalRegime::SubsurfaceEquatorial {
                                                ui.small(
                                                    RichText::new("Subsurface core (91-200 m): excluded from surface drifter velocity.")
                                                        .color(Color32::from_rgb(125, 211, 252)),
                                                );
                                            }
                                            ui.small(format!("ENSO Response: {}", cur.enso_response));
                                            ui.hyperlink_to(cur.source_title, cur.source_url)
                                                .on_hover_text(cur.source_url);
                                        });
                                }
                            }
                        });
                    });
            });

        egui::TopBottomPanel::bottom("thermocline_footer")
            .min_height(150.0)
            .show(ctx, |ui| {
                ui.columns(2, |cols| {
                    let l_ui = &mut cols[0];
                    l_ui.horizontal_wrapped(|ui| {
                        ui.label(
                            RichText::new("Historical Scrubber (1950–2026)")
                                .strong()
                                .color(Color32::from_rgb(186, 230, 253)),
                        );
                        if ui
                            .button(if self.timeline_playing { "Pause [Space]" } else { "Play 1950->2026 [Space]" })
                            .on_hover_text("Play or pause chronological 1950-2026 ENSO playback (Hotkey: Space)")
                            .clicked()
                        {
                            self.timeline_playing = !self.timeline_playing;
                            self.timeline_accum_sec = 0.0;
                        }
                        if ui
                            .button("< Prev")
                            .on_hover_text("Step back one 3-month ONI season (Hotkey: Left Arrow)")
                            .clicked()
                        {
                            self.step_season(-1);
                        }
                        if ui
                            .button("Next >")
                            .on_hover_text("Step forward one 3-month ONI season (Hotkey: Right Arrow)")
                            .clicked()
                        {
                            self.step_season(1);
                        }
                    });

                    // Interactive 1950–2026 ONI Sparkline Chart
                    let total_slots = ((END_YEAR - START_YEAR) * 12 + 8) as usize;
                    let chart_w = l_ui.available_width().max(220.0);
                    let (chart_rect, chart_resp) =
                        l_ui.allocate_exact_size(Vec2::new(chart_w, 44.0), Sense::click_and_drag());
                    if chart_resp.clicked() || chart_resp.dragged() {
                        if let Some(pos) = chart_resp.interact_pointer_pos() {
                            self.timeline_playing = false;
                            let frac = ((pos.x - chart_rect.left()) / chart_rect.width().max(1.0))
                                .clamp(0.0, 1.0);
                            let slot = (frac * total_slots as f32).round() as usize;
                            let y = START_YEAR + (slot / 12) as i32;
                            let s = slot % 12;
                            self.apply_historical_date(y, s, false, None);
                        }
                    }

                    let cp = l_ui.painter_at(chart_rect);
                    cp.rect_filled(chart_rect, CornerRadius::same(4), Color32::from_rgb(8, 17, 32));
                    let zero_y = chart_rect.center().y;
                    let scale_y = (chart_rect.height() * 0.44) / 2.7_f32;
                    cp.line_segment(
                        [Pos2::new(chart_rect.left(), zero_y), Pos2::new(chart_rect.right(), zero_y)],
                        Stroke::new(0.8_f32, Color32::from_rgba_unmultiplied(148, 163, 184, 90)),
                    );
                    cp.line_segment(
                        [
                            Pos2::new(chart_rect.left(), zero_y - 0.5 * scale_y),
                            Pos2::new(chart_rect.right(), zero_y - 0.5 * scale_y),
                        ],
                        Stroke::new(0.6_f32, Color32::from_rgba_unmultiplied(251, 146, 60, 75)),
                    );
                    cp.line_segment(
                        [
                            Pos2::new(chart_rect.left(), zero_y + 0.5 * scale_y),
                            Pos2::new(chart_rect.right(), zero_y + 0.5 * scale_y),
                        ],
                        Stroke::new(0.6_f32, Color32::from_rgba_unmultiplied(56, 189, 248, 75)),
                    );

                    for slot in (0..=total_slots).step_by(2) {
                        let y = START_YEAR + (slot / 12) as i32;
                        let s = slot % 12;
                        if let Some(oni) = lookup_oni(y, s) {
                            let x = chart_rect.left()
                                + (slot as f32 / total_slots.max(1) as f32) * chart_rect.width();
                            let bar_y = zero_y - (oni as f32) * scale_y;
                            let col = if oni >= 1.5 {
                                Color32::from_rgb(239, 68, 68)
                            } else if oni >= 0.5 {
                                Color32::from_rgb(251, 146, 60)
                            } else if oni <= -1.5 {
                                Color32::from_rgb(37, 99, 235)
                            } else if oni <= -0.5 {
                                Color32::from_rgb(56, 189, 248)
                            } else {
                                Color32::from_rgb(100, 116, 139)
                            };
                            cp.line_segment(
                                [Pos2::new(x, zero_y), Pos2::new(x, bar_y)],
                                Stroke::new(1.3_f32, col),
                            );
                        }
                    }

                    let active_slot = (((self.state.year - START_YEAR) * 12)
                        + self.state.season_index as i32)
                        .clamp(0, total_slots as i32) as f32;
                    let cur_x = chart_rect.left()
                        + (active_slot / total_slots.max(1) as f32) * chart_rect.width();
                    cp.line_segment(
                        [Pos2::new(cur_x, chart_rect.top()), Pos2::new(cur_x, chart_rect.bottom())],
                        Stroke::new(1.8_f32, Color32::from_rgb(250, 204, 21)),
                    );

                    l_ui.horizontal(|ui| {
                        let mut yr = self.state.year;
                        let max_s = latest_valid_season_for_year(yr);
                        let mut sidx = self.state.season_index.min(max_s);
                        let changed_y = ui
                            .add(egui::Slider::new(&mut yr, START_YEAR..=END_YEAR).text("Year"))
                            .changed();
                        let changed_s = ui
                            .add(
                                egui::Slider::new(&mut sidx, 0..=max_s)
                                    .text("Season")
                                    .custom_formatter(|v, _| {
                                        SEASONS[(v.round() as usize).min(11)].to_string()
                                    }),
                            )
                            .changed();
                        if changed_y || changed_s {
                            self.timeline_playing = false;
                            self.apply_historical_date(yr, sidx, false, None);
                        }
                    });

                    // Right column: Visual 20°C Thermocline Cross-Section (120°E -> 80°W) & Metrics
                    let r_ui = &mut cols[1];
                    let d_west = equatorial_20c_isotherm_depth_historical_m(
                        140.0,
                        self.state.year,
                        self.state.oni_celsius,
                    );
                    let d_east = equatorial_20c_isotherm_depth_historical_m(
                        275.0,
                        self.state.year,
                        self.state.oni_celsius,
                    );
                    r_ui.label(
                        RichText::new("Equatorial Pacific 20°C Thermocline Isotherm (120°E -> 80°W)")
                            .strong()
                            .color(Color32::from_rgb(186, 230, 253)),
                    );

                    let th_w = r_ui.available_width().max(200.0);
                    let (th_rect, _) =
                        r_ui.allocate_exact_size(Vec2::new(th_w, 56.0), Sense::hover());
                    let tp = r_ui.painter_at(th_rect);
                    tp.rect_filled(th_rect, CornerRadius::same(4), Color32::from_rgb(8, 22, 46));
                    let pts_n = 32_usize;
                    let mut iso_pts = Vec::with_capacity(pts_n + 1);
                    for i in 0..=pts_n {
                        let t = i as f64 / pts_n as f64;
                        let lon360 = 120.0 + t * 160.0;
                        let depth_m = equatorial_20c_isotherm_depth_historical_m(
                            lon360,
                            self.state.year,
                            self.state.oni_celsius,
                        );
                        let px = th_rect.left() + (t as f32) * th_rect.width();
                        let py = th_rect.top()
                            + ((depth_m as f32 / 240.0).clamp(0.05, 0.96)) * th_rect.height();
                        iso_pts.push(Pos2::new(px, py));
                    }
                    let mut warm_mesh = egui::epaint::Mesh::default();
                    let warm_col = if self.state.oni_celsius >= 0.5 {
                        Color32::from_rgba_unmultiplied(249, 115, 22, 130)
                    } else if self.state.oni_celsius <= -0.5 {
                        Color32::from_rgba_unmultiplied(56, 189, 248, 120)
                    } else {
                        Color32::from_rgba_unmultiplied(20, 184, 166, 115)
                    };
                    for &pt in &iso_pts {
                        warm_mesh.colored_vertex(Pos2::new(pt.x, th_rect.top()), warm_col);
                        warm_mesh.colored_vertex(pt, warm_col);
                    }
                    for i in 0..pts_n {
                        let a = (i * 2) as u32;
                        warm_mesh.add_triangle(a, a + 2, a + 3);
                        warm_mesh.add_triangle(a, a + 3, a + 1);
                    }
                    tp.add(egui::Shape::mesh(warm_mesh));
                    for pair in iso_pts.windows(2) {
                        tp.line_segment(
                            [pair[0], pair[1]],
                            Stroke::new(2.2_f32, Color32::from_rgb(250, 204, 21)),
                        );
                    }
                    let west_x = th_rect.left() + ((140.0 - 120.0) / 160.0) * th_rect.width();
                    let east_x = th_rect.left() + ((275.0 - 120.0) / 160.0) * th_rect.width();
                    tp.line_segment(
                        [Pos2::new(west_x, th_rect.top()), Pos2::new(west_x, th_rect.bottom())],
                        Stroke::new(0.8_f32, Color32::from_rgba_unmultiplied(226, 232, 240, 90)),
                    );
                    tp.line_segment(
                        [Pos2::new(east_x, th_rect.top()), Pos2::new(east_x, th_rect.bottom())],
                        Stroke::new(0.8_f32, Color32::from_rgba_unmultiplied(226, 232, 240, 90)),
                    );
                    tp.rect_stroke(
                        th_rect,
                        CornerRadius::same(4),
                        Stroke::new(0.9_f32, Color32::from_rgb(51, 85, 120)),
                        egui::StrokeKind::Inside,
                    );

                    r_ui.label(format!(
                        "West (140°E): {:.0} m  |  East (85°W): {:.0} m  |  Tilt ΔZ: {:+.0} m",
                        d_west,
                        d_east,
                        d_west - d_east
                    ));
                });

                ui.label(
                    RichText::new(&self.status_note).color(Color32::from_rgb(167, 243, 208)),
                );
                ui.separator();
                ui.horizontal(|ui| {
                    ui.small("Copyright CRV Science Network  |  GPL-3.0-or-later  |  PelagicGlobe v1.1.0");
                    ui.separator();
                    if ui.small_button("About PelagicGlobe").clicked() {
                        self.show_about_window = true;
                    }
                });
            });

        egui::CentralPanel::default().show(ctx, |ui| {
            let (rect, response) = ui.allocate_exact_size(ui.available_size(), Sense::drag());
            if response.hovered() {
                let scroll = ui.input(|i| i.smooth_scroll_delta.y);
                if scroll.abs() > 0.1 {
                    self.state.zoom = (self.state.zoom * (1.0 + scroll * 0.0015)).clamp(0.55, 2.8);
                }
            }
            if response.double_clicked() {
                self.reset_camera();
            }
            if response.dragged() {
                let delta = response.drag_delta();
                if self.view_mode == MapViewMode::WorldMap2D {
                    let lim_x = rect.width() * 0.65;
                    let lim_y = rect.height() * 0.55;
                    self.pan_2d.x = (self.pan_2d.x + delta.x).clamp(-lim_x, lim_x);
                    self.pan_2d.y = (self.pan_2d.y + delta.y).clamp(-lim_y, lim_y);
                } else {
                    self.state.yaw_rad += delta.x * 0.008;
                    self.state.pitch_rad =
                        (self.state.pitch_rad + delta.y * 0.008).clamp(-1.35, 1.35);
                    self.state.auto_rotate = false;
                }
            }

            let painter = ui.painter_at(rect);
            painter.rect_filled(rect, CornerRadius::ZERO, Color32::from_rgb(5, 11, 22));

            let morph = self.effective_morph();
            let center = rect.center();
            let radius = rect.width().min(rect.height()) * 0.42 * self.state.zoom;
            let jump_limit = if self.view_mode == MapViewMode::WorldMap2D {
                self.world_map_rect(rect).width() * 0.40
            } else {
                radius * 0.28
            };

            if self.view_mode == MapViewMode::WorldMap2D {
                self.render_surface_mesh_2d(&painter, rect);
            } else {
                if morph < 0.95 {
                    painter.circle_filled(
                        center,
                        radius * (1.045 - 0.12 * morph),
                        Color32::from_rgba_unmultiplied(36, 112, 198, ((1.0 - morph) * 42.0) as u8),
                    );
                    painter.circle_filled(
                        center,
                        radius * (1.02 - 0.12 * morph),
                        Color32::from_rgba_unmultiplied(16, 48, 86, 105),
                    );
                }
                self.render_surface_mesh_3d_or_spilhaus(&painter, rect);
            }

            self.render_graticule(&painter, rect, jump_limit);
            self.render_coastlines(&painter, rect, jump_limit);

            if self.view_mode != MapViewMode::WorldMap2D && self.state.show_spilhaus_overlay {
                for pair in self.seam_cache.windows(2) {
                    if (pair[0].0 - pair[1].0).abs() > 180.0 {
                        continue;
                    }
                    if let (Some((a, _)), Some((b, _))) = (
                        self.project_lonlat(pair[0].0, pair[0].1, rect),
                        self.project_lonlat(pair[1].0, pair[1].1, rect),
                    ) {
                        if a.distance(b) <= jump_limit {
                            painter.line_segment(
                                [a, b],
                                Stroke::new(1.5_f32, Color32::from_rgb(250, 204, 21)),
                            );
                        }
                    }
                }
            }

            if self.state.show_currents {
                for &(lon, lat) in &self.drifters {
                    if let Some((p, _)) = self.project_lonlat(lon, lat, rect) {
                        let (_sst, anom) = self.sample_effective_sst_anom(lon, lat);
                        let col = if anom > 0.35 {
                            Color32::from_rgba_unmultiplied(251, 146, 60, 210)
                        } else if anom < -0.35 {
                            Color32::from_rgba_unmultiplied(56, 189, 248, 210)
                        } else {
                            Color32::from_rgba_unmultiplied(167, 243, 208, 175)
                        };
                        let (vx, vy) = self.sample_effective_velocity(lon, lat);
                        let inv_cos = 1.0 / lat.to_radians().cos().max(0.15);
                        let next_lon = (lon + vx * 2.2 * inv_cos + 180.0).rem_euclid(360.0) - 180.0;
                        let next_lat = (lat + vy * 2.2).clamp(-84.0, 84.0);
                        if let Some((p_tip, _)) = self.project_lonlat(next_lon, next_lat, rect) {
                            let d = p.distance(p_tip);
                            if d > 0.6 && d <= 14.0 {
                                painter.line_segment([p, p_tip], Stroke::new(1.25_f32, col));
                                continue;
                            }
                        }
                        painter.circle_filled(p, 1.3, col);
                    }
                }
            }

            self.render_selected_current_track(&painter, rect, jump_limit);
            self.render_nino_boxes(&painter, rect, jump_limit);
            self.render_geographic_labels(&painter, rect);
            self.render_docked_spilhaus_inset(&painter, rect);
            self.render_hud_badge_and_colorbar(&painter, rect);

            if let Some(hover_pos) = response.hover_pos() {
                if let Some((hlon, hlat)) = self.screen_pos_to_lonlat(hover_pos, rect) {
                    let telem = self.sample_ocean_telemetry(hlon, hlat);
                    painter.circle_stroke(
                        hover_pos,
                        6.0,
                        Stroke::new(1.2_f32, Color32::from_rgb(250, 204, 21)),
                    );
                    let card = Rect::from_min_size(
                        Pos2::new(rect.left() + 12.0, rect.top() + 92.0),
                        Vec2::new(294.0, 74.0),
                    );
                    if rect.height() > 220.0 {
                        painter.rect_filled(
                            card,
                            CornerRadius::same(6),
                            Color32::from_rgba_unmultiplied(8, 19, 35, 222),
                        );
                        painter.rect_stroke(
                            card,
                            CornerRadius::same(6),
                            Stroke::new(1.0_f32, Color32::from_rgba_unmultiplied(250, 204, 21, 145)),
                            egui::StrokeKind::Inside,
                        );
                        let lon_dir = if telem.lon_deg >= 0.0 { "E" } else { "W" };
                        let lat_dir = if telem.lat_deg >= 0.0 { "N" } else { "S" };
                        let surface_tag = if telem.is_ice {
                            "Polar Ice"
                        } else if telem.is_land {
                            "Land"
                        } else {
                            "Ocean"
                        };
                        painter.text(
                            Pos2::new(card.left() + 10.0, card.top() + 14.0),
                            egui::Align2::LEFT_CENTER,
                            format!(
                                "Cursor: {:.2}°{lat_dir}, {:.2}°{lon_dir} [{surface_tag}]",
                                telem.lat_deg.abs(),
                                telem.lon_deg.abs()
                            ),
                            egui::FontId::proportional(10.5),
                            Color32::from_rgb(254, 240, 138),
                        );
                        let speed = telem.u_ms.hypot(telem.v_ms);
                        let thermo_txt = match telem.thermo_m {
                            Some(d) => format!("{d:.0} m"),
                            None => "Masked (outside Trop. Pac.)".to_string(),
                        };
                        painter.text(
                            Pos2::new(card.left() + 10.0, card.top() + 32.0),
                            egui::Align2::LEFT_CENTER,
                            format!(
                                "SST: {:.1}°C (Δ {:+.2}°C) | |V|: {:.2} m/s",
                                telem.sst_c, telem.anom_c, speed
                            ),
                            egui::FontId::proportional(10.0),
                            Color32::from_rgb(226, 232, 240),
                        );
                        painter.text(
                            Pos2::new(card.left() + 10.0, card.top() + 48.0),
                            egui::Align2::LEFT_CENTER,
                            format!("Flow u/v: {:+.2} / {:+.2} m/s", telem.u_ms, telem.v_ms),
                            egui::FontId::proportional(10.0),
                            Color32::from_rgb(167, 243, 208),
                        );
                        painter.text(
                            Pos2::new(card.left() + 10.0, card.top() + 63.0),
                            egui::Align2::LEFT_CENTER,
                            format!("20°C Thermocline: {thermo_txt}"),
                            egui::FontId::proportional(10.0),
                            Color32::from_rgb(125, 211, 252),
                        );
                    }
                }
            }
        });

        self.render_command_palette(ctx);
        if self.show_about_window {
            let mut open = self.show_about_window;
            egui::Window::new("About PelagicGlobe")
                .open(&mut open)
                .collapsible(false)
                .resizable(false)
                .show(ctx, |ui| {
                    ui.heading("PelagicGlobe v1.1.0");
                    ui.label("3D Ocean Currents Globe, Conformal Spilhaus World Ocean Layer & NOAA ENSO Workstation");
                    ui.separator();
                    ui.label(RichText::new("Copyright CRV Science Network").strong());
                    ui.label("License: GPL-3.0-or-later (GNU General Public License v3.0 or later)");
                    ui.add_space(4.0);
                    ui.horizontal(|ui| {
                        ui.label("CRV Science Network:");
                        ui.hyperlink_to("https://crvscience.com", "https://crvscience.com");
                    });
                    ui.horizontal(|ui| {
                        ui.label("Open-Source Apps:");
                        ui.hyperlink_to("https://crvscience.com/apps", "https://crvscience.com/apps");
                    });
                    ui.horizontal(|ui| {
                        ui.label("Source Repository:");
                        ui.hyperlink_to(
                            "bpwhite/crvscience-pelagic-globe",
                            "https://github.com/bpwhite/crvscience-pelagic-globe",
                        );
                    });
                    ui.add_space(6.0);
                    if ui.button("Close").clicked() {
                        self.show_about_window = false;
                    }
                });
            if !open {
                self.show_about_window = false;
            }
        }
    }

    fn export_enso_digest_file(&mut self, as_csv: bool) {
        let prefix = if self.export_prefix.trim().is_empty() {
            "pelagic_globe"
        } else {
            self.export_prefix.trim()
        };
        let digest = build_enso_comparative_digest();
        if as_csv {
            let path = format!("{prefix}_enso_digest.csv");
            match std::fs::write(&path, digest.to_csv()) {
                Ok(()) => self.status_note = format!("Exported ENSO comparative digest CSV -> {path}"),
                Err(e) => self.status_note = format!("Digest CSV export failed: {e}"),
            }
        } else {
            let path = format!("{prefix}_enso_digest.json");
            match digest.to_json_pretty().and_then(|j| {
                std::fs::write(&path, j).map_err(serde_json::Error::io)
            }) {
                Ok(()) => self.status_note = format!("Exported ENSO comparative digest JSON -> {path}"),
                Err(e) => self.status_note = format!("Digest JSON export failed: {e}"),
            }
        }
    }

    fn build_palette_entries(&self) -> Vec<PelagicPaletteEntry> {
        let mut entries = vec![
            PelagicPaletteEntry {
                category: "Morph",
                title: "3D Ocean Globe (0% Spilhaus Morph)".into(),
                subtitle: "Switch to interactive rotating 3D sphere (Hotkey: 1)".into(),
                action: PelagicPaletteAction::SetMorph(MapViewMode::Globe3D, 0.0),
            },
            PelagicPaletteEntry {
                category: "Morph",
                title: "2D Equirectangular World Map".into(),
                subtitle: "Switch to global equirectangular projection (Hotkey: 2)".into(),
                action: PelagicPaletteAction::SetMorph(MapViewMode::WorldMap2D, 0.0),
            },
            PelagicPaletteEntry {
                category: "Morph",
                title: "Hybrid 50% Globe <-> Spilhaus Morph".into(),
                subtitle: "Set mid-morph transition between 3D sphere and Spilhaus square".into(),
                action: PelagicPaletteAction::SetMorph(MapViewMode::Globe3D, 0.5),
            },
            PelagicPaletteEntry {
                category: "Morph",
                title: "Conformal Spilhaus World Ocean Square (100% Morph)".into(),
                subtitle: "Full ESRI:54099 Spilhaus World Ocean projection (Hotkey: 3)".into(),
                action: PelagicPaletteAction::SetMorph(MapViewMode::SpilhausSquare, 1.0),
            },
            PelagicPaletteEntry {
                category: "Inset",
                title: "Spilhaus Inset: Docked Square".into(),
                subtitle: "Render synchronized ESRI:54099 square inset in bottom-right".into(),
                action: PelagicPaletteAction::SetInset(SpilhausInsetMode::DockedSquare),
            },
            PelagicPaletteEntry {
                category: "Inset",
                title: "Spilhaus Inset: 3x3 Periodic Mosaic".into(),
                subtitle: "Render 3x3 tiled Adams/Spilhaus conformal mosaic inset".into(),
                action: PelagicPaletteAction::SetInset(SpilhausInsetMode::Mosaic3x3),
            },
            PelagicPaletteEntry {
                category: "Inset",
                title: "Spilhaus Inset: Hidden".into(),
                subtitle: "Hide the bottom-right Spilhaus inset panel".into(),
                action: PelagicPaletteAction::SetInset(SpilhausInsetMode::Hidden),
            },
            PelagicPaletteEntry {
                category: "Layer",
                title: "Layer: 20°C Thermocline Depth (m)".into(),
                subtitle: "Tropical Pacific 20°C isotherm overlay (120°E–80°W, |lat|<=20°)".into(),
                action: PelagicPaletteAction::SetLayer(DataLayer::ThermoclineDepth),
            },
            PelagicPaletteEntry {
                category: "Layer",
                title: "Layer: ENSO SST Anomaly (Δ°C)".into(),
                subtitle: "Diverging sea-surface temperature anomaly vs climatology".into(),
                action: PelagicPaletteAction::SetLayer(DataLayer::EnsoSstAnomaly),
            },
            PelagicPaletteEntry {
                category: "Layer",
                title: "Layer: Surface Current Velocity (m/s)".into(),
                subtitle: "Mixed-layer geostrophic & wind-driven flow speed".into(),
                action: PelagicPaletteAction::SetLayer(DataLayer::CurrentVelocity),
            },
            PelagicPaletteEntry {
                category: "Layer",
                title: "Layer: Absolute Sea Surface Temperature (°C)".into(),
                subtitle: "Seasonal climatology + historical ENSO/heatwave SST".into(),
                action: PelagicPaletteAction::SetLayer(DataLayer::AbsoluteSst),
            },
            PelagicPaletteEntry {
                category: "Overlay",
                title: "Toggle NOAA Niño 1+2 / 3 / 3.4 / 4 Monitoring Boxes".into(),
                subtitle: "Show or hide official equatorial Pacific monitoring boxes".into(),
                action: PelagicPaletteAction::ToggleNinoBoxes,
            },
            PelagicPaletteEntry {
                category: "Overlay",
                title: "Toggle Spilhaus Boundary Cut Seam & Grid".into(),
                subtitle: "Show or hide East Asia – Bering Strait – South America cut seam".into(),
                action: PelagicPaletteAction::ToggleSpilhausSeam,
            },
            PelagicPaletteEntry {
                category: "Overlay",
                title: "Toggle Lagrangian Current Drifters (4,200 RK2)".into(),
                subtitle: "Show or hide animated surface current ribbons".into(),
                action: PelagicPaletteAction::ToggleCurrents,
            },
            PelagicPaletteEntry {
                category: "Export",
                title: "Export ENSO Niño-3.4 / Thermocline Digest CSV".into(),
                subtitle: "Write comparative CSV across all historical El Niño / La Niña presets".into(),
                action: PelagicPaletteAction::ExportDigestCsv,
            },
            PelagicPaletteEntry {
                category: "Export",
                title: "Export ENSO Niño-3.4 / Thermocline Digest JSON".into(),
                subtitle: "Write structured comparative JSON across all historical ENSO presets".into(),
                action: PelagicPaletteAction::ExportDigestJson,
            },
            PelagicPaletteEntry {
                category: "Export",
                title: "Export 1950–2026 NOAA ONI CSV".into(),
                subtitle: "Write full ERSSTv6 3-month ONI historical table".into(),
                action: PelagicPaletteAction::ExportOniCsv,
            },
            PelagicPaletteEntry {
                category: "Export",
                title: "Export Conformal Spilhaus SVG".into(),
                subtitle: "Write Jacobian-rotated Spilhaus World Ocean SVG".into(),
                action: PelagicPaletteAction::ExportSpilhausSvg,
            },
        ];

        for preset in ENSO_PRESETS {
            let flavor = classify_historical_flavor(preset.year, preset.oni).label();
            entries.push(PelagicPaletteEntry {
                category: "ENSO",
                title: format!("{} [{}]", preset.title, preset.slug),
                subtitle: format!("{flavor} — {}", preset.summary),
                action: PelagicPaletteAction::SelectPreset(preset.slug),
            });
        }
        entries
    }

    fn apply_palette_action(&mut self, action: &PelagicPaletteAction) {
        match action {
            PelagicPaletteAction::SetMorph(mode, morph) => {
                self.view_mode = *mode;
                self.state.spilhaus_morph = *morph;
                if *mode == MapViewMode::WorldMap2D {
                    self.pan_2d = Vec2::ZERO;
                }
            }
            PelagicPaletteAction::SetInset(inset) => {
                self.state.spilhaus_inset_mode = *inset;
            }
            PelagicPaletteAction::SetLayer(layer) => {
                self.state.data_layer = *layer;
            }
            PelagicPaletteAction::ToggleNinoBoxes => {
                self.state.show_nino_boxes = !self.state.show_nino_boxes;
            }
            PelagicPaletteAction::ToggleSpilhausSeam => {
                self.state.show_spilhaus_overlay = !self.state.show_spilhaus_overlay;
            }
            PelagicPaletteAction::ToggleCurrents => {
                self.state.show_currents = !self.state.show_currents;
            }
            PelagicPaletteAction::SelectPreset(slug) => {
                if let Some(p) = ENSO_PRESETS.iter().find(|p| p.slug == *slug) {
                    self.timeline_playing = false;
                    self.apply_historical_date(
                        p.year,
                        p.season_index,
                        true,
                        Some(p.summary.to_string()),
                    );
                }
            }
            PelagicPaletteAction::ExportDigestCsv => {
                self.export_enso_digest_file(true);
            }
            PelagicPaletteAction::ExportDigestJson => {
                self.export_enso_digest_file(false);
            }
            PelagicPaletteAction::ExportOniCsv => {
                let prefix = if self.export_prefix.trim().is_empty() {
                    "pelagic_globe"
                } else {
                    self.export_prefix.trim()
                };
                let path = format!("{prefix}_oni_1950_2026.csv");
                match std::fs::write(&path, crate::enso::export_oni_csv()) {
                    Ok(()) => self.status_note = format!("Exported NOAA ONI CSV -> {path}"),
                    Err(e) => self.status_note = format!("CSV export failed: {e}"),
                }
            }
            PelagicPaletteAction::ExportSpilhausSvg => {
                let prefix = if self.export_prefix.trim().is_empty() {
                    "pelagic_globe"
                } else {
                    self.export_prefix.trim()
                };
                let season = SEASONS[self.state.season_index.min(11)];
                let path = format!("{prefix}_spilhaus_{}_{}.svg", self.state.year, season);
                let svg = export_spilhaus_svg(
                    self.state.year,
                    self.state.season_index,
                    self.state.oni_celsius,
                );
                match std::fs::write(&path, svg) {
                    Ok(()) => self.status_note = format!("Exported Spilhaus SVG -> {path}"),
                    Err(e) => self.status_note = format!("SVG export failed: {e}"),
                }
            }
        }
    }

    fn render_command_palette(&mut self, ctx: &egui::Context) {
        if !self.command_palette_open {
            return;
        }
        let mut keep_open = true;
        let enter_pressed = ctx.input(|i| i.key_pressed(egui::Key::Enter));
        let entries = self.build_palette_entries();
        let q = self.command_palette_query.trim();
        let mut scored: Vec<(i32, usize)> = entries
            .iter()
            .enumerate()
            .filter_map(|(idx, e)| {
                let hay = format!("{} {} {}", e.category, e.title, e.subtitle);
                fuzzy_match_score(q, &hay).map(|s| (s, idx))
            })
            .collect();
        scored.sort_by(|a, b| b.0.cmp(&a.0).then_with(|| a.1.cmp(&b.1)));

        let mut activated: Option<PelagicPaletteAction> = None;
        egui::Window::new("Quick-Jump Ocean & ENSO Palette (Ctrl+K)")
            .open(&mut keep_open)
            .collapsible(false)
            .resizable(false)
            .default_width(540.0)
            .anchor(egui::Align2::CENTER_TOP, Vec2::new(0.0, 48.0))
            .show(ctx, |ui| {
                ui.horizontal(|ui| {
                    ui.label(RichText::new("Search:").strong());
                    let resp = ui.add(
                        egui::TextEdit::singleline(&mut self.command_palette_query)
                            .id(egui::Id::new("pelagic_command_palette_input"))
                            .hint_text("Type Spilhaus morph, ENSO flavor, thermocline, or export...")
                            .desired_width(420.0),
                    );
                    if self.command_palette_just_opened {
                        resp.request_focus();
                        self.command_palette_just_opened = false;
                    }
                    if ui.small_button("Clear").clicked() {
                        self.command_palette_query.clear();
                        resp.request_focus();
                    }
                });
                ui.separator();

                let first_match = scored.first().map(|&(_, idx)| entries[idx].action.clone());
                egui::ScrollArea::vertical()
                    .max_height(340.0)
                    .show(ui, |ui| {
                        if scored.is_empty() {
                            ui.small("No matching Spilhaus modes, ENSO episodes, overlays, or exports.");
                        } else {
                            for &(_, idx) in scored.iter().take(16) {
                                let entry = &entries[idx];
                                let resp = ui
                                    .selectable_label(
                                        false,
                                        format!("[{}] {}", entry.category, entry.title),
                                    )
                                    .on_hover_text(&entry.subtitle);
                                ui.small(
                                    RichText::new(&entry.subtitle)
                                        .color(Color32::from_rgb(148, 163, 184)),
                                );
                                ui.add_space(2.0);
                                if resp.clicked() {
                                    activated = Some(entry.action.clone());
                                }
                            }
                        }
                    });

                if activated.is_none() && enter_pressed {
                    activated = first_match;
                }
            });

        if !keep_open {
            self.command_palette_open = false;
            self.command_palette_just_opened = false;
        }
        if let Some(action) = activated {
            self.apply_palette_action(&action);
            self.command_palette_open = false;
            self.command_palette_just_opened = false;
        }
    }
}

#[derive(Debug, Clone)]
enum PelagicPaletteAction {
    SetMorph(MapViewMode, f32),
    SetInset(SpilhausInsetMode),
    SetLayer(DataLayer),
    ToggleNinoBoxes,
    ToggleSpilhausSeam,
    ToggleCurrents,
    SelectPreset(&'static str),
    ExportDigestCsv,
    ExportDigestJson,
    ExportOniCsv,
    ExportSpilhausSvg,
}

#[derive(Debug, Clone)]
struct PelagicPaletteEntry {
    category: &'static str,
    title: String,
    subtitle: String,
    action: PelagicPaletteAction,
}

fn fuzzy_match_score(query: &str, text: &str) -> Option<i32> {
    let q = query.trim().to_ascii_lowercase();
    if q.is_empty() {
        return Some(0);
    }
    let t = text.to_ascii_lowercase();
    if let Some(pos) = t.find(&q) {
        return Some(1000 - pos as i32);
    }
    let mut t_chars = t.char_indices();
    let mut score = 200i32;
    let mut prev_idx: Option<usize> = None;
    for qc in q.chars() {
        if qc.is_whitespace() {
            continue;
        }
        let found = t_chars.find(|&(_, tc)| tc == qc);
        match found {
            Some((idx, _)) => {
                if let Some(prev) = prev_idx {
                    score -= (idx.saturating_sub(prev + 1)) as i32;
                }
                prev_idx = Some(idx);
            }
            None => return None,
        }
    }
    Some(score)
}

impl eframe::App for PelagicGlobeApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.ui_step(ctx, true);
    }
}

struct NoopWaker;

impl std::task::Wake for NoopWaker {
    fn wake(self: std::sync::Arc<Self>) {}
}

fn block_on_future<F: std::future::Future>(fut: F) -> F::Output {
    let mut fut = std::pin::pin!(fut);
    let waker = std::task::Waker::from(std::sync::Arc::new(NoopWaker));
    let mut cx = std::task::Context::from_waker(&waker);
    loop {
        match fut.as_mut().poll(&mut cx) {
            std::task::Poll::Ready(val) => return val,
            std::task::Poll::Pending => std::thread::yield_now(),
        }
    }
}

pub fn render_headless_screenshot(
    out_path: &std::path::Path,
    view: u8,
) -> Result<(), Box<dyn std::error::Error>> {
    use eframe::egui_wgpu::{wgpu, Renderer, ScreenDescriptor};

    if let Some(parent) = out_path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }

    let width: u32 = 1600;
    let height: u32 = 1000;
    let instance = wgpu::Instance::new(&wgpu::InstanceDescriptor::default());
    let adapter = block_on_future(instance.request_adapter(&wgpu::RequestAdapterOptions {
        power_preference: wgpu::PowerPreference::HighPerformance,
        compatible_surface: None,
        force_fallback_adapter: false,
    }))
    .ok_or("No GPU adapter available for offscreen rendering")?;
    let (device, queue) = block_on_future(adapter.request_device(
        &wgpu::DeviceDescriptor {
            label: Some("pelagic-globe-offscreen-device"),
            ..Default::default()
        },
        None,
    ))?;

    let format = wgpu::TextureFormat::Rgba8Unorm;
    let texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("pelagic-globe-offscreen-target"),
        size: wgpu::Extent3d {
            width,
            height,
            depth_or_array_layers: 1,
        },
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format,
        usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC,
        view_formats: &[],
    });
    let texture_view = texture.create_view(&wgpu::TextureViewDescriptor::default());
    let bytes_per_row = width * 4;
    let staging = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("pelagic-globe-offscreen-staging"),
        size: (bytes_per_row * height) as u64,
        usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });

    let mut renderer = Renderer::new(&device, format, None, 1, false);
    let egui_ctx = egui::Context::default();
    egui_ctx.set_pixels_per_point(1.0);
    egui_ctx.set_visuals(egui::Visuals::dark());

    let mut snap = SimulationSnapshot::default();
    snap.auto_rotate = false;
    snap.show_nino_boxes = true;
    snap.show_currents = true;

    match view {
        2 => {
            snap.year = 2015;
            snap.season_index = 11; // NDJ 2015
            snap.oni_celsius = 2.6;
            snap.spilhaus_morph = 0.0;
            snap.spilhaus_inset_mode = SpilhausInsetMode::Hidden;
        }
        3 => {
            snap.year = 1973;
            snap.season_index = 11; // NDJ 1973
            snap.oni_celsius = -2.1;
            snap.spilhaus_morph = 1.0;
            snap.show_spilhaus_overlay = true;
            snap.spilhaus_inset_mode = SpilhausInsetMode::Hidden;
        }
        _ => {
            snap.year = 1997;
            snap.season_index = 10; // OND 1997
            snap.oni_celsius = 2.4;
            snap.yaw_rad = 2.42;
            snap.pitch_rad = 0.08;
            snap.zoom = 1.0;
            snap.spilhaus_morph = 0.0;
            snap.show_spilhaus_overlay = true;
            snap.spilhaus_inset_mode = SpilhausInsetMode::DockedSquare;
            snap.selected_current_id = Some("necc_pacific".to_string());
        }
    }

    let mut app = PelagicGlobeApp::new_headless(snap, None);
    match view {
        2 => {
            app.view_mode = MapViewMode::WorldMap2D;
            app.show_world_map = true;
            app.show_graticule = true;
            app.apply_historical_date(2015, 11, false, None);
        }
        3 => {
            app.view_mode = MapViewMode::SpilhausSquare;
            app.state.spilhaus_morph = 1.0;
            app.apply_historical_date(1973, 11, false, None);
        }
        _ => {
            app.view_mode = MapViewMode::Globe3D;
            app.apply_historical_date(1997, 10, true, None);
            app.state.yaw_rad = 2.42;
            app.state.pitch_rad = 0.08;
        }
    }

    for _ in 0..24 {
        app.step_drifters(0.05);
    }

    let warmup_delays_ms = [40u64, 20, 20, 20];
    for (frame_idx, delay_ms) in warmup_delays_ms.into_iter().enumerate() {
        std::thread::sleep(std::time::Duration::from_millis(delay_ms));
        let raw_input = egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(width as f32, height as f32),
            )),
            viewports: std::iter::once((
                egui::ViewportId::ROOT,
                egui::ViewportInfo {
                    native_pixels_per_point: Some(1.0),
                    ..Default::default()
                },
            ))
            .collect(),
            time: Some(frame_idx as f64 / 60.0),
            predicted_dt: 1.0 / 60.0,
            ..Default::default()
        };
        let full_output = egui_ctx.run(raw_input, |ctx| {
            app.ui_step(ctx, false);
        });
        for (id, image_delta) in &full_output.textures_delta.set {
            renderer.update_texture(&device, &queue, *id, image_delta);
        }
        let paint_jobs = egui_ctx.tessellate(full_output.shapes, full_output.pixels_per_point);
        let screen_descriptor = ScreenDescriptor {
            size_in_pixels: [width, height],
            pixels_per_point: 1.0,
        };
        let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("pelagic-globe-offscreen-encoder"),
        });
        let user_cmd_bufs = renderer.update_buffers(
            &device,
            &queue,
            &mut encoder,
            &paint_jobs,
            &screen_descriptor,
        );
        {
            let render_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("pelagic-globe-offscreen-pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &texture_view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color::BLACK),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: None,
                occlusion_query_set: None,
                timestamp_writes: None,
            });
            let mut render_pass = render_pass.forget_lifetime();
            renderer.render(&mut render_pass, &paint_jobs, &screen_descriptor);
        }
        if frame_idx == 3 {
            encoder.copy_texture_to_buffer(
                wgpu::TexelCopyTextureInfo {
                    texture: &texture,
                    mip_level: 0,
                    origin: wgpu::Origin3d::ZERO,
                    aspect: wgpu::TextureAspect::All,
                },
                wgpu::TexelCopyBufferInfo {
                    buffer: &staging,
                    layout: wgpu::TexelCopyBufferLayout {
                        offset: 0,
                        bytes_per_row: Some(bytes_per_row),
                        rows_per_image: Some(height),
                    },
                },
                wgpu::Extent3d {
                    width,
                    height,
                    depth_or_array_layers: 1,
                },
            );
        }
        queue.submit(user_cmd_bufs.into_iter().chain(std::iter::once(encoder.finish())));
        device.poll(wgpu::Maintain::Wait);
        for id in &full_output.textures_delta.free {
            renderer.free_texture(id);
        }
    }

    let slice = staging.slice(..);
    let (tx, rx) = std::sync::mpsc::channel();
    slice.map_async(wgpu::MapMode::Read, move |res| {
        let _ = tx.send(res);
    });
    device.poll(wgpu::Maintain::Wait);
    rx.recv()??;
    let mapped = slice.get_mapped_range();
    image::save_buffer_with_format(
        out_path,
        &mapped,
        width,
        height,
        image::ColorType::Rgba8,
        image::ImageFormat::Png,
    )?;
    drop(mapped);
    staging.unmap();
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run_app_key(app: &mut PelagicGlobeApp, ctx: &egui::Context, key: egui::Key, modifiers: egui::Modifiers) {
        let _ = ctx.run(
            egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(
                    Pos2::ZERO,
                    Vec2::new(1280.0, 800.0),
                )),
                events: vec![egui::Event::Key {
                    key,
                    physical_key: None,
                    pressed: true,
                    repeat: false,
                    modifiers,
                }],
                ..Default::default()
            },
            |c| app.ui_step(c, true),
        );
    }

    #[test]
    fn command_palette_ctrl_k_opens_and_switches_spilhaus_morph() {
        let mut app = PelagicGlobeApp::new_headless(SimulationSnapshot::default(), None);
        let ctx = egui::Context::default();
        assert!(!app.command_palette_open);

        run_app_key(&mut app, &ctx, egui::Key::K, egui::Modifiers::COMMAND);
        assert!(app.command_palette_open, "Ctrl+K must open command palette");

        let entries = app.build_palette_entries();
        let spilhaus_entry = entries
            .iter()
            .filter_map(|e| {
                let hay = format!("{} {} {}", e.category, e.title, e.subtitle);
                fuzzy_match_score("spilhaus 100%", &hay).map(|s| (s, e))
            })
            .max_by_key(|x| x.0)
            .expect("matching Spilhaus morph entry")
            .1
            .clone();
        app.apply_palette_action(&spilhaus_entry.action);
        assert_eq!(app.view_mode, MapViewMode::SpilhausSquare);
        assert!((app.state.spilhaus_morph - 1.0).abs() < 1e-6);
    }

    #[test]
    fn command_palette_switches_enso_flavor_and_thermocline_layer() {
        let mut app = PelagicGlobeApp::new_headless(SimulationSnapshot::default(), None);
        let entries = app.build_palette_entries();

        let modoki = entries
            .iter()
            .filter_map(|e| {
                let hay = format!("{} {} {}", e.category, e.title, e.subtitle);
                fuzzy_match_score("modoki-2009", &hay).map(|s| (s, e))
            })
            .max_by_key(|x| x.0)
            .expect("modoki entry")
            .1
            .clone();
        app.apply_palette_action(&modoki.action);
        assert_eq!(app.state.year, 2009);
        assert!((app.state.oni_celsius - 1.5).abs() < 1e-6);

        let thermo = entries
            .iter()
            .filter_map(|e| {
                let hay = format!("{} {} {}", e.category, e.title, e.subtitle);
                fuzzy_match_score("thermocline depth", &hay).map(|s| (s, e))
            })
            .max_by_key(|x| x.0)
            .expect("thermocline entry")
            .1
            .clone();
        app.apply_palette_action(&thermo.action);
        assert_eq!(app.state.data_layer, DataLayer::ThermoclineDepth);
    }

    #[test]
    fn command_palette_escape_and_about_f1_behave_cleanly() {
        let mut app = PelagicGlobeApp::new_headless(SimulationSnapshot::default(), None);
        let ctx = egui::Context::default();

        run_app_key(&mut app, &ctx, egui::Key::K, egui::Modifiers::COMMAND);
        assert!(app.command_palette_open);
        run_app_key(&mut app, &ctx, egui::Key::Escape, egui::Modifiers::NONE);
        assert!(!app.command_palette_open);

        run_app_key(&mut app, &ctx, egui::Key::F1, egui::Modifiers::NONE);
        assert!(app.show_about_window);
        run_app_key(&mut app, &ctx, egui::Key::Escape, egui::Modifiers::NONE);
        assert!(!app.show_about_window);
    }

    #[test]
    fn fuzzy_match_scores_subsequences_and_rejects_missing_queries() {
        assert!(fuzzy_match_score("", "Anything").is_some());
        let sub = fuzzy_match_score("spilhaus", "Conformal Spilhaus World Ocean Square")
            .expect("substring match");
        let seq = fuzzy_match_score("cswo", "Conformal Spilhaus World Ocean Square")
            .expect("subsequence match");
        assert!(sub > seq);
        assert!(fuzzy_match_score("xqzv", "Conformal Spilhaus World Ocean Square").is_none());
    }
}
