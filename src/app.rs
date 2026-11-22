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
    classify_historical_flavor, climatological_sst_c, equatorial_20c_isotherm_depth_historical_m,
    latest_valid_season_for_year, lookup_oni, lookup_oni_clamped, spatial_sst_anomaly_historical,
    EnsoPhase, END_YEAR, ENSO_PRESETS, NINO_BOXES, SEASONS, START_YEAR,
};
use crate::raster::{
    erddap_oisst_supported, fetch_live_noaa_erddap, FieldSource, SatelliteRasterGrid,
};
use crate::spilhaus::{
    lonlat_to_spilhaus_norm, sample_boundary_seam, spilhaus_norm_to_lonlat,
};
use crate::storage::BookmarkStore;
use crate::{DataLayer, SimulationSnapshot, SpilhausInsetMode};

const SPILHAUS_GRID_N: usize = 96;
const EQUIRECT_COLS: usize = 120;
const EQUIRECT_ROWS: usize = 60;
const INSET_GRID_N: usize = 48;

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
    valid: bool,
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
    coastline_rings: Vec<Vec<ProjectedCoastPt>>,
    seam_cache: Vec<(f64, f64)>,
    ocean_seeds: Vec<(f64, f64)>,
    drifters: Vec<(f64, f64)>,
    active_raster: Option<SatelliteRasterGrid>,
    erddap_rx: Option<std::sync::mpsc::Receiver<(String, Result<SatelliteRasterGrid, String>)>>,
    erddap_in_flight: Option<String>,
    timeline_playing: bool,
    timeline_accum_sec: f32,
    status_note: String,
    pub show_about_window: bool,
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
            coastline_rings,
            seam_cache,
            ocean_seeds,
            drifters,
            active_raster: raster,
            erddap_rx: None,
            erddap_in_flight: None,
            timeline_playing: false,
            timeline_accum_sec: 0.0,
            status_note,
            show_about_window: false,
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

        if let Some(cached) = self
            .store
            .as_ref()
            .and_then(|s| s.get_cached_raster(&iso).ok().flatten())
        {
            let cells = cached.valid_cell_count();
            let has_vel = cached.has_velocity;
            self.active_raster = Some(cached);
            self.state.field_source = FieldSource::SatelliteRaster;
            self.status_note = if has_vel {
                format!("Loaded cached NOAA ERDDAP raster ({iso}, {cells} cells, OISST + NESDIS currents).")
            } else {
                format!("Loaded cached NOAA OISST v2.1 raster ({iso}, {cells} cells; analytical currents pre-2017).")
            };
            return;
        }

        let stride = self.state.erddap_stride;
        let iso_thread = iso.clone();
        let (tx, rx) = std::sync::mpsc::channel();
        self.erddap_in_flight = Some(iso.clone());
        self.erddap_rx = Some(rx);
        self.status_note = format!("Fetching NOAA ERDDAP {iso} in background...");
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

    fn effective_morph(&self) -> f32 {
        match self.view_mode {
            MapViewMode::SpilhausSquare => 1.0,
            MapViewMode::WorldMap2D => 0.0,
            MapViewMode::Globe3D => self.state.spilhaus_morph.clamp(0.0, 1.0),
        }
    }

    fn ocean_color(&self, lon_deg: f64, lat_deg: f64) -> [f32; 3] {
        let mut sst_override = None;
        let mut anom_override = None;
        let mut vel_override = None;
        if self.state.field_source == FieldSource::SatelliteRaster {
            if let Some(grid) = &self.active_raster {
                if let Some((sst, anom, u, v, _sla)) = grid.sample_bilinear(lon_deg, lat_deg) {
                    sst_override = Some(sst);
                    anom_override = Some(anom);
                    if grid.has_velocity {
                        vel_override = Some((u, v));
                    }
                }
            }
        }

        match self.state.data_layer {
            DataLayer::EnsoSstAnomaly => {
                let anom = anom_override.unwrap_or_else(|| {
                    spatial_sst_anomaly_historical(
                        lon_deg,
                        lat_deg,
                        self.state.year,
                        self.state.season_index,
                        self.state.oni_celsius,
                    )
                }) as f32;
                diverged_anom_rgb(anom)
            }
            DataLayer::CurrentVelocity => {
                let (vx, vy) = vel_override
                    .unwrap_or_else(|| sample_velocity_ms(lon_deg, lat_deg, self.state.oni_celsius));
                let t = ((vx * vx + vy * vy).sqrt() as f32 / 1.4).clamp(0.0, 1.0);
                if t < 0.5 {
                    lerp_rgb([8.0, 24.0, 54.0], [24.0, 154.0, 188.0], t * 2.0)
                } else {
                    lerp_rgb([24.0, 154.0, 188.0], [250.0, 204.0, 21.0], (t - 0.5) * 2.0)
                }
            }
            DataLayer::AbsoluteSst => {
                let sst = sst_override.unwrap_or_else(|| {
                    let base = climatological_sst_c(lon_deg, lat_deg, self.state.season_index);
                    let anom = spatial_sst_anomaly_historical(
                        lon_deg,
                        lat_deg,
                        self.state.year,
                        self.state.season_index,
                        self.state.oni_celsius,
                    );
                    (base + anom).clamp(-2.0, 33.5)
                }) as f32;
                let t = ((sst + 2.0) / 34.0).clamp(0.0, 1.0);
                if t < 0.5 {
                    lerp_rgb([16.0, 32.0, 98.0], [20.0, 148.0, 162.0], t * 2.0)
                } else {
                    lerp_rgb([20.0, 148.0, 162.0], [246.0, 108.0, 62.0], (t - 0.5) * 2.0)
                }
            }
            DataLayer::ThermoclineDepth => {
                let lon360 = lon_deg.rem_euclid(360.0);
                let depth = equatorial_20c_isotherm_depth_historical_m(
                    lon360,
                    self.state.year,
                    self.state.oni_celsius,
                ) as f32;
                let t = ((depth - 20.0) / 170.0).clamp(0.0, 1.0);
                if t < 0.5 {
                    lerp_rgb([28.0, 34.0, 108.0], [18.0, 158.0, 126.0], t * 2.0)
                } else {
                    lerp_rgb([18.0, 158.0, 126.0], [244.0, 164.0, 44.0], (t - 0.5) * 2.0)
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

        let raw_ndotl = (x3 * 0.48 + y3 * 0.48 + z3 * 0.74).clamp(0.0, 1.0);
        let lit = (0.34 + 0.66 * raw_ndotl) * (1.0 - morph) + morph;

        ProjectedScreenVert {
            pos: Pos2::new(center.x + bx * radius, center.y + by * radius),
            z3,
            lit,
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
            let mid_lon = p.0 + 0.5 * scale * u1;
            let mid_lat = (p.1 + 0.5 * scale * v1).clamp(-82.0, 82.0);
            let (u2, v2) = self.sample_effective_velocity(mid_lon, mid_lat);
            let next_lon = (p.0 + scale * u2 + 180.0).rem_euclid(360.0) - 180.0;
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
            let base_rgb = if v.is_land && self.show_world_map {
                land_color(v.lon_deg, v.lat_deg, v.is_ice)
            } else {
                self.ocean_color(v.lon_deg, v.lat_deg)
            };
            colors.push(to_color32(base_rgb, pv.lit));
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
            if (pair[0].0 - pair[1].0).abs() > 180.0 {
                continue;
            }
            if let (Some((a, _)), Some((b, _))) = (
                self.project_lonlat(pair[0].0, pair[0].1, rect),
                self.project_lonlat(pair[1].0, pair[1].1, rect),
            ) {
                if a.distance(b) <= jump_limit {
                    painter.line_segment([a, b], track_stroke);
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
        if self.state.spilhaus_inset_mode == SpilhausInsetMode::Hidden {
            return;
        }
        let side = 184.0_f32;
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
        let n = INSET_GRID_N;
        let stride = n + 1;
        let mut mesh = egui::epaint::Mesh::default();
        for v in &self.inset_grid {
            let px = inner.center().x + v.su * (inner.width() * 0.5);
            let py = inner.center().y - v.sv * (inner.height() * 0.5);
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
            let inset_shore = Stroke::new(0.8_f32, Color32::from_rgba_unmultiplied(226, 232, 240, 215));
            let inset_jump = inner.width() * 0.25;
            for ring in &self.coastline_rings {
                for seg in ring.windows(2) {
                    if (seg[0].lon_deg - seg[1].lon_deg).abs() > 180.0 {
                        continue;
                    }
                    let ax = inner.center().x + seg[0].su * (inner.width() * 0.5);
                    let ay = inner.center().y - seg[0].sv * (inner.height() * 0.5);
                    let bx = inner.center().x + seg[1].su * (inner.width() * 0.5);
                    let by = inner.center().y - seg[1].sv * (inner.height() * 0.5);
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
                let px = inner.center().x + (u as f32) * (inner.width() * 0.5);
                let py = inner.center().y - (v as f32) * (inner.height() * 0.5);
                painter.circle_filled(
                    Pos2::new(px, py),
                    1.0,
                    Color32::from_rgba_unmultiplied(125, 211, 252, 210),
                );
            }
        }

        painter.rect_stroke(
            inner,
            CornerRadius::ZERO,
            Stroke::new(1.1_f32, Color32::from_rgb(250, 204, 21)),
            egui::StrokeKind::Inside,
        );
        painter.rect_stroke(
            inset,
            CornerRadius::same(6),
            Stroke::new(1.2_f32, Color32::from_rgb(56, 189, 248)),
            egui::StrokeKind::Outside,
        );
        painter.text(
            Pos2::new(inset.center().x, inset.top() + 9.0),
            egui::Align2::CENTER_CENTER,
            "Spilhaus Square (ESRI:54099)",
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
                if let Ok((iso, res)) = rx.try_recv() {
                    self.erddap_rx = None;
                    self.erddap_in_flight = None;
                    match res {
                        Ok(grid) => {
                            let cells = grid.valid_cell_count();
                            let has_vel = grid.has_velocity;
                            if let Some(store) = &self.store {
                                let _ = store.put_cached_raster(&iso, &grid);
                            }
                            self.active_raster = Some(grid);
                            self.state.field_source = FieldSource::SatelliteRaster;
                            self.status_note = if has_vel {
                                format!("Loaded NOAA ERDDAP raster ({iso}, {cells} cells, OISST + NESDIS currents).")
                            } else {
                                format!("Loaded NOAA OISST v2.1 raster ({iso}, {cells} cells; analytical currents pre-2017).")
                            };
                        }
                        Err(err) => {
                            self.status_note = err;
                        }
                    }
                }
            }

            ctx.input(|i| {
                if i.key_pressed(egui::Key::F1) {
                    self.show_about_window = !self.show_about_window;
                }
                if i.key_pressed(egui::Key::Escape) && self.show_about_window {
                    self.show_about_window = false;
                }
            });

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
                    .selectable_label(self.view_mode == MapViewMode::Globe3D && self.state.spilhaus_morph < 0.99, "3D Globe")
                    .clicked()
                {
                    self.view_mode = MapViewMode::Globe3D;
                    self.state.spilhaus_morph = 0.0;
                }
                if ui
                    .selectable_label(self.view_mode == MapViewMode::WorldMap2D, "World Map (2D)")
                    .clicked()
                {
                    self.view_mode = MapViewMode::WorldMap2D;
                    self.pan_2d = Vec2::ZERO;
                }
                if ui
                    .selectable_label(self.view_mode == MapViewMode::SpilhausSquare, "Spilhaus Square")
                    .clicked()
                {
                    self.view_mode = MapViewMode::SpilhausSquare;
                    self.state.spilhaus_morph = 1.0;
                }

                ui.separator();
                if ui.button("Save SQLite Bookmark").clicked() {
                    if let Some(store) = &self.store {
                        let label = format!(
                            "{} {} ({:+.1}C)",
                            self.state.year,
                            SEASONS[self.state.season_index.min(11)],
                            self.state.oni_celsius
                        );
                        if store.save_bookmark(&label, &self.state).is_ok() {
                            self.saved_bookmarks = store.list_bookmarks().unwrap_or_default();
                            self.status_note =
                                format!("Saved bookmark '{label}' to WAL SQLite store.");
                        }
                    }
                }
                ui.separator();
                if ui.button("About").clicked() {
                    self.show_about_window = !self.show_about_window;
                }
            });
        });

        egui::SidePanel::left("controls_panel")
            .default_width(308.0)
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
                                .changed()
                            {
                                self.view_mode = if self.state.spilhaus_morph >= 0.99 {
                                    MapViewMode::SpilhausSquare
                                } else {
                                    MapViewMode::Globe3D
                                };
                            }
                            ui.checkbox(&mut self.show_world_map, "World Map & Coastlines");
                            ui.checkbox(&mut self.show_graticule, "Lat/Lon Graticule");
                            ui.checkbox(
                                &mut self.state.show_spilhaus_overlay,
                                "Spilhaus Seam & Conformal Grid",
                            );
                            ui.checkbox(&mut self.state.auto_rotate, "Auto-Rotate 3D Globe");
                            ui.checkbox(&mut self.state.show_currents, "Lagrangian Current Ribbons");
                            ui.checkbox(
                                &mut self.state.show_nino_boxes,
                                "NOAA Niño 1+2 / 3 / 3.4 / 4 Boxes",
                            );
                            ui.add(
                                egui::Slider::new(&mut self.state.flow_speed, 0.25..=2.5)
                                    .text("Advection Speed"),
                            );

                            ui.separator();
                            ui.heading("Field Source (Model vs. Satellite Raster)");
                            ui.selectable_value(
                                &mut self.state.field_source,
                                FieldSource::AnalyticalModel,
                                FieldSource::AnalyticalModel.label(),
                            );
                            if ui
                                .selectable_label(
                                    self.state.field_source == FieldSource::SatelliteRaster,
                                    FieldSource::SatelliteRaster.label(),
                                )
                                .clicked()
                            {
                                self.state.field_source = FieldSource::SatelliteRaster;
                                if self.active_raster.is_none() {
                                    self.status_note = "No satellite raster loaded yet — click 'Fetch NOAA Satellite Raster' below.".to_string();
                                }
                            }
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
                                ui.selectable_value(
                                    &mut self.state.data_layer,
                                    layer,
                                    layer.label(),
                                );
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
                                    if ui.selectable_label(is_active, preset.title).clicked() {
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
                                    if self.saved_bookmarks.is_empty() {
                                        ui.small("No saved bookmarks yet — click Save SQLite Bookmark above.");
                                    } else {
                                        let mut restored: Option<(String, SimulationSnapshot)> = None;
                                        for (id, name, snap) in &self.saved_bookmarks {
                                            if ui.button(format!("#{id} {name}")).clicked() {
                                                restored = Some((name.clone(), snap.clone()));
                                            }
                                        }
                                        if let Some((name, snap)) = restored {
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
                            ui.heading("Major Ocean Currents (14)");
                            for cur in MAJOR_CURRENTS {
                                let selected =
                                    self.state.selected_current_id.as_deref() == Some(cur.id);
                                if ui
                                    .selectable_label(
                                        selected,
                                        format!("{} ({})", cur.name, cur.basin),
                                    )
                                    .clicked()
                                {
                                    self.state.selected_current_id = Some(cur.id.to_string());
                                    self.state.yaw_rad = -(cur.focus_lon_deg as f32).to_radians();
                                    self.state.pitch_rad =
                                        (cur.focus_lat_deg as f32).to_radians() * 0.55;
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
                            .button(if self.timeline_playing { "Pause" } else { "Play 1950->2026" })
                            .clicked()
                        {
                            self.timeline_playing = !self.timeline_playing;
                            self.timeline_accum_sec = 0.0;
                        }
                        if ui.button("< Prev").clicked() {
                            self.timeline_playing = false;
                            let (py, ps) = if self.state.season_index > 0 {
                                (self.state.year, self.state.season_index - 1)
                            } else if self.state.year > START_YEAR {
                                (self.state.year - 1, latest_valid_season_for_year(self.state.year - 1))
                            } else {
                                (END_YEAR, latest_valid_season_for_year(END_YEAR))
                            };
                            self.apply_historical_date(py, ps, false, None);
                        }
                        if ui.button("Next >").clicked() {
                            self.timeline_playing = false;
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
                self.pan_2d = Vec2::ZERO;
                self.state.zoom = 1.0;
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
                        radius * (1.02 - 0.12 * morph),
                        Color32::from_rgba_unmultiplied(16, 48, 86, 95),
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
                        let anom = if self.state.field_source == FieldSource::SatelliteRaster {
                            self.active_raster
                                .as_ref()
                                .and_then(|g| g.sample_bilinear(lon, lat).map(|s| s.1))
                                .unwrap_or_else(|| {
                                    spatial_sst_anomaly_historical(
                                        lon,
                                        lat,
                                        self.state.year,
                                        self.state.season_index,
                                        self.state.oni_celsius,
                                    )
                                })
                        } else {
                            spatial_sst_anomaly_historical(
                                lon,
                                lat,
                                self.state.year,
                                self.state.season_index,
                                self.state.oni_celsius,
                            )
                        };
                        let col = if anom > 0.35 {
                            Color32::from_rgba_unmultiplied(251, 146, 60, 210)
                        } else if anom < -0.35 {
                            Color32::from_rgba_unmultiplied(56, 189, 248, 210)
                        } else {
                            Color32::from_rgba_unmultiplied(167, 243, 208, 175)
                        };
                        let (vx, vy) = self.sample_effective_velocity(lon, lat);
                        let next_lon = (lon + vx * 2.2).clamp(-179.9, 179.9);
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
        });

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
            snap.selected_current_id = Some("equatorial_counter".to_string());
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
