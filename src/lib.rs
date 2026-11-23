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

pub mod app;
pub mod coastlines;
pub mod currents;
pub mod enso;
pub mod globe_wgsl;
pub mod icon;
pub mod raster;
pub mod spilhaus;
pub mod storage;

pub use crate::coastlines::{WorldLandMask, CONTINENT_POLYLINES};
use crate::currents::sample_velocity_ms;
use crate::enso::{spatial_sst_anomaly_historical, EnsoPhase, SEASONS};
use crate::raster::FieldSource;
use crate::spilhaus::{
    lonlat_to_spilhaus_norm, spilhaus_local_jacobian_velocity, spilhaus_norm_to_lonlat,
};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DataLayer {
    EnsoSstAnomaly,
    CurrentVelocity,
    AbsoluteSst,
    ThermoclineDepth,
}

impl DataLayer {
    pub fn label(self) -> &'static str {
        match self {
            Self::EnsoSstAnomaly => "ENSO SST Anomaly (Δ°C)",
            Self::CurrentVelocity => "Surface Velocity (m/s)",
            Self::AbsoluteSst => "Sea Surface Temp (°C)",
            Self::ThermoclineDepth => "20°C Thermocline Depth (m)",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SpilhausInsetMode {
    DockedSquare,
    Mosaic3x3,
    Hidden,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationSnapshot {
    pub schema_version: String,
    pub year: i32,
    pub season_index: usize,
    pub oni_celsius: f64,
    pub spilhaus_morph: f32,
    pub show_spilhaus_overlay: bool,
    pub spilhaus_inset_mode: SpilhausInsetMode,
    pub field_source: FieldSource,
    pub erddap_stride: usize,
    pub data_layer: DataLayer,
    pub show_nino_boxes: bool,
    pub show_currents: bool,
    pub particle_count: usize,
    pub flow_speed: f32,
    pub auto_rotate: bool,
    pub yaw_rad: f32,
    pub pitch_rad: f32,
    pub zoom: f32,
    pub selected_current_id: Option<String>,
}

impl Default for SimulationSnapshot {
    fn default() -> Self {
        Self {
            schema_version: "pelagic_globe.v1".to_string(),
            year: 2026,
            season_index: 8, // ASO 2026 (+2.5 C per NOAA CPC ONI v6)
            oni_celsius: 2.5,
            spilhaus_morph: 0.0,
            show_spilhaus_overlay: true,
            spilhaus_inset_mode: SpilhausInsetMode::DockedSquare,
            field_source: FieldSource::AnalyticalModel,
            erddap_stride: 8,
            data_layer: DataLayer::EnsoSstAnomaly,
            show_nino_boxes: true,
            show_currents: true,
            particle_count: 4200,
            flow_speed: 1.0,
            auto_rotate: true,
            yaw_rad: 2.55, // Centered over Equatorial Pacific & Americas
            pitch_rad: 0.14,
            zoom: 1.0,
            selected_current_id: Some("sec_pacific".to_string()),
        }
    }
}

/// Renders a standalone conformal Spilhaus World Ocean SVG (`ESRI:54099`) with Jacobian-rotated
/// surface velocity vectors and Natural-Earth continental coastlines.
pub fn export_spilhaus_svg(year: i32, season_index: usize, oni: f64) -> String {
    let size = 720.0_f64;
    let pad = 36.0_f64;
    let inner = size - 2.0 * pad;
    let mask = WorldLandMask::build();
    let mut svg = format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {size} {size}\" width=\"{size}\" height=\"{size}\">\n         <rect width=\"100%\" height=\"100%\" fill=\"#061325\"/>\n         <rect x=\"{pad}\" y=\"{pad}\" width=\"{inner}\" height=\"{inner}\" fill=\"#0a2240\" stroke=\"#facc15\" stroke-width=\"2\"/>\n"
    );

    let grid = 40;
    for r in 0..grid {
        let v = -0.96 + 1.92 * (r as f64 / (grid - 1) as f64);
        for c in 0..grid {
            let u = -0.96 + 1.92 * (c as f64 / (grid - 1) as f64);
            let (lon, lat) = spilhaus_norm_to_lonlat(u, v);
            if mask.is_land(lon, lat) {
                continue;
            }
            let anom = spatial_sst_anomaly_historical(lon, lat, year, season_index, oni);
            let (vx, vy) = sample_velocity_ms(lon, lat, oni);
            let (vu, vv) = spilhaus_local_jacobian_velocity(lon, lat, vx, vy);
            let px = pad + (u + 1.0) * 0.5 * inner;
            let py = pad + (1.0 - v) * 0.5 * inner;
            let stroke = if anom >= 0.4 {
                "#fb923c"
            } else if anom <= -0.4 {
                "#38bdf8"
            } else {
                "#34d399"
            };
            let dx = vu * 7.5;
            let dy = -vv * 7.5;
            let tx = px + dx;
            let ty = py + dy;
            svg.push_str(&format!(
                "<line x1=\"{px:.1}\" y1=\"{py:.1}\" x2=\"{tx:.1}\" y2=\"{ty:.1}\" stroke=\"{stroke}\" stroke-width=\"1.4\" stroke-opacity=\"0.82\"/>\n"
            ));
            let mag = dx.hypot(dy);
            if mag > 1.8 {
                let ux = dx / mag;
                let uy = dy / mag;
                let ax = tx - ux * 2.6 + uy * 1.4;
                let ay = ty - uy * 2.6 - ux * 1.4;
                let bx = tx - ux * 2.6 - uy * 1.4;
                let by = ty - uy * 2.6 + ux * 1.4;
                svg.push_str(&format!(
                    "<polygon points=\"{tx:.1},{ty:.1} {ax:.1},{ay:.1} {bx:.1},{by:.1}\" fill=\"{stroke}\" fill-opacity=\"0.85\"/>\n"
                ));
            }
        }
    }

    let max_jump = inner * 0.25;
    for poly in CONTINENT_POLYLINES {
        let mut pts_str = String::new();
        let mut pts_count = 0_usize;
        let mut prev: Option<(f64, f64, f64)> = None;
        for &(lat, lon) in *poly {
            let (su, sv) = lonlat_to_spilhaus_norm(lon, lat);
            let px = pad + (su + 1.0) * 0.5 * inner;
            let py = pad + (1.0 - sv) * 0.5 * inner;
            let should_break = match prev {
                Some((plon, ppx, ppy)) => {
                    (lon - plon).abs() > 180.0 || (px - ppx).hypot(py - ppy) > max_jump
                }
                None => false,
            };
            if should_break {
                if pts_count >= 2 {
                    svg.push_str(&format!(
                        "<polyline points=\"{pts_str}\" fill=\"none\" stroke=\"#dbe4ee\" stroke-width=\"1.2\" stroke-opacity=\"0.9\"/>\n"
                    ));
                }
                pts_str.clear();
                pts_count = 0;
            }
            if pts_count > 0 {
                pts_str.push(' ');
            }
            pts_str.push_str(&format!("{px:.1},{py:.1}"));
            pts_count += 1;
            prev = Some((lon, px, py));
        }
        if pts_count >= 2 {
            svg.push_str(&format!(
                "<polyline points=\"{pts_str}\" fill=\"none\" stroke=\"#dbe4ee\" stroke-width=\"1.2\" stroke-opacity=\"0.9\"/>\n"
            ));
        }
    }

    svg.push_str(&format!(
        "<text x=\"{pad}\" y=\"24\" fill=\"#e2e8f0\" font-family=\"monospace\" font-size=\"14\">Spilhaus World Ocean Square (ESRI:54099) | {year} {} | ONI {oni:+.1} C ({})</text>\n</svg>\n",
        SEASONS[season_index.min(11)],
        EnsoPhase::from_oni(oni).label()
    ));
    svg
}
