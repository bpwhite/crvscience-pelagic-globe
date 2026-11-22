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

use crate::raster::FieldSource;
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
