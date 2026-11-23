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

#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]
#![deny(unsafe_code)]

use std::fs;
use std::path::PathBuf;

use clap::Parser;
use pelagic_globe::app::{render_headless_screenshot, PelagicGlobeApp};
use pelagic_globe::currents::MAJOR_CURRENTS;
use pelagic_globe::enso::{
    equatorial_20c_isotherm_depth_historical_m, export_oni_csv, EnsoPhase, ENSO_PRESETS, SEASONS,
};
use pelagic_globe::globe_wgsl::build_spilhaus_globe_mesh;
use pelagic_globe::icon::app_egui_icon;
use pelagic_globe::raster::{
    fetch_live_noaa_erddap, load_raster_from_path, FieldSource, SatelliteRasterGrid,
};
use pelagic_globe::spilhaus::{
    lonlat_to_spilhaus_norm, sample_boundary_seam, spilhaus_norm_to_lonlat, SPILHAUS_LAT0_DEG,
    SPILHAUS_LON0_DEG,
};
use pelagic_globe::{export_spilhaus_svg, SimulationSnapshot};

#[cfg(windows)]
#[allow(unsafe_code)]
fn attach_parent_console_if_cli_args() {
    if std::env::args_os().len() > 1 {
        unsafe {
            let _ = windows_sys::Win32::System::Console::AttachConsole(u32::MAX);
        }
    }
}

#[derive(Parser, Debug)]
#[command(name = "pelagic-globe", version, about = "3D Ocean Currents, Spilhaus World Ocean Layer & NOAA ENSO Workstation")]
struct Cli {
    /// Select an ENSO historical episode slug (e.g., super-2026, super-2015, super-1997, lanina-1973)
    #[arg(long)]
    episode: Option<String>,

    /// Override Oceanic Niño Index (ONI) value in degrees Celsius
    #[arg(long, allow_hyphen_values = true)]
    oni: Option<f64>,

    /// Run headless mathematical & simulation self-tests and exit
    #[arg(long)]
    self_test: bool,

    /// Project a single `lon,lat` coordinate into Spilhaus normalized & raw space
    #[arg(long, allow_hyphen_values = true)]
    project_coord: Option<String>,

    /// Export embedded 1950–2026 NOAA CPC ERSSTv6 ONI table to CSV and exit
    #[arg(long)]
    dump_oni_csv: Option<PathBuf>,

    /// Export a standalone conformal Spilhaus World Ocean SVG layer for the active ONI state
    #[arg(long)]
    export_spilhaus_svg: Option<PathBuf>,

    /// Fetch live NOAA CoastWatch ERDDAP OISST v2.1 + NESDIS currents for YYYY-MM-DD or 'last'
    #[arg(long)]
    fetch_erddap: Option<String>,

    /// Grid stride for ERDDAP requests (4 = 1.0 deg, 8 = 2.0 deg)
    #[arg(long, default_value_t = 8)]
    erddap_stride: usize,

    /// Import a local ERDDAP .json, .csv, or equirectangular .png/.jpg satellite raster
    #[arg(long)]
    import_raster: Option<PathBuf>,

    /// Render an offscreen 1600x1000 RGBA PNG screenshot without opening a GUI window
    #[arg(long)]
    headless_screenshot: Option<PathBuf>,

    /// Select view mode (1=3D Globe 1997 El Nino, 2=2D World Map 2015 El Nino, 3=Spilhaus 1973 La Nina)
    #[arg(long, default_value_t = 1)]
    screenshot_view: u8,
}

fn main() -> eframe::Result<()> {
    #[cfg(windows)]
    attach_parent_console_if_cli_args();

    let cli = Cli::parse();

    if let Some(out_path) = &cli.headless_screenshot {
        if let Err(err) = render_headless_screenshot(out_path, cli.screenshot_view) {
            eprintln!("Headless screenshot error: {err}");
            std::process::exit(1);
        }
        println!(
            "Rendered offscreen 1600x1000 screenshot (view {}) -> {}",
            cli.screenshot_view,
            out_path.display()
        );
        return Ok(());
    }

    let mut snapshot = SimulationSnapshot::default();
    let mut loaded_raster: Option<SatelliteRasterGrid> = None;

    if let Some(slug) = &cli.episode {
        if let Some(preset) = ENSO_PRESETS.iter().find(|p| p.slug == slug) {
            snapshot.year = preset.year;
            snapshot.season_index = preset.season_index;
            snapshot.oni_celsius = preset.oni;
        }
    }
    if let Some(oni) = cli.oni {
        snapshot.oni_celsius = oni;
    }

    if let Some(coord) = &cli.project_coord {
        let parts: Vec<&str> = coord.split(',').collect();
        if parts.len() == 2 {
            let lon: f64 = parts[0].trim().parse().unwrap_or(0.0);
            let lat: f64 = parts[1].trim().parse().unwrap_or(0.0);
            let (u, v) = lonlat_to_spilhaus_norm(lon, lat);
            let (rlon, rlat) = spilhaus_norm_to_lonlat(u, v);
            println!(
                "Input ({lon:.4}°, {lat:.4}°) -> Spilhaus norm ({u:+.6}, {v:+.6}) -> Inverse ({rlon:.4}°, {rlat:.4}°)"
            );
        }
        return Ok(());
    }

    if let Some(path) = &cli.dump_oni_csv {
        fs::write(path, export_oni_csv()).expect("Failed to write ONI CSV");
        println!("Wrote 1950-2026 NOAA CPC ONI CSV to {}", path.display());
        return Ok(());
    }

    if let Some(path) = &cli.export_spilhaus_svg {
        fs::write(
            path,
            export_spilhaus_svg(snapshot.year, snapshot.season_index, snapshot.oni_celsius),
        )
        .expect("Failed to write Spilhaus SVG");
        println!("Wrote conformal Spilhaus SVG layer to {}", path.display());
        return Ok(());
    }

    if let Some(date_str) = &cli.fetch_erddap {
        let grid = fetch_live_noaa_erddap(date_str, cli.erddap_stride)
            .expect("NOAA ERDDAP fetch failed");
        println!(
            "Fetched {} ({}x{}, {} valid ocean cells, timestamp={})",
            grid.source_label,
            grid.width,
            grid.height,
            grid.valid_cell_count(),
            grid.timestamp_utc
        );
        snapshot.field_source = FieldSource::SatelliteRaster;
        snapshot.erddap_stride = cli.erddap_stride;
        loaded_raster = Some(grid);
    }

    if let Some(raster_path) = &cli.import_raster {
        let grid = load_raster_from_path(raster_path).expect("Failed to load local raster");
        println!(
            "Imported local raster {} ({}x{}, {} valid cells)",
            grid.source_label,
            grid.width,
            grid.height,
            grid.valid_cell_count()
        );
        snapshot.field_source = FieldSource::SatelliteRaster;
        loaded_raster = Some(grid);
    }

    if cli.self_test {
        let (u0, v0) = lonlat_to_spilhaus_norm(SPILHAUS_LON0_DEG, SPILHAUS_LAT0_DEG);
        let seam = sample_boundary_seam(32);
        let (verts, idx) = build_spilhaus_globe_mesh(64, snapshot.oni_celsius);
        let d_east = equatorial_20c_isotherm_depth_historical_m(275.0, snapshot.year, snapshot.oni_celsius);
        println!("=== PelagicGlobe v1.1.0 Self-Test ===");
        println!("Origin (66.95E, 49.56S) -> Spilhaus ({u0:+.8}, {v0:+.8})");
        println!("Seam vertices: {}, 3D morph mesh: {} verts / {} tris", seam.len(), verts.len(), idx.len() / 3);
        println!(
            "Episode {} {} | Active ONI: {:+.1} C ({}) | Eastern Pacific 20C isotherm: {d_east:.1} m",
            snapshot.year,
            SEASONS[snapshot.season_index.min(11)],
            snapshot.oni_celsius,
            EnsoPhase::from_oni(snapshot.oni_celsius).label()
        );
        println!("Verified currents loaded: {}", MAJOR_CURRENTS.len());
        return Ok(());
    }

    let native_options = eframe::NativeOptions {
        viewport: eframe::egui::ViewportBuilder::default()
            .with_inner_size([1360.0, 840.0])
            .with_icon(app_egui_icon())
            .with_title("PelagicGlobe v1.1.0 — 3D Ocean Currents, Spilhaus Layer & NOAA ENSO Workstation"),
        ..Default::default()
    };

    eframe::run_native(
        "PelagicGlobe",
        native_options,
        Box::new(move |_cc| Ok(Box::new(PelagicGlobeApp::new_with_raster(snapshot, loaded_raster)))),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cli_parses_headless_screenshot_and_view() {
        let parsed = Cli::try_parse_from([
            "pelagic-globe",
            "--headless-screenshot",
            "validation/headless_screenshot_3.png",
            "--screenshot-view",
            "3",
        ])
        .expect("clap parses headless screenshot arguments");
        assert_eq!(
            parsed.headless_screenshot,
            Some(PathBuf::from("validation/headless_screenshot_3.png"))
        );
        assert_eq!(parsed.screenshot_view, 3);
    }
}
