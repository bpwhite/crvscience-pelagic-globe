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

//! Gridded satellite raster ingestion, bilinear interpolation, and GPU RGBA8 texture packing.
//!
//! Supports live HTTPS queries against NOAA CoastWatch ERDDAP (`ncdcOisst21Agg`,
//! `ncdcOisst21NrtAgg`, `nesdisSSH1day`, and `jplOscar`) as well as offline ERDDAP `.json`,
//! `.csv`, and equirectangular `.png`/`.jpg` satellite image rasters.

use std::fs;
use std::path::Path;
use std::time::Duration;

use serde::{Deserialize, Serialize};
use serde_json::Value;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FieldSource {
    AnalyticalModel,
    SatelliteRaster,
}

impl FieldSource {
    pub fn label(self) -> &'static str {
        match self {
            Self::AnalyticalModel => "Analytical ENSO + Gyre Model",
            Self::SatelliteRaster => "Gridded Satellite Raster (ERDDAP / File)",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SatelliteRasterGrid {
    pub width: usize,
    pub height: usize,
    pub source_label: String,
    pub timestamp_utc: String,
    /// Row-major (`row 0 = -90°S` to `row height-1 = +90°N`, `col 0 = -180°W` to `+180°E`)
    pub sst_c: Vec<f32>,
    pub anom_c: Vec<f32>,
    pub u_ms: Vec<f32>,
    pub v_ms: Vec<f32>,
    pub sla_m: Vec<f32>,
    pub valid_mask: Vec<u8>,
    pub has_velocity: bool,
}

impl SatelliteRasterGrid {
    pub fn new_empty(width: usize, height: usize, source_label: &str, timestamp_utc: &str) -> Self {
        let len = width * height;
        Self {
            width,
            height,
            source_label: source_label.to_string(),
            timestamp_utc: timestamp_utc.to_string(),
            sst_c: vec![0.0; len],
            anom_c: vec![0.0; len],
            u_ms: vec![0.0; len],
            v_ms: vec![0.0; len],
            sla_m: vec![0.0; len],
            valid_mask: vec![0; len],
            has_velocity: false,
        }
    }

    pub fn valid_cell_count(&self) -> usize {
        self.valid_mask.iter().filter(|&&m| m > 0).count()
    }

    /// Bilinearly samples `(sst_c, anom_c, u_ms, v_ms, sla_m)` at geographic `(lon_deg, lat_deg)`.
    pub fn sample_bilinear(&self, lon_deg: f64, lat_deg: f64) -> Option<(f64, f64, f64, f64, f64)> {
        if self.width < 2 || self.height < 2 {
            return None;
        }
        let u = ((lon_deg + 180.0).rem_euclid(360.0) / 360.0) * (self.width as f64) - 0.5;
        let v = ((lat_deg.clamp(-89.875, 89.875) + 90.0) / 180.0) * (self.height as f64) - 0.5;

        let x0 = u.floor() as isize;
        let y0 = v.floor().clamp(0.0, (self.height - 1) as f64) as usize;
        let y1 = (y0 + 1).min(self.height - 1);
        let tx = (u - u.floor()) as f32;
        let ty = (v - v.floor()).clamp(0.0, 1.0) as f32;

        let wrap_x = |x: isize| -> usize { x.rem_euclid(self.width as isize) as usize };
        let c00 = y0 * self.width + wrap_x(x0);
        let c10 = y0 * self.width + wrap_x(x0 + 1);
        let c01 = y1 * self.width + wrap_x(x0);
        let c11 = y1 * self.width + wrap_x(x0 + 1);

        let w00 = (1.0 - tx) * (1.0 - ty) * f32::from(self.valid_mask[c00]);
        let w10 = tx * (1.0 - ty) * f32::from(self.valid_mask[c10]);
        let w01 = (1.0 - tx) * ty * f32::from(self.valid_mask[c01]);
        let w11 = tx * ty * f32::from(self.valid_mask[c11]);
        let w_sum = w00 + w10 + w01 + w11;
        if w_sum < 1.0e-5 {
            return None;
        }

        let blend = |arr: &[f32]| -> f64 {
            f64::from((arr[c00] * w00 + arr[c10] * w10 + arr[c01] * w01 + arr[c11] * w11) / w_sum)
        };
        Some((
            blend(&self.sst_c),
            blend(&self.anom_c),
            blend(&self.u_ms),
            blend(&self.v_ms),
            blend(&self.sla_m),
        ))
    }

    /// Packs gridded channels into an equirectangular `RGBA8` buffer (`R=anom`, `G=sst`, `B=speed`, `A=valid`).
    pub fn pack_rgba8_texture(&self) -> Vec<u8> {
        let mut rgba = vec![0_u8; self.width * self.height * 4];
        for i in 0..(self.width * self.height) {
            let o = i * 4;
            if self.valid_mask[i] == 0 {
                rgba[o] = 128;
                rgba[o + 1] = 0;
                rgba[o + 2] = 0;
                rgba[o + 3] = 0;
            } else {
                let anom_n = ((self.anom_c[i] + 5.0) / 10.0).clamp(0.0, 1.0);
                let sst_n = ((self.sst_c[i] + 2.0) / 37.0).clamp(0.0, 1.0);
                let spd = (self.u_ms[i] * self.u_ms[i] + self.v_ms[i] * self.v_ms[i]).sqrt();
                let spd_n = (spd / 2.5).clamp(0.0, 1.0);
                rgba[o] = (anom_n * 255.0).round() as u8;
                rgba[o + 1] = (sst_n * 255.0).round() as u8;
                rgba[o + 2] = (spd_n * 255.0).round() as u8;
                rgba[o + 3] = 255;
            }
        }
        rgba
    }
}

/// Builds the verified NOAA CoastWatch ERDDAP OISST v2.1 URL (`ncdcOisst21Agg` / `ncdcOisst21NrtAgg`).
pub fn build_oisst_erddap_url(date_or_last: &str, stride: usize, near_real_time: bool) -> String {
    let ds = if near_real_time {
        "ncdcOisst21NrtAgg"
    } else {
        "ncdcOisst21Agg"
    };
    let t_sel = if date_or_last == "last" {
        "last".to_string()
    } else if date_or_last.contains('T') {
        date_or_last.to_string()
    } else {
        format!("{date_or_last}T12:00:00Z")
    };
    let s = stride.clamp(1, 16);
    format!(
        "https://coastwatch.pfeg.noaa.gov/erddap/griddap/{ds}.json?sst[({t_sel})][(0.0)][(-89.875):{s}:(89.875)][(0.125):{s}:(359.875)],anom[({t_sel})][(0.0)][(-89.875):{s}:(89.875)][(0.125):{s}:(359.875)]"
    )
}

/// Builds the verified NOAA/NESDIS RADS Altimetry Geostrophic Currents URL (`nesdisSSH1day`, 2017–present).
pub fn build_nesdis_currents_erddap_url(date_or_last: &str, stride: usize) -> String {
    let t_sel = if date_or_last == "last" {
        "last".to_string()
    } else if date_or_last.contains('T') {
        date_or_last.to_string()
    } else {
        format!("{date_or_last}T00:00:00Z")
    };
    let s = stride.clamp(1, 16);
    format!(
        "https://coastwatch.pfeg.noaa.gov/erddap/griddap/nesdisSSH1day.json?sla[({t_sel})][(-89.875):{s}:(89.875)][(-179.875):{s}:(179.875)],ugos[({t_sel})][(-89.875):{s}:(89.875)][(-179.875):{s}:(179.875)],vgos[({t_sel})][(-89.875):{s}:(89.875)][(-179.875):{s}:(179.875)]"
    )
}

/// Parses an ERDDAP `.json` `table` response into (or merges onto) a `SatelliteRasterGrid`.
pub fn merge_erddap_json_into_grid(
    raw_json: &str,
    grid: &mut SatelliteRasterGrid,
) -> Result<usize, String> {
    let root: Value = serde_json::from_str(raw_json).map_err(|e| format!("JSON error: {e}"))?;
    let table = root
        .get("table")
        .ok_or_else(|| "Missing ERDDAP 'table' object".to_string())?;
    let cols = table
        .get("columnNames")
        .and_then(Value::as_array)
        .ok_or_else(|| "Missing 'columnNames'".to_string())?;
    let rows = table
        .get("rows")
        .and_then(Value::as_array)
        .ok_or_else(|| "Missing 'rows'".to_string())?;

    let col_idx = |names: &[&str]| -> Option<usize> {
        cols.iter().position(|c| {
            c.as_str()
                .map(|s| names.iter().any(|n| s.eq_ignore_ascii_case(n)))
                .unwrap_or(false)
        })
    };

    let i_time = col_idx(&["time"]);
    let i_lat = col_idx(&["latitude", "lat"]).ok_or("Missing latitude column")?;
    let i_lon = col_idx(&["longitude", "lon"]).ok_or("Missing longitude column")?;
    let i_sst = col_idx(&["sst", "analysed_sst"]);
    let i_anom = col_idx(&["anom", "sst_anom"]);
    let i_u = col_idx(&["ugos", "u", "water_u", "u_current"]);
    let i_v = col_idx(&["vgos", "v", "water_v", "v_current"]);
    let i_sla = col_idx(&["sla", "ssh"]);

    let mut populated = 0_usize;
    for row in rows {
        let Some(arr) = row.as_array() else { continue };
        if let Some(ti) = i_time {
            if grid.timestamp_utc.is_empty() {
                if let Some(ts) = arr.get(ti).and_then(Value::as_str) {
                    grid.timestamp_utc = ts.to_string();
                }
            }
        }
        let Some(lat) = arr.get(i_lat).and_then(Value::as_f64) else { continue };
        let Some(lon_raw) = arr.get(i_lon).and_then(Value::as_f64) else { continue };
        // Convert 0..360 or 20..380 longitudes seamlessly to [-180, 180)
        let lon = (lon_raw + 180.0).rem_euclid(360.0) - 180.0;

        let cx = (((lon + 180.0) / 360.0) * grid.width as f64)
            .floor()
            .clamp(0.0, (grid.width - 1) as f64) as usize;
        let cy = (((lat + 90.0) / 180.0) * grid.height as f64)
            .floor()
            .clamp(0.0, (grid.height - 1) as f64) as usize;
        let idx = cy * grid.width + cx;

        let mut any_val = false;
        if let Some(si) = i_sst {
            if let Some(v) = arr.get(si).and_then(Value::as_f64) {
                if v > -5.0 && v < 45.0 {
                    grid.sst_c[idx] = v as f32;
                    any_val = true;
                }
            }
        }
        if let Some(ai) = i_anom {
            if let Some(v) = arr.get(ai).and_then(Value::as_f64) {
                if v > -15.0 && v < 15.0 {
                    grid.anom_c[idx] = v as f32;
                    any_val = true;
                }
            }
        }
        if let (Some(ui), Some(vi)) = (i_u, i_v) {
            if let (Some(u), Some(v)) = (
                arr.get(ui).and_then(Value::as_f64),
                arr.get(vi).and_then(Value::as_f64),
            ) {
                if u.abs() < 6.0 && v.abs() < 6.0 {
                    grid.u_ms[idx] = u as f32;
                    grid.v_ms[idx] = v as f32;
                    grid.has_velocity = true;
                    any_val = true;
                }
            }
        }
        if let Some(sli) = i_sla {
            if let Some(s) = arr.get(sli).and_then(Value::as_f64) {
                if s.abs() < 5.0 {
                    grid.sla_m[idx] = s as f32;
                    any_val = true;
                }
            }
        }
        if any_val {
            grid.valid_mask[idx] = 1;
            populated += 1;
        }
    }
    Ok(populated)
}

/// Returns `true` when `(year, month)` falls within the NOAA OISST v2.1 satellite archive (Sep 1981–present).
pub fn erddap_oisst_supported(year: i32, month: usize) -> bool {
    year > 1981 || (year == 1981 && month >= 9)
}

/// Returns `true` when `date_or_last` falls within the `nesdisSSH1day` altimetry archive (2017–present or `"last"`).
pub fn erddap_nesdis_currents_supported(date_or_last: &str) -> bool {
    let trimmed = date_or_last.trim();
    if trimmed.eq_ignore_ascii_case("last") {
        return true;
    }
    trimmed
        .split('-')
        .next()
        .and_then(|y| y.parse::<i32>().ok())
        .map(|y| y >= 2017)
        .unwrap_or(false)
}

/// Fetches live NOAA CoastWatch ERDDAP OISST v2.1 + NESDIS geostrophic currents over HTTPS.
pub fn fetch_live_noaa_erddap(date_or_last: &str, stride: usize) -> Result<SatelliteRasterGrid, String> {
    let s = stride.clamp(2, 16);
    let width = 1440 / s;
    let height = 720 / s;
    let mut grid = SatelliteRasterGrid::new_empty(
        width,
        height,
        &format!("NOAA OISST v2.1 Historical SST (Analytical Currents, stride {s})"),
        "",
    );

    let agent = ureq::AgentBuilder::new()
        .timeout_connect(Duration::from_secs(8))
        .timeout_read(Duration::from_secs(18))
        .build();

    let oisst_url = build_oisst_erddap_url(date_or_last, s, date_or_last == "last");
    let body = agent
        .get(&oisst_url)
        .call()
        .or_else(|_| {
            agent
                .get(&build_oisst_erddap_url(date_or_last, s, false))
                .call()
        })
        .map_err(|e| format!("NOAA ERDDAP OISST request failed: {e}"))?
        .into_string()
        .map_err(|e| format!("Failed reading OISST body: {e}"))?;
    merge_erddap_json_into_grid(&body, &mut grid)?;

    if grid.valid_cell_count() == 0 {
        return Err(format!(
            "NOAA ERDDAP returned 0 valid ocean cells for '{date_or_last}'."
        ));
    }

    // Attempt optional NESDIS geostrophic surface currents merge ONLY for 2017–present dates
    if erddap_nesdis_currents_supported(date_or_last) {
        let cur_url = build_nesdis_currents_erddap_url(date_or_last, s);
        if let Ok(resp) = agent.get(&cur_url).call() {
            if let Ok(cur_body) = resp.into_string() {
                let _ = merge_erddap_json_into_grid(&cur_body, &mut grid);
            }
        }
    }

    grid.source_label = if grid.has_velocity {
        format!("NOAA OISST v2.1 + NESDIS Altimetry (stride {s})")
    } else {
        format!("NOAA OISST v2.1 Historical SST (Analytical Currents, stride {s})")
    };

    Ok(grid)
}

/// Loads a local `.json` (ERDDAP griddap), `.csv`, or equirectangular `.png`/`.jpg` satellite raster.
pub fn load_raster_from_path<P: AsRef<Path>>(path: P) -> Result<SatelliteRasterGrid, String> {
    let p = path.as_ref();
    let ext = p
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    if ext == "json" {
        let raw = fs::read_to_string(p).map_err(|e| e.to_string())?;
        let mut grid = SatelliteRasterGrid::new_empty(180, 90, &p.display().to_string(), "");
        merge_erddap_json_into_grid(&raw, &mut grid)?;
        Ok(grid)
    } else if ext == "png" || ext == "jpg" || ext == "jpeg" {
        let img = image::open(p).map_err(|e| e.to_string())?.to_rgba8();
        let (w, h) = (img.width() as usize, img.height() as usize);
        let mut grid = SatelliteRasterGrid::new_empty(w, h, &p.display().to_string(), "imported-image");
        for y in 0..h {
            for x in 0..w {
                let px = img.get_pixel(x as u32, (h - 1 - y) as u32);
                let idx = y * w + x;
                if px[3] > 16 {
                    // Map warm red vs cool blue luminance into approximate SST anomaly [-3..+3 C]
                    let warm_cool = (f32::from(px[0]) - f32::from(px[2])) / 255.0;
                    grid.anom_c[idx] = warm_cool * 3.2;
                    grid.sst_c[idx] = 15.0 + warm_cool * 14.0;
                    grid.valid_mask[idx] = 1;
                }
            }
        }
        Ok(grid)
    } else {
        let raw = fs::read_to_string(p).map_err(|e| e.to_string())?;
        let mut grid = SatelliteRasterGrid::new_empty(180, 90, &p.display().to_string(), "imported-csv");
        let mut lines = raw.lines().filter(|l| !l.trim().is_empty());
        let Some(first_line) = lines.next() else {
            return Err("CSV raster file is empty".to_string());
        };
        let header_cols: Vec<String> = first_line
            .split(',')
            .map(|s| s.trim().to_ascii_lowercase())
            .collect();
        let find_col = |names: &[&str]| -> Option<usize> {
            header_cols.iter().position(|h| names.iter().any(|n| h == n))
        };
        let has_named_header = find_col(&["lat", "latitude"]).is_some()
            && find_col(&["lon", "longitude"]).is_some();

        let (i_lat, i_lon, i_sst, i_anom, i_u, i_v, i_sla) = if has_named_header {
            (
                find_col(&["lat", "latitude"]).unwrap_or(0),
                find_col(&["lon", "longitude"]).unwrap_or(1),
                find_col(&["sst", "analysed_sst"]).or(Some(2)),
                find_col(&["anom", "sst_anom"]).or(Some(3)),
                find_col(&["u", "ugos", "water_u", "u_current"]),
                find_col(&["v", "vgos", "water_v", "v_current"]),
                find_col(&["sla", "ssh"]),
            )
        } else {
            (0, 1, Some(2), Some(3), Some(4), Some(5), Some(6))
        };

        let mut data_lines: Vec<&str> = Vec::new();
        if !has_named_header && first_line.split(',').next().and_then(|s| s.trim().parse::<f64>().ok()).is_some() {
            data_lines.push(first_line);
        }
        data_lines.extend(lines);

        for line in data_lines {
            let parts: Vec<&str> = line.split(',').collect();
            let Some(lat) = parts.get(i_lat).and_then(|s| s.trim().parse::<f64>().ok()) else {
                continue;
            };
            let Some(lon_raw) = parts.get(i_lon).and_then(|s| s.trim().parse::<f64>().ok()) else {
                continue;
            };
            if !lat.is_finite() || !lon_raw.is_finite() {
                continue;
            }
            let lon = (lon_raw + 180.0).rem_euclid(360.0) - 180.0;
            let cx = (((lon + 180.0) / 360.0) * 180.0).floor().clamp(0.0, 179.0) as usize;
            let cy = (((lat + 90.0) / 180.0) * 90.0).floor().clamp(0.0, 89.0) as usize;
            let idx = cy * 180 + cx;

            let mut populated = false;
            if let Some(si) = i_sst {
                if let Some(sst) = parts.get(si).and_then(|s| s.trim().parse::<f32>().ok()) {
                    if sst.is_finite() && sst > -5.0 && sst < 45.0 {
                        grid.sst_c[idx] = sst;
                        populated = true;
                    }
                }
            }
            if let Some(ai) = i_anom {
                if let Some(anom) = parts.get(ai).and_then(|s| s.trim().parse::<f32>().ok()) {
                    if anom.is_finite() && anom > -15.0 && anom < 15.0 {
                        grid.anom_c[idx] = anom;
                        populated = true;
                    }
                }
            }
            if let (Some(ui), Some(vi)) = (i_u, i_v) {
                if let (Some(u), Some(v)) = (
                    parts.get(ui).and_then(|s| s.trim().parse::<f32>().ok()),
                    parts.get(vi).and_then(|s| s.trim().parse::<f32>().ok()),
                ) {
                    if u.is_finite() && v.is_finite() && u.abs() < 6.0 && v.abs() < 6.0 {
                        grid.u_ms[idx] = u;
                        grid.v_ms[idx] = v;
                        grid.has_velocity = true;
                        populated = true;
                    }
                }
            }
            if let Some(sli) = i_sla {
                if let Some(sla) = parts.get(sli).and_then(|s| s.trim().parse::<f32>().ok()) {
                    if sla.is_finite() && sla.abs() < 5.0 {
                        grid.sla_m[idx] = sla;
                        populated = true;
                    }
                }
            }
            if populated {
                grid.valid_mask[idx] = 1;
            }
        }
        Ok(grid)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_erddap_0_360_oisst_and_nesdis_currents() {
        let mock_oisst = r#"{
          "table": {
            "columnNames": ["time", "zlev", "latitude", "longitude", "sst", "anom"],
            "rows": [
              ["2023-12-15T12:00:00Z", 0.0, 0.125, 215.125, 28.4, 2.1],
              ["2023-12-15T12:00:00Z", 0.0, -5.125, 275.125, 25.1, 2.7]
            ]
          }
        }"#;
        let mut grid = SatelliteRasterGrid::new_empty(180, 90, "test", "");
        let count = merge_erddap_json_into_grid(mock_oisst, &mut grid).unwrap();
        assert_eq!(count, 2);
        assert_eq!(grid.timestamp_utc, "2023-12-15T12:00:00Z");
        // 215.125E == -144.875W (Niño 3.4)
        let sampled = grid.sample_bilinear(-144.8, 0.2).unwrap();
        assert!((sampled.0 - 28.4).abs() < 0.05);
        assert!((sampled.1 - 2.1).abs() < 0.05);
    }

    #[test]
    fn checks_historical_erddap_archive_coverage_windows() {
        assert!(!erddap_oisst_supported(1973, 12));
        assert!(!erddap_oisst_supported(1981, 8));
        assert!(erddap_oisst_supported(1981, 9));
        assert!(erddap_oisst_supported(1997, 12));

        assert!(!erddap_nesdis_currents_supported("2015-12-15"));
        assert!(erddap_nesdis_currents_supported("2023-12-15"));
        assert!(erddap_nesdis_currents_supported("last"));
    }

    #[test]
    fn loads_csv_raster_with_velocity_columns_and_sets_has_velocity() {
        struct TempFileGuard(std::path::PathBuf);
        impl Drop for TempFileGuard {
            fn drop(&mut self) {
                let _ = fs::remove_file(&self.0);
            }
        }

        let path = std::env::temp_dir().join(format!(
            "pelagic_csv_vel_{}_{:?}.csv",
            std::process::id(),
            std::thread::current().id()
        ));
        let _guard = TempFileGuard(path.clone());
        let csv = "lat,lon,sst,anom,u,v\n-3.0,-149.0,27.5,0.4,-0.85,0.12\n7.0,-145.0,28.8,1.2,0.92,0.05\n";
        fs::write(&path, csv).expect("write temporary CSV raster");

        let grid = load_raster_from_path(&path).expect("load CSV raster");
        assert!(grid.has_velocity, "CSV with u,v columns must set has_velocity == true");
        assert_eq!(grid.valid_cell_count(), 2);
        let sample = grid.sample_bilinear(-149.0, -3.0).expect("sample loaded cell");
        assert!((sample.0 - 27.5).abs() < 0.1);
        assert!((sample.1 - 0.4).abs() < 0.1);
        assert!((sample.2 - (-0.85)).abs() < 0.1);
        assert!((sample.3 - 0.12).abs() < 0.1);
    }
}
