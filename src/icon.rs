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

//! Embedded procedural RGBA window icon (`ViewportBuilder::with_icon`) matching `assets/pelagic_globe.ico`.

use eframe::egui;

fn blend_px(buf: &mut [u8], idx: usize, r: f32, g: f32, b: f32, a: f32) {
    let alpha = a.clamp(0.0, 1.0);
    if alpha <= 0.001 {
        return;
    }
    let da = f32::from(buf[idx + 3]) / 255.0;
    let out_a = alpha + da * (1.0 - alpha);
    if out_a <= 1e-5 {
        return;
    }
    let br = f32::from(buf[idx]);
    let bg = f32::from(buf[idx + 1]);
    let bb = f32::from(buf[idx + 2]);
    buf[idx] = ((r * alpha + br * da * (1.0 - alpha)) / out_a)
        .round()
        .clamp(0.0, 255.0) as u8;
    buf[idx + 1] = ((g * alpha + bg * da * (1.0 - alpha)) / out_a)
        .round()
        .clamp(0.0, 255.0) as u8;
    buf[idx + 2] = ((b * alpha + bb * da * (1.0 - alpha)) / out_a)
        .round()
        .clamp(0.0, 255.0) as u8;
    buf[idx + 3] = (out_a * 255.0).round().clamp(0.0, 255.0) as u8;
}

fn sd_rounded_box(px: f32, py: f32, bx: f32, by: f32, r: f32) -> f32 {
    let qx = px.abs() - bx + r;
    let qy = py.abs() - by + r;
    qx.max(0.0).hypot(qy.max(0.0)) + qx.max(qy).min(0.0) - r
}

/// Renders an anti-aliased `size x size` RGBA PelagicGlobe badge matching `assets/pelagic_globe.ico`.
pub fn render_pelagic_globe_rgba(size: usize) -> Vec<u8> {
    let mut buf = vec![0u8; size * size * 4];
    let aa = 2.2 / (size as f32).max(1.0);
    let border_w = 0.038;

    for py in 0..size {
        let y = ((py as f32 + 0.5) / size as f32) * 2.0 - 1.0;
        for px in 0..size {
            let x = ((px as f32 + 0.5) / size as f32) * 2.0 - 1.0;
            let idx = (py * size + px) * 4;
            let box_d = sd_rounded_box(x, y, 0.88, 0.88, 0.24);
            if box_d < aa {
                let bg_a = (0.5 - box_d / aa).clamp(0.0, 1.0);
                let grad = (y + 1.0) * 0.5;
                blend_px(
                    &mut buf,
                    idx,
                    8.0 + 6.0 * (1.0 - grad),
                    18.0 + 18.0 * (1.0 - grad),
                    36.0 + 26.0 * grad,
                    bg_a,
                );
                if box_d.abs() < border_w {
                    blend_px(
                        &mut buf,
                        idx,
                        56.0,
                        189.0,
                        248.0,
                        (1.0 - box_d.abs() / border_w) * 0.85 * bg_a,
                    );
                }
            }

            let dia_d = (x.abs() + y.abs() - 0.84).abs();
            if dia_d < 0.032 + aa {
                let da = (0.5 + (0.022 - dia_d) / aa).clamp(0.0, 1.0);
                blend_px(&mut buf, idx, 250.0, 204.0, 21.0, da * 0.82);
            }

            let r = x.hypot(y);
            let globe_r = 0.58;
            if r < globe_r + aa {
                let ga = (0.5 + (globe_r - r) / aa).clamp(0.0, 1.0);
                let nx = x / globe_r;
                let ny = y / globe_r;
                let nz = (1.0 - nx * nx - ny * ny).max(0.0).sqrt();
                let shade = (0.45 + 0.55 * (nx * -0.4 + ny * -0.45 + nz * 0.78)).max(0.25);
                blend_px(
                    &mut buf,
                    idx,
                    14.0 + 36.0 * shade,
                    78.0 + 115.0 * shade,
                    145.0 + 95.0 * shade,
                    ga,
                );

                let eq_wave = 0.14 * (nx * 3.8 + 0.4).sin() - 0.03;
                let d_eq = (ny - eq_wave).abs();
                if d_eq < 0.055 + aa {
                    let ea = (0.5 + (0.045 - d_eq) / aa).clamp(0.0, 1.0) * ga;
                    let warm_t = ((nx + 0.8) * 0.6).clamp(0.0, 1.0);
                    blend_px(
                        &mut buf,
                        idx,
                        251.0 * warm_t + 56.0 * (1.0 - warm_t),
                        146.0 * warm_t + 220.0 * (1.0 - warm_t),
                        60.0 * warm_t + 248.0 * (1.0 - warm_t),
                        ea * 0.92,
                    );
                }

                let gyre_n = (ny - (-0.42 + 0.14 * nx * nx)).abs();
                if gyre_n < 0.038 + aa {
                    let ca = (0.5 + (0.03 - gyre_n) / aa).clamp(0.0, 1.0) * ga;
                    blend_px(&mut buf, idx, 125.0, 211.0, 252.0, ca * 0.88);
                }

                let mer_d = (nx.abs() - 0.45 * (1.0 - ny * ny).max(0.0).sqrt()).abs();
                if mer_d < 0.025 {
                    blend_px(&mut buf, idx, 224.0, 242.0, 254.0, 0.45 * ga);
                }

                let rim_d = (r - globe_r).abs();
                if rim_d < 0.028 {
                    blend_px(
                        &mut buf,
                        idx,
                        125.0,
                        211.0,
                        252.0,
                        (1.0 - rim_d / 0.028) * 0.75 * ga,
                    );
                }
            }
        }
    }

    buf
}

/// Returns the 64x64 RGBA application icon for `ViewportBuilder::with_icon`.
pub fn app_egui_icon() -> egui::IconData {
    egui::IconData {
        rgba: render_pelagic_globe_rgba(64),
        width: 64,
        height: 64,
    }
}
