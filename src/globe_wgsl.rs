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

//! `wgpu` 24.0 WGSL shader source and conformal Spilhaus-parameterized 3D globe mesh builder.

use bytemuck::{Pod, Zeroable};
use glam::Vec3;

use crate::currents::sample_velocity_ms;
use crate::enso::{equatorial_20c_isotherm_depth_m, spatial_sst_anomaly_c};
use crate::spilhaus::spilhaus_norm_to_lonlat;

#[repr(C)]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
pub struct GlobeVertex {
    pub sphere_pos: [f32; 3],
    pub spilhaus_uv: [f32; 2],
    pub lon_lat_deg: [f32; 2],
    pub sst_anom_and_speed: [f32; 2],
    pub thermocline_m: f32,
}

/// Builds a tear-free 3D globe / Spilhaus morph mesh parameterized directly in Spilhaus `(u, v) ∈ [-1, 1]²`.
pub fn build_spilhaus_globe_mesh(subdivisions: usize, oni: f64) -> (Vec<GlobeVertex>, Vec<u32>) {
    let n = subdivisions.max(32);
    let cols = n + 1;
    let mut vertices = Vec::with_capacity(cols * cols);
    let mut indices = Vec::with_capacity(n * n * 6);

    for row in 0..=n {
        let v_norm = -0.999 + 1.998 * (row as f64 / n as f64);
        for col in 0..=n {
            let u_norm = -0.999 + 1.998 * (col as f64 / n as f64);
            let (lon_deg, lat_deg) = spilhaus_norm_to_lonlat(u_norm, v_norm);
            let lon_r = (lon_deg as f32).to_radians();
            let lat_r = (lat_deg as f32).to_radians();
            let pos = Vec3::new(
                lat_r.cos() * lon_r.cos(),
                lat_r.sin(),
                -lat_r.cos() * lon_r.sin(),
            );
            let anom = spatial_sst_anomaly_c(lon_deg, lat_deg, oni) as f32;
            let (u_ms, v_ms) = sample_velocity_ms(lon_deg, lat_deg, oni);
            let speed = (u_ms * u_ms + v_ms * v_ms).sqrt() as f32;
            let thermo = equatorial_20c_isotherm_depth_m(lon_deg.rem_euclid(360.0), oni) as f32;

            vertices.push(GlobeVertex {
                sphere_pos: pos.to_array(),
                spilhaus_uv: [u_norm as f32, v_norm as f32],
                lon_lat_deg: [lon_deg as f32, lat_deg as f32],
                sst_anom_and_speed: [anom, speed],
                thermocline_m: thermo,
            });
        }
    }

    for row in 0..n {
        for col in 0..n {
            let i0 = (row * cols + col) as u32;
            let i1 = i0 + 1;
            let i2 = ((row + 1) * cols + col) as u32;
            let i3 = i2 + 1;
            indices.extend_from_slice(&[i0, i2, i1, i1, i2, i3]);
        }
    }

    (vertices, indices)
}

pub const GLOBE_SHADER_WGSL: &str = r#"
struct Uniforms {
    view_proj: mat4x4<f32>,
    camera_pos: vec4<f32>,
    morph_overlay_layer_time: vec4<f32>, // x = spilhaus_morph, y = show_spilhaus_overlay, z = data_layer, w = time
};

@group(0) @binding(0)
var<uniform> u_scene: Uniforms;

struct VertexInput {
    @location(0) sphere_pos: vec3<f32>,
    @location(1) spilhaus_uv: vec2<f32>,
    @location(2) lon_lat_deg: vec2<f32>,
    @location(3) sst_anom_and_speed: vec2<f32>,
    @location(4) thermocline_m: f32,
};

struct VertexOutput {
    @builtin(position) clip_pos: vec4<f32>,
    @location(0) world_pos: vec3<f32>,
    @location(1) normal: vec3<f32>,
    @location(2) spilhaus_uv: vec2<f32>,
    @location(3) sst_speed: vec2<f32>,
    @location(4) lat_deg: f32,
};

@vertex
fn vs_main(in: VertexInput) -> VertexOutput {
    let morph = clamp(u_scene.morph_overlay_layer_time.x, 0.0, 1.0);
    let plane_pos = vec3<f32>(in.spilhaus_uv.x * 1.45, in.spilhaus_uv.y * 1.45, 0.0);
    let pos = mix(in.sphere_pos, plane_pos, morph);
    let nrm = normalize(mix(in.sphere_pos, vec3<f32>(0.0, 0.0, 1.0), morph));

    var out: VertexOutput;
    out.clip_pos = u_scene.view_proj * vec4<f32>(pos, 1.0);
    out.world_pos = pos;
    out.normal = nrm;
    out.spilhaus_uv = in.spilhaus_uv;
    out.sst_speed = in.sst_anom_and_speed;
    out.lat_deg = in.lon_lat_deg.y;
    return out;
}

@fragment
fn fs_main(in: VertexOutput) -> @location(0) vec4<f32> {
    let light_dir = normalize(vec3<f32>(0.65, 0.45, 0.62));
    let view_dir = normalize(u_scene.camera_pos.xyz - in.world_pos);
    let ndotl = max(dot(in.normal, light_dir), 0.14);
    let half_vec = normalize(light_dir + view_dir);
    let spec = pow(max(dot(in.normal, half_vec), 0.0), 42.0) * 0.28;
    let fresnel = pow(1.0 - max(dot(in.normal, view_dir), 0.0), 3.4) * 0.24;

    var base = vec3<f32>(0.04, 0.16, 0.30);
    let anom = in.sst_speed.x;
    if (anom >= 0.0) {
        base = mix(base, vec3<f32>(0.96, 0.34, 0.16), clamp(anom / 2.8, 0.0, 1.0));
    } else {
        base = mix(base, vec3<f32>(0.12, 0.68, 0.96), clamp(-anom / 2.2, 0.0, 1.0));
    }

    var color = base * (0.35 + 0.65 * ndotl) + vec3<f32>(spec + fresnel * 0.55);

    // Optional Spilhaus conformal square isolines & outer continental seam highlight
    if (u_scene.morph_overlay_layer_time.y > 0.5) {
        let gx = abs(fract(in.spilhaus_uv.x * 4.0 + 0.5) - 0.5);
        let gy = abs(fract(in.spilhaus_uv.y * 4.0 + 0.5) - 0.5);
        let grid = 1.0 - smoothstep(0.0, 0.018, min(gx, gy));
        let edge = smoothstep(0.965, 0.998, max(abs(in.spilhaus_uv.x), abs(in.spilhaus_uv.y)));
        color = mix(color, vec3<f32>(0.38, 0.88, 0.95), grid * 0.28);
        color = mix(color, vec3<f32>(0.98, 0.78, 0.26), edge * 0.82);
    }

    return vec4<f32>(color, 1.0);
}
"#;
