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

//! Exact conformal Spilhaus World Ocean Map in a Square (`ESRI:54099` / `PROJ +proj=spilhaus`).
//!
//! Formulated as an oblique aspect of Oscar S. Adams's World in a Square II (1929)
//! with projection center `(66.94970198°E, 49.56371678°S)` and azimuth `40.17823482°`
//! (Šavrič, Burrows & Kennedy; PROJ 9.6+; Fukushima 2009 / Bulirsch 1965 elliptic algorithms).

use std::f64::consts::PI;

/// Complete elliptic integral of the first kind `K(1/2) = F(π/2 | m = 0.5)`.
pub const K05: f64 = 1.854_074_677_301_372;

/// Default Spilhaus center longitude in degrees (`66.94970198°E`, Southern Indian Ocean).
pub const SPILHAUS_LON0_DEG: f64 = 66.949_701_98;
/// Default Spilhaus center latitude in degrees (`49.56371678°S`).
pub const SPILHAUS_LAT0_DEG: f64 = -49.563_716_78;
/// Default Spilhaus azimuth rotation in degrees (`40.17823482°`).
pub const SPILHAUS_AZIMUTH_DEG: f64 = 40.178_234_82;

const ARCSN_COEFFS: [f64; 9] = [
    0.25,
    0.11875,
    0.0703125,
    0.047_037_760_416_666_664,
    0.034_057_617_187_5,
    0.026_056_142_953_725_96,
    0.020_749_664_306_640_624,
    0.017_029_054_024_640_25,
    0.014_304_488_897_323_608,
];

const AGM_ARITH: [f64; 4] = [
    1.0,
    0.853_553_390_593_273_7,
    0.847_224_902_923_494_2,
    0.847_213_084_835_192_9,
];

const AGM_GEOM: [f64; 4] = [
    0.707_106_781_186_547_6,
    0.840_896_415_253_714_6,
    0.847_201_266_746_891_6,
    0.847_213_084_752_765_4,
];

const AGM_LIMIT: f64 = 0.847_213_084_793_979_2;

/// Applies the forward D3/PROJ oblique Euler rotation `[-66.94970198°, +49.56371678°, +40.17823482°]`.
pub fn oblique_rotate_forward(lon_rad: f64, lat_rad: f64) -> (f64, f64) {
    let dl = (-SPILHAUS_LON0_DEG).to_radians();
    let dp = (-SPILHAUS_LAT0_DEG).to_radians();
    let dg = SPILHAUS_AZIMUTH_DEG.to_radians();

    let l = lon_rad + dl;
    let cos_lat = lat_rad.cos();
    let x = cos_lat * l.cos();
    let y = cos_lat * l.sin();
    let z = lat_rad.sin();

    let (sin_p, cos_p) = dp.sin_cos();
    let (sin_g, cos_g) = dg.sin_cos();

    let k = z * cos_p + x * sin_p;
    let x2 = x * cos_p - z * sin_p;
    let y2 = y * cos_g - k * sin_g;
    let z2 = (k * cos_g + y * sin_g).clamp(-1.0, 1.0);

    (y2.atan2(x2), z2.asin())
}

/// Applies the exact inverse of `oblique_rotate_forward`.
pub fn oblique_rotate_inverse(lam_rot: f64, phi_rot: f64) -> (f64, f64) {
    let dl = (-SPILHAUS_LON0_DEG).to_radians();
    let dp = (-SPILHAUS_LAT0_DEG).to_radians();
    let dg = SPILHAUS_AZIMUTH_DEG.to_radians();

    let cos_lat = phi_rot.cos();
    let x2 = cos_lat * lam_rot.cos();
    let y2 = cos_lat * lam_rot.sin();
    let z2 = phi_rot.sin();

    let (sin_p, cos_p) = dp.sin_cos();
    let (sin_g, cos_g) = dg.sin_cos();

    let k = z2 * cos_g - y2 * sin_g;
    let y = y2 * cos_g + z2 * sin_g;
    let x = x2 * cos_p + k * sin_p;
    let z = (k * cos_p - x2 * sin_p).clamp(-1.0, 1.0);

    let mut lon = y.atan2(x) - dl;
    lon = (lon + PI).rem_euclid(2.0 * PI) - PI;
    (lon, z.asin())
}

fn arcsn_series(mut sin2: f64, mut cos2: f64) -> f64 {
    let mut scale = 1.0;
    while sin2 >= 0.03768 {
        let denom = (1.0 + cos2.max(0.0).sqrt()) * (1.0 + (1.0 - 0.5 * sin2).max(0.0).sqrt());
        sin2 /= denom;
        cos2 = (1.0 - sin2).max(0.0);
        scale *= 2.0;
    }
    let s2 = sin2.max(0.0);
    let u = &ARCSN_COEFFS;
    let poly = 1.0
        + s2 * (u[0]
            + s2 * (u[1]
                + s2 * (u[2]
                    + s2 * (u[3]
                        + s2 * (u[4] + s2 * (u[5] + s2 * (u[6] + s2 * (u[7] + s2 * u[8]))))))));
    scale * s2.sqrt() * poly
}

/// Incomplete elliptic integral of the first kind `F(φ | 1/2)` given `sin²φ` and `cos²φ`.
pub fn elliptic_f05(sin2_phi: f64, cos2_phi: f64) -> f64 {
    if sin2_phi > 0.9 {
        let denom = 1.0 + cos2_phi;
        K05 - arcsn_series(2.0 * cos2_phi / denom, (1.0 - cos2_phi) / denom)
    } else {
        arcsn_series(sin2_phi, cos2_phi)
    }
}

/// Jacobi elliptic functions `(sn, cn, dn)` for parameter `m = 1/2` via descending AGM.
pub fn elliptic_j05(u: f64) -> (f64, f64, f64) {
    let mut dn = 1.0;
    let mut c = AGM_LIMIT;
    let v = u * c;
    let (mut sn, mut cn) = v.sin_cos();

    if sn.abs() > 1.0e-15 {
        let mut a = cn / sn;
        c *= a;
        for i in (0..4).rev() {
            let b = AGM_ARITH[i];
            a *= c;
            c *= dn;
            dn = (AGM_GEOM[i] + a) / (b + a);
            a = c / b;
        }
        let inv = 1.0 / (c * c + 1.0).sqrt();
        sn = if sn >= 0.0 { inv } else { -inv };
        cn = c * sn;
    }
    (sn, cn, dn)
}

fn conformal_square_forward(cosa: f64, cosb: f64, cosc: f64) -> (f64, f64) {
    let p = cosa * cosb;
    let cos2c = cosc * cosc;
    let q = ((cosb * cosb + cos2c) * (cosa * cosa + cos2c)).max(0.0).sqrt();
    let sum = cosa + cosb;
    let diff = cosb - cosa;

    let (cos2m, cos2n) = if p >= 0.0 {
        let cn = q + p;
        let cm = if cn > 0.0 { cos2c / cn } else { 0.0 };
        (cm, cn)
    } else {
        let cm = q - p;
        let cn = if cm > 0.0 { cos2c / cm } else { 0.0 };
        (cm, cn)
    };

    let sin2m = if 1.0 + p + q > 0.0 {
        (sum * sum) / (1.0 + p + q)
    } else {
        0.0
    };
    let sin2n = if 1.0 - p + q > 0.0 {
        (diff * diff) / (1.0 - p + q)
    } else {
        0.0
    };

    let x = elliptic_f05(sin2m, cos2m);
    let y = elliptic_f05(sin2n, cos2n);
    (
        if sum < 0.0 { -x } else { x },
        if diff < 0.0 { -y } else { y },
    )
}

fn conformal_square_inverse(x: f64, y: f64) -> (f64, f64, f64) {
    let (snu, cnu, dnu) = elliptic_j05(K05 - 0.5 * (x + y));
    let (snv, cnv, dnv) = elliptic_j05(0.5 * (x - y));
    let denom = cnv * cnv + 0.5 * snu * snu * snv * snv;
    if denom.abs() < 1.0e-18 {
        return (0.0, 0.0, -1.0);
    }
    let re = cnu * cnv / denom;
    let im = -snu * dnu * snv * dnv / denom;
    let d2 = re * re + im * im;
    let s = 2.0 / (1.0 + d2);
    (-s * im, s * re, (1.0 - d2) / (1.0 + d2))
}

/// Projects geographic `(lon_deg, lat_deg)` into normalized Spilhaus square coordinates `(u, v) ∈ [-1, 1]²`.
pub fn lonlat_to_spilhaus_norm(lon_deg: f64, lat_deg: f64) -> (f64, f64) {
    let (x, y) = lonlat_to_spilhaus_raw(lon_deg, lat_deg);
    ((x / K05).clamp(-1.0, 1.0), (y / K05).clamp(-1.0, 1.0))
}

/// Projects geographic `(lon_deg, lat_deg)` into raw conformal Spilhaus square coordinates `(x, y) ∈ [-K05, K05]²`.
pub fn lonlat_to_spilhaus_raw(lon_deg: f64, lat_deg: f64) -> (f64, f64) {
    let (lam, phi) = oblique_rotate_forward(lon_deg.to_radians(), lat_deg.to_radians());
    let v = (0.5 * phi).tan().clamp(-1.0, 1.0);
    let cos_phi1 = (1.0 - v * v).max(0.0).sqrt();
    conformal_square_forward(cos_phi1 * (0.5 * lam).sin(), v, cos_phi1 * (0.5 * lam).cos())
}

/// Inverse-projects normalized Spilhaus square coordinates `(u, v) ∈ [-1, 1]²` back to `(lon_deg, lat_deg)`.
pub fn spilhaus_norm_to_lonlat(u_norm: f64, v_norm: f64) -> (f64, f64) {
    spilhaus_raw_to_lonlat(u_norm.clamp(-1.0, 1.0) * K05, v_norm.clamp(-1.0, 1.0) * K05)
}

/// Inverse-projects raw Spilhaus square coordinates `(x, y) ∈ [-K05, K05]²` back to `(lon_deg, lat_deg)`.
pub fn spilhaus_raw_to_lonlat(x: f64, y: f64) -> (f64, f64) {
    spilhaus_raw_to_lonlat_unclamped(x.clamp(-K05, K05), y.clamp(-K05, K05))
}

/// Inverse-projects unclamped raw Spilhaus plane coordinates `(x, y)` (for periodic `3×3` mosaics)
/// back to `(lon_deg, lat_deg)`.
pub fn spilhaus_raw_to_lonlat_unclamped(x: f64, y: f64) -> (f64, f64) {
    let (u, v, w) = conformal_square_inverse(x, y);
    let lam_rot = 2.0 * u.atan2(w);
    let phi_rot = 2.0 * v.atan();
    let (lon_rad, lat_rad) = oblique_rotate_inverse(lam_rot, phi_rot);
    (lon_rad.to_degrees(), lat_rad.to_degrees())
}

/// Inverse-projects unclamped normalized Spilhaus plane coordinates `(u_norm, v_norm)` (e.g. `[-3, 3]²`)
/// back to `(lon_deg, lat_deg)`.
pub fn spilhaus_norm_to_lonlat_unclamped(u_norm: f64, v_norm: f64) -> (f64, f64) {
    spilhaus_raw_to_lonlat_unclamped(u_norm * K05, v_norm * K05)
}

/// Rotates geographic surface velocity `(u_east_ms, v_north_ms)` at `(lon_deg, lat_deg)` into
/// local conformal Spilhaus square coordinates `(vu_spilhaus, vv_spilhaus)` using the
/// forward projection Jacobian `∂(u, v) / ∂(East, North)`.
pub fn spilhaus_local_jacobian_velocity(
    lon_deg: f64,
    lat_deg: f64,
    u_east_ms: f64,
    v_north_ms: f64,
) -> (f64, f64) {
    let speed = u_east_ms.hypot(v_north_ms);
    if speed < 1.0e-6 {
        return (0.0, 0.0);
    }
    let (u0, v0) = lonlat_to_spilhaus_norm(lon_deg, lat_deg);
    let eps = 0.25_f64;
    let cos_lat = lat_deg.to_radians().cos().max(0.15);
    let dlon = eps / cos_lat;

    let (ue_p, ve_p) = lonlat_to_spilhaus_norm(lon_deg + dlon, lat_deg);
    let (ue_m, ve_m) = lonlat_to_spilhaus_norm(lon_deg - dlon, lat_deg);
    let de_p = ((ue_p - u0).hypot(ve_p - v0), (ue_p - u0, ve_p - v0));
    let de_m = ((u0 - ue_m).hypot(v0 - ve_m), (u0 - ue_m, v0 - ve_m));
    let (du_e, dv_e) = if de_p.0 <= de_m.0 { de_p.1 } else { de_m.1 };

    let lat_p = (lat_deg + eps).clamp(-89.5, 89.5);
    let lat_m = (lat_deg - eps).clamp(-89.5, 89.5);
    let (un_p, vn_p) = lonlat_to_spilhaus_norm(lon_deg, lat_p);
    let (un_m, vn_m) = lonlat_to_spilhaus_norm(lon_deg, lat_m);
    let dn_p = ((un_p - u0).hypot(vn_p - v0), (un_p - u0, vn_p - v0));
    let dn_m = ((u0 - un_m).hypot(v0 - vn_m), (u0 - un_m, v0 - vn_m));
    let (du_n, dv_n) = if dn_p.0 <= dn_m.0 { dn_p.1 } else { dn_m.1 };

    let mag_e = du_e.hypot(dv_e);
    let mag_n = du_n.hypot(dv_n);
    if mag_e < 1.0e-9 || mag_n < 1.0e-9 {
        return (u_east_ms, v_north_ms);
    }
    let e_hat = (du_e / mag_e, dv_e / mag_e);
    let n_hat = (du_n / mag_n, dv_n / mag_n);
    (
        u_east_ms * e_hat.0 + v_north_ms * n_hat.0,
        u_east_ms * e_hat.1 + v_north_ms * n_hat.1,
    )
}

/// Samples the Spilhaus boundary cut curve in `(lon_deg, lat_deg)` around the square perimeter.
pub fn sample_boundary_seam(points_per_edge: usize) -> Vec<(f64, f64)> {
    let n = points_per_edge.max(8);
    let mut pts = Vec::with_capacity(n * 4 + 1);
    let eps = 0.9985;
    for i in 0..n {
        let t = -eps + 2.0 * eps * (i as f64 / n as f64);
        pts.push(spilhaus_norm_to_lonlat(t, eps));
    }
    for i in 0..n {
        let t = eps - 2.0 * eps * (i as f64 / n as f64);
        pts.push(spilhaus_norm_to_lonlat(eps, t));
    }
    for i in 0..n {
        let t = eps - 2.0 * eps * (i as f64 / n as f64);
        pts.push(spilhaus_norm_to_lonlat(t, -eps));
    }
    for i in 0..=n {
        let t = -eps + 2.0 * eps * (i as f64 / n as f64);
        pts.push(spilhaus_norm_to_lonlat(-eps, t));
    }
    pts
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kerguelen_origin_maps_to_center() {
        let (u, v) = lonlat_to_spilhaus_norm(SPILHAUS_LON0_DEG, SPILHAUS_LAT0_DEG);
        assert!(u.abs() < 1.0e-10, "u = {u}");
        assert!(v.abs() < 1.0e-10, "v = {v}");
    }

    #[test]
    fn bering_strait_lies_on_spilhaus_outer_seam() {
        let (_u, v) = lonlat_to_spilhaus_norm(-169.0, 65.3);
        assert!(
            (v - 1.0).abs() < 2.0e-3,
            "Expected Bering Strait on top edge (+1), got v = {v}"
        );
    }

    #[test]
    fn round_trip_accuracy_across_world_oceans() {
        let coords = [
            (SPILHAUS_LON0_DEG, SPILHAUS_LAT0_DEG),
            (-145.0, 0.0),   // Nino 3.4 Equatorial Pacific
            (-125.0, 44.6),  // Offshore Oregon / California Current
            (142.0, 35.0),   // Kuroshio Extension
            (-70.0, 36.0),   // Gulf Stream
            (25.0, -37.0),   // Agulhas Retroflection
            (-65.0, -58.0),  // Drake Passage ACC
            (118.0, -8.5),   // Lombok Strait ITF
        ];
        for (lon, lat) in coords {
            let (x, y) = lonlat_to_spilhaus_raw(lon, lat);
            let (lon2, lat2) = spilhaus_raw_to_lonlat(x, y);
            let dlon = ((lon - lon2 + 180.0).rem_euclid(360.0) - 180.0).abs();
            let dlat = (lat - lat2).abs();
            assert!(dlon < 1.0e-7 && dlat < 1.0e-7, "Mismatch for ({lon}, {lat}) -> ({lon2}, {lat2})");
        }
    }

    #[test]
    fn jacobian_rotates_velocity_vectors_and_unclamped_mosaic_is_continuous() {
        let (vu, vv) = spilhaus_local_jacobian_velocity(-145.0, 0.0, 1.0, 0.0);
        // In Equatorial Pacific (Nino 3.4), Eastward flow rotates obliquely in Spilhaus (u, v) space
        assert!(vv.abs() > 0.3, "Expected Spilhaus Jacobian to rotate Eastward flow in Equatorial Pacific (got vu={vu:.3}, vv={vv:.3})");
        assert!((vu.hypot(vv) - 1.0).abs() < 0.05, "Jacobian preserves speed magnitude");

        let (lon_in, lat_in) = spilhaus_norm_to_lonlat_unclamped(0.998, 0.25);
        let (lon_out, lat_out) = spilhaus_norm_to_lonlat_unclamped(1.002, 0.25);
        let dlon = ((lon_in - lon_out + 180.0).rem_euclid(360.0) - 180.0).abs();
        let dlat = (lat_in - lat_out).abs();
        assert!(dlon < 1.0 && dlat < 1.0, "Unclamped mosaic edge transition must be continuous");
    }
}
