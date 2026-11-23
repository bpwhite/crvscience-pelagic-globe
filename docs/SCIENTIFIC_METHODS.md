# PelagicGlobe Scientific Methods & Mathematical Reference

*Version 1.1.0 — Copyright CRV Science Network — Licensed under `GPL-3.0-or-later`.*

## 1. Conformal Spilhaus World Ocean Projection (`ESRI:54099`)

PelagicGlobe implements the oblique Adams World in a Square II formulation (`src/spilhaus.rs`) parameterized by:

```
lon_0 = +66.94970198 deg    (Southern Indian Ocean near Kerguelen)
lat_0 = -49.56371678 deg
alpha = +40.17823482 deg    (oblique azimuth rotation)
K(1/2) = 1.854074677301372  (complete elliptic integral of the first kind, m = 0.5)
```

### Forward & Inverse Spherical Rotation

Geographic coordinates `(lambda, phi)` are rotated onto the oblique sphere `(lambda_rot, phi_rot)` via three-axis Euler rotation `[-lon_0, -lat_0, +alpha]` in Cartesian `(x, y, z) = (cos phi cos l, cos phi sin l, sin phi)` where `l = lambda - lon_0`:

```
k  = z cos(-lat_0) + x sin(-lat_0)
x2 = x cos(-lat_0) - z sin(-lat_0)
y2 = y cos(alpha)  - k sin(alpha)
z2 = k cos(alpha)  + y sin(alpha)
lambda_rot = atan2(y2, x2),    phi_rot = asin(clamp(z2, -1, 1))
```

### Elliptic Conformal Square Mapping & Periodic `3×3` Mosaic

Using `v = tan(phi_rot / 2)`, `cos_phi1 = sqrt(1 - v^2)`, `a = cos_phi1 sin(lambda_rot / 2)`, `b = v`, `c = cos_phi1 cos(lambda_rot / 2)`, the incomplete elliptic integral `F(phi | 1/2)` maps the oblique sphere onto `(x, y) in [-K(1/2), +K(1/2)]^2`, normalized as `(u, v) = (x / K(1/2), y / K(1/2))`.

For the inverse (and unclamped periodic `3×3` mosaic over `(u, v) in [-3, +3]^2`), Bulirsch descending arithmetic-geometric mean (AGM) evaluates the Jacobi elliptic functions `(sn, cn, dn)` at complex argument `z = [K(1/2) - (x + y)/2] + i [(x - y)/2]`:

```
cn(z | 1/2) = (cn_u cn_v - i sn_u dn_u sn_v dn_v) / (cn_v^2 + 0.5 sn_u^2 sn_v^2)
```

which projects inversely from the Riemann sphere back to `(lambda_rot, phi_rot)` and geographic `(lon, lat)` with round-trip error `< 10^-7 deg`.

### Local Projection Jacobian for Surface Velocity Vectors

Because the Spilhaus map rotates geographic North and East differently across the square, physical East-North velocity vectors `(u_east, v_north)` in `export_spilhaus_svg` are rotated via the local conformal Jacobian (`spilhaus_local_jacobian_velocity`):

```
dlon = eps / max(cos(lat), 0.15),    eps = 0.25 deg
J_e  = d(u, v) / d(East),            J_n = d(u, v) / d(North)   (seam-safe one-sided differences)
e_hat = J_e / ||J_e||,               n_hat = J_n / ||J_n||
(v_u, v_v) = u_east * e_hat + v_north * n_hat
```

## 2. Ocean Circulation & Spherical RK2 Advection (`src/currents.rs`, `src/app.rs`)

### Subsurface Exclusion & Metric-Weighted Jet Filaments

Surface velocity `(u_east, v_north)` combines planetary zonal bands (ACC Gaussian jet in `[-66°, -42°]` and trade/westerly drift `u_bg = -0.28 cos(3.2 lat) cos(lat)` in `|lat| < 44°`) with 13 surface major current polylines. The subsurface Cromwell Equatorial Undercurrent (`ThermalRegime::SubsurfaceEquatorial`, `91–200 m`) is excluded from surface sampling (`include_subsurface = false`).

To eliminate shared-vertex Gaussian double-counting across adjacent polyline segments `windows(2)` and prevent polar directional distortion:

1. Each polyline segment `(x0, y0) -> (x1, y1)` is unwrapped across `±180°` and projected in physical spherical surface coordinates `(dx_m = dx * cos(lat_mid), dy_m = dy)`.
2. For each current `c`, only the single closest segment (`min dist^2`) contributes:
   ```
   dist^2 = ((lon - lon_c) cos(lat))^2 + (lat - lat_c)^2
   w      = exp(-0.5 * dist^2 / sigma^2),    sigma^2 = 34.0 deg^2
   speed  = 0.62 * v_peak * scale_ENSO(c, ONI) * w
   (du, dv) = speed * (dx_m, dy_m) / sqrt(dx_m^2 + dy_m^2)
   ```

### Spherical RK2 Lagrangian Drifter Integration

Drifters integrate geographic `(lon, lat)` using a 2nd-order Runge-Kutta (midpoint) step with spherical `1 / cos(lat)` zonal metric scaling (`inv_cos = 1 / max(cos(lat), 0.15)`):

```
(u1, v1) = V(lon_n, lat_n)
lon_mid  = lon_n + 0.5 * h * u1 / cos(lat_n)
lat_mid  = clamp(lat_n + 0.5 * h * v1, -82, 82)
(u2, v2) = V(lon_mid, lat_mid)
lon_{n+1} = wrap_180(lon_n + h * u2 / cos(lat_mid))
lat_{n+1} = clamp(lat_n + h * v2, -82, 82)
```

## 3. Historical ENSO SST Flavors & Tropical Pacific `20°C` Thermocline (`src/enso.rs`)

- **NOAA CPC ERSSTv6 ONI Table**: Covers `1950 DJF` through `2026 ASO (+2.5°C)` across 77 years (`ONI_TABLE_1950_2026`) with seasonal clamping for future 2026 windows.
- **Spatial SST Anomaly**: Sums equatorial warm/cold tongue (`core_center = 198°E` for Modoki vs `246°E` for canonical Eastern Pacific), Peruvian upwelling (`Niño 1+2`), North American coastally trapped Kelvin wave corridor, Western Pacific warm-pool horseshoe counter-anomaly, NE Pacific Warm Blob (`2014–2016`), North Atlantic heatwave (`2023–2024`), PDO multi-decadal shift, and secular 1950–2026 ocean warming trend.
- **Masked Tropical Pacific `20°C` Isotherm Depth**:
  ```
  t = clamp((lon_360 - 120) / 160, 0, 1)         for |lat| <= 20 deg and 120E <= lon_360 <= 280E
  Z_20(lon) = clamp(172 - 124 * t^1.15 + shift_ENSO(t, year, ONI), 16, 225) m
  ```
  Outside `|lat| <= 20°` or `[120°E, 80°W]`, `tropical_pacific_20c_isotherm_depth_m` returns `None` and both the CPU and WGSL (`thermocline_m < 0.0`) renderers mask non-Pacific basins cleanly.

## 4. Primary Data Sources

1. NOAA CPC ERSSTv6 Oceanic Niño Index (ONI): <https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/enso/oni/v6/>
2. NOAA CoastWatch ERDDAP OISST v2.1 (`ncdcOisst21Agg` / `ncdcOisst21NrtAgg`): <https://coastwatch.pfeg.noaa.gov/erddap/griddap/ncdcOisst21Agg.html>
3. NOAA/NESDIS Altimetric Geostrophic Currents (`nesdisSSH1day`): <https://coastwatch.pfeg.noaa.gov/erddap/griddap/nesdisSSH1day.graph>
4. PROJ 9.9 Spilhaus Projection (`+proj=spilhaus`, `ESRI:54099`): <https://proj.org/en/stable/operations/projections/spilhaus.html>
5. NOAA PMEL TAO/TRITON Equatorial Pacific Thermocline & Kelvin Wave Observations: <https://www.pmel.noaa.gov/elnino/el-nino-and-sea-level/>
