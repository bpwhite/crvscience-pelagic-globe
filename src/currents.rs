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

//! Global surface & equatorial subsurface ocean circulation catalog and calm RK2 Lagrangian ribbon engine.

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ThermalRegime {
    Warm,
    Cold,
    SubsurfaceEquatorial,
    Circumpolar,
}

#[derive(Debug, Clone)]
pub struct MajorCurrent {
    pub id: &'static str,
    pub name: &'static str,
    pub basin: &'static str,
    pub regime: ThermalRegime,
    pub transport_sv: &'static str,
    pub peak_speed_ms: f64,
    pub focus_lon_deg: f64,
    pub focus_lat_deg: f64,
    pub enso_response: &'static str,
    pub source_title: &'static str,
    pub source_url: &'static str,
    /// Polyline spine `(lon_deg, lat_deg)` driving ribbon jet intensification.
    pub waypoints: &'static [(f64, f64)],
}

pub const MAJOR_CURRENTS: &[MajorCurrent] = &[
    MajorCurrent {
        id: "acc",
        name: "Antarctic Circumpolar Current (ACC)",
        basin: "Southern Ocean",
        regime: ThermalRegime::Circumpolar,
        transport_sv: "135 Sv (Drake Passage) – 147 Sv (S. of Tasmania)",
        peak_speed_ms: 1.05,
        focus_lon_deg: -65.0,
        focus_lat_deg: -58.0,
        enso_response: "Southern Annular Mode & Rossby wave teleconnection shifts frontal meanders.",
        source_title: "Antarctic Circumpolar Current — Orsi et al. / Wikipedia",
        source_url: "https://en.wikipedia.org/wiki/Antarctic_Circumpolar_Current",
        waypoints: &[
            (-180.0, -56.0), (-120.0, -58.0), (-65.0, -58.5), (-20.0, -53.0),
            (30.0, -50.0), (90.0, -52.0), (147.0, -54.0), (180.0, -56.0),
        ],
    },
    MajorCurrent {
        id: "gulf_stream",
        name: "Gulf Stream & North Atlantic Drift",
        basin: "North Atlantic",
        regime: ThermalRegime::Warm,
        transport_sv: "30 Sv (Florida Strait) – 150 Sv (55°W)",
        peak_speed_ms: 2.50,
        focus_lon_deg: -72.0,
        focus_lat_deg: 35.5,
        enso_response: "Subtle winter jet-stream & NAO modulation following strong Pacific El Niño winters.",
        source_title: "NOAA National Ocean Service — Boundary Currents",
        source_url: "https://oceanservice.noaa.gov/education/tutorial_currents/04currents3.html",
        waypoints: &[
            (-81.0, 25.0), (-79.0, 30.0), (-74.0, 35.5), (-63.0, 40.0), (-45.0, 44.5), (-20.0, 52.0), (4.0, 62.0),
        ],
    },
    MajorCurrent {
        id: "kuroshio",
        name: "Kuroshio Current & Extension",
        basin: "North Pacific",
        regime: ThermalRegime::Warm,
        transport_sv: "25 Sv (Luzon) – 65 Sv (SE of Japan)",
        peak_speed_ms: 2.20,
        focus_lon_deg: 138.0,
        focus_lat_deg: 33.0,
        enso_response: "NEC bifurcation latitude shifts northward during El Niño, modulating Kuroshio intrusion.",
        source_title: "Kuroshio Current — Wikipedia / JMA",
        source_url: "https://en.wikipedia.org/wiki/Kuroshio_Current",
        waypoints: &[
            (123.0, 18.0), (124.5, 24.0), (131.0, 30.5), (139.5, 34.8), (152.0, 36.0), (172.0, 38.0),
        ],
    },
    MajorCurrent {
        id: "agulhas",
        name: "Agulhas Current & Retroflection",
        basin: "Southwest Indian",
        regime: ThermalRegime::Warm,
        transport_sv: "70 Sv net (95–136 Sv peak at Agulhas Bank)",
        peak_speed_ms: 2.45,
        focus_lon_deg: 27.0,
        focus_lat_deg: -35.0,
        enso_response: "Coupled to Indian Ocean Dipole & Indonesian Throughflow heat supply.",
        source_title: "Agulhas Current — Bryden et al. / Wikipedia",
        source_url: "https://en.wikipedia.org/wiki/Agulhas_Current",
        waypoints: &[
            (38.0, -20.0), (33.0, -28.0), (27.5, -34.0), (20.5, -38.5), (26.0, -40.5), (42.0, -41.0),
        ],
    },
    MajorCurrent {
        id: "necc_pacific",
        name: "Pacific North Equatorial Countercurrent (NECC)",
        basin: "Equatorial Pacific",
        regime: ThermalRegime::Warm,
        transport_sv: "20–25 Sv (21 Sv West, 14.2 Sv Central, 12 Sv East)",
        peak_speed_ms: 1.15,
        focus_lon_deg: -145.0,
        focus_lat_deg: 6.5,
        enso_response: "STRONGLY INTENSIFIES (+45% to +80%) during El Niño, surging warm pool water east.",
        source_title: "Equatorial Counter Current — Wikipedia",
        source_url: "https://en.wikipedia.org/wiki/Equatorial_Counter_Current",
        waypoints: &[
            (132.0, 6.5), (160.0, 6.8), (-170.0, 7.0), (-140.0, 7.2), (-110.0, 7.5), (-85.0, 7.2),
        ],
    },
    MajorCurrent {
        id: "cromwell_euc",
        name: "Cromwell Current (Equatorial Undercurrent)",
        basin: "Equatorial Pacific (91–200 m depth)",
        regime: ThermalRegime::SubsurfaceEquatorial,
        transport_sv: "Up to 30 Sv along 20°C thermocline core",
        peak_speed_ms: 1.50,
        focus_lon_deg: -125.0,
        focus_lat_deg: 0.0,
        enso_response: "Weakens sharply when eastern Pacific sea-level slope flattens during peak El Niño.",
        source_title: "Cromwell Current — Wikipedia",
        source_url: "https://en.wikipedia.org/wiki/Cromwell_Current",
        waypoints: &[
            (152.0, 0.0), (180.0, 0.0), (-150.0, 0.0), (-120.0, 0.0), (-92.0, -0.3),
        ],
    },
    MajorCurrent {
        id: "sec_pacific",
        name: "Pacific South Equatorial Current (SEC)",
        basin: "Tropical Pacific",
        regime: ThermalRegime::Warm,
        transport_sv: "25–55 Sv westward trade-wind limb",
        peak_speed_ms: 1.10,
        focus_lon_deg: -140.0,
        focus_lat_deg: -3.5,
        enso_response: "WEAKENS & REVERSES in central-western Pacific under El Niño westerly wind bursts; accelerates in La Niña.",
        source_title: "NOAA PMEL — ENSO & South Equatorial Current",
        source_url: "https://www.pmel.noaa.gov/elnino/featured-publication/enso-and-short-term-variability-south-equatorial-current-entering-coral-sea/",
        waypoints: &[
            (-86.0, -4.0), (-115.0, -3.5), (-145.0, -3.0), (-175.0, -3.5), (158.0, -6.0),
        ],
    },
    MajorCurrent {
        id: "nec_pacific",
        name: "Pacific North Equatorial Current (NEC)",
        basin: "North Pacific",
        regime: ThermalRegime::Warm,
        transport_sv: "30–45 Sv westward (>12,000 km span)",
        peak_speed_ms: 0.75,
        focus_lon_deg: -155.0,
        focus_lat_deg: 14.0,
        enso_response: "Bifurcates off Luzon into Kuroshio and Mindanao (13–39 Sv); shifts north during El Niño.",
        source_title: "Mindanao & North Equatorial Current — Wikipedia",
        source_url: "https://en.wikipedia.org/wiki/Mindanao_Current",
        waypoints: &[
            (-112.0, 15.0), (-145.0, 14.5), (-178.0, 14.0), (152.0, 13.8), (127.0, 14.2),
        ],
    },
    MajorCurrent {
        id: "humboldt",
        name: "Humboldt (Peru) Current",
        basin: "Southeast Pacific",
        regime: ThermalRegime::Cold,
        transport_sv: "15–20 Sv broad coastal upwelling (~900 km wide)",
        peak_speed_ms: 0.65,
        focus_lon_deg: -78.0,
        focus_lat_deg: -18.0,
        enso_response: "UPWELLING COLLAPSES during strong El Niño as warm water deepens the 20°C isotherm off Peru.",
        source_title: "NOAA National Ocean Service — Boundary Currents",
        source_url: "https://oceanservice.noaa.gov/education/tutorial_currents/04currents3.html",
        waypoints: &[
            (-75.0, -42.0), (-73.5, -30.0), (-77.0, -18.0), (-81.5, -6.0),
        ],
    },
    MajorCurrent {
        id: "california",
        name: "California Current",
        basin: "Northeast Pacific (Oregon to Baja)",
        regime: ThermalRegime::Cold,
        transport_sv: "10–14 Sv eastern boundary limb",
        peak_speed_ms: 0.55,
        focus_lon_deg: -125.5,
        focus_lat_deg: 41.0,
        enso_response: "Coastally trapped El Niño Kelvin waves raise sea level (~15 cm) and suppress cold PNW/California upwelling.",
        source_title: "NOAA Fisheries — West Coast El Niño Ocean Indicators",
        source_url: "https://www.fisheries.noaa.gov/feature-story/ocean-indicators-were-tracking-el-ninos-arrival-west-coast",
        waypoints: &[
            (-130.0, 48.0), (-125.8, 44.0), (-123.5, 37.0), (-119.5, 31.0), (-114.0, 23.0),
        ],
    },
    MajorCurrent {
        id: "itf",
        name: "Indonesian Throughflow (ITF)",
        basin: "Maritime Continent (Pacific -> Indian)",
        regime: ThermalRegime::Warm,
        transport_sv: "~15 Sv (12.7 Sv Makassar + Lifamatola per INSTANT)",
        peak_speed_ms: 1.20,
        focus_lon_deg: 118.5,
        focus_lat_deg: -8.5,
        enso_response: "Decreases during El Niño as Western Pacific sea surface height drops; strengthens in La Niña.",
        source_title: "Indonesian Throughflow — INSTANT / Wikipedia",
        source_url: "https://en.wikipedia.org/wiki/Indonesian_Throughflow",
        waypoints: &[
            (126.0, 4.0), (119.0, -2.0), (116.5, -9.0), (105.0, -13.0), (80.0, -15.0),
        ],
    },
    MajorCurrent {
        id: "eac",
        name: "East Australian Current (EAC)",
        basin: "Southwest Pacific",
        regime: ThermalRegime::Warm,
        transport_sv: "22–30 Sv (vigorous Tasman eddy shedding)",
        peak_speed_ms: 1.80,
        focus_lon_deg: 154.0,
        focus_lat_deg: -29.0,
        enso_response: "Separation latitude (~31°S–33°S) and Tasman Sea heatwaves covary with ENSO.",
        source_title: "Ocean Gyres & Boundary Currents — Wikipedia",
        source_url: "https://en.wikipedia.org/wiki/Ocean_gyre",
        waypoints: &[
            (153.0, -16.0), (154.2, -25.0), (152.5, -33.0), (158.0, -37.0), (170.0, -38.0),
        ],
    },
    MajorCurrent {
        id: "brazil",
        name: "Brazil Current",
        basin: "South Atlantic",
        regime: ThermalRegime::Warm,
        transport_sv: "16.2 Sv (meets Malvinas Current at 38°S)",
        peak_speed_ms: 1.10,
        focus_lon_deg: -44.0,
        focus_lat_deg: -26.0,
        enso_response: "Modulated by South Atlantic subtropical high anomalies during extreme ENSO events.",
        source_title: "Western Boundary Transport Comparison — Wikipedia",
        source_url: "https://en.wikipedia.org/wiki/Agulhas_Current",
        waypoints: &[
            (-34.5, -8.0), (-38.0, -18.0), (-45.0, -27.0), (-53.0, -37.5), (-35.0, -40.0),
        ],
    },
    MajorCurrent {
        id: "benguela",
        name: "Benguela Current",
        basin: "Southeast Atlantic",
        regime: ThermalRegime::Cold,
        transport_sv: "15–20 Sv upwelling + Agulhas Rings (3–9 Sv)",
        peak_speed_ms: 0.65,
        focus_lon_deg: 12.0,
        focus_lat_deg: -24.0,
        enso_response: "Receives warm-core Agulhas Rings and feeds the Atlantic South Equatorial Current.",
        source_title: "South Atlantic Gyre — Wikipedia",
        source_url: "https://en.wikipedia.org/wiki/South_Atlantic_Gyre",
        waypoints: &[
            (17.0, -34.0), (13.0, -25.0), (10.0, -15.0), (-5.0, -8.0), (-25.0, -5.0),
        ],
    },
];

/// Evaluates the ENSO-coupled surface velocity `(u_east_ms, v_north_ms)` at `(lon_deg, lat_deg)`.
pub fn sample_velocity_ms(lon_deg: f64, lat_deg: f64, oni: f64) -> (f64, f64) {
    let mut u = 0.0;
    let mut v = 0.0;

    // Broad planetary wind-driven background bands (trades, westerlies, ACC)
    let lat_rad = lat_deg.to_radians();
    if lat_deg < -42.0 && lat_deg > -66.0 {
        let acc = (-0.5 * ((lat_deg + 54.0) / 6.5).powi(2)).exp();
        u += 0.72 * acc;
    } else if lat_deg.abs() < 44.0 {
        // Trade winds & mid-latitude westerly drift
        u += -0.28 * (3.2 * lat_rad).cos() * lat_rad.cos();
    }

    // Add Gaussian jet filaments along all 14 verified major currents
    for cur in MAJOR_CURRENTS {
        let mut scale = 1.0;
        if cur.id == "sec_pacific" {
            // Weakens & reverses in western/central Pacific during strong El Niño
            scale = (1.0 - 0.52 * oni).clamp(-0.65, 1.85);
        } else if cur.id == "necc_pacific" {
            scale = (1.0 + 0.32 * oni).clamp(0.35, 2.1);
        } else if cur.id == "humboldt" || cur.id == "california" {
            scale = (1.0 - 0.22 * oni).clamp(0.35, 1.5);
        } else if cur.id == "itf" {
            scale = (1.0 - 0.18 * oni).clamp(0.45, 1.45);
        }

        for seg in cur.waypoints.windows(2) {
            let (mut x0, y0) = seg[0];
            let (mut x1, y1) = seg[1];
            // Unwrap dateline
            if (x1 - x0).abs() > 180.0 {
                if x1 < x0 {
                    x1 += 360.0;
                } else {
                    x0 += 360.0;
                }
            }
            let mut px = lon_deg;
            if (px - 0.5 * (x0 + x1)) > 180.0 {
                px -= 360.0;
            } else if (px - 0.5 * (x0 + x1)) < -180.0 {
                px += 360.0;
            }

            let dx = x1 - x0;
            let dy = y1 - y0;
            let len2 = dx * dx + dy * dy + 1.0e-6;
            let t = (((px - x0) * dx + (lat_deg - y0) * dy) / len2).clamp(0.0, 1.0);
            let cx = x0 + t * dx;
            let cy = y0 + t * dy;
            let dist2 = ((px - cx) * lat_rad.cos()).powi(2) + (lat_deg - cy).powi(2);
            let sigma2 = 34.0; // ~5.8 deg jet half-width
            if dist2 < 4.5 * sigma2 {
                let w = (-0.5 * dist2 / sigma2).exp();
                let norm = (dx * dx + dy * dy).sqrt().max(1.0e-4);
                let speed = 0.62 * cur.peak_speed_ms * scale * w;
                u += speed * (dx / norm);
                v += speed * (dy / norm);
            }
        }
    }

    (u.clamp(-2.6, 2.6), v.clamp(-2.6, 2.6))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sec_reverses_and_necc_accelerates_under_super_el_nino() {
        let (u_sec_lanina, _) = sample_velocity_ms(-150.0, -3.2, -1.8);
        let (u_sec_super, _) = sample_velocity_ms(-150.0, -3.2, 2.5);
        assert!(u_sec_lanina < -0.5, "SEC should flow strongly west in La Niña");
        assert!(u_sec_super > u_sec_lanina + 0.6, "SEC should weaken/reverse in Super El Niño");

        let (u_necc_neutral, _) = sample_velocity_ms(-145.0, 7.0, 0.0);
        let (u_necc_super, _) = sample_velocity_ms(-145.0, 7.0, 2.5);
        assert!(u_necc_super > u_necc_neutral, "NECC must intensify eastward in El Niño");
    }
}
