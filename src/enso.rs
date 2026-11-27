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

//! NOAA Climate Prediction Center ERSSTv6 Oceanic Niño Index (ONI) dataset (1950–2026)
//! and coupled equatorial Pacific sea-surface temperature (SST) & 20°C thermocline model.
//!
//! Data source: <https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/enso/oni/v6/>

use serde::{Deserialize, Serialize};

pub const SEASONS: [&str; 12] = [
    "DJF", "JFM", "FMA", "MAM", "AMJ", "MJJ", "JJA", "JAS", "ASO", "SON", "OND", "NDJ",
];

pub const START_YEAR: i32 = 1950;
pub const END_YEAR: i32 = 2026;

/// NOAA CPC ERSSTv6 Oceanic Niño Index (ONI) values (`°C`) for 1950–2026 (`f32::NAN` for future 2026 seasons).
pub const ONI_TABLE_1950_2026: [[f32; 12]; 77] = [
    [-1.3, -1.2, -1.1, -1.1, -1.1, -0.9, -0.6, -0.6, -0.6, -0.6, -0.7, -0.8], // 1950
    [-0.7, -0.4,  0.0,  0.2,  0.3,  0.5,  0.7,  0.9,  1.0,  1.0,  1.0,  0.7], // 1951
    [ 0.5,  0.3,  0.3,  0.2, -0.1, -0.4, -0.4, -0.2,  0.0,  0.0, -0.1,  0.1], // 1952
    [ 0.3,  0.5,  0.5,  0.6,  0.8,  0.7,  0.6,  0.7,  0.6,  0.8,  0.5,  0.6], // 1953
    [ 0.4,  0.4,  0.0, -0.3, -0.5, -0.5, -0.7, -0.8, -0.7, -0.5, -0.4, -0.4], // 1954
    [-0.5, -0.5, -0.6, -0.8, -0.8, -0.9, -0.7, -0.9, -1.1, -1.4, -1.6, -1.5], // 1955
    [-1.2, -0.8, -0.6, -0.5, -0.5, -0.5, -0.5, -0.5, -0.4, -0.4, -0.4, -0.4], // 1956
    [-0.3, -0.1,  0.3,  0.6,  0.7,  0.8,  1.0,  1.0,  1.1,  1.2,  1.4,  1.6], // 1957
    [ 1.7,  1.6,  1.2,  0.8,  0.6,  0.5,  0.4,  0.3,  0.3,  0.3,  0.4,  0.6], // 1958
    [ 0.6,  0.6,  0.5,  0.4,  0.2,  0.0, -0.2, -0.2,  0.0,  0.1,  0.2,  0.1], // 1959
    [ 0.0,  0.0,  0.1,  0.2,  0.0,  0.0,  0.0,  0.2,  0.2,  0.1, -0.1, -0.1], // 1960
    [ 0.0,  0.0,  0.0,  0.1,  0.2,  0.2,  0.1, -0.2, -0.3, -0.2, -0.1,  0.0], // 1961
    [-0.2, -0.2, -0.1, -0.2, -0.2, -0.1,  0.0, -0.1, -0.2, -0.2, -0.2, -0.3], // 1962
    [-0.3, -0.2,  0.0,  0.1,  0.1,  0.5,  0.8,  1.1,  1.1,  1.1,  1.2,  1.1], // 1963
    [ 0.9,  0.5,  0.1, -0.3, -0.6, -0.7, -0.7, -0.7, -0.7, -0.8, -0.8, -0.7], // 1964
    [-0.5, -0.2,  0.0,  0.2,  0.5,  0.8,  1.1,  1.3,  1.6,  1.7,  1.8,  1.6], // 1965
    [ 1.2,  1.0,  0.8,  0.4,  0.2,  0.2,  0.3,  0.1, -0.1, -0.1, -0.2, -0.3], // 1966
    [-0.5, -0.5, -0.5, -0.4, -0.1,  0.1,  0.0, -0.3, -0.4, -0.4, -0.4, -0.4], // 1967
    [-0.6, -0.7, -0.6, -0.5, -0.1,  0.2,  0.4,  0.4,  0.3,  0.5,  0.7,  0.9], // 1968
    [ 1.0,  1.0,  0.8,  0.7,  0.6,  0.5,  0.4,  0.4,  0.6,  0.7,  0.8,  0.7], // 1969
    [ 0.5,  0.2,  0.2,  0.3,  0.1, -0.3, -0.5, -0.7, -0.7, -0.7, -0.8, -1.1], // 1970
    [-1.3, -1.4, -1.1, -0.8, -0.7, -0.7, -0.7, -0.6, -0.7, -0.8, -0.9, -0.8], // 1971
    [-0.6, -0.3,  0.0,  0.3,  0.5,  0.7,  1.0,  1.3,  1.5,  1.7,  1.8,  1.8], // 1972
    [ 1.6,  1.1,  0.6,  0.0, -0.5, -0.8, -1.0, -1.1, -1.3, -1.6, -1.9, -2.0], // 1973
    [-1.9, -1.6, -1.2, -1.0, -0.7, -0.6, -0.5, -0.4, -0.5, -0.7, -0.9, -0.7], // 1974
    [-0.6, -0.5, -0.6, -0.7, -0.8, -0.9, -1.1, -1.1, -1.3, -1.3, -1.5, -1.6], // 1975
    [-1.5, -1.1, -0.6, -0.4, -0.2,  0.0,  0.2,  0.4,  0.6,  0.8,  0.8,  0.8], // 1976
    [ 0.7,  0.7,  0.3,  0.3,  0.3,  0.5,  0.4,  0.4,  0.5,  0.7,  0.8,  0.8], // 1977
    [ 0.7,  0.5,  0.2, -0.1, -0.3, -0.2, -0.4, -0.4, -0.4, -0.3, -0.1,  0.1], // 1978
    [ 0.2,  0.2,  0.3,  0.3,  0.3,  0.1,  0.1,  0.2,  0.3,  0.5,  0.5,  0.6], // 1979
    [ 0.6,  0.5,  0.4,  0.3,  0.5,  0.4,  0.2, -0.1, -0.2, -0.1,  0.1,  0.0], // 1980
    [-0.3, -0.4, -0.4, -0.3, -0.3, -0.3, -0.3, -0.2, -0.1, -0.2, -0.2, -0.1], // 1981
    [ 0.0,  0.0,  0.1,  0.3,  0.7,  0.7,  0.8,  1.0,  1.5,  1.8,  2.0,  2.1], // 1982
    [ 2.1,  1.9,  1.5,  1.2,  1.0,  0.6,  0.3, -0.1, -0.4, -0.8, -1.0, -0.9], // 1983
    [-0.6, -0.5, -0.4, -0.5, -0.5, -0.5, -0.4, -0.2, -0.3, -0.6, -0.9, -1.1], // 1984
    [-1.0, -0.9, -0.8, -0.8, -0.8, -0.6, -0.4, -0.3, -0.3, -0.3, -0.3, -0.4], // 1985
    [-0.4, -0.4, -0.2, -0.2, -0.1,  0.0,  0.2,  0.5,  0.6,  0.9,  1.1,  1.2], // 1986
    [ 1.2,  1.1,  1.0,  0.9,  1.0,  1.1,  1.3,  1.5,  1.5,  1.4,  1.2,  1.0], // 1987
    [ 0.6,  0.4,  0.0, -0.4, -1.0, -1.3, -1.3, -1.1, -1.3, -1.5, -1.8, -1.8], // 1988
    [-1.6, -1.4, -1.1, -0.8, -0.6, -0.4, -0.4, -0.4, -0.3, -0.2, -0.1,  0.0], // 1989
    [ 0.2,  0.2,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3], // 1990
    [ 0.4,  0.3,  0.2,  0.3,  0.5,  0.6,  0.8,  0.7,  0.7,  0.8,  1.2,  1.4], // 1991
    [ 1.5,  1.5,  1.4,  1.3,  1.1,  0.7,  0.4,  0.1, -0.1, -0.2, -0.2,  0.0], // 1992
    [ 0.1,  0.3,  0.5,  0.7,  0.8,  0.6,  0.3,  0.2,  0.2,  0.2,  0.2,  0.2], // 1993
    [ 0.1,  0.0,  0.0,  0.2,  0.3,  0.2,  0.3,  0.3,  0.5,  0.7,  0.9,  1.0], // 1994
    [ 0.9,  0.7,  0.5,  0.2,  0.0, -0.1, -0.3, -0.5, -0.8, -1.0, -1.0, -1.0], // 1995
    [-0.9, -0.7, -0.5, -0.4, -0.2, -0.2, -0.2, -0.3, -0.3, -0.3, -0.4, -0.4], // 1996
    [-0.4, -0.3,  0.0,  0.3,  0.7,  1.1,  1.5,  1.8,  2.0,  2.2,  2.3,  2.4], // 1997
    [ 2.2,  1.9,  1.4,  1.0,  0.4, -0.1, -0.6, -0.8, -1.0, -1.1, -1.3, -1.5], // 1998
    [-1.5, -1.2, -0.9, -0.8, -0.8, -0.8, -0.9, -0.9, -1.0, -1.1, -1.3, -1.5], // 1999
    [-1.5, -1.3, -1.0, -0.8, -0.7, -0.6, -0.5, -0.4, -0.4, -0.5, -0.8, -0.8], // 2000
    [-0.8, -0.5, -0.4, -0.3, -0.2, -0.1,  0.0, -0.1, -0.1, -0.3, -0.3, -0.3], // 2001
    [-0.1,  0.1,  0.1,  0.2,  0.3,  0.5,  0.6,  0.6,  0.8,  1.0,  1.2,  1.1], // 2002
    [ 0.9,  0.6,  0.4,  0.0, -0.2, -0.1,  0.1,  0.1,  0.2,  0.2,  0.3,  0.3], // 2003
    [ 0.3,  0.3,  0.2,  0.1,  0.1,  0.2,  0.4,  0.5,  0.6,  0.6,  0.6,  0.7], // 2004
    [ 0.6,  0.6,  0.4,  0.4,  0.3,  0.1,  0.0, -0.1, -0.1, -0.3, -0.6, -0.8], // 2005
    [-0.9, -0.8, -0.6, -0.4, -0.2,  0.0,  0.1,  0.3,  0.5,  0.7,  0.9,  0.9], // 2006
    [ 0.6,  0.2, -0.1, -0.3, -0.3, -0.4, -0.5, -0.8, -1.0, -1.4, -1.6, -1.7], // 2007
    [-1.8, -1.6, -1.4, -1.0, -0.8, -0.6, -0.3, -0.2, -0.2, -0.4, -0.6, -0.8], // 2008
    [-0.9, -0.8, -0.6, -0.3,  0.0,  0.3,  0.5,  0.6,  0.7,  0.9,  1.3,  1.5], // 2009
    [ 1.5,  1.2,  0.9,  0.4, -0.2, -0.7, -1.0, -1.2, -1.4, -1.5, -1.6, -1.5], // 2010
    [-1.3, -1.1, -0.8, -0.7, -0.5, -0.3, -0.4, -0.5, -0.8, -0.9, -1.0, -0.9], // 2011
    [-0.7, -0.6, -0.4, -0.3, -0.2,  0.0,  0.3,  0.4,  0.4,  0.3,  0.1,  0.0], // 2012
    [-0.2, -0.3, -0.2, -0.2, -0.3, -0.4, -0.4, -0.3, -0.3, -0.1, -0.1, -0.1], // 2013
    [-0.2, -0.2,  0.0,  0.3,  0.4,  0.2,  0.1,  0.1,  0.2,  0.5,  0.7,  0.8], // 2014
    [ 0.7,  0.7,  0.7,  0.9,  1.0,  1.2,  1.4,  1.7,  2.0,  2.3,  2.5,  2.6], // 2015
    [ 2.5,  2.2,  1.7,  1.1,  0.6,  0.1, -0.2, -0.3, -0.4, -0.5, -0.5, -0.4], // 2016
    [-0.1,  0.1,  0.3,  0.3,  0.4,  0.3,  0.1, -0.1, -0.2, -0.4, -0.6, -0.8], // 2017
    [-0.7, -0.7, -0.5, -0.3, -0.1,  0.1,  0.2,  0.3,  0.5,  0.8,  1.0,  1.1], // 2018
    [ 1.0,  0.9,  0.9,  0.8,  0.7,  0.6,  0.4,  0.2,  0.3,  0.5,  0.7,  0.8], // 2019
    [ 0.7,  0.7,  0.6,  0.3,  0.0, -0.2, -0.3, -0.4, -0.8, -1.0, -1.1, -1.1], // 2020
    [-1.0, -0.9, -0.7, -0.6, -0.4, -0.3, -0.3, -0.5, -0.6, -0.8, -0.9, -0.8], // 2021
    [-0.8, -0.7, -0.8, -0.9, -0.8, -0.7, -0.7, -0.8, -0.9, -0.9, -0.8, -0.7], // 2022
    [-0.5, -0.3, -0.1,  0.2,  0.5,  0.7,  1.0,  1.3,  1.5,  1.7,  1.9,  2.0], // 2023
    [ 1.8,  1.5,  1.2,  0.8,  0.4,  0.2,  0.1,  0.0, -0.1, -0.2, -0.3, -0.4], // 2024
    [-0.5, -0.2, -0.1,  0.0,  0.0,  0.0, -0.1, -0.3, -0.4, -0.6, -0.6, -0.6], // 2025
    [-0.4, -0.2,  0.1,  0.5,  0.9,  1.4,  1.8,  2.2,  2.5, f32::NAN, f32::NAN, f32::NAN], // 2026
];

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum EnsoPhase {
    SuperElNino,
    StrongElNino,
    ModerateElNino,
    WeakElNino,
    Neutral,
    WeakLaNina,
    ModerateLaNina,
    StrongLaNina,
}

impl EnsoPhase {
    pub fn from_oni(oni: f64) -> Self {
        if oni >= 2.0 {
            Self::SuperElNino
        } else if oni >= 1.5 {
            Self::StrongElNino
        } else if oni >= 1.0 {
            Self::ModerateElNino
        } else if oni >= 0.5 {
            Self::WeakElNino
        } else if oni <= -1.5 {
            Self::StrongLaNina
        } else if oni <= -1.0 {
            Self::ModerateLaNina
        } else if oni <= -0.5 {
            Self::WeakLaNina
        } else {
            Self::Neutral
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::SuperElNino => "Very Strong ('Super') El Niño",
            Self::StrongElNino => "Strong El Niño",
            Self::ModerateElNino => "Moderate El Niño",
            Self::WeakElNino => "Weak El Niño",
            Self::Neutral => "ENSO Neutral",
            Self::WeakLaNina => "Weak La Niña",
            Self::ModerateLaNina => "Moderate La Niña",
            Self::StrongLaNina => "Strong La Niña",
        }
    }
}

#[derive(Debug, Clone)]
pub struct EnsoPreset {
    pub slug: &'static str,
    pub title: &'static str,
    pub year: i32,
    pub season_index: usize,
    pub oni: f64,
    pub summary: &'static str,
}

pub const ENSO_PRESETS: &[EnsoPreset] = &[
    EnsoPreset {
        slug: "super-2026",
        title: "2026 ASO (+2.5°C) — Developing Super El Niño",
        year: 2026,
        season_index: 8,
        oni: 2.5,
        summary: "Latest NOAA CPC ERSSTv6 reading (ASO 2026 +2.5°C ONI / +2.0°C RONI); intense subsurface Kelvin wave warming across Niño 3.4 and Niño 1+2.",
    },
    EnsoPreset {
        slug: "strong-2023",
        title: "2023–24 NDJ (+2.0°C) — Very Strong El Niño",
        year: 2023,
        season_index: 11,
        oni: 2.0,
        summary: "Peaked at +2.0°C ONI in NDJ 2023; drove coastally trapped Kelvin waves along the U.S. West Coast and global coral bleaching.",
    },
    EnsoPreset {
        slug: "lanina-2020",
        title: "2020–23 OND (-1.1°C) — Triple-Dip La Niña",
        year: 2020,
        season_index: 10,
        oni: -1.1,
        summary: "Multi-year cold-tongue intensification with strengthened Pacific trade winds, rapid South Equatorial Current, and deep Western Pacific warm pool.",
    },
    EnsoPreset {
        slug: "super-2015",
        title: "2015–16 NDJ (+2.6°C) — Record Super El Niño",
        year: 2015,
        season_index: 11,
        oni: 2.6,
        summary: "Highest 3-month ONI (+2.6°C in NDJ 2015) in the 1950–present NOAA CPC ERSSTv6 table; pronounced central-eastern Pacific warm anomaly + NE Pacific Warm Blob.",
    },
    EnsoPreset {
        slug: "modoki-2009",
        title: "2009–10 NDJ (+1.5°C) — Central Pacific (Modoki) El Niño",
        year: 2009,
        season_index: 11,
        oni: 1.5,
        summary: "Canonical Central Pacific (Modoki) warm-pool displacement peaking near the International Date Line (+1.5°C ONI) with muted Peru coastal warming.",
    },
    EnsoPreset {
        slug: "lanina-1999",
        title: "1998–2000 NDJ (-1.5°C) — Strong La Niña",
        year: 1999,
        season_index: 11,
        oni: -1.5,
        summary: "Sustained post-1998 cold phase with vigorous Humboldt/Peru and California Current upwelling and steepened equatorial thermocline.",
    },
    EnsoPreset {
        slug: "super-1997",
        title: "1997–98 NDJ (+2.4°C) — Classic Super El Niño",
        year: 1997,
        season_index: 11,
        oni: 2.4,
        summary: "Archetypal eastern-Pacific El Niño (+2.4°C ONI); 20°C isotherm flattened to ~150 m in the eastern Pacific and collapsed Peruvian upwelling.",
    },
    EnsoPreset {
        slug: "lanina-1988",
        title: "1988–89 NDJ (-1.8°C) — Strong La Niña",
        year: 1988,
        season_index: 11,
        oni: -1.8,
        summary: "Powerful equatorial cold tongue (-1.8°C ONI) driven by enhanced easterly trade winds and shoaled eastern Pacific thermocline.",
    },
    EnsoPreset {
        slug: "super-1982",
        title: "1982–83 DJF (+2.1°C) — Historic Super El Niño",
        year: 1983,
        season_index: 0,
        oni: 2.1,
        summary: "Landmark 1982–83 event (+2.1°C ONI, +2.4°C RONI) that motivated deployment of NOAA PMEL's TAO/TRITON equatorial mooring array.",
    },
    EnsoPreset {
        slug: "lanina-1973",
        title: "1973–74 NDJ (-2.0°C) — Record La Niña",
        year: 1973,
        season_index: 11,
        oni: -2.0,
        summary: "Coldest 3-month ONI (-2.0°C in NDJ 1973) in the NOAA CPC ERSSTv6 instrumental series.",
    },
    EnsoPreset {
        slug: "neutral-1990",
        title: "1990 DJF (+0.2°C) — ENSO Neutral Baseline",
        year: 1990,
        season_index: 0,
        oni: 0.2,
        summary: "ENSO Neutral reference state (1990 DJF +0.2°C ONI): normal easterly trades, Western Pacific warm pool (~30°C), and sloped 20°C thermocline.",
    },
];

#[derive(Debug, Clone, Copy)]
pub struct NinoBox {
    pub id: &'static str,
    pub name: &'static str,
    pub lat_min: f64,
    pub lat_max: f64,
    pub lon_west: f64,
    pub lon_east: f64,
    pub color_rgb: [u8; 3],
}

pub const NINO_BOXES: &[NinoBox] = &[
    NinoBox {
        id: "nino12",
        name: "Niño 1+2 (0°–10°S, 90°W–80°W)",
        lat_min: -10.0,
        lat_max: 0.0,
        lon_west: -90.0,
        lon_east: -80.0,
        color_rgb: [251, 146, 60],
    },
    NinoBox {
        id: "nino3",
        name: "Niño 3 (5°N–5°S, 150°W–90°W)",
        lat_min: -5.0,
        lat_max: 5.0,
        lon_west: -150.0,
        lon_east: -90.0,
        color_rgb: [248, 113, 113],
    },
    NinoBox {
        id: "nino34",
        name: "Niño 3.4 ONI Core (5°N–5°S, 170°W–120°W)",
        lat_min: -5.0,
        lat_max: 5.0,
        lon_west: -170.0,
        lon_east: -120.0,
        color_rgb: [250, 204, 21],
    },
    NinoBox {
        id: "nino4",
        name: "Niño 4 (5°N–5°S, 160°E–150°W)",
        lat_min: -5.0,
        lat_max: 5.0,
        lon_west: 160.0,
        lon_east: -150.0,
        color_rgb: [56, 189, 248],
    },
];

pub fn lookup_oni(year: i32, season_index: usize) -> Option<f64> {
    if !(START_YEAR..=END_YEAR).contains(&year) || season_index >= 12 {
        return None;
    }
    let val = ONI_TABLE_1950_2026[(year - START_YEAR) as usize][season_index];
    if val.is_nan() {
        None
    } else {
        Some((f64::from(val) * 10.0).round() / 10.0)
    }
}

/// Returns the latest non-NaN 3-month season index (`0..=11`) available for `year`.
pub fn latest_valid_season_for_year(year: i32) -> usize {
    let clamped_y = year.clamp(START_YEAR, END_YEAR);
    let row = &ONI_TABLE_1950_2026[(clamped_y - START_YEAR) as usize];
    row.iter().rposition(|v| !v.is_nan()).unwrap_or(0)
}

/// Clamps `season_index` to the latest published ONI season for `year` and returns `(season_index, oni_celsius)`.
pub fn lookup_oni_clamped(year: i32, season_index: usize) -> (usize, f64) {
    let clamped_y = year.clamp(START_YEAR, END_YEAR);
    let max_s = latest_valid_season_for_year(clamped_y);
    let s = season_index.min(max_s);
    let val = lookup_oni(clamped_y, s).unwrap_or(0.0);
    (s, val)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum HistoricalEnsoFlavor {
    CanonicalEasternPacific,
    CentralPacificModoki,
    MixedSuperWithWarmBlob,
    GlobalHeatwaveElNino,
    DevelopingSuperElNino,
    LaNinaColdTongue,
    NeutralBaseline,
}

impl HistoricalEnsoFlavor {
    pub fn label(self) -> &'static str {
        match self {
            Self::CanonicalEasternPacific => "Canonical Eastern Pacific (EP)",
            Self::CentralPacificModoki => "Central Pacific (Modoki)",
            Self::MixedSuperWithWarmBlob => "Mixed Super + NE Pacific Warm Blob",
            Self::GlobalHeatwaveElNino => "Global Heatwave El Niño",
            Self::DevelopingSuperElNino => "Developing 2026 Super El Niño",
            Self::LaNinaColdTongue => "La Niña Cold Tongue",
            Self::NeutralBaseline => "ENSO Neutral Baseline",
        }
    }
}

pub fn classify_historical_flavor(year: i32, oni: f64) -> HistoricalEnsoFlavor {
    if year >= 2026 && oni >= 1.0 {
        HistoricalEnsoFlavor::DevelopingSuperElNino
    } else if (2023..=2024).contains(&year) && oni >= 0.8 {
        HistoricalEnsoFlavor::GlobalHeatwaveElNino
    } else if (2015..=2016).contains(&year) && oni >= 1.0 {
        HistoricalEnsoFlavor::MixedSuperWithWarmBlob
    } else if matches!(
        year,
        1991 | 1992 | 1994 | 1995 | 2002..=2005 | 2009 | 2010 | 2018 | 2019
    ) && oni >= 0.5
    {
        HistoricalEnsoFlavor::CentralPacificModoki
    } else if oni >= 0.5 {
        HistoricalEnsoFlavor::CanonicalEasternPacific
    } else if oni <= -0.5 {
        HistoricalEnsoFlavor::LaNinaColdTongue
    } else {
        HistoricalEnsoFlavor::NeutralBaseline
    }
}

/// Evaluates the spatial ENSO SST anomaly (`Δ°C`) at `(lon_deg, lat_deg)` for a given Niño 3.4 ONI (`°C`).
///
/// Captures the Equatorial Pacific cold-tongue / warm-tongue structure (`160°E` to `80°W`),
/// South American (Peru/Ecuador) upwelling anomaly, coastally trapped Kelvin wave warming along
/// Baja/California/Oregon (`NOAA Fisheries`), Western Pacific boomerang cooling, and Indian Ocean dipole coupling.
pub fn spatial_sst_anomaly_c(lon_deg: f64, lat_deg: f64, oni: f64) -> f64 {
    if oni.abs() < 1.0e-4 {
        return 0.0;
    }
    // Convert longitude to continuous Pacific domain [0..360) where 120E=120, 180=180, 120W=240, 80W=280
    let lon360 = lon_deg.rem_euclid(360.0);

    // 1. Primary Equatorial Pacific tongue (165E..280E, centered at 225E / 135W)
    let eq_lat = (-0.5 * (lat_deg / 6.8).powi(2)).exp();
    let eq_lon = if (145.0..=285.0).contains(&lon360) {
        let center = 232.0; // 128W
        let width = if lon360 < center { 42.0 } else { 36.0 };
        (-0.5 * ((lon360 - center) / width).powi(2)).exp()
    } else {
        0.0
    };
    let core_tongue = 1.18 * oni * eq_lat * eq_lon;

    // 2. Peruvian / Ecuadorian coastal upwelling wedge (Niño 1+2, 270E..283E, 18S..2N)
    let peru_coast = (-0.5 * ((lon360 - 277.0) / 7.5).powi(2)).exp()
        * (-0.5 * ((lat_deg + 6.5) / 8.5).powi(2)).exp()
        * (0.95 * oni);

    // 3. Coastally trapped NorthEast Pacific Kelvin wave corridor (Baja -> California -> Oregon/Alaska)
    let ne_pacific_kelvin = if (12.0..=56.0).contains(&lat_deg) {
        let coast_lon = 255.0 - 0.58 * (lat_deg - 18.0);
        let d_lon = lon360 - coast_lon;
        if (-18.0..=4.0).contains(&d_lon) {
            (-0.5 * (d_lon / 6.5).powi(2)).exp()
                * (-0.5 * ((lat_deg - 34.0) / 15.0).powi(2)).exp()
                * (0.58 * oni)
        } else {
            0.0
        }
    } else {
        0.0
    };

    // 4. Western Pacific Warm Pool / Horseshoe counter-anomaly (115E..160E & subtropical 18N/18S)
    let wp_cool = (-0.5 * ((lon360 - 138.0) / 20.0).powi(2)).exp()
        * (-0.5 * (lat_deg / 14.0).powi(2)).exp()
        * (-0.38 * oni);
    let horseshoe = (-0.5 * ((lon360 - 185.0) / 34.0).powi(2)).exp()
        * ((-0.5 * ((lat_deg - 24.0) / 9.0).powi(2)).exp()
            + (-0.5 * ((lat_deg + 24.0) / 9.0).powi(2)).exp())
        * (-0.34 * oni);

    // 5. Western Indian Ocean warming teleconnection (45E..85E)
    let indian = (-0.5 * ((lon360 - 62.0) / 22.0).powi(2)).exp()
        * (-0.5 * (lat_deg / 13.0).powi(2)).exp()
        * (0.24 * oni);

    (core_tongue + peru_coast + ne_pacific_kelvin + wp_cool + horseshoe + indian).clamp(-3.8, 4.2)
}

fn wrap_lon_dist(a_deg: f64, b_deg: f64) -> f64 {
    let mut d = (a_deg - b_deg).rem_euclid(360.0);
    if d > 180.0 {
        d -= 360.0;
    }
    d
}

/// Evaluates the historically differentiated spatial SST anomaly (`Δ°C`) for a specific
/// `(year, season_index, oni)` combining EP vs. CP/Modoki flavor, Peru coastal upwelling vs.
/// NE Pacific Warm Blob (2014–2016), North Atlantic heatwave (2023–2024), PDO regime, and warming trend.
pub fn spatial_sst_anomaly_historical(
    lon_deg: f64,
    lat_deg: f64,
    year: i32,
    season_index: usize,
    oni: f64,
) -> f64 {
    let lon360 = lon_deg.rem_euclid(360.0);
    let flavor = classify_historical_flavor(year, oni);
    let season_angle = ((season_index.min(11) as f64) / 12.0) * std::f64::consts::TAU;
    let season_factor = 1.0 + 0.08 * season_angle.cos();

    let (core_center, peru_scale, wp_scale, blob_amp, natl_amp) = match flavor {
        HistoricalEnsoFlavor::CanonicalEasternPacific => (246.0, 1.28, -0.40, 0.0, 0.0),
        HistoricalEnsoFlavor::CentralPacificModoki => (198.0, 0.22, -0.30, 0.0, 0.0),
        HistoricalEnsoFlavor::MixedSuperWithWarmBlob => (224.0, 0.85, -0.42, 2.2, 0.0),
        HistoricalEnsoFlavor::GlobalHeatwaveElNino => (238.0, 1.12, -0.34, 0.45, 1.45),
        HistoricalEnsoFlavor::DevelopingSuperElNino => (236.0, 1.08, -0.38, 0.35, 0.55),
        HistoricalEnsoFlavor::LaNinaColdTongue => (228.0, 0.96, -0.42, 0.0, 0.0),
        HistoricalEnsoFlavor::NeutralBaseline => (232.0, 0.92, -0.36, 0.0, 0.0),
    };

    let mut total = 0.0;

    if oni.abs() >= 1.0e-4 {
        let eq_lat = (-0.5 * (lat_deg / 6.8).powi(2)).exp();
        let eq_lon = if (140.0..=286.0).contains(&lon360) {
            let width = if lon360 < core_center { 42.0 } else { 34.0 };
            (-0.5 * ((lon360 - core_center) / width).powi(2)).exp()
        } else {
            0.0
        };
        total += 1.16 * oni * season_factor * eq_lat * eq_lon;

        // Peruvian / Ecuadorian coastal upwelling anomaly (canonical EP strong vs Modoki muted)
        total += (-0.5 * ((lon360 - 277.0) / 7.5).powi(2)).exp()
            * (-0.5 * ((lat_deg + 6.5) / 8.5).powi(2)).exp()
            * (peru_scale * oni);

        // Coastally trapped Northeast Pacific Kelvin wave corridor
        if (12.0..=56.0).contains(&lat_deg) {
            let coast_lon = 255.0 - 0.58 * (lat_deg - 18.0);
            let d_lon = lon360 - coast_lon;
            if (-18.0..=4.0).contains(&d_lon) {
                total += (-0.5 * (d_lon / 6.5).powi(2)).exp()
                    * (-0.5 * ((lat_deg - 34.0) / 15.0).powi(2)).exp()
                    * (0.56 * oni);
            }
        }

        // Western Pacific warm-pool / Modoki flanking counter-anomaly & subtropical horseshoe
        total += (-0.5 * ((lon360 - 138.0) / 20.0).powi(2)).exp()
            * (-0.5 * (lat_deg / 14.0).powi(2)).exp()
            * (wp_scale * oni);

        let hs_lon = if flavor == HistoricalEnsoFlavor::CentralPacificModoki {
            174.0
        } else {
            188.0
        };
        total += (-0.5 * ((lon360 - hs_lon) / 34.0).powi(2)).exp()
            * ((-0.5 * ((lat_deg - 24.0) / 9.0).powi(2)).exp()
                + (-0.5 * ((lat_deg + 24.0) / 9.0).powi(2)).exp())
            * (-0.32 * oni);

        // Indian Ocean dipole / basin-wide El Niño lag warming
        total += (-0.5 * ((lon360 - 62.0) / 22.0).powi(2)).exp()
            * (-0.5 * (lat_deg / 13.0).powi(2)).exp()
            * (0.25 * oni);
    }

    // NE Pacific Marine Heatwave ("The Blob" in 2014–2016)
    let effective_blob = if blob_amp > 0.0 {
        blob_amp
    } else if (2014..=2016).contains(&year) {
        1.85
    } else {
        0.0
    };
    if effective_blob > 0.0 {
        let dl = wrap_lon_dist(lon360, 218.0); // 142W
        total += effective_blob
            * (-0.5 * (dl / 22.0).powi(2)).exp()
            * (-0.5 * ((lat_deg - 44.0) / 9.5).powi(2)).exp();
    }

    // North Atlantic Record Marine Heatwave (2023–2024)
    let effective_natl = if natl_amp > 0.0 {
        natl_amp
    } else if (2023..=2024).contains(&year) {
        1.35
    } else {
        0.0
    };
    if effective_natl > 0.0 {
        let dl = wrap_lon_dist(lon360, 318.0); // 42W
        total += effective_natl
            * (-0.5 * (dl / 26.0).powi(2)).exp()
            * (-0.5 * ((lat_deg - 33.0) / 12.0).powi(2)).exp();
    }

    // Pacific Decadal Oscillation (PDO) multi-decadal North Pacific background shift
    let pdo_sign = if year <= 1976 || (1999..=2013).contains(&year) || (2020..=2022).contains(&year)
    {
        -0.42
    } else if (1977..=1998).contains(&year) {
        0.42
    } else {
        0.18
    };
    let pdo_dl = wrap_lon_dist(lon360, 195.0);
    total += pdo_sign
        * (-0.5 * (pdo_dl / 28.0).powi(2)).exp()
        * (-0.5 * ((lat_deg - 38.0) / 10.0).powi(2)).exp();

    // Secular background ocean warming trend (1950–2026) relative to mid-century baseline
    let t_decade = ((f64::from(year.clamp(START_YEAR, END_YEAR)) - 1950.0) / 76.0).clamp(0.0, 1.0);
    let secular = (-0.28 + 0.76 * t_decade) * (1.0 - 0.35 * (lat_deg.abs() / 90.0));
    total += secular;

    total.clamp(-4.2, 4.6)
}

/// Evaluates seasonal climatological SST (°C) with Western Pacific Warm Pool & eastern-boundary cold upwelling.
pub fn climatological_sst_c(lon_deg: f64, lat_deg: f64, season_index: usize) -> f64 {
    let lon360 = lon_deg.rem_euclid(360.0);
    // Thermal equator lags astronomical solstice by ~1 season index (peaks JFM south, JAS north)
    let decl_angle = (((season_index.min(11) as f64) - 1.0) / 12.0) * std::f64::consts::TAU;
    let lat_shift = -6.5 * decl_angle.cos();
    let eff_lat = (lat_deg - lat_shift).clamp(-89.0, 89.0);
    let mut sst = 28.5 * eff_lat.to_radians().cos().max(0.0).powf(1.35) - 1.6;

    // Western Pacific Warm Pool & Indo-Pacific Warm Pool
    let wp_dl = wrap_lon_dist(lon360, 140.0);
    sst += 2.4 * (-0.5 * (wp_dl / 32.0).powi(2)).exp() * (-0.5 * (lat_deg / 16.0).powi(2)).exp();

    // Cold eastern-boundary upwelling tongues (Humboldt, Benguela, California)
    let humboldt_dl = wrap_lon_dist(lon360, 278.0);
    sst -= 2.8
        * (-0.5 * (humboldt_dl / 14.0).powi(2)).exp()
        * (-0.5 * ((lat_deg + 10.0) / 14.0).powi(2)).exp();

    let benguela_dl = wrap_lon_dist(lon360, 12.0);
    sst -= 2.4
        * (-0.5 * (benguela_dl / 12.0).powi(2)).exp()
        * (-0.5 * ((lat_deg + 18.0) / 13.0).powi(2)).exp();

    let cal_dl = wrap_lon_dist(lon360, 236.0);
    sst -= 2.2
        * (-0.5 * (cal_dl / 14.0).powi(2)).exp()
        * (-0.5 * ((lat_deg - 32.0) / 12.0).powi(2)).exp();

    sst.clamp(-1.9, 33.5)
}

/// Computes the equatorial Pacific `20°C` thermocline isotherm depth (`m`, positive downward)
/// from `120°E` (`lon360 = 120`) to `80°W` (`lon360 = 280`).
///
/// Grounded in NOAA PMEL TAO/TRITON observations:
/// - Neutral (`ONI = 0`): slopes from `~165 m` in the Western Pacific (`140°E`) to `~48 m` off Ecuador/Peru (`85°W`).
/// - Super El Niño (`ONI ~ +2.4..+2.6°C`): flattens across the basin (`~115 m` west to `~150 m` east).
/// - Strong La Niña (`ONI ~ -1.8°C`): steepens (`~195 m` west to `~22 m` shoaled upwelling in the east).
pub fn equatorial_20c_isotherm_depth_m(lon360_deg: f64, oni: f64) -> f64 {
    let t = ((lon360_deg - 120.0) / 160.0).clamp(0.0, 1.0); // 0 at 120E, 1 at 80W
    let neutral_depth = 172.0 - 124.0 * t.powf(1.15);
    let enso_tilt_shift = oni * (-22.0 * (1.0 - t) + 41.0 * t.powf(1.25));
    (neutral_depth + enso_tilt_shift).clamp(16.0, 225.0)
}

/// Computes the historical episode-differentiated equatorial Pacific `20°C` isotherm depth (`m`).
pub fn equatorial_20c_isotherm_depth_historical_m(lon360_deg: f64, year: i32, oni: f64) -> f64 {
    let t = ((lon360_deg - 120.0) / 160.0).clamp(0.0, 1.0);
    let neutral_depth = 172.0 - 124.0 * t.powf(1.15);
    let shift = match classify_historical_flavor(year, oni) {
        HistoricalEnsoFlavor::CanonicalEasternPacific => {
            oni * (-24.0 * (1.0 - t) + 46.0 * t.powf(1.25))
        }
        HistoricalEnsoFlavor::CentralPacificModoki => {
            let mid_bell = (-0.5 * ((t - 0.45) / 0.22).powi(2)).exp();
            oni * (-14.0 * (1.0 - t) + 22.0 * t.powf(1.25) + 20.0 * mid_bell)
        }
        _ => oni * (-22.0 * (1.0 - t) + 41.0 * t.powf(1.25)),
    };
    (neutral_depth + shift).clamp(16.0, 225.0)
}

/// Evaluates the tropical Pacific `20°C` thermocline isotherm depth (`m`) inside the
/// tropical Pacific basin (`|lat| <= 20°`, `120°E <= lon360 <= 280°E` / `80°W`).
///
/// Returns `None` outside the tropical Pacific (polar latitudes, Atlantic, Indian Ocean)
/// so global maps and cursor telemetry do not extrapolate equatorial Pacific thermocline tilts
/// onto extratropical or non-Pacific basins.
pub fn tropical_pacific_20c_isotherm_depth_m(
    lon_deg: f64,
    lat_deg: f64,
    year: i32,
    oni: f64,
) -> Option<f64> {
    if lat_deg.abs() > 20.0 {
        return None;
    }
    let lon360 = lon_deg.rem_euclid(360.0);
    if !(120.0..=280.0).contains(&lon360) {
        return None;
    }
    Some(equatorial_20c_isotherm_depth_historical_m(lon360, year, oni))
}

pub fn export_oni_csv() -> String {
    let mut out = String::from("year,season,oni_celsius,phase\n");
    for (yi, row) in ONI_TABLE_1950_2026.iter().enumerate() {
        let year = START_YEAR + yi as i32;
        for (si, &val) in row.iter().enumerate() {
            if !val.is_nan() {
                let phase = EnsoPhase::from_oni(f64::from(val)).label();
                out.push_str(&format!(
                    "{year},{},{val:.1},\"{phase}\"\n",
                    SEASONS[si]
                ));
            }
        }
    }
    out
}

/// Comparative ENSO Niño-3.4, eastern-Pacific upwelling, and 20°C thermocline anomaly
/// metrics for a historical El Niño / La Niña / Neutral preset episode.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EnsoDigestRow {
    pub slug: String,
    pub title: String,
    pub year: i32,
    pub season: String,
    pub oni_celsius: f64,
    pub phase: String,
    pub flavor: String,
    pub nino34_sst_anom_c: f64,
    pub nino12_peru_sst_anom_c: f64,
    pub dateline_sst_anom_c: f64,
    pub west_pac_thermocline_m: f64,
    pub east_pac_thermocline_m: f64,
    pub thermocline_slope_m: f64,
    pub delta_east_thermocline_vs_neutral_m: f64,
}

/// Comparative digest across all landmark historical ENSO episodes in `ENSO_PRESETS`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EnsoComparativeDigest {
    pub baseline_slug: String,
    pub baseline_east_thermocline_m: f64,
    pub rows: Vec<EnsoDigestRow>,
}

/// Builds a comparative Niño-3.4 / thermocline anomaly digest across all landmark
/// historical ENSO episodes (`ENSO_PRESETS`), anchored against the `neutral-1990` baseline.
pub fn build_enso_comparative_digest() -> EnsoComparativeDigest {
    let neutral = ENSO_PRESETS
        .iter()
        .find(|p| p.slug == "neutral-1990")
        .unwrap_or(&ENSO_PRESETS[ENSO_PRESETS.len() - 1]);
    let neutral_east_m =
        equatorial_20c_isotherm_depth_historical_m(275.0, neutral.year, neutral.oni);

    let rows = ENSO_PRESETS
        .iter()
        .map(|p| {
            let nino34 =
                spatial_sst_anomaly_historical(-145.0, 0.0, p.year, p.season_index, p.oni);
            let nino12 =
                spatial_sst_anomaly_historical(-82.0, -5.0, p.year, p.season_index, p.oni);
            let dateline =
                spatial_sst_anomaly_historical(175.0, 0.0, p.year, p.season_index, p.oni);
            let west_m = equatorial_20c_isotherm_depth_historical_m(150.0, p.year, p.oni);
            let east_m = equatorial_20c_isotherm_depth_historical_m(275.0, p.year, p.oni);
            EnsoDigestRow {
                slug: p.slug.to_string(),
                title: p.title.to_string(),
                year: p.year,
                season: SEASONS[p.season_index.min(11)].to_string(),
                oni_celsius: p.oni,
                phase: EnsoPhase::from_oni(p.oni).label().to_string(),
                flavor: classify_historical_flavor(p.year, p.oni).label().to_string(),
                nino34_sst_anom_c: nino34,
                nino12_peru_sst_anom_c: nino12,
                dateline_sst_anom_c: dateline,
                west_pac_thermocline_m: west_m,
                east_pac_thermocline_m: east_m,
                thermocline_slope_m: west_m - east_m,
                delta_east_thermocline_vs_neutral_m: east_m - neutral_east_m,
            }
        })
        .collect();

    EnsoComparativeDigest {
        baseline_slug: neutral.slug.to_string(),
        baseline_east_thermocline_m: neutral_east_m,
        rows,
    }
}

impl EnsoComparativeDigest {
    pub fn to_csv(&self) -> String {
        let mut out = String::from(
            "slug,year,season,oni_celsius,phase,flavor,nino34_sst_anom_c,nino12_peru_sst_anom_c,dateline_sst_anom_c,west_pac_thermocline_m,east_pac_thermocline_m,thermocline_slope_m,delta_east_thermocline_vs_neutral_m\n",
        );
        for r in &self.rows {
            out.push_str(&format!(
                "{},{},{},{:+.1},\"{}\",\"{}\",{:.3},{:.3},{:.3},{:.2},{:.2},{:.2},{:+.2}\n",
                r.slug,
                r.year,
                r.season,
                r.oni_celsius,
                r.phase,
                r.flavor,
                r.nino34_sst_anom_c,
                r.nino12_peru_sst_anom_c,
                r.dateline_sst_anom_c,
                r.west_pac_thermocline_m,
                r.east_pac_thermocline_m,
                r.thermocline_slope_m,
                r.delta_east_thermocline_vs_neutral_m,
            ));
        }
        out
    }

    pub fn to_json_pretty(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn verifies_landmark_noaa_cpc_oni_peaks() {
        assert_eq!(lookup_oni(2015, 11), Some(2.6));
        assert_eq!(lookup_oni(1997, 11), Some(2.4));
        assert_eq!(lookup_oni(2023, 11), Some(2.0));
        assert_eq!(lookup_oni(2026, 8), Some(2.5));
        assert_eq!(lookup_oni(1973, 11), Some(-2.0));
        assert_eq!(lookup_oni(2026, 11), None);
    }

    #[test]
    fn el_nino_flattens_eastern_pacific_thermocline_toward_150m() {
        let neutral_east = equatorial_20c_isotherm_depth_m(275.0, 0.0);
        let super_el_nino_east = equatorial_20c_isotherm_depth_m(275.0, 2.4);
        assert!(neutral_east < 60.0);
        assert!((super_el_nino_east - 145.0).abs() < 15.0);
    }

    #[test]
    fn clamps_2026_season_to_latest_published_aso() {
        assert_eq!(latest_valid_season_for_year(2026), 8);
        assert_eq!(latest_valid_season_for_year(2015), 11);
        assert_eq!(lookup_oni_clamped(2026, 11), (8, 2.5));
        assert_eq!(lookup_oni_clamped(1973, 11), (11, -2.0));
    }

    #[test]
    fn differentiates_historical_enso_flavors_and_regional_heatwaves() {
        // 1997-98 Canonical EP El Nino is substantially warmer off Peru than 2009-10 Modoki
        let anom_peru_1997 = spatial_sst_anomaly_historical(-82.0, -5.0, 1997, 11, 2.4);
        let anom_peru_2009 = spatial_sst_anomaly_historical(-82.0, -5.0, 2009, 11, 1.5);
        assert!(
            anom_peru_1997 > anom_peru_2009 + 1.5,
            "Expected 1997 EP off Peru ({anom_peru_1997:.2}) >> 2009 Modoki ({anom_peru_2009:.2})"
        );

        // 2009-10 Modoki peaks further west near the Dateline than 1997-98 EP
        let anom_dl_2009 = spatial_sst_anomaly_historical(175.0, 0.0, 2009, 11, 1.5);
        let anom_dl_1997 = spatial_sst_anomaly_historical(175.0, 0.0, 1997, 11, 2.4);
        assert!(
            anom_dl_2009 > anom_dl_1997 + 0.3,
            "Expected 2009 Modoki at 175E ({anom_dl_2009:.2}) > 1997 EP ({anom_dl_1997:.2})"
        );

        // 2015-16 exhibits a strong NE Pacific Warm Blob compared to 1997-98
        let blob_2015 = spatial_sst_anomaly_historical(-142.0, 44.0, 2015, 11, 2.6);
        let blob_1997 = spatial_sst_anomaly_historical(-142.0, 44.0, 1997, 11, 2.4);
        assert!(
            blob_2015 > blob_1997 + 1.4,
            "Expected 2015 NE Pacific Blob ({blob_2015:.2}) > 1997 ({blob_1997:.2})"
        );

        // 2023-24 exhibits North Atlantic heatwave relative to 1982-83
        let natl_2023 = spatial_sst_anomaly_historical(-42.0, 33.0, 2023, 11, 2.0);
        let natl_1982 = spatial_sst_anomaly_historical(-42.0, 33.0, 1982, 11, 2.1);
        assert!(
            natl_2023 > natl_1982 + 1.0,
            "Expected 2023 N. Atlantic heatwave ({natl_2023:.2}) > 1982 ({natl_1982:.2})"
        );

        // 1973-74 Record La Nina produces strong negative anomaly in Nino 3.4
        let nino34_1973 = spatial_sst_anomaly_historical(-145.0, 0.0, 1973, 11, -2.0);
        assert!(
            nino34_1973 < -1.8,
            "Expected 1973 La Nina Nino3.4 ({nino34_1973:.2}) < -1.8"
        );
    }

    #[test]
    fn masks_thermocline_depth_outside_tropical_pacific_and_validates_neutral_1990_preset() {
        assert_eq!(tropical_pacific_20c_isotherm_depth_m(-140.0, 0.0, 1997, 2.4).is_some(), true);
        assert_eq!(tropical_pacific_20c_isotherm_depth_m(-140.0, 65.0, 1997, 2.4), None);
        assert_eq!(tropical_pacific_20c_isotherm_depth_m(-140.0, -70.0, 1997, 2.4), None);
        assert_eq!(tropical_pacific_20c_isotherm_depth_m(-40.0, 0.0, 1997, 2.4), None);
        assert_eq!(tropical_pacific_20c_isotherm_depth_m(75.0, 0.0, 1997, 2.4), None);

        let neutral = ENSO_PRESETS
            .iter()
            .find(|p| p.slug == "neutral-1990")
            .expect("neutral-1990 preset exists");
        assert_eq!(neutral.year, 1990);
        assert_eq!(lookup_oni(neutral.year, neutral.season_index), Some(neutral.oni));
    }

    #[test]
    fn comparative_digest_covers_all_presets_with_thermocline_flattening() {
        let digest = build_enso_comparative_digest();
        assert_eq!(digest.baseline_slug, "neutral-1990");
        assert_eq!(digest.rows.len(), ENSO_PRESETS.len());

        let ep1997 = digest
            .rows
            .iter()
            .find(|r| r.slug == "super-1997")
            .expect("super-1997 row");
        let ln1973 = digest
            .rows
            .iter()
            .find(|r| r.slug == "lanina-1973")
            .expect("lanina-1973 row");

        // Strong El Niño deepens eastern Pacific 20°C isotherm & flattens west-east slope
        assert!(ep1997.delta_east_thermocline_vs_neutral_m > 50.0);
        assert!(ln1973.delta_east_thermocline_vs_neutral_m < -5.0);
        assert!(ln1973.thermocline_slope_m > ep1997.thermocline_slope_m);
        assert!(ep1997.nino34_sst_anom_c > 2.0);
        assert!(ln1973.nino34_sst_anom_c < -1.8);
    }

    #[test]
    fn comparative_digest_exports_valid_csv_and_json() {
        let digest = build_enso_comparative_digest();
        let csv = digest.to_csv();
        assert!(csv.starts_with("slug,year,season,oni_celsius,"));
        assert_eq!(csv.lines().count(), ENSO_PRESETS.len() + 1);
        assert!(csv.contains("modoki-2009"));

        let json_raw = digest.to_json_pretty().expect("serialize ENSO digest");
        let roundtrip: EnsoComparativeDigest =
            serde_json::from_str(&json_raw).expect("deserialize ENSO digest");
        assert_eq!(roundtrip.rows.len(), ENSO_PRESETS.len());
    }

    #[test]
    fn enso_phase_classification_boundaries_are_consistent() {
        assert_eq!(EnsoPhase::from_oni(2.1), EnsoPhase::SuperElNino);
        assert_eq!(EnsoPhase::from_oni(1.6), EnsoPhase::StrongElNino);
        assert_eq!(EnsoPhase::from_oni(1.1), EnsoPhase::ModerateElNino);
        assert_eq!(EnsoPhase::from_oni(0.6), EnsoPhase::WeakElNino);
        assert_eq!(EnsoPhase::from_oni(0.0), EnsoPhase::Neutral);
        assert_eq!(EnsoPhase::from_oni(-0.7), EnsoPhase::WeakLaNina);
        assert_eq!(EnsoPhase::from_oni(-1.2), EnsoPhase::ModerateLaNina);
        assert_eq!(EnsoPhase::from_oni(-1.9), EnsoPhase::StrongLaNina);
    }
}
