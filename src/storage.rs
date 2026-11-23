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

//! Local SQLite (`rusqlite` WAL) bookmark persistence and JSON snapshot exchange.

use std::path::Path;

use rusqlite::{params, Connection, Result};

use crate::SimulationSnapshot;

pub struct BookmarkStore {
    conn: Connection,
}

impl BookmarkStore {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let conn = Connection::open(path)?;
        conn.pragma_update(None, "journal_mode", "WAL")?;
        conn.execute_batch(
            "CREATE TABLE IF NOT EXISTS saved_views (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                name TEXT NOT NULL,
                created_utc TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
                snapshot_json TEXT NOT NULL
            );",
        )?;
        Ok(Self { conn })
    }

    pub fn save_bookmark(&self, name: &str, snapshot: &SimulationSnapshot) -> Result<()> {
        let json = serde_json::to_string(snapshot).unwrap_or_default();
        self.conn.execute(
            "INSERT INTO saved_views (name, snapshot_json) VALUES (?1, ?2)",
            params![name, json],
        )?;
        Ok(())
    }

    pub fn list_bookmarks(&self) -> Result<Vec<(i64, String, SimulationSnapshot)>> {
        let mut stmt = self
            .conn
            .prepare("SELECT id, name, snapshot_json FROM saved_views ORDER BY id DESC LIMIT 32")?;
        let rows = stmt.query_map([], |row| {
            let id: i64 = row.get(0)?;
            let name: String = row.get(1)?;
            let raw: String = row.get(2)?;
            let snap: SimulationSnapshot = serde_json::from_str(&raw).unwrap_or_default();
            Ok((id, name, snap))
        })?;
        rows.collect()
    }

    pub fn delete_bookmark(&self, id: i64) -> Result<usize> {
        self.conn
            .execute("DELETE FROM saved_views WHERE id = ?1", params![id])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn saves_lists_and_deletes_bookmarks() {
        let store = BookmarkStore::open(":memory:").expect("open in-memory SQLite store");
        let snap = SimulationSnapshot::default();
        store.save_bookmark("2026 ASO Test", &snap).unwrap();
        let list = store.list_bookmarks().unwrap();
        assert_eq!(list.len(), 1);
        assert_eq!(list[0].1, "2026 ASO Test");
        let deleted = store.delete_bookmark(list[0].0).unwrap();
        assert_eq!(deleted, 1);
        assert!(store.list_bookmarks().unwrap().is_empty());
    }
}

impl BookmarkStore {
    pub fn ensure_raster_cache_table(&self) -> Result<()> {
        self.conn.execute_batch(
            "CREATE TABLE IF NOT EXISTS raster_tile_cache (
                cache_key TEXT PRIMARY KEY,
                fetched_utc TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
                raster_json TEXT NOT NULL
            );",
        )
    }

    pub fn put_cached_raster(&self, key: &str, grid: &crate::raster::SatelliteRasterGrid) -> Result<()> {
        self.ensure_raster_cache_table()?;
        let payload = serde_json::to_string(grid).unwrap_or_default();
        self.conn.execute(
            "INSERT OR REPLACE INTO raster_tile_cache (cache_key, raster_json) VALUES (?1, ?2)",
            params![key, payload],
        )?;
        Ok(())
    }

    pub fn get_cached_raster(&self, key: &str) -> Result<Option<crate::raster::SatelliteRasterGrid>> {
        self.ensure_raster_cache_table()?;
        let mut stmt = self
            .conn
            .prepare("SELECT raster_json FROM raster_tile_cache WHERE cache_key = ?1")?;
        let mut rows = stmt.query(params![key])?;
        if let Some(row) = rows.next()? {
            let raw: String = row.get(0)?;
            Ok(serde_json::from_str(&raw).ok())
        } else {
            Ok(None)
        }
    }
}
