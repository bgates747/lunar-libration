#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Cache a longer JPL Horizons interval for a complete orbit starting at perigee.

Run this explicitly to refresh the cache; the renderer never accesses the
network. The original 28-day model-validation reference is left unchanged.
"""
import csv
from datetime import datetime, timezone
import io
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError
import time
import sys

# Resolve the sibling source modules from this repository, independent of cwd.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from project_paths import DATA_DIR, PROJECT_ROOT


def fetch(kind, parameters):
    url = "https://ssd.jpl.nasa.gov/api/horizons.api?" + urlencode(parameters)
    for attempt in range(3):
        try:
            with urlopen(Request(url, headers={"User-Agent": "lunar-libration/1.0"}), timeout=45) as response:
                payload = json.load(response)
            break
        except (HTTPError, URLError) as exc:
            if isinstance(exc, HTTPError) and exc.code not in (429, 500, 502, 503, 504):
                raise
            if attempt == 2:
                raise
            print(f"Retrying {kind} after a temporary Horizons response", flush=True)
            time.sleep(2*(attempt+1))
    if payload.get("error"):
        raise RuntimeError(payload["error"])
    result = payload["result"]
    if "$$SOE" not in result or "$$EOE" not in result:
        raise RuntimeError(f"No ephemeris table returned for {kind}")
    header, table = result.split("$$SOE", 1)
    table = table.split("$$EOE", 1)[0].strip()
    rows = list(csv.reader(io.StringIO(table)))
    return kind, {"request_parameters": parameters, "request_url": url,
                  "signature": payload.get("signature"), "header": header, "rows": rows}


def main():
    original = json.loads((DATA_DIR / "moon-model-reference.json").read_text())
    jobs = []
    for kind in ("moon", "sun", "orientation"):
        parameters = dict(original[kind]["request_parameters"])
        parameters["START_TIME"] = "'2026-09-01 00:00'"
        parameters["STOP_TIME"] = "'2026-10-10 00:00'"
        jobs.append((kind, parameters))
    # Keep API requests sequential so we do not compete with ourselves for
    # Horizons service capacity. Temporary service failures get bounded retries.
    tables = {}
    for job in jobs:
        kind, table = fetch(*job)
        tables[kind] = table
        print(f"Fetched {kind}: {len(table['rows'])} epochs", flush=True)
    counts = {kind: len(table["rows"]) for kind, table in tables.items()}
    if len(set(counts.values())) != 1 or counts["moon"] < 900:
        raise RuntimeError(f"Unexpected ephemeris row counts: {counts}")
    record = {
        "source": "NASA/JPL Horizons",
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "Animation of a complete geometric revolution starting at the first available perigee",
        "coordinate_frame": "ICRF equatorial axes; geometric Earth-centred position and velocity vectors",
        "orientation_time_note": original["orientation_time_note"],
        "sampling": "Hourly September 1–October 10, 2026; renderer interpolates and locates radial extrema",
        **tables,
    }
    path = (DATA_DIR / "animation-ephemeris.json")
    path.write_text(json.dumps(record, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps({"saved": path.relative_to(PROJECT_ROOT).as_posix(), "row_counts": counts,
                      "retrieved_at_utc": record["retrieved_at_utc"]}))


if __name__ == "__main__":
    main()
