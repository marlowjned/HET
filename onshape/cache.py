"""
Local cache of Onshape STEP pulls, keyed by emag_height.

Only emag_height ever needs to go through Onshape: it's the one sweep
parameter that cascades through Onshape's mates (repositions the plates/
chamber_spacer/etc. -- confirmed empirically, see onshape/cad_cache/
validation in this module's git history). inner_coil_id/od and
outer_coil_id/od are left at their Onshape defaults on every pull and
resized locally in Gmsh instead (purely local resizes, no repositioning),
so they're never part of this cache's key -- see gmsh_getdp/README.md for
that half of the pipeline once it exists.

Caches by document microversion too: if the Onshape document is edited
again later, a previously-cached height is treated as stale and re-pulled,
rather than silently serving geometry from before the edit.
"""
import sqlite3
import time
from pathlib import Path

from .client import get_client

CACHE_DIR = Path(__file__).resolve().parent / "cad_cache"
DB_PATH = CACHE_DIR / "index.db"

DID = "c41ea30140eac94bca56baf9"
WID = "f9e3107b7398b0d9cf6622d0"
EID = "82b7a1bc6b1233ce1a780568"  # BPL-700 Assembly element

POLL_INTERVAL_S = 2.0
POLL_TIMEOUT_S = 120.0


def _init_db(conn: sqlite3.Connection) -> None:
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS pulls (
            configuration TEXT PRIMARY KEY,
            document_microversion TEXT NOT NULL,
            step_file TEXT NOT NULL,
            fetched_at REAL NOT NULL,
            translation_id TEXT
        )
        """
    )
    conn.commit()


def _current_microversion(client) -> str:
    resp = client.get(f"/elements/d/{DID}/w/{WID}/e/{EID}/configuration")
    resp.raise_for_status()
    return resp.json()["sourceMicroversion"]


def _pull_from_onshape(client, configuration: str) -> tuple[Path, str]:
    body = {
        "formatName": "STEP",
        "configuration": configuration,
        "destinationName": f"sweep_{configuration}",
        "storeInDocument": False,
        "notifyUser": False,
        "grouping": True,
    }
    resp = client.post(f"/assemblies/d/{DID}/w/{WID}/e/{EID}/translations", json=body)
    resp.raise_for_status()
    data = resp.json()
    translation_id = data["id"]

    elapsed = 0.0
    while data["requestState"] not in ("DONE", "FAILED"):
        if elapsed > POLL_TIMEOUT_S:
            raise TimeoutError(f"Onshape translation {translation_id} did not finish within {POLL_TIMEOUT_S}s")
        time.sleep(POLL_INTERVAL_S)
        elapsed += POLL_INTERVAL_S
        data = client.get(f"/translations/{translation_id}").json()

    if data["requestState"] == "FAILED":
        raise RuntimeError(f"Onshape translation failed: {data.get('failureReason')}")

    fid = data["resultExternalDataIds"][0]
    dl = client.get(f"/documents/d/{DID}/externaldata/{fid}")
    dl.raise_for_status()

    safe_name = configuration.replace("/", "_").replace(" ", "")
    step_path = CACHE_DIR / f"{safe_name}.step"
    step_path.write_bytes(dl.content)
    return step_path, translation_id


def get_step_for_height(emag_height_in: float, force_refresh: bool = False) -> Path:
    """
    Local path to a STEP export of the assembly at the given emag_height
    (inches), pulling from Onshape only if not already cached for the
    CURRENT document microversion.
    """
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    configuration = f"emag_height={emag_height_in:g}in"

    conn = sqlite3.connect(DB_PATH)
    _init_db(conn)

    with get_client() as client:
        microversion = _current_microversion(client)

        if not force_refresh:
            row = conn.execute(
                "SELECT step_file, document_microversion FROM pulls WHERE configuration = ?",
                (configuration,),
            ).fetchone()
            if row is not None:
                step_file, cached_microversion = row
                step_path = Path(step_file)
                if cached_microversion == microversion and step_path.exists():
                    conn.close()
                    return step_path
                # stale (document changed since this was cached) -- fall through and re-pull

        step_path, translation_id = _pull_from_onshape(client, configuration)

        conn.execute(
            "INSERT OR REPLACE INTO pulls (configuration, document_microversion, step_file, fetched_at, translation_id) "
            "VALUES (?, ?, ?, ?, ?)",
            (configuration, microversion, str(step_path), time.time(), translation_id),
        )
        conn.commit()
        conn.close()
        return step_path
