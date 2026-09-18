"""
Local cache of Onshape STEP pulls, keyed by the full 5-parameter
configuration (emag_height, inner_coil_id, inner_coil_od, outer_coil_id,
outer_coil_od).

All 5 are real Onshape Configuration parameters that cascade through the
assembly's mates -- confirmed via a live GET on .../configuration, which
lists exactly these 5 and no others. Earlier versions of this module sent
only emag_height and left the 4 coil ID/OD params at their Onshape
defaults, assuming they were purely local resizes with nothing to gain
from going through the API -- that assumption was wrong: Onshape models a
real coil solid per winding, and inner_coil_id/outer_coil_id also size the
soft-iron core each coil wraps (there's no separate core-diameter
parameter -- the core IS sized by the coil ID it sits inside). So every
geometry axis in the sweep now goes through Onshape, and the coil-sleeve
geometry Gmsh used to build itself (occ.addCylinder/occ.cut) comes from
the STEP import instead -- see gmsh_getdp/assembly/build_assembly.py.

The multi-parameter configuration string is built by Onshape's own
`configurationencodings` endpoint (POST .../elements/d/{did}/e/{eid}/
configurationencodings, body {"parameters": [{"parameterId":...,
"parameterValue":...}, ...]}) rather than hand-joining "key=value" pairs --
confirmed live that multi-param `encodedId`s are ";"-separated
"parameterId=value" pairs (e.g. "emag_height=1.5in;inner_coil_id=0.9in"),
but there's no reason to duplicate Onshape's own encoding logic once the
endpoint exists to do it (and no `wid` needed for it -- it's a stateless
string-encoding call, unlike the workspace-scoped translation below).

Caches by document microversion too: if the Onshape document is edited
again later, a previously-cached configuration is treated as stale and
re-pulled, rather than silently serving geometry from before the edit.
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

# The 5 real Onshape Configuration parameters (all LENGTH/inch quantities),
# confirmed via GET /elements/d/{DID}/w/{WID}/e/{EID}/configuration.
PARAM_IDS = ["emag_height", "inner_coil_id", "inner_coil_od", "outer_coil_id", "outer_coil_od"]

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


def _encode_configuration(client, params_in: dict) -> str:
    """params_in: {parameter_id: value_in_inches}. Returns Onshape's own
    ";"-joined encodedId string for the full param set (see module
    docstring for why this goes through Onshape's endpoint rather than
    being hand-built)."""
    body = {
        "parameters": [
            {"parameterId": pid, "parameterValue": f"{params_in[pid]:g}in"}
            for pid in PARAM_IDS
        ]
    }
    resp = client.post(f"/elements/d/{DID}/e/{EID}/configurationencodings", json=body)
    resp.raise_for_status()
    return resp.json()["encodedId"]


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


def get_step_for_config(params: dict, force_refresh: bool = False) -> Path:
    """
    Local path to a STEP export of the assembly at the given geometry
    configuration, pulling from Onshape only if not already cached for the
    CURRENT document microversion.

    params: dict covering all 5 of PARAM_IDS (inches). Every key must be
    present -- there's no "leave at Onshape default" case anymore now that
    all 5 params are meaningful (see module docstring).
    """
    missing = set(PARAM_IDS) - set(params)
    if missing:
        raise ValueError(f"get_step_for_config: missing params {sorted(missing)}")

    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(DB_PATH)
    _init_db(conn)

    with get_client() as client:
        configuration = _encode_configuration(client, params)
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
