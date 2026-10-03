"""
Local cache of Onshape STEP pulls, keyed by document + element + the full
6-parameter configuration.

RETARGETED 2026-10-02 to the BPL-1000 document. The previous target was a
different thruster (BPL-700, DID c41ea301.../WID f9e3107b.../EID
82b7a1bc...), whose assembly exposed a DIFFERENT set of 5 configuration
parameters: emag_height, inner_coil_id, inner_coil_od, outer_coil_id,
outer_coil_od. Only `emag_height` survives by name, and even its range
moved (BPL-700 swept 1.125-2.000 in; BPL-1000 declares [1.625, 10.0] with
default 1.625). The BPL-700 ids are kept in OLD_BPL700 below so an old
cached pull can still be identified, not so it can be mixed in.

The BPL-1000 parameters separate the core diameter from the winding build,
which the BPL-700 ones did not: there, the core was sized implicitly by the
coil ID it sat inside (no separate core-diameter parameter existed). Here
`iem_coreod`/`oem_coreod` set the core directly and `iem_coilt`/`oem_coilt`
set the winding's radial thickness, so coil OD ~ coreod + 2*coilt. That is
a parametrization convenience, NOT an extra degree of freedom -- the
chamber bore still bounds what the inner electromagnet's OD can be, so a
swept (coreod, coilt) pair has to be checked for fit, not assumed valid.
`iem_platet` (inner pole-plate thickness) has no BPL-700 equivalent at all.

Every geometry axis goes through Onshape, and the coil-sleeve geometry Gmsh
used to build itself (occ.addCylinder/occ.cut) comes from the STEP import
instead -- see gmsh_getdp/assembly/build_assembly.py.

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

The cache key includes the document and element ids as well. It did not
until 2026-10-02, when retargeting to BPL-1000 exposed the hole: the key
was (configuration string, microversion), so a DIFFERENT document whose
encoded configuration string happened to match would have been served the
old thruster's STEP out of cad_cache/, and every downstream solve would
have silently modelled the wrong hardware.
"""
import sqlite3
import time
from pathlib import Path

from .client import get_client

CACHE_DIR = Path(__file__).resolve().parent / "cad_cache"
DB_PATH = CACHE_DIR / "index.db"

DID = "16412e6567201641a88630bb"
WID = "50bbbc3be2fb38fd925a2589"
EID = "6bfc5ca353d5636ad82fb365"  # BPL-1000 Assembly element

# Previous target, for provenance only -- see module docstring.
OLD_BPL700 = dict(DID="c41ea30140eac94bca56baf9",
                  WID="f9e3107b7398b0d9cf6622d0",
                  EID="82b7a1bc6b1233ce1a780568")

# The 6 real Onshape Configuration parameters (all LENGTH/inch quantities),
# confirmed 2026-10-02 via GET /elements/d/{DID}/w/{WID}/e/{EID}/configuration.
# iem_ = inner electromagnet, oem_ = outer electromagnet.
PARAM_IDS = ["emag_height", "iem_coreod", "iem_coilt", "iem_platet",
             "oem_coreod", "oem_coilt"]

# Declared min/max/default per parameter, from the same live call. Onshape
# clamps out-of-range values silently rather than erroring, which would turn
# a sweep point into a duplicate of its neighbour without saying so, so
# get_step_for_config range-checks against this instead of trusting the
# caller. Kept here rather than re-fetched per call: it only changes when
# the CAD author edits the configuration, which also moves the microversion.
PARAM_RANGE_IN = {
    "emag_height": (1.625, 10.0, 1.625),
    "iem_coreod":  (0.0,   2.0,  1.0),
    "iem_coilt":   (0.0,   0.375, 0.25),
    "iem_platet":  (0.0,   0.5,  0.375),
    "oem_coreod":  (0.375, 1.25, 1.0),
    "oem_coilt":   (0.0,   1.0,  0.5),
}

DEFAULTS_IN = {k: v[2] for k, v in PARAM_RANGE_IN.items()}

POLL_INTERVAL_S = 2.0
POLL_TIMEOUT_S = 120.0


def _init_db(conn: sqlite3.Connection) -> None:
    # pulls_v2, not pulls: the v1 table was keyed on the configuration string
    # alone and its rows are all BPL-700 pulls. Rather than migrate rows whose
    # document is only knowable by inference, v2 starts clean with the document
    # and element in the key. v1 is left in place, unread, as a record.
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS pulls_v2 (
            document_id TEXT NOT NULL,
            element_id TEXT NOT NULL,
            configuration TEXT NOT NULL,
            document_microversion TEXT NOT NULL,
            step_file TEXT NOT NULL,
            fetched_at REAL NOT NULL,
            translation_id TEXT,
            PRIMARY KEY (document_id, element_id, configuration)
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

    # Document id in the filename as well as the DB key: cad_cache/ is
    # gitignored and holds BPL-700 files from before the retarget, and two
    # documents can encode the same configuration string.
    safe_name = configuration.replace("/", "_").replace(" ", "")
    step_path = CACHE_DIR / f"{DID[:8]}_{safe_name}.step"
    step_path.write_bytes(dl.content)
    return step_path, translation_id


def get_step_for_config(params: dict, force_refresh: bool = False) -> Path:
    """
    Local path to a STEP export of the assembly at the given geometry
    configuration, pulling from Onshape only if not already cached for the
    CURRENT document microversion.

    params: dict covering all 6 of PARAM_IDS (inches). Every key must be
    present -- there's no "leave at Onshape default" case, since all 6 are
    meaningful (see module docstring). Use DEFAULTS_IN to fill the axes a
    given sweep is not varying, explicitly.
    """
    missing = set(PARAM_IDS) - set(params)
    extra = set(params) - set(PARAM_IDS)
    if missing or extra:
        raise ValueError(f"get_step_for_config: missing params {sorted(missing)}, "
                         f"unknown params {sorted(extra)}")
    for pid in PARAM_IDS:
        lo, hi, _ = PARAM_RANGE_IN[pid]
        if not (lo - 1e-9 <= params[pid] <= hi + 1e-9):
            raise ValueError(
                f"get_step_for_config: {pid}={params[pid]} in is outside Onshape's "
                f"declared range [{lo}, {hi}] in. Onshape clamps silently, so this "
                f"would have returned a duplicate of the range endpoint instead.")

    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(DB_PATH)
    _init_db(conn)

    with get_client() as client:
        configuration = _encode_configuration(client, params)
        microversion = _current_microversion(client)

        if not force_refresh:
            row = conn.execute(
                "SELECT step_file, document_microversion FROM pulls_v2 "
                "WHERE document_id = ? AND element_id = ? AND configuration = ?",
                (DID, EID, configuration),
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
            "INSERT OR REPLACE INTO pulls_v2 (document_id, element_id, configuration, "
            "document_microversion, step_file, fetched_at, translation_id) "
            "VALUES (?, ?, ?, ?, ?, ?, ?)",
            (DID, EID, configuration, microversion, str(step_path), time.time(), translation_id),
        )
        conn.commit()
        conn.close()
        return step_path
