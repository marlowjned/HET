# Onshape API access

Auth layer + caching STEP-pull for scripting against the parametrized
BPL-700 CAD in Onshape (document `c41ea30140eac94bca56baf9`, workspace
`f9e3107b7398b0d9cf6622d0`, "BPL-700 Assembly" element
`82b7a1bc6b1233ce1a780568`) -- generating a real, correctly-assembled STEP
export for a given `emag_height` without going through the Onshape UI by
hand.

The assembly exposes 5 Configuration parameters: `emag_height`,
`inner_coil_id`, `inner_coil_od`, `outer_coil_id`, `outer_coil_od` --
confirmed via a live `GET .../configuration` call, not assumed. **All 5
are sent to Onshape on every pull** (via `get_step_for_config`, keyed on
the full 5-param dict). An earlier version of this module sent only
`emag_height` and left the 4 coil ID/OD params at their Onshape defaults,
on the assumption they were purely local resizes not worth an API call --
that assumption was wrong: Onshape models a real coil solid per winding
(`inner_coil`/`outer_coil` in the STEP export), and the coil ID params
also size the soft-iron core each coil wraps (there's no separate
core-diameter parameter -- the core IS sized by the coil ID it sits
inside). See `cache.py`'s docstring for the multi-param encoding details.

## Setup

1. Generate an API key pair at the
   [Onshape Developer Portal](https://dev-portal.onshape.com) (the secret
   is shown once -- copy both values immediately).
2. Copy `.env.example` (repo root) to `.env` and fill in the real
   access/secret key. `.env` is gitignored -- never commit it.
3. `pip install onshape-rest-api-client python-dotenv httpx` (httpx comes
   in as a transitive dependency of the first, but is used directly here).

```python
from onshape.client import get_client

with get_client() as client:
    resp = client.get("/users/sessioninfo")
    resp.raise_for_status()
    print(resp.json()["name"])
```

## Why a hand-written `auth.py` instead of the installed client library

`onshape-rest-api-client` (PyPI, OpenAPI-generated straight from Onshape's
own spec) only ships a Bearer-token `AuthenticatedClient` -- meant for
OAuth2 access tokens, not Onshape's access-key/secret-key **HMAC** scheme.
There's no built-in support for signing requests with an API key pair.
`auth.py`'s `OnshapeApiKeyAuth` is a plain `httpx.Auth` implementing that
signing, ported from Onshape's own reference sample
(`onshape-public/apikey`, `python/apikey/onshape.py`) to modern Python 3 --
fetched and read directly before writing this, not reconstructed from
memory.

`client.py` builds a plain `httpx.Client` with that auth attached, rather
than using `onshape_rest_api_client.Client`/`.request()` calls. Two
reasons, both hit while validating this:

- **Base URL must include `/api`.** The generated client's per-endpoint
  code requests relative paths like `/users/sessioninfo` (no `/api`
  prefix), so `ONSHAPE_BASE_URL` has to be `https://cad.onshape.com/api`,
  not just `https://cad.onshape.com` -- easy to get wrong once and see a
  confusing "Expecting value" JSON-decode error (the server serves an HTML
  page for the un-prefixed path, not JSON).
- **The generated response models have real bugs.** `session_info`'s
  typed model (`BTUserOAuth2SummaryInfo`) crashes deserializing a `null`
  nested `company` field (`TypeError: 'NoneType' object is not
  iterable`) -- confirmed the raw HTTP response was valid, correctly
  authenticated JSON; the crash is in the generated model's `from_dict`,
  not in this project's code. Given that, it's not worth taking on the
  generated client's model layer for every endpoint. The installed
  package is still useful as a **reference for exact endpoint paths and
  request shapes** (browse `onshape_rest_api_client/api/<tag>/<operation>.py`
  -- generated straight from Onshape's live spec) -- just parse responses
  as plain dicts via `response.json()` instead of the typed `.sync()`
  helpers.

## Getting a STEP export: `cache.py`

```python
from onshape.cache import get_step_for_height

step_path = get_step_for_height(1.5)  # emag_height, inches
```

Caches by `emag_height` in a local SQLite index
(`onshape/cad_cache/index.db`, gitignored along with the `.step` files
themselves -- large, regenerable, machine-local). A cache entry also
records the document's `sourceMicroversion` at fetch time; if the Onshape
document has been edited again since, the cached entry is treated as
stale and re-pulled automatically, so a mid-sweep CAD fix can't silently
leave old geometry in the cache. A cache hit costs one cheap
`get_configuration` read (to check the microversion) and no translation
call; a miss costs a translation request + poll + download (~10s here).

### Two non-obvious things that cost real debugging time getting here

- **The `configuration` parameter has to be a JSON body field on the
  right endpoint, not a query string.** `POST .../assemblies/d/{did}/w/
  {wid}/e/{eid}/export/step` (the dedicated "export as STEP" endpoint,
  and what the installed client's `create_assembly_export_step` wraps)
  silently **ignores** a `configuration` query parameter -- it exports at
  whatever the document's currently-saved configuration state is, no
  error, no warning. The one that actually honors it is the generic
  `POST .../assemblies/d/{did}/w/{wid}/e/{eid}/translations` endpoint
  (`translate_format` in the installed client), with `"configuration":
  "emag_height=1.5in"` (the `encodedId` form from
  `encodeConfigurationMap`, not the URL-encoded `queryParam` form) as a
  plain body field alongside `"formatName": "STEP"`. Confirmed against
  Onshape's own docs
  ([api-adv/configs](https://onshape-public.github.io/docs/api-adv/configs/))
  after silently getting the *default* geometry back twice.
- **`grouping: true` must be passed explicitly.** Without it (even though
  the installed client's stub documents `True` as the default), the
  export comes back as a **zip of flat, per-part STEP files at each
  part's own local origin** -- the same un-assembled, duplicate-position
  shape the old MATLAB pipeline had to work around with
  `apply_transforms.m` (see `../matlab/README.md`). With `grouping: true`,
  the result is one real STEP file with correct assembled placement
  (verified: the 4 `outer_emag_core` copies land at 4 distinct spread-out
  positions, not stacked on top of each other).

## Next

`get_step_for_config` and `gmsh_getdp/assembly/build_assembly.py` (whole-
assembly STEP import, real coil solids, incremental iron fuse) are wired
together and validated end-to-end with a real solve. Nonlinear iron is done,
and `gmsh_getdp/assembly/emag_sweep.py` is the first sweep driver built on
this cache (emag_height x inner_coil_id). At `emag_height=1.125in` the
spacer has zero height and Onshape omits the `chamber_spacer` body from the
export rather than erroring.
