"""
optimise.py

Optimisation of the field-shaping features picked on 2026-10-20 (no cap lip):
  - front-plate taper   fp_t (bore-rim thickness), fp_w (taper radial width)
  - pole-cap taper      cap_r (radius), cap_t (rim thickness), cap_w (taper width)
  - inner screen        is_t (wall), is_gap (tip gap to the cap back face plane),
                        is_flare (tip widened outward over the last 3 mm)
  - outer screen        os_gap, os_t, os_r

Every run is scored at the spec outer/inner ratio (q = 0.635, EMAG_SPEC.md)
with the channel-averaged peak B_r scaled to 300 G, and records the channel
shape, the coil drive, and the field in every iron part (mean / p95 / max |B|,
saturation ratio = p95 / KNEE_T, and the flux-mean field in the inner core and
in the real outer posts, undoing the 2-D ring smear).

    python optimise.py baseline|sweeps|search|tolerance|figures|all

Solves run in parallel and are cached in runs/opt/<hash>.
"""
import csv
import hashlib
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from explore import front_plate, pole_cap  # noqa: E402
from geometry import CHANNEL, DEFAULTS, R_IRON, R_RING, R_SCREEN, Y_COIL_TOP, Y_EXIT, Y_FP_IN, ring_geometry  # noqa: E402
from metrics import combine, shape_metrics  # noqa: E402
from shaping import coil_power, map_fig, profile_fig, _mpl, PALETTE  # noqa: E402
from solve import RUNS, solve_case  # noqa: E402

FIG = HERE / "figures" / "optimise"
RES = HERE / "results" / "optimise"
Q_SPEC = 0.635
N_IN, N_OUT = 230, 184           # spec turns (EMAG_SPEC.md)
KNEE_T = 1.5                     # A36 proxy B-H knee, same as ../assembly/emag_sweep.py
SAT_FLAG = 0.8                   # saturation ratio above which a run is flagged / rejected
WORKERS = 12

R_KEEP_IN = 20.2                 # inner screen / cap outer limit (channel keep-out 20.62)
R_KEEP_OS = 49.0                 # outer screen outer limit (outer winding sheet at 49.3)
R_CAP_MAX = 20.62                # pole cap radius limit
OS_R_MIN = CHANNEL["r_out"] + 3  # outer screen inner radius limit
from geometry import CONFIG  # noqa: E402
if CONFIG == "bn":
    from geometry import BN  # noqa: E402
    R_KEEP_IN = BN["bore"] / 2 - 0.1                    # inner screen slips inside the BN bore
    R_CAP_MAX = BN["bore"] / 2 - 0.1                    # pole cap inside the BN bore
    OS_R_MIN = BN["od"] / 2 + 0.0                       # outer screen slips over the BN body
    R_KEEP_OS = 2.5 * 25.4 - DEFAULTS["ocore_r"] - DEFAULTS["ocoil_build"] - 0.5   # outer coils at 50.0

BASE = dict(fp_t=1.5, fp_w=15.0, cap_r=20.0, cap_t=1.5, cap_w=6.0,
            is_t=3.0, is_gap=6.0, is_flare=0.0, os_gap=6.0, os_t=3.0, os_r=40.0)

SWEEPS = {
    "fp_t": [1.5, 2.5, 3.5, 5.0, 7.0, 9.525],
    "fp_w": [5.0, 10.0, 15.0, 20.0, 25.0],
    "cap_r": [17.5, 18.5, 19.5, 20.0, 20.6],
    "cap_t": [1.5, 2.5, 3.5, 5.0, 9.525],
    "cap_w": [2.0, 4.0, 6.0, 8.0, 10.0],
    "is_t": [1.5, 2.0, 3.0, 4.0, 5.0],
    "is_gap": [1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 15.0],
    "is_flare": [0.0, 1.0, 2.0, 2.5],
    "os_gap": [1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 15.0],
    "os_t": [1.5, 3.0, 5.0],
    "os_r": [35.0, 40.0, 45.0],
}
LABEL = {"fp_t": "front-plate rim thickness (mm)", "fp_w": "front-plate taper width (mm)",
         "cap_r": "pole-cap radius (mm)", "cap_t": "pole-cap rim thickness (mm)",
         "cap_w": "pole-cap taper width (mm)", "is_t": "inner screen wall (mm)",
         "is_gap": "inner screen tip gap (mm)", "is_flare": "inner screen tip flare (mm)",
         "os_gap": "outer screen tip gap (mm)", "os_t": "outer screen wall (mm)",
         "os_r": "outer screen inner radius (mm)"}
IS_R0 = DEFAULTS["icoil_ro"] + 0.4


# --------------------------------------------------------------------------
# geometry
def feasible(x):
    """Return None if x is buildable, else the reason."""
    if min(x["fp_t"], x["cap_t"], x["is_t"], x["os_t"]) < 1.5 - 1e-9:
        return "wall < 1.5 mm"
    if x["fp_t"] > 9.525 + 1e-9 or x["cap_t"] > 9.525 + 1e-9:
        return "rim thicker than plate"
    if x["cap_r"] > R_CAP_MAX + 1e-9:
        return "cap past keep-out"
    if x["cap_w"] > x["cap_r"] - x.get("core_r", DEFAULTS["core_r"]) + 1e-9:
        return "cap taper into the core"
    if x.get("is_r0", IS_R0) + x["is_t"] + x["is_flare"] > R_KEEP_IN + 1e-9:
        return "inner screen past keep-out"
    if x["os_r"] < OS_R_MIN - 1e-9 or x["os_r"] + x["os_t"] > R_KEEP_OS + 1e-9:
        return "outer screen clashes"
    if min(x["is_gap"], x["os_gap"]) < 0.5 or x["is_flare"] < 0:
        return "gap/flare out of range"
    return None


def inner_screen(t, gap, flare, r0=None):
    y_tip = Y_FP_IN + gap
    r0 = IS_R0 if r0 is None else r0
    r1 = r0 + t
    if flare <= 0:
        return [(r0, Y_COIL_TOP), (r1, Y_COIL_TOP), (r1, y_tip), (r0, y_tip)]
    # outward flare: full width r1 + flare over the last 3 mm, 45-deg-ish ramp back to r1
    return [(r0, Y_COIL_TOP), (r1, Y_COIL_TOP), (r1, y_tip + 3 + flare), (r1 + flare, y_tip + 3),
            (r1 + flare, y_tip), (r0, y_tip)]


def outer_screen(r0, t, gap):
    y_tip = Y_FP_IN + gap
    return [(r0, Y_COIL_TOP), (r0 + t, Y_COIL_TOP), (r0 + t, y_tip), (r0, y_tip)]


WIRE = 0.89                     # AWG 20 heavy polyimide, mm over insulation (spec windings: 46 x 5, 46 x 4)
N_TURNS_IN, N_TURNS_OUT = 230, 184


def layered(n_per_layer, turns):
    """(length, build) of a square-wound layered coil: n turns per layer, enough layers for `turns`."""
    layers = -(-turns // n_per_layer)
    return n_per_layer * WIRE, layers * WIRE


def cap_poly(R, T, t, w):
    """Pole cap of radius R and thickness T, front face on the exit plane; back-side chamfer from the
    edge face (length t) to radius R - w on the back face."""
    if t >= T:
        return [(0, Y_EXIT), (R, Y_EXIT), (R, Y_EXIT + T), (0, Y_EXIT + T)]
    return [(0, Y_EXIT), (R, Y_EXIT), (R, Y_EXIT + t), (R - w, Y_EXIT + T), (0, Y_EXIT + T)]


def to_params(x):
    fp_t = None if x["fp_t"] >= 9.525 else x["fp_t"]
    cap_t = None if x["cap_t"] >= 9.525 else x["cap_t"]
    p = dict(front_plate_poly=front_plate(t_tip=fp_t, w=x["fp_w"], bore=x.get("fp_bore")),
             pole_cap_poly=pole_cap(r=x["cap_r"], t_tip=cap_t, w=x["cap_w"]),
             extra_iron=[])
    if "cap_T" in x:
        # pole cap thickness other than the 9.525 plate: the inner coil starts right behind it
        p["pole_cap_poly"] = cap_poly(x["cap_r"], x["cap_T"], x["cap_t"], x["cap_w"])
        p["pole_cap_t"] = x["cap_T"]
    # optional keys (absent from older runs, so their cache keys are unchanged)
    if not x.get("no_is"):
        p["extra_iron"].append(inner_screen(x["is_t"], x["is_gap"], x["is_flare"], x.get("is_r0")))
    if not x.get("no_os"):
        p["extra_iron"].append(outer_screen(x["os_r"], x["os_t"], x["os_gap"]))
    if x.get("spacer"):
        # back plate moved L upstream. spacer_mode 1 "taller": cores, posts and coils all longer (same turns);
        # 2 "iron spacer": cores/posts longer, coils unchanged; 3 "inert spacer": cores, posts, coils
        # unchanged, non-magnetic gap of L between their rear ends and the back plate
        L, mode = x["spacer"], int(x.get("spacer_mode", 1))
        p["spacer"] = L
        # 4 "inert spacer, inner core only": as 2, but the inner core keeps its length with a
        # non-magnetic gap of L behind it (the outer posts are extended in iron)
        if mode in (2, 3, 4):
            p["icoil_y"] = (Y_FP_IN, Y_COIL_TOP)
            p["ocoil_y"] = (Y_FP_IN, Y_COIL_TOP)
        if mode == 3:
            p["core_gap"] = p["ring_gap"] = L
        if mode == 4:
            p["core_gap"] = L
    if x.get("core_gap"):
        # inner core shortened at its rear end by this much, non-magnetic insert in its place;
        # back plate and posts unchanged. coil_short: the inner coil ends at the same plane (the insert
        # is a disc as wide as the pole cap, under both core and coil), same ampere-turns
        p["core_gap"] = x["core_gap"]
        if x.get("coil_short"):
            p["icoil_y"] = (Y_FP_IN, Y_COIL_TOP - x["core_gap"])
    if x.get("f_in", 1.0) != 1.0 or x.get("f_out", 1.0) != 1.0:
        # coils wound over only the front fraction f of the core / posts (same turns, build grows as 1/f)
        from shaping import coils
        p.update(coils(f_in=x.get("f_in", 1.0), f_out=x.get("f_out", 1.0)))
    if "n_in" in x:
        # layered windings (fixed turns, real wire pitch), seated against the pole cap / front plate
        y0 = Y_EXIT + x.get("cap_T", 9.525) + x.get("in_off", 0.0)      # in_off: coil not fully seated
        li, bi = layered(int(x["n_in"]), N_TURNS_IN)
        p["icoil_y"] = (y0, y0 + li)
        p["icoil_ro"] = x.get("core_r", DEFAULTS["core_r"]) + bi
    if "n_out" in x:
        lo, bo = layered(int(x["n_out"]), N_TURNS_OUT)
        p["ocoil_y"] = (Y_FP_IN + x.get("out_off", 0.0), Y_FP_IN + x.get("out_off", 0.0) + lo)
        p["ocoil_build"] = bo
    if "core_r" in x:
        p["core_r"] = x["core_r"]
    if "core_r" in x and "n_in" not in x:
        # inner core radius changed: the inner winding keeps the same turns of the same wire (same window
        # area), so its build stays B_IN0 / f_in on top of the new core
        from shaping import B_IN0
        p["icoil_ro"] = x["core_r"] + B_IN0 / x.get("f_in", 1.0)
    if x.get("split_g"):
        # non-magnetic insert of thickness split_g cutting the inner core (and the pole cap if it reaches
        # it), centred split_z mm upstream of the pole cap back face; radius split_r (default core radius)
        yc = Y_FP_IN + x["split_z"]
        rr = x.get("split_r", DEFAULTS["core_r"])
        g = x["split_g"]
        p["air_cut"] = [[(0, yc - g / 2), (rr, yc - g / 2), (rr, yc + g / 2), (0, yc + g / 2)]]
    if x.get("ring_gap"):
        # outer posts shortened at their rear end by this much (non-magnetic in its place)
        p["ring_gap"] = x["ring_gap"]
    if x.get("no_bp"):
        p["back_plate"] = False
    if "mu" in x:
        p["mu_r"] = x["mu"]
    if "mu_screen" in x:
        p["mu_screen"] = x["mu_screen"]
    return p


def key(x):
    s = json.dumps({k: round(float(x[k]), 3) for k in sorted(x)})
    return "opt/" + hashlib.sha1(s.encode()).hexdigest()[:12]


def parts(x):
    """Labelled iron parts, in the geometry's priority order (first wins)."""
    p = to_params(x)
    r0, r1, _ = ring_geometry(DEFAULTS)
    yb = Y_COIL_TOP + p.get("spacer", 0.0)                 # back plate inner face
    ye = yb - p.get("core_gap", 0.0)                        # rear end of the inner core
    yr = yb - p.get("ring_gap", 0.0)                        # rear end of the outer posts
    rc = p.get("core_r", DEFAULTS["core_r"])
    return [("inner core", [(0, Y_EXIT), (rc, Y_EXIT), (rc, ye), (0, ye)]),
            ("pole cap", p["pole_cap_poly"]),
            ("front plate", p["front_plate_poly"]),
            ("back plate", [(0, yb), (200, yb), (200, yb + 9.525), (0, yb + 9.525)]),
            *([("inner screen", inner_screen(x["is_t"], x["is_gap"], x["is_flare"], x.get("is_r0")))]
              if not x.get("no_is") else []),
            *([("outer screen", outer_screen(x["os_r"], x["os_t"], x["os_gap"]))] if not x.get("no_os") else []),
            ("outer posts", [(r0, Y_FP_IN), (r1, Y_FP_IN), (r1, yr), (r0, yr)])]


# --------------------------------------------------------------------------
# solving (parallel, cached) and per-element weights
def _elem_weights(run_dir):
    """Element centroid (mm) and volume weight (2 pi r dA, mm^3) of every iron
    and ring element, from the mesh; saved next to the basis fields."""
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.open(str(run_dir / "axi.msh"))
    tags, coords, _ = gmsh.model.mesh.getNodes()
    xy = np.zeros((int(tags.max()) + 1, 2))
    xy[tags] = coords.reshape(-1, 3)[:, :2] * 1e3
    cents, vols = [], []
    groups = {t for d, t in gmsh.model.getPhysicalGroups(2)}
    for reg in (R_IRON, R_RING, R_SCREEN):
        if reg not in groups:          # R_SCREEN exists only when the screens get their own mu
            continue
        for ent in gmsh.model.getEntitiesForPhysicalGroup(2, reg):
            types, _, nodes = gmsh.model.mesh.getElements(2, ent)
            for t, nd in zip(types, nodes):
                nn = gmsh.model.mesh.getElementProperties(t)[3]
                tri = xy[np.asarray(nd, dtype=np.int64)].reshape(-1, nn, 2)[:, :3]
                c = tri.mean(axis=1)
                a = 0.5 * np.abs((tri[:, 1, 0] - tri[:, 0, 0]) * (tri[:, 2, 1] - tri[:, 0, 1])
                                 - (tri[:, 2, 0] - tri[:, 0, 0]) * (tri[:, 1, 1] - tri[:, 0, 1]))
                cents.append(c)
                vols.append(2 * np.pi * c[:, 0] * a)
    gmsh.finalize()
    np.savez_compressed(run_dir / "elems.npz", c=np.vstack(cents), v=np.concatenate(vols))


def _solve_one(args):
    name, x, map_ = args
    t0 = time.time()
    res = solve_case(name, to_params(x), map_=map_)
    d = res["dir"]
    if not (d / "elems.npz").exists():
        _elem_weights(d)
    meta = d / "meta.json"
    if not meta.exists():
        meta.write_text(json.dumps(dict(x=x, solve_s=time.time() - t0,
                                        n_elem=res.get("derived", {}).get("n_elem"))))
    for f in ("axisym.res", "axisym.pre"):     # GetDP's own outputs, ~10 MB each
        (d / f).unlink(missing_ok=True)
    return name


def solve_all(xs, map_=False):
    todo = [(key(x), x, map_) for x in xs if feasible(x) is None]
    new = [a for a in todo if not (RUNS / a[0] / "elems.npz").exists()
           or (map_ and not (RUNS / a[0] / "map.npz").exists())]
    if new:
        print(f"solving {len(new)} new geometries on {WORKERS} workers", flush=True)
        with ProcessPoolExecutor(WORKERS) as ex:
            list(ex.map(_solve_one, new))


# --------------------------------------------------------------------------
# metrics
def iron_metrics(x, res, ni_in, ni_out):
    d = res["dir"]
    E = np.load(d / "elems.npz")
    bi, bo = res["iron"]["inner"], res["iron"]["outer"]
    pts = bi[:, :2] * 1e3
    from scipy.spatial import cKDTree
    dist, idx = cKDTree(E["c"]).query(pts)
    assert np.max(dist) < 1e-3, f"element match failed ({np.max(dist):.2e} mm) in {d}"
    v = E["v"][idx]
    B = ni_in * bi[:, 2:4] + ni_out * bo[:, 2:4]                  # T, (Br, By)
    Bm = np.linalg.norm(B, axis=1)
    from matplotlib.path import Path as MPath
    label = np.full(len(pts), -1)
    for i, (_, poly) in enumerate(parts(x)):
        inside = MPath(poly).contains_points(pts, radius=1e-6) & (label < 0)
        label[inside] = i
    out = {}
    worst = 0.0
    for i, (nm, _) in enumerate(parts(x)):
        m = label == i
        if not m.any():
            continue
        s = nm.replace(" ", "_")
        b, w = Bm[m], v[m]
        o = np.argsort(b)
        cw = np.cumsum(w[o]) / w.sum()
        p95 = float(b[o][np.searchsorted(cw, 0.95)])
        if nm == "outer posts":
            # the ring smears 4 posts; real post field = ring flux / real post area.
            # B_y flux through each axial station of the ring, 1 mm bins.
            yb = np.round(pts[m, 1])
            A_post = 4 * np.pi * DEFAULTS["ocore_r"] ** 2
            phi = np.array([np.sum(B[m][yb == yy, 1] * w[yb == yy]) for yy in np.unique(yb)])
            dy = np.array([1.0 for _ in np.unique(yb)])
            b_post = np.abs(phi / dy) / A_post                       # mean B_y over the real posts
            out[f"{s}_Bflux_max_T"] = float(b_post.max())
            out[f"{s}_sat"] = out[f"{s}_Bflux_max_T"] / KNEE_T
        else:
            out[f"{s}_Bmean_T"] = float(np.sum(b * w) / w.sum())
            out[f"{s}_Bp95_T"] = p95
            out[f"{s}_Bmax_T"] = float(b.max())
            out[f"{s}_sat"] = p95 / KNEE_T
        if nm == "inner core":
            # flux-mean B_y along the shaft under the coil (exact in 2-D axisym)
            sel = m & (pts[:, 1] > Y_FP_IN) & (pts[:, 1] < Y_COIL_TOP + x.get("spacer", 0.0))
            yb = np.round(pts[sel, 1])
            A = np.pi * x.get("core_r", DEFAULTS["core_r"]) ** 2
            phi = np.array([np.sum(B[sel][yb == yy, 1] * v[sel][yb == yy]) for yy in np.unique(yb)])
            out["inner_core_Bflux_max_T"] = float(np.max(np.abs(phi)) / A)
            out[f"{s}_sat"] = max(p95, out["inner_core_Bflux_max_T"]) / KNEE_T
        worst = max(worst, out[f"{s}_sat"])
    out["sat_worst"] = worst
    out["sat_worst_part"] = max((k for k in out if k.endswith("_sat") and k != "sat_worst"), key=out.get)[:-4]
    return out


def channel_extras(res, ni_in, ni_out, scale_pk):
    """scale_pk: +-1, the sign that makes the channel B_r positive."""
    """Wall profiles, plume leakage, field-line angle at the walls near the anode."""
    ch = (ni_in * res["channel"]["inner"][:, 2:4] + ni_out * res["channel"]["outer"][:, 2:4]) * 1e4 * scale_pk
    r = res["channel"]["inner"][:, 0] * 1e3
    y = res["channel"]["inner"][:, 1] * 1e3
    rin, rout = r.min(), r.max()
    ya, ye = CHANNEL["y_anode"], CHANNEL["y_exit"]

    def wall(rw):
        m = np.abs(r - rw) < 1e-6
        o = np.argsort(y[m])
        return y[m][o], ch[m][o]
    out = {}
    yi, bi = wall(rin)
    yo, bo = wall(rout)
    out["exit_Br_inner_wall_G"] = float(np.interp(ye, yi, bi[:, 0]))
    out["exit_Br_outer_wall_G"] = float(np.interp(ye, yo, bo[:, 0]))
    out["exit_Br_in_over_out"] = out["exit_Br_inner_wall_G"] / out["exit_Br_outer_wall_G"]
    # angle of B from the wall (axial) direction, |B|-weighted over the last 5 mm before the anode;
    # 0 deg = field lines graze the wall (shielding-like), 90 deg = perpendicular
    for nm, (yy, bb) in (("inner", (yi, bi)), ("outer", (yo, bo))):
        m = (yy >= ya - 5) & (yy <= ya)
        ang = np.degrees(np.arctan2(np.abs(bb[m, 0]), np.abs(bb[m, 1])))
        w = np.hypot(bb[m, 0], bb[m, 1])
        out[f"anode_wall_angle_{nm}_deg"] = float(np.sum(ang * w) / w.sum())
    return out


def evaluate_x(x):
    name = key(x)
    res = solve_case(name, to_params(x))
    y, Br, By = combine(res, 1.0, Q_SPEC)
    m, s = shape_metrics(y, Br, By)
    ni_in, ni_out = abs(s), abs(s) * Q_SPEC
    params = to_params(x)
    m.update(case=name.split("/")[1], NI_inner=ni_in, NI_outer=ni_out,
             I_inner_A=ni_in / N_IN, I_outer_A=ni_out / N_OUT,
             P_W=coil_power(params, ni_in, ni_out),
             G_per_kAt_inner=300.0 / ni_in * 1e3)
    yy = y
    Brs, Bys = Br * s, By * s
    Bmag = np.hypot(Brs, Bys)
    for dmm in (10, 20):
        m[f"plume_B_{dmm}mm_pct"] = float(np.interp(Y_EXIT - dmm, yy, Bmag) / 300 * 100)
    m.update(channel_extras(res, ni_in, ni_out, np.sign(s)))
    m.update(iron_metrics(x, res, ni_in, ni_out))
    meta = json.loads((res["dir"] / "meta.json").read_text()) if (res["dir"] / "meta.json").exists() else {}
    m["solve_s"] = meta.get("solve_s")
    m["n_elem"] = meta.get("n_elem")
    m.update(x)
    m["objective"] = objective(m)
    return m, (yy, Brs, Bys), res


def objective(m, P0=None):
    """Lower is better. Shape terms in % of peak (the exit term uses B_r AT the
    exit, not the peak position: the peak is flat-topped and its argmax hops
    between grid points by 0.5-1 mm on sub-mm geometry changes), a mild power term, and a
    hard penalty past the saturation flag. Power is signed: 5 points per 100 %
    of the baseline's power, saved or spent. sigma is NOT scored (user, 2026-10-20)."""
    P0 = P0 or OBJ_P0
    j = (100 - m["exit_Br_pct"]) + abs(m["anode_B_pct"]) + m["gauss_rms_pct"] + 5 * (m["P_W"] / P0 - 1)
    j += 100 * max(0.0, m["sat_worst"] - SAT_FLAG)
    return j


OBJ_P0 = 2.5        # W, the baseline's power is about this; reset in baseline()

COLS = (["study", "case", "objective", "peak_offset_mm", "exit_Br_pct", "anode_Br_pct", "anode_B_pct",
         "half_depth_mm", "sigma_mm", "gauss_rms_pct", "plume_B_10mm_pct", "plume_B_20mm_pct",
         "exit_Br_inner_wall_G", "exit_Br_outer_wall_G", "exit_Br_in_over_out",
         "anode_wall_angle_inner_deg", "anode_wall_angle_outer_deg",
         "q", "NI_inner", "NI_outer", "I_inner_A", "I_outer_A", "P_W", "G_per_kAt_inner",
         "sat_worst", "sat_worst_part", "inner_core_Bflux_max_T", "outer_posts_Bflux_max_T"]
        + [f"{p}_{k}" for p in ("inner_core", "pole_cap", "front_plate", "back_plate", "inner_screen",
                                "outer_screen") for k in ("Bmean_T", "Bp95_T", "Bmax_T", "sat")]
        + ["outer_posts_sat"] + list(BASE) + ["n_elem", "solve_s"])


def write(study, rows, mode="w"):
    RES.mkdir(parents=True, exist_ok=True)
    f = RES / f"{study}.csv"
    with open(f, mode, newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, extrasaction="ignore")
        if mode == "w":
            w.writeheader()
        for r in rows:
            w.writerow({k: (round(float(v), 4) if isinstance(v, (float, np.floating)) else v)
                        for k, v in dict(r, study=study, q=Q_SPEC).items()})


def show(rows, title, var=None):
    print(f"\n{title}")
    print(f"{'var':>8s} {'obj':>6s} {'pk':>5s} {'anBr%':>6s} {'an|B|%':>6s} {'sig':>5s} {'gRMS':>5s} "
          f"{'P_W':>5s} {'sat':>5s} {'worst':>13s} {'coreB':>5s} {'postB':>5s} {'angI':>5s} {'angO':>5s}")
    for r in rows:
        v = f"{r[var]:8.2f}" if var else f"{r['case'][:8]:>8s}"
        print(f"{v} {r['objective']:6.2f} {r['peak_offset_mm']:5.1f} {r['anode_Br_pct']:6.1f} "
              f"{r['anode_B_pct']:6.1f} {r['sigma_mm']:5.1f} {r['gauss_rms_pct']:5.1f} {r['P_W']:5.2f} "
              f"{r['sat_worst']:5.2f} {r['sat_worst_part']:>13s} {r['inner_core_Bflux_max_T']:5.2f} "
              f"{r['outer_posts_Bflux_max_T']:5.2f} {r['anode_wall_angle_inner_deg']:5.0f} "
              f"{r['anode_wall_angle_outer_deg']:5.0f}")


# --------------------------------------------------------------------------
# studies
def baseline():
    solve_all([BASE])
    m, _, _ = evaluate_x(BASE)
    global OBJ_P0
    OBJ_P0 = m["P_W"]
    m["objective"] = objective(m)
    write("baseline", [m])
    show([m], "baseline: 10-15 tapers + inner/outer screens (gap 6), no lip")
    return m


def sweeps():
    baseline()
    xs = {v: [dict(BASE, **{v: val}) for val in vals] for v, vals in SWEEPS.items()}
    solve_all([x for l in xs.values() for x in l])
    allrows = []
    for v, l in xs.items():
        rows = []
        for x in l:
            if (why := feasible(x)) is not None:
                print(f"  skip {v}={x[v]}: {why}")
                continue
            m, _, _ = evaluate_x(x)
            m["swept"] = v
            rows.append(m)
        show(rows, f"sweep {v}", v)
        allrows += [dict(r, case=f"{v}={r[v]:g}") for r in rows]
    write("sweeps", allrows)
    sweep_fig(allrows)
    return allrows


def pattern_search(vars_, x0=None, steps=(2.0, 1.0, 0.5), max_iter=60):
    """Compass search on the cached grid: try +-step on each variable, move to
    the best improving point, halve the step when none improves."""
    baseline()
    x = dict(x0 or BASE)
    best, _, _ = evaluate_x(x)
    trace = [dict(best, case=f"it0")]
    it = 0
    for st in steps:
        while it < max_iter:
            it += 1
            cands = []
            for v in vars_:
                for sgn in (-1, 1):
                    c = dict(x, **{v: round(x[v] + sgn * st, 3)})
                    if feasible(c) is None:
                        cands.append(c)
            solve_all(cands)
            ms = [evaluate_x(c)[0] for c in cands]
            mb = min(ms, key=lambda m: m["objective"])
            if mb["objective"] < best["objective"] - 1e-3:
                best = mb
                x = {k: mb[k] for k in BASE}
                trace.append(dict(best, case=f"it{it} step{st:g}"))
                print(f"it {it} step {st}: obj {best['objective']:.3f}  " +
                      " ".join(f"{v}={x[v]:g}" for v in vars_), flush=True)
            else:
                break
    write("search_trace", trace)
    (RES / "optimum.json").write_text(json.dumps({k: best[k] for k in BASE}, indent=1))
    show(trace, "search trace")
    return x, best


def tolerance(x=None, d=0.5):
    x = x or json.loads((RES / "optimum.json").read_text())
    baseline()
    xs = [x] + [dict(x, **{v: x[v] + s * d}) for v in BASE for s in (-1, 1) if v not in ("os_r",)]
    xs = [c for c in xs if feasible(c) is None]
    solve_all(xs)
    rows = []
    for c in xs:
        m, _, _ = evaluate_x(c)
        dv = [v for v in BASE if abs(c[v] - x[v]) > 1e-9]
        m["case"] = "optimum" if not dv else f"{dv[0]} {c[dv[0]] - x[dv[0]]:+g}"
        rows.append(m)
    write("tolerance", rows)
    show(rows, "tolerance +-0.5 mm")
    return rows


# --------------------------------------------------------------------------
# figures
def sweep_fig(rows):
    plt = _mpl()
    FIG.mkdir(parents=True, exist_ok=True)
    vs = list(SWEEPS)
    keys = [("peak_offset_mm", "peak pos (mm)"), ("anode_B_pct", "|B| anode (%)"),
            ("sigma_mm", "sigma (mm)"), ("gauss_rms_pct", "Gauss misfit (%)"),
            ("P_W", "power (W)"), ("sat_worst", "worst sat. ratio")]
    fig, axs = plt.subplots(len(keys), len(vs), figsize=(2.1 * len(vs), 1.7 * len(keys)),
                            constrained_layout=True, sharey="row")
    for j, v in enumerate(vs):
        sub = sorted([r for r in rows if r["swept"] == v], key=lambda r: r[v])
        xv = [r[v] for r in sub]
        for i, (k, lab) in enumerate(keys):
            ax = axs[i, j]
            ax.plot(xv, [r[k] for r in sub], "o-", color=PALETTE[0], ms=3, lw=1.4)
            ax.axvline(BASE[v], color="0.6", lw=0.8, ls=":")
            if k == "sat_worst":
                ax.axhline(SAT_FLAG, color=PALETTE[7], lw=0.8, ls="--")
            if k == "peak_offset_mm":
                ax.axhline(0, color="0.3", lw=0.6)
            if j == 0:
                ax.set_ylabel(lab)
            if i == len(keys) - 1:
                ax.set_xlabel(LABEL[v], fontsize=7)
    fig.suptitle("One-at-a-time sweeps around the baseline (dotted), spec ratio, 300 G peak", x=0.01, ha="left")
    fig.savefig(FIG / "sweeps.png", dpi=120)
    plt.close(fig)


def sat_fig(items):
    """items: [(label, metrics)]. Saturation ratio per part."""
    plt = _mpl()
    names = ["inner_core", "pole_cap", "front_plate", "back_plate", "inner_screen", "outer_screen", "outer_posts"]
    fig, ax = plt.subplots(figsize=(8, 3.8), constrained_layout=True)
    w = 0.8 / len(items)
    for i, (lab, m) in enumerate(items):
        ax.bar(np.arange(len(names)) + i * w, [m.get(f"{n}_sat", np.nan) for n in names], w,
               color=PALETTE[i], label=lab)
    ax.axhline(SAT_FLAG, color=PALETTE[7], lw=1, ls="--", label=f"flag ({SAT_FLAG:g})")
    ax.set_xticks(np.arange(len(names)) + w * (len(items) - 1) / 2, [n.replace("_", " ") for n in names])
    ax.set_ylabel(f"saturation ratio (p95 |B| / {KNEE_T:g} T;\nposts: flux-mean B / {KNEE_T:g} T)")
    ax.set_title("Iron loading per part at 300 G peak", loc="left")
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(FIG / "saturation.png", dpi=130)
    plt.close(fig)


def figures():
    import shaping
    shaping.FIG = FIG
    xo = json.loads((RES / "optimum.json").read_text())
    picks = [("baseline", BASE), ("optimum", xo)]
    solve_all([x for _, x in picks], map_=True)
    out = []
    for lab, x in picks:
        m, prof, _ = evaluate_x(x)
        res = solve_case(key(x), to_params(x), map_=True)
        out.append((lab, m, prof, res, x))
    profile_fig([(lab, prof, m) for lab, m, prof, _, _ in out], "before_after", "Baseline vs optimum")
    map_fig([(lab, res, m, to_params(x)) for lab, m, _, res, x in out], "before_after", ncol=2)
    sat_fig([(lab, m) for lab, m, _, _, _ in out])


if __name__ == "__main__":
    what = sys.argv[1:] or ["all"]
    for w in what:
        if w in ("baseline", "all"):
            baseline()
        if w in ("sweeps", "all"):
            sweeps()
        if w == "search":
            pattern_search(sys.argv[sys.argv.index("search") + 1].split(","))
            break
        if w == "tolerance":
            tolerance()
        if w == "figures":
            figures()
