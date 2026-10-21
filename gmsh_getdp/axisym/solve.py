"""
solve.py

Build one axisymmetric geometry and solve its two basis fields (1 A-turn on
the inner coil; 1 A-turn per outer coil). The iron is linear, so any
excitation is NI_in * b_inner + NI_out * b_outer, and both are kept.

    from solve import solve_case
    res = solve_case("baseline", {}, map_=True)
    res["channel"]["inner"]   # (N, 4): r, y, B_r, B_y on the channel grid (m, T)
"""
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from geometry import build  # noqa: E402

GETDP = HERE.parents[1] / "tools" / "getdp-3.5.0-Windows64" / "getdp.exe"
RUNS = HERE / "runs"


def _table(path):
    """GetDP SimpleTable (x y z v... per line) -> (N, 4) array of x, y, B_x,
    B_y for a vector, (N, 3) x, y, value for a scalar."""
    a = np.loadtxt(path, ndmin=2)
    if a.shape[1] >= 6:
        return a[:, [0, 1, 3, 4]]
    return a[:, [0, 1, 3]]


def _getdp(run_dir, ni_in, ni_out, map_, ni_trim=0):
    cmd = [str(GETDP), "axisym.pro", "-msh", "axi.msh", "-solve", "Res", "-pos", "Post",
           "-setnumber", "NI_INNER", str(ni_in), "-setnumber", "NI_OUTER", str(ni_out),
           "-setnumber", "NI_TRIM", str(ni_trim),
           "-setnumber", "Map", "1" if map_ else "0", "-v", "2"]
    r = subprocess.run(cmd, cwd=run_dir, capture_output=True, text=True)
    if r.returncode != 0 or not (run_dir / "b_channel.txt").exists():
        raise RuntimeError(f"getdp failed in {run_dir}:\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}")


def solve_case(name, params, map_=False, force=False):
    run_dir = RUNS / name
    cache = run_dir / "basis.npz"
    if cache.exists() and not force and (not map_ or (run_dir / "map.npz").exists()):
        d = dict(np.load(cache, allow_pickle=True))
        names = [k[3:] for k in d if k.startswith("ch_")]
        res = {"channel": {k: d[f"ch_{k}"] for k in names},
               "iron": {k: d[f"ir_{k}"] for k in names}, "dir": run_dir}
        if map_:
            m = dict(np.load(run_dir / "map.npz"))
            res["map"] = {k: (m[f"b_{k}"], m[f"raz_{k}"]) for k in names}
        return res

    run_dir.mkdir(parents=True, exist_ok=True)
    for f in run_dir.iterdir():
        f.unlink()
    t0 = time.time()
    derived = build(params, run_dir)
    shutil.copy(HERE / "axisym.pro", run_dir / "axisym.pro")
    t_mesh = time.time() - t0

    res = {"channel": {}, "iron": {}, "map": {}, "dir": run_dir, "derived": derived}
    bases = [("inner", (1, 0, 0)), ("outer", (0, 1, 0))]
    if params.get("trim"):
        bases.append(("trim", (0, 0, 1)))
    for basis, (ni_in, ni_out, ni_tr) in bases:
        _getdp(run_dir, ni_in, ni_out, map_, ni_tr)
        res["channel"][basis] = _table(run_dir / "b_channel.txt")
        res["iron"][basis] = _table(run_dir / "b_iron.txt")
        if map_:
            res["map"][basis] = (_table(run_dir / "b_map.txt"), _table(run_dir / "raz_map.txt"))
        for f in ("b_channel.txt", "b_iron.txt", "b_map.txt", "raz_map.txt"):
            (run_dir / f).unlink(missing_ok=True)
    for f in ("axi.pre", "axi.res"):
        (run_dir / f).unlink(missing_ok=True)

    np.savez_compressed(cache, **{f"ch_{k}": v for k, v in res["channel"].items()},
                        **{f"ir_{k}": v for k, v in res["iron"].items()})
    if map_:
        np.savez_compressed(run_dir / "map.npz", **{f"b_{k}": v[0] for k, v in res["map"].items()},
                            **{f"raz_{k}": v[1] for k, v in res["map"].items()})
    print(f"[{name}] {derived['n_elem']} elements, mesh {t_mesh:.1f}s, total {time.time() - t0:.1f}s",
          flush=True)
    return res


if __name__ == "__main__":
    r = solve_case("smoke", {}, force=True)
    for k, v in r["channel"].items():
        print(k, v.shape, v[:3])
    print(r["iron"]["inner"].shape, r["iron"]["inner"][:2])
