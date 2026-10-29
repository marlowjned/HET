"""Figures and a summary for a het2d.py run.

    python analyse.py runs/example [--avg-last 2e-6]

Reads:
- meta.json and monitor.csv: currents, particle counts, probe potentials;
- profiles.npz: particle moments binned in z, from the late part of the run;
- diags/fields: WarpX openPMD output (phi, E, rho per species; cell-centred, axes z, r).
Writes figures and summary.json into the run folder.
"""
import argparse
import json
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from geometry import design_b as g  # noqa: E402

Q_E, M_E = 1.602176634e-19, 9.1093837139e-31
M_AR = 39.948 * 1.66053906660e-27


def read_fields(run, t_from):
    """Average the openPMD field outputs with t >= t_from. Returns (z, r, {name: (nz, nr) array}, n, times)."""
    import openpmd_api as io
    s = io.Series(str(run / "diags/fields/openpmd_%06T.h5"), io.Access.read_only)
    acc, n, times, z, r = {}, 0, [], None, None
    for i in s.iterations:
        it = s.iterations[i]
        t = it.time * it.time_unit_SI
        if t < t_from:
            continue
        data = {}
        for name in ("phi", "rho_electrons", "rho_ions"):
            rc = it.meshes[name][io.Mesh_Record_Component.SCALAR]
            data[name] = rc.load_chunk()
        transverse = "x" if "x" in list(it.meshes["E"]) else "r"
        for comp, key in ((transverse, "Er"), ("z", "Ez")):
            data[key] = it.meshes["E"][comp].load_chunk()
        s.flush()
        m = it.meshes["phi"]
        lab = list(m.axis_labels)
        if z is None:
            assert lab[-2:] in (["z", "r"], ["z", "x"]), lab
            dz, dr = m.grid_spacing
            z0, r0 = m.grid_global_offset
            shape = data["phi"].shape[-2:]
            z = z0 + (np.arange(shape[0]) + 0.5) * dz
            r = r0 + (np.arange(shape[1]) + 0.5) * dr
        for k, v in data.items():
            acc[k] = acc.get(k, 0) + (v[0] if v.ndim == 3 else v)
        n += 1
        times.append(t)
    s.close()
    return z, r, {k: v / max(n, 1) for k, v in acc.items()}, n, times


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run", type=Path)
    ap.add_argument("--avg-last", type=float, default=None, help="average fields over the last [s] (default: last 25 %%)")
    a = ap.parse_args()
    run = a.run.resolve()
    meta = json.loads((run / "meta.json").read_text())
    mon = pd.read_csv(run / "monitor.csv")
    t_end = mon["t_us"].iloc[-1] * 1e-6
    avg = a.avg_last or 0.25 * t_end
    t_from = t_end - avg
    late = mon[mon["t_us"] * 1e-6 >= t_from]
    v_d = meta["args"]["v_anode"]
    mdot = meta["args"]["mdot"]
    i_mass = meta["mass_equivalent_current_A"]

    # ---------------------------------------------------------------- currents and counts
    fig, axs = plt.subplots(2, 2, figsize=(13, 8))
    t = mon["t_us"]
    roll = lambda s: s.rolling(8, min_periods=1).mean()                 # noqa: E731
    ax = axs[0, 0]
    ax.plot(t, roll(mon["I_discharge"]), label="discharge current (electrons to the anode)")
    ax.plot(t, roll(mon["I_i_z_hi"]), label="ion current through the cathode plane (beam)")
    ax.plot(t, roll(mon["I_i_x_lo"] + mon["I_i_x_hi"]), label="ion current to the radial walls")
    ax.plot(t, roll(-(mon["I_e_x_lo"] + mon["I_e_x_hi"])), label="electron current to the radial walls")
    ax.axhline(i_mass, color="k", ls=":", lw=0.8, label=f"mass-flow equivalent {i_mass:.2f} A")
    ax.set_ylabel("A")
    ax.set_xlabel("t [us]")
    ax.set_yscale("symlog", linthresh=0.1)
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)
    ax.set_title("currents (rolling mean of 8 records)", fontsize=10)
    ax = axs[0, 1]
    ax.plot(t, mon["n_e"], label="electron macroparticles")
    ax.plot(t, mon["n_i"], label="ion macroparticles")
    ax.set_xlabel("t [us]")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax.set_title("macroparticle counts", fontsize=10)
    ax = axs[1, 0]
    for c, lab in (("phi_mid_anode5", "5 mm from the anode face"), ("phi_mid_exit", "exit plane"),
                   ("phi_mid_plume15", "15 mm downstream")):
        ax.plot(t, roll(mon[c]), label=lab)
    ax.set_xlabel("t [us]")
    ax.set_ylabel("V")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax.set_title("potential at the channel mid-radius", fontsize=10)

    # ---------------------------------------------------------------- fields
    z, r, f, nf, times = read_fields(run, t_from)
    ne = -f["rho_electrons"] / Q_E
    ni = f["rho_ions"] / Q_E
    w = (r / r.sum()) if meta.get("geometry", "rz") == "rz" else np.full(r.size, 1 / r.size)  # radial mean
    rm = lambda a2: (a2 * w[None, :]).sum(1)                               # noqa: E731
    zmm = z * 1e3
    ax = axs[1, 1]
    ax.plot(zmm, rm(f["phi"]), label="potential [V]")
    ax2 = ax.twinx()
    ax2.plot(zmm, rm(f["Ez"]) / 1e3, color="C3", label="E_z [kV/m]")
    ax.set_xlabel("z [mm] (anode face at -41.43, exit plane 0)")
    ax.set_ylabel("V")
    ax2.set_ylabel("kV/m", color="C3")
    ax.axvline(0, color="k", lw=0.6, ls="--")
    ax.grid(alpha=0.3)
    ax.set_title(f"radial mean over the last {avg * 1e6:.2f} us ({nf} outputs)", fontsize=10)
    fig.tight_layout()
    fig.savefig(run / "currents.png", dpi=140)
    plt.close(fig)

    # maps
    fig, axs = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
    ext = [zmm[0] - 0.5 * (zmm[1] - zmm[0]), zmm[-1] + 0.5 * (zmm[1] - zmm[0]), r[0] * 1e3, r[-1] * 1e3]
    for ax, (arr, title, cmap) in zip(axs, ((f["phi"], "potential [V]", "viridis"),
                                            (ne, "electron density [m^-3]", "magma"),
                                            (ni, "ion density [m^-3]", "magma"))):
        im = ax.imshow(arr.T, origin="lower", extent=ext, aspect="auto", cmap=cmap)
        plt.colorbar(im, ax=ax, label=title)
        ax.axvline(0, color="w", lw=0.6, ls="--")
        ax.set_ylabel("r [mm]")
        ax.set_title(title, fontsize=10)
    axs[-1].set_xlabel("z [mm]")
    fig.tight_layout()
    fig.savefig(run / "maps.png", dpi=140)
    plt.close(fig)

    # ---------------------------------------------------------------- particle moments
    summary = {"t_end_us": t_end * 1e6, "avg_window_us": avg * 1e6, "field_outputs_averaged": nf,
               "stop_reason": meta.get("stop_reason"), "ms_per_step": meta.get("ms_per_step")}
    pf = run / "profiles.npz"
    fig, axs = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
    axs[0, 0].plot(zmm, rm(ne), label="electrons")
    axs[0, 0].plot(zmm, rm(ni), label="ions")
    axs[0, 0].set_ylabel("m^-3")
    axs[0, 0].set_title("density (radial mean)", fontsize=10)
    axs[0, 0].legend(fontsize=8)
    axs[0, 1].plot(zmm, rm(f["phi"]))
    axs[0, 1].set_ylabel("V")
    axs[0, 1].set_title("potential (radial mean)", fontsize=10)
    if pf.exists():
        p = np.load(pf)
        zc = 0.5 * (p["z_edges"][1:] + p["z_edges"][:-1]) * 1e3
        with np.errstate(invalid="ignore", divide="ignore"):
            vz, vr, vt = (p[k] / p["we"] for k in ("e_vz", "e_vr", "e_vt"))
            v2 = p["e_v2"] / p["we"]
            te = M_E * (v2 - vz ** 2 - vr ** 2 - vt ** 2) / (3 * Q_E)      # thermal part, drifts removed
            ivz = p["i_vz"] / p["wi"]
        axs[1, 0].plot(zc, te, label="electron temperature [eV] (drifts removed)")
        axs[1, 0].plot(zc, 0.5 * M_E * vt ** 2 / Q_E, label="azimuthal (Hall) drift energy [eV]")
        axs[1, 0].set_ylabel("eV")
        axs[1, 0].legend(fontsize=8)
        axs[1, 0].set_title(f"electrons ({int(p['samples'])} samples)", fontsize=10)
        axs[1, 1].plot(zc, ivz / 1e3)
        axs[1, 1].set_ylabel("km/s")
        axs[1, 1].set_title("mean ion axial velocity", fontsize=10)
        i_hi = np.isfinite(ivz)
        v_exit = float(np.nanmean(ivz[-5:])) if i_hi.any() else float("nan")
        summary.update(Te_peak_eV=float(np.nanmax(te)), Te_peak_z_mm=float(zc[np.nanargmax(te)]),
                       ion_vz_at_cathode_plane_km_s=v_exit / 1e3)
    for ax in axs.ravel():
        ax.axvline(0, color="k", lw=0.6, ls="--")
        ax.grid(alpha=0.3)
    for ax in axs[1]:
        ax.set_xlabel("z [mm]")
    fig.tight_layout()
    fig.savefig(run / "profiles.png", dpi=140)
    plt.close(fig)

    i_d = float(late["I_discharge"].mean())
    i_b = float(late["I_i_z_hi"].mean())
    summary.update(
        discharge_current_A=i_d, beam_current_A=i_b, discharge_power_W=v_d * i_d,
        current_utilisation=i_b / i_d if i_d else None, mass_utilisation=i_b / i_mass,
        ion_wall_current_A=float((late["I_i_x_lo"] + late["I_i_x_hi"]).mean()),
        peak_ne_m3=float(np.max(rm(ne))), peak_ne_z_mm=float(zmm[np.argmax(rm(ne))]),
        phi_exit_V=float(np.interp(0.0, z, rm(f["phi"]))))
    if "ion_vz_at_cathode_plane_km_s" in summary and np.isfinite(summary["ion_vz_at_cathode_plane_km_s"]):
        thrust = i_b / Q_E * M_AR * summary["ion_vz_at_cathode_plane_km_s"] * 1e3
        summary.update(thrust_mN_estimate=thrust * 1e3, isp_s_estimate=thrust / (mdot * 9.80665),
                       anode_efficiency_estimate=thrust ** 2 / (2 * mdot * v_d * i_d) if i_d else None)
    (run / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
