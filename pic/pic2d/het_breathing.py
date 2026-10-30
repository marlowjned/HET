"""2-D PIC-MCC/DSMC model of Design B for breathing-mode studies (WarpX; no bplsim).

What it adds to het2d.py (the first example):
- **Kinetic neutrals.** Injected at the anode face with the propellant mass flow. Excitation and ionization are
  particle-particle collisions (DSMC), so the propellant depletes and refills, which is what breathing is.
  Elastic electron and ion collisions use a fixed neutral background (MCC) of the undepleted profile. Particle
  elastic collisions would split a new, roughly 1000x lighter neutral off every event for a negligible recoil:
  about 8000 extra neutral macroparticles per 100 steps in the argon test.
- **A quasi-neutral cathode.** Electrons are added in a thin slab upstream of the cathode plane wherever ions
  outnumber electrons there, so no virtual cathode forms.
- **A propellant table.** Krypton is the set point (xs/Kr, fitted to the rate tables of the 1-D reference by
  xs/make_kr.py; no ion-neutral collisions, as in the 1-D reference); argon (warpx-data) is there for tests.
- **Anomalous transport as in the 1-D reference** (HallThruster.jl TwoZoneBohm): nu = c(z) omega_ce, c going
  linearly from --alpha-in to --alpha-out over --transition (0.1 L_ch) centred on the exit plane. The defaults are
  the 1-D reference case (c1 = 1/800, c2 = 1/16; ../hallthruster).

- **Sheath walls** (--walls sheath, the default). The 1-D reference loses electron energy to the channel walls
  (HallThruster.jl WallSheath, boron nitride) but keeps the electrons. Here an electron that reaches a radial wall is
  reflected if its normal energy is below a cut-off, otherwise absorbed and re-emitted at the same place as a cold
  secondary (--see-te). The cut-off is set so that Maxwellian electrons at --wall-te lose the power of the 1-D
  model (sheath_cut). Ions and neutrals reflect, as in the 1-D reference (no ion wall losses).

- **Start from the 1-D reference** (--init-1d, default the HallThruster.jl run `ref`). From a cold start the 1-D
  model needs about 400 us (an ignition spike, then relaxation swings) before the breathing settles; the PIC starts
  from its time-averaged state instead: neutral density, plasma density, Te and ion velocity against z. The same
  neutral profile is the fixed background of the elastic electron collisions (MCC), so the classical cross-field
  mobility sees the depleted neutrals (the undepleted profile overstates it about 250 times at the exit plane).
  Its breathing modulation is still missing; ionization and excitation do see the particle neutrals.

Geometry (set-point channel), magnet map (165 G) and permittivity scaling are as in het2d.py.

    source ../activate.sh
    python het_breathing.py --propellant Ar --t-end 0.2e-6 --tag kn_test      # machinery test on argon data
    python het_breathing.py --t-end 300e-6 --tag breathing                     # krypton (needs xs/Kr; GPU)
"""
import argparse
import csv
import json
import os
import sys
import time
from math import pi, sqrt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
from geometry import design_b as g  # noqa: E402
from het2d import BFILES, BSUM, KB, M_E, Q_E, subst, write_anomalous_xs  # noqa: E402

AMU = 1.66053906660e-27
WD = Path.home() / "tools/warpx-data/MCC_cross_sections"
# electron processes: name -> (file, energy loss [eV] or None); ion processes: name -> (file, angle model)
PROPELLANTS = {
    "Ar": {"mass_u": 39.948, "dir": WD / "Ar",
           "electron": {"elastic": ("electron_scattering.dat", None), "excitation1": ("excitation_1.dat", 11.5),
                        "ionization": ("ionization.dat", 15.7596112)},
           "ion": {"elastic": ("ion_scattering.dat", "isotropic"), "elastic_back": ("ion_back_scatter.dat", "backward")}},
    # krypton: files and energies come from xs/Kr/processes.json, written by the LXCat converter
    "Kr": {"mass_u": 83.798, "dir": HERE / "xs" / "Kr"},
}


def anomalous_expr(args):
    """nu_anom(z) = c(z) e |B|(z) / m_e, |B| from the fit to the magnet map at the channel mid-radius; c goes
    linearly from alpha_in to alpha_out over the transition length centred on the exit plane (z = 0), as
    HallThruster.jl's TwoZoneBohm."""
    fit = json.loads(BSUM.read_text())["bmag_fit_mid_radius"]
    bexpr = ("if(z < zpk, bfu + (bpk - bfu)*exp(-0.5*((z - zpk)/sgu)^2), "
             "bfd + (bpk - bfd)*exp(-0.5*((z - zpk)/sgd)^2))")
    alpha = "(a_in + (a_out - a_in)*min(max(z/ltr + 0.5, 0), 1))"
    consts = dict(bpk=fit["B_peak_T"], zpk=fit["z_peak_m"], sgu=fit["sigma_up_m"], sgd=fit["sigma_down_m"],
                  bfu=fit["B_floor_up_T"], bfd=fit["B_floor_down_T"], a_in=args.alpha_in, a_out=args.alpha_out,
                  ltr=args.transition, qm=Q_E / M_E)
    nu_max = max(args.alpha_in, args.alpha_out) * Q_E / M_E * fit["B_peak_T"] * 1.05
    return f"{alpha}*qm*{bexpr}", consts, nu_max


def profiles_1d(tag):
    """Time-averaged axial profiles of a HallThruster.jl run (../hallthruster/runs/<tag>), z in PIC coordinates
    (exit plane at 0): neutral density, plasma density, Te [eV], singly charged ion velocity."""
    o = json.loads((HERE.parent / "hallthruster" / "runs" / tag / "output.json").read_text())
    a = o["output"]["average"]
    gas = next(iter(a["neutrals"]))
    z = np.array(a["z"]) - o["input"]["config"]["thruster"]["geometry"]["channel_length"]
    return {"z": z, "nn": np.array(a["neutrals"][gas]["n"]), "ne": np.array(a["ne"]), "te": np.array(a["Tev"]),
            "ui": np.array(a["ions"][gas][0]["u"]), "tag": tag}


def pwl(z, v, knots):
    """Parser expression of v(z), piecewise linear through the knots, as a balanced tree of if() (depth log2 of
    the number of segments, so it stays cheap when evaluated per particle)."""
    vk = np.interp(knots, z, v).tolist()                     # plain floats: the parser cannot read np.float64(...)
    knots = [float(k) for k in knots]

    def seg(i):
        s = (vk[i + 1] - vk[i]) / (knots[i + 1] - knots[i])
        return f"({vk[i]!r} + {s!r}*(z - {knots[i]!r}))"

    def tree(lo, hi):
        if hi - lo == 1:
            return seg(lo)
        mid = (lo + hi) // 2
        return f"if(z < {knots[mid]!r}, {tree(lo, mid)}, {tree(mid, hi)})"
    return tree(0, len(knots) - 1)


def sheath_cut(te, M, sigma0=0.24, eps_star=45.0, e_see=4.0):
    """Electron reflection cut-off [V] at the radial walls that loses the electron power of HallThruster.jl's
    WallSheath model (boron nitride: sigma0 0.24, eps* 45 eV; loss scale 1) at electron temperature te [eV].

    1-D: P = nu_ew (2 Te + (1 - g) phi_s), nu_ew = h u_B / (dr (1 - g)), h = 0.86/sqrt(3), u_B = sqrt(e Te / M),
         g = min(g_max, sigma0 + 1.5 Te / eps* (1 - sigma0)), phi_s = Te ln((1 - g) sqrt(M / (2 pi m_e))).
    PIC: Maxwellian electrons reach the two walls at (2 / dr) (vbar / 4) exp(-x) per electron, each lost one carries
         (x + 2) Te and is replaced by a secondary of mean energy e_see; x = cut / Te. dr cancels."""
    g_max = 1.0 - 8.3 * sqrt(M_E / M)
    gam = min(g_max, sigma0 + 1.5 * te / eps_star * (1.0 - sigma0))
    phi_s = te * np.log((1.0 - gam) * sqrt(M / (2 * pi * M_E)))
    u_b, vbar = sqrt(Q_E * te / M), sqrt(8 * Q_E * te / (pi * M_E))
    target = 0.86 / sqrt(3) * u_b / (1.0 - gam) * (2 * te + (1.0 - gam) * phi_s)
    f = lambda x: 0.5 * vbar * np.exp(-x) * ((x + 2) * te - e_see) - target  # noqa: E731
    lo, hi = 0.5, 40.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
    return 0.5 * (lo + hi) * te


def propellant(name):
    p = dict(PROPELLANTS[name])
    pj = p["dir"] / "processes.json"
    if pj.exists():
        p.update(json.loads(pj.read_text()))
    if "electron" not in p:
        sys.exit(f"no cross sections for {name}: expected {pj} (see REPORT.md, krypton cross sections)")
    for proc in list(p["electron"].values()) + list(p["ion"].values()):
        if not (p["dir"] / proc[0]).exists():
            sys.exit(f"missing cross-section file {p['dir'] / proc[0]}")
    p["mass"] = p["mass_u"] * AMU
    return p


def parse():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    a = ap.add_argument
    a("--tag", default="run")
    a("--propellant", choices=sorted(PROPELLANTS), default="Kr")
    a("--mdot", type=float, default=3.0e-6, help="propellant mass flow [kg/s] (set point 3.0 mg/s)")
    a("--v-anode", type=float, default=300.0)
    a("--t-gas", type=float, default=600.0, help="neutral temperature at injection [K]")
    a("--steps", type=int, default=None)
    a("--t-end", type=float, default=1e-6)
    a("--dt", type=float, default=5e-11)
    a("--cell", type=float, default=0.3e-3)
    a("--gamma2", type=float, default=400.0, help="vacuum permittivity scaling factor")
    a("--z-cathode", type=float, default=30e-3, help="cathode plane, downstream of the exit [m]")
    a("--init-1d", default="ref", help="HallThruster.jl run (../hallthruster/runs/<tag>) to start from, or none")
    a("--n0", type=float, default=5e16, help="initial plasma density [m^-3] without --init-1d; cathode weights")
    a("--te0", type=float, default=10.0)
    a("--ppc", type=float, default=20, help="initial plasma macroparticles per cell and species")
    a("--ppc-neutral", type=float, default=30, help="neutral macroparticles per cell at the channel density")
    a("--resample-ppc", type=int, default=80)
    a("--cathode-te", type=float, default=10.0, help="temperature of the cathode electrons [eV]")
    a("--inject-depth", type=float, default=5e-3, help="cathode slab centre, upstream of the cathode plane [m]")
    a("--slab", type=float, default=1e-3, help="cathode slab half-width [m]")
    a("--inject-every", type=int, default=20, help="steps between cathode corrections")
    a("--walls", choices=["sheath", "reflect", "absorb"], default="sheath",
      help="radial walls: sheath (energy loss, see the docstring), reflect (no loss), absorb (all lost)")
    a("--wall-te", type=float, default=20.0, help="Te [eV] at which the sheath walls match the 1-D wall power")
    a("--sheath-cut", type=float, default=None, help="electron reflection cut-off [V] (default: from --wall-te)")
    a("--see-te", type=float, default=2.0, help="temperature of the re-emitted wall electrons [eV]")
    a("--alpha-in", type=float, default=1 / 800, help="anomalous c in the channel (1-D reference c1)")
    a("--alpha-out", type=float, default=1 / 16, help="anomalous c in the plume (1-D reference c2)")
    a("--transition", type=float, default=0.1 * g.DEPTH * 1e-3, help="c1 -> c2 ramp length [m] (0.1 L_ch)")
    a("--dsmc-e-every", type=int, default=2)
    a("--dsmc-i-every", type=int, default=20)
    a("--threads", type=int, default=12)
    a("--diag-every", type=float, default=0.5e-6)
    a("--monitor-every", type=int, default=200, help="steps between records (10 ns at the default dt)")
    a("--max-particles", type=float, default=8e6)
    a("--max-wall", type=float, default=6600.0)
    a("--chunk", type=int, default=2000)
    a("--checkpoint-every", type=float, default=20e-6, help="simulated time between checkpoints [s] (0: none)")
    a("--restart", default=None, help="checkpoint directory to continue from (same arguments otherwise)")
    a("--seed", type=int, default=1)
    return ap.parse_args()


def main():
    args = parse()
    os.environ.setdefault("OMP_NUM_THREADS", str(args.threads))
    from pywarpx import callbacks, particle_containers, picmi
    prop = propellant(args.propellant)
    M = prop["mass"]
    out = HERE / "runs" / args.tag
    out.mkdir(parents=True, exist_ok=True)
    os.chdir(out)

    # ------------------------------------------------------------------ domain (planar x-z, x = r)
    r_lo, r_hi = g.R_IN * 1e-3, g.R_OUT * 1e-3
    z_lo, z_hi = -g.DEPTH * 1e-3, args.z_cathode
    nr = max(8, 8 * round((r_hi - r_lo) / args.cell / 8))
    nz = max(8, 8 * round((z_hi - z_lo) / args.cell / 8))
    dx, dz = (r_hi - r_lo) / nr, (z_hi - z_lo) / nz
    l_y = 2 * pi * 0.5 * (r_lo + r_hi)                         # per-metre-of-y quantities -> totals
    area = pi * (r_hi ** 2 - r_lo ** 2)                        # = l_y * (r_hi - r_lo)
    nsteps = args.steps or int(round(args.t_end / args.dt))
    pbc = "reflecting" if args.walls == "reflect" else "absorbing"
    grid = picmi.Cartesian2DGrid(
        number_of_cells=[nr, nz], lower_bound=[r_lo, z_lo], upper_bound=[r_hi, z_hi],
        lower_boundary_conditions=["neumann", "dirichlet"], upper_boundary_conditions=["neumann", "dirichlet"],
        lower_boundary_conditions_particles=[pbc, "absorbing"], upper_boundary_conditions_particles=[pbc, "absorbing"],
        warpx_potential_lo_z=args.v_anode, warpx_potential_hi_z=0.0, warpx_blocking_factor=8, warpx_max_grid_size=32)
    solver = picmi.ElectrostaticSolver(grid=grid, method="Multigrid", required_precision=1e-6,
                                       maximum_iterations=200, warpx_self_fields_verbosity=0)

    # ------------------------------------------------------------------ species
    vth_e = sqrt(args.te0 * Q_E / M_E)
    vth_n = sqrt(KB * args.t_gas / M)
    vbar = sqrt(8 * KB * args.t_gas / (pi * M))
    flux_n = args.mdot / M / area                              # neutral flux at the anode face [m^-2 s^-1]
    n_ch = 2 * flux_n / vbar                                   # free-molecular channel density (half-Maxwellian)
    w_slot = r_hi - r_lo
    prof = None if args.init_1d == "none" else profiles_1d(args.init_1d)
    knots = np.linspace(z_lo, z_hi, int(round((z_hi - z_lo) / 1e-3)) + 1)       # 1 mm
    if prof:
        nn_expr = f"exp({pwl(prof['z'], np.log(prof['nn']), knots)})"
        nn_max = 1.05 * float(np.interp(knots, prof["z"], prof["nn"]).max())
    else:
        nn_expr = subst("n_ch*if(z <= 0, 1, (2/pi)*atan(w_slot/(2*z)))", dict(n_ch=n_ch, w_slot=w_slot, pi=pi))
        nn_max = n_ch
    cell_v = dx * dz                                           # per metre of y
    w_n = n_ch * cell_v / args.ppc_neutral                     # neutral macroparticle weight (per metre of y)
    ppc_inj_n = flux_n * dx * args.dt / w_n                    # injected neutral macroparticles per cell per step
    resample = dict(warpx_do_resampling=True, warpx_resampling_algorithm="leveling_thinning",
                    warpx_resampling_trigger_max_avg_ppc=args.resample_ppc, warpx_resampling_min_ppc=10)
    save = {f"warpx_save_particles_at_{b}": True for b in ("zlo", "zhi")}
    if args.walls != "reflect":
        save.update({f"warpx_save_particles_at_{b}": True for b in ("xlo", "xhi")})
    always = {f"warpx_reflection_model_{b}": "1" for b in ("xlo", "xhi")}
    refl_n = always if args.walls != "reflect" else {}          # neutrals always bounce off the walls
    refl_i, refl_e, cut = {}, {}, None
    if args.walls == "sheath":
        cut = args.sheath_cut or sheath_cut(args.wall_te, M)
        v_cut = sqrt(2 * Q_E * cut / M_E)
        refl_i = always
        refl_e = {f"warpx_reflection_model_{b}": f"if(v < {v_cut!r}, 1, 0)" for b in ("xlo", "xhi")}
    if prof:
        ne_expr = f"exp({pwl(prof['z'], np.log(np.maximum(prof['ne'], 1e12)), knots)})"
        se_expr = f"sqrt({Q_E / M_E!r}*max({pwl(prof['z'], prof['te'], knots)}, 0.1))"
        init_e = picmi.AnalyticDistribution(density_expression=ne_expr, warpx_momentum_spread_expressions=[se_expr] * 3)
        init_i = picmi.AnalyticDistribution(density_expression=ne_expr,
                                            momentum_expressions=[None, None, pwl(prof["z"], prof["ui"], knots)],
                                            warpx_momentum_spread_expressions=[repr(vth_n)] * 3)
    else:
        init_e = picmi.UniformDistribution(density=args.n0, rms_velocity=[vth_e] * 3)
        init_i = picmi.UniformDistribution(density=args.n0, rms_velocity=[vth_n] * 3)
    electrons = picmi.Species(particle_type="electron", name="electrons", initial_distribution=init_e,
                              **save, **resample, **refl_e)
    ions = picmi.Species(particle_type=args.propellant, name="ions", charge="q_e", mass=M, initial_distribution=init_i,
                         **save, **resample, **refl_i)
    v_drift = vbar / 2                                         # mean axial speed of the injected half-Maxwellian
    init_n = picmi.AnalyticDistribution(density_expression=nn_expr, momentum_expressions=[None, None, repr(v_drift)],
                                        warpx_momentum_spread_expressions=[repr(vth_n)] * 3)
    inj_n = picmi.UniformFluxDistribution(flux=flux_n, flux_normal_axis="z", surface_flux_position=z_lo + 0.05 * dz,
                                          flux_direction=+1, rms_velocity=[vth_n] * 3,
                                          gaussian_flux_momentum_distribution=True)
    neutrals = picmi.Species(particle_type=args.propellant, name="neutrals", charge=0, mass=M,
                             initial_distribution=[init_n, inj_n], warpx_do_not_deposit=True,
                             warpx_do_not_gather=True, warpx_reflection_model_zlo="1",   # the anode face re-emits
                             warpx_save_particles_at_zhi=True, warpx_do_resampling=True,
                             warpx_resampling_algorithm="leveling_thinning",
                             warpx_resampling_trigger_max_avg_ppc=int(3 * args.ppc_neutral),
                             warpx_resampling_min_ppc=10, **refl_n)

    # ------------------------------------------------------------------ collisions
    xs = prop["dir"]
    e_dsmc, e_mcc = {}, {}
    for name, (f, energy) in prop["electron"].items():
        d = {"cross_section": str(xs / f)}
        if energy is not None:
            d["energy"] = energy
        if name.startswith("elastic"):
            e_mcc[name] = d                                  # fixed background: no neutral splitting
        else:
            e_dsmc[name] = d                                 # excitation, ionization: depleting neutrals
            if name == "ionization":
                d["target_species"] = neutrals
    coll_en = picmi.DSMCCollisions(name="coll_en", species=[electrons, neutrals], product_species=[electrons, ions],
                                   ndt_supercycle=args.dsmc_e_every, scattering_processes=e_dsmc)
    coll_e_el = picmi.MCCCollisions(name="coll_e_el", species=electrons, background_density=nn_expr,
                                    background_temperature=args.t_gas, background_mass=M, max_background_density=nn_max,
                                    ndt_supercycle=args.dsmc_e_every, scattering_processes=e_mcc)
    i_proc = {name: {"cross_section": str(xs / f), "scattering_angle_model": model}
              for name, (f, model) in prop["ion"].items()}
    coll_in = (picmi.MCCCollisions(name="coll_in", species=ions, background_density=nn_expr,
                                   background_temperature=args.t_gas, background_mass=M, max_background_density=nn_max,
                                   ndt_supercycle=args.dsmc_i_every, scattering_processes=i_proc) if i_proc else None)
    k_anom = 1e-13
    write_anomalous_xs(out / "anomalous_xs.dat", k_anom)
    nu_expr, nu_consts, nu_max = anomalous_expr(args)
    coll_anom = picmi.MCCCollisions(
        name="coll_anom", species=electrons, background_density=subst(f"({nu_expr})/k_anom", {**nu_consts, "k_anom": k_anom}),
        background_temperature=300.0, background_mass=1e-15, max_background_density=nu_max / k_anom,
        ndt_supercycle=args.dsmc_e_every, scattering_processes={"elastic": {"cross_section": str(out / "anomalous_xs.dat")}})

    # ------------------------------------------------------------------ simulation and diagnostics
    n_diag = max(1, int(round(args.diag_every / args.dt)))
    diag = picmi.FieldDiagnostic(name="fields", grid=grid, period=n_diag, write_dir="diags",
                                 data_list=["phi", "Ex", "Ez", "rho_electrons", "rho_ions"],
                                 warpx_format="openpmd", warpx_openpmd_backend="h5",
                                 warpx_particle_fields_to_plot=[picmi.ParticleFieldDiagnostic(name="w", func="1",
                                                                                              do_average=False)],
                                 warpx_particle_fields_species=["neutrals"])
    if args.restart:
        args.restart = str(Path(args.restart).resolve(strict=True))
    sim = picmi.Simulation(solver=solver, time_step_size=args.dt, max_steps=nsteps, verbose=1, particle_shape=1,
                           warpx_amr_restart=args.restart,
                           warpx_collisions=[c for c in (coll_en, coll_e_el, coll_in, coll_anom) if c],
                           warpx_use_filter=False,
                           warpx_random_seed=args.seed)
    sim.add_applied_field(picmi.LoadAppliedField(read_fields_from_path=str(BFILES["xz"]), load_E=False, load_B=True))
    sim.add_species(electrons, layout=picmi.PseudoRandomLayout(n_macroparticles_per_cell=args.ppc, grid=grid))
    sim.add_species(ions, layout=picmi.PseudoRandomLayout(n_macroparticles_per_cell=args.ppc, grid=grid))
    sim.add_species(neutrals, layout=[picmi.PseudoRandomLayout(n_macroparticles_per_cell=args.ppc_neutral, grid=grid),
                                      picmi.PseudoRandomLayout(n_macroparticles_per_cell=ppc_inj_n, grid=grid)])
    sim.add_diagnostic(diag)
    if args.checkpoint_every > 0:
        sim.add_diagnostic(picmi.Checkpoint(name="chk", period=max(1, int(round(args.checkpoint_every / args.dt))),
                                            write_dir="chk", warpx_file_prefix="chk"))

    i_mass = args.mdot / M * Q_E
    meta = {"args": vars(args), "propellant": {k: v for k, v in prop.items() if k != "dir"} | {"dir": str(prop["dir"])},
            "nr": nr, "nz": nz, "dx": dx, "dz": dz, "r_lo": r_lo, "r_hi": r_hi, "z_lo": z_lo, "z_hi": z_hi,
            "l_y": l_y, "area_m2": area, "steps": nsteps, "flux_neutral": flux_n, "n_neutral_channel": n_ch,
            "w_neutral": w_n, "ppc_inject_neutral": ppc_inj_n, "mass_equivalent_current_A": i_mass,
            "nu_anom_expr": nu_expr, "nu_anom_consts": nu_consts, "sheath_cut_V": cut,
            "init_1d": prof["tag"] if prof else None, "neutral_background_max": nn_max, "restart": args.restart}
    (out / "meta.json").write_text(json.dumps(meta, indent=2, default=str))
    print(f"{args.propellant}: {nr} x {nz} cells ({dx * 1e3:.3f} x {dz * 1e3:.3f} mm), dt {args.dt:g} s, {nsteps} steps; "
          f"neutrals {n_ch:.3g} m^-3, mass-equivalent current {i_mass:.3f} A; walls {args.walls}"
          + (f" (cut-off {cut:.1f} V)" if cut else ""), flush=True)
    if args.monitor_every % args.inject_every:
        sys.exit("--monitor-every must be a multiple of --inject-every (the boundary tallies run with the cathode)")

    # ------------------------------------------------------------------ callbacks
    import warnings
    warnings.filterwarnings("ignore")
    inv_g2 = 1.0 / args.gamma2

    def scale_rho():
        mf = sim.fields.get("rho_fp", level=0)
        ng = mf.n_grow_vect
        mf.mult(inv_g2, 0, mf.n_comp, max(ng[0], ng[1]))

    callbacks.installafterdeposition(scale_rho)
    pcs = {s: particle_containers.ParticleContainerWrapper(s) for s in ("electrons", "ions", "neutrals")}
    buf = particle_containers.ParticleBoundaryBufferWrapper()
    rng = np.random.default_rng(args.seed)
    z_c = z_hi - args.inject_depth
    w_inj = args.n0 * cell_v / args.ppc                        # cathode electron weight (per metre of y)
    vth_c = sqrt(args.cathode_te * Q_E / M_E)
    st = {"inj_w": 0.0, "t0": time.time(), "last": time.time()}

    def slab_weight(name):
        z = pcs[name].get_particle_z(copy_to_host=True)
        if not z:
            return 0.0
        w = pcs[name].get_particle_weight(copy_to_host=True)
        return float(sum(np.sum(wi[np.abs(zi - z_c) < args.slab]) for zi, wi in zip(z, w)))

    def cathode():
        """Quasi-neutral cathode: add electrons in the slab where ions outnumber them."""
        if sim.extension.warpx.getistep(0) % args.inject_every:
            return
        deficit = slab_weight("ions") - slab_weight("electrons")
        if deficit <= 0:
            return
        n = int(np.ceil(deficit / w_inj))
        x = rng.uniform(r_lo, r_hi, n)
        z = rng.uniform(z_c - args.slab, z_c + args.slab, n)
        u = rng.normal(0.0, vth_c, (3, n))
        pcs["electrons"].add_particles(x=x, y=np.zeros(n), z=z, ux=u[0], uy=u[1], uz=u[2], w=np.full(n, deficit / n))
        st["inj_w"] += deficit

    mon = out / "monitor.csv"
    cols = ["step", "t_us", "n_e", "n_i", "n_n", "neutral_mass_ug", "I_e_anode", "I_i_anode", "I_e_cath", "I_i_beam",
            "I_n_out_A_equiv", "I_cathode_inj", "I_discharge", "I_e_wall", "P_wall_W", "phi_mid_exit", "wall_s"]
    if not (args.restart and mon.exists()):
        with open(mon, "w", newline="") as fh:
            csv.writer(fh).writerow(cols)

    def host(a):
        return a.get() if hasattr(a, "get") else np.asarray(a)     # cupy on GPU builds

    def bufcat(name, b, comp):
        if not buf.get_particle_boundary_buffer_size(name, b):
            return np.zeros(0)
        return np.concatenate([host(a).ravel() for a in buf.get_particle_boundary_buffer(name, b, comp, 0)])

    acc = dict.fromkeys(("e_zlo", "i_zlo", "e_zhi", "i_zhi", "n_zhi", "e_wall", "p_wall"), 0.0)
    tally = {"e_zlo": ("electrons", "z_lo"), "i_zlo": ("ions", "z_lo"), "e_zhi": ("electrons", "z_hi"),
             "i_zhi": ("ions", "z_hi"), "n_zhi": ("neutrals", "z_hi")}
    vth_see = sqrt(args.see_te * Q_E / M_E)

    def boundaries():
        """Tally what left through the anode and cathode planes. With sheath walls, re-emit the electrons the walls
        absorbed as cold secondaries at the same place: their energy is lost, the electrons are kept."""
        if sim.extension.warpx.getistep(0) % args.inject_every:
            return
        for k, (name, b) in tally.items():
            acc[k] += float(bufcat(name, b, "w").sum())
        if args.walls == "sheath":
            for b, x0, sgn in (("x_lo", r_lo, 1.0), ("x_hi", r_hi, -1.0)):
                w = bufcat("electrons", b, "w")
                if not w.size:
                    continue
                z = bufcat("electrons", b, "z")
                u2 = sum(bufcat("electrons", b, c) ** 2 for c in ("ux", "uy", "uz"))
                acc["e_wall"] += float(w.sum())
                acc["p_wall"] += float(np.sum(w * 0.5 * M_E * u2))
                n = w.size
                un = sgn * vth_see * np.sqrt(-2.0 * np.log(1.0 - rng.random(n)))      # flux-weighted normal speed
                ut = rng.normal(0.0, vth_see, (2, n))
                pcs["electrons"].add_particles(x=np.full(n, x0 + sgn * 1e-3 * dx), y=np.zeros(n), z=z,
                                               ux=un, uy=ut[0], uz=ut[1], w=w)
        buf.clear_buffer()

    def monitor():
        step = sim.extension.warpx.getistep(0)
        if step % args.monitor_every:
            return
        dtw = args.monitor_every * args.dt
        c = l_y * Q_E / dtw
        ie_a, ii_a, ie_c, ii_c, in_out = (c * acc[k] for k in ("e_zlo", "i_zlo", "e_zhi", "i_zhi", "n_zhi"))
        i_wall, p_wall = c * acc["e_wall"], l_y * acc["p_wall"] / dtw
        for k in acc:
            acc[k] = 0.0
        i_inj = c * st["inj_w"]
        st["inj_w"] = 0.0
        wn = pcs["neutrals"].get_particle_weight(copy_to_host=True)
        mass_ug = (float(sum(np.sum(x) for x in wn)) if wn else 0.0) * l_y * M * 1e9
        from pywarpx import fields
        phi = host(fields.PhiFPWrapper(0)[...])
        zn = np.linspace(z_lo, z_hi, phi.shape[1])
        row = [step, step * args.dt * 1e6, pcs["electrons"].get_particle_count(), pcs["ions"].get_particle_count(),
               pcs["neutrals"].get_particle_count(), mass_ug, ie_a, ii_a, ie_c, ii_c, in_out, i_inj, ie_a - ii_a,
               i_wall, p_wall, float(np.interp(0.0, zn, phi[phi.shape[0] // 2])), time.time() - st["t0"]]
        with open(mon, "a", newline="") as fh:
            csv.writer(fh).writerow(row)
        if time.time() - st["last"] > 30 or step == nsteps:
            st["last"] = time.time()
            s0 = st.setdefault("step0", step)
            print(f"[{step}/{nsteps}] t {row[1]:.3f} us, {(time.time() - st['t0']) / max(step - s0, 1) * 1e3:.2f} ms/step, "
                  f"N_e {row[2]}, N_i {row[3]}, N_n {row[4]}, neutral mass {mass_ug:.3f} ug, I_d {row[12]:.3f} A, "
                  f"I_beam {ii_c:.3f} A, cathode {i_inj:.3f} A, wall {i_wall:.3f} A / {p_wall:.0f} W", flush=True)

    callbacks.installafterstep(cathode)
    callbacks.installafterstep(boundaries)
    callbacks.installafterstep(monitor)

    t0 = time.time()
    done, stop, first = 0, "", None
    while done < nsteps:
        n = min(args.chunk, nsteps - done)
        sim.step(n)
        done = sim.extension.warpx.getistep(0)
        first = done - n if first is None else first             # the step this run started from (restarts)
        tot = sum(pcs[s].get_particle_count() for s in pcs)
        if tot > args.max_particles:
            stop = f"stopped early: {tot} macroparticles"
        elif time.time() - t0 > args.max_wall:
            stop = f"stopped early: wall time > {args.max_wall:g} s"
        if stop:
            print(stop, flush=True)
            break
    wall = time.time() - t0
    meta.update(wall_s=wall, ms_per_step=wall / max(done - (first or 0), 1) * 1e3, steps_done=done,
                stop_reason=stop or "completed")
    (out / "meta.json").write_text(json.dumps(meta, indent=2, default=str))
    print(f"done: {done} steps ({done * args.dt * 1e6:.3f} us) in {wall:.0f} s; output in {out}")


if __name__ == "__main__":
    main()
