"""2-D axisymmetric (r-z) PIC-MCC model of the Design B channel and near plume: argon, 300 V (WarpX; no bplsim).

Geometry: planar x-z by default. It is the channel annulus unrolled: x = r, z axial, y = azimuth with
circumference 2 pi r_mid, which keeps the channel area exactly. WarpX 26.10 leaves the radial Poisson boundary
conditions of an RZ grid unset when r_lo > 0 (REPORT I5); --geometry rz is kept for when that is fixed.

Domain: the channel annulus, r from 22.12 to 31.45 mm. It runs from the anode face (z = -41.43 mm, 300 V) to a
cathode plane 30 mm downstream of the exit (0 V), where electrons are injected.
- **Radial walls:** zero-normal-field condition. They absorb ions, and reflect electrons below a wall-normal energy
  cutoff (absorbing those above it), as a placeholder for the sheath of a floating BN wall. WarpX's electrostatic
  solver has no dielectric.
- **Magnetic field:** the Design B map from ../magnetics (165 G design point).
- **Neutrals:** a fixed argon background sized to 2 mg/s, with argon cross sections (warpx-data):
  - electrons: elastic, excitation, ionization;
  - ions: isotropic and backward elastic.
- **Anomalous transport:** Bohm-type, as a fictitious collision process with frequency alpha * omega_ce.
- **Vacuum permittivity:** scaled by gamma2. The charge density is divided by gamma2 right after deposition, which
  makes Debye lengths and plasma periods resolvable. Densities stay physical, so the collision rates are real.

Every placeholder and assumption is listed in ../REPORT.md.

    source ../activate.sh
    python het2d.py --steps 2000 --tag smoke        # smoke test
    python het2d.py --t-end 10e-6 --tag example     # the example run (background, ~1-2 h on 16 threads)
    python analyse.py runs/example                   # figures and summary
"""
import argparse
import csv
import json
import os
import sys
import time
from math import atan, pi, sqrt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from geometry import design_b as g  # noqa: E402

XS = Path.home() / "tools/warpx-data/MCC_cross_sections/Ar"
BFILES = {"rz": HERE.parent / "magnetics/out/bfield_rz.h5", "xz": HERE.parent / "magnetics/out/bfield_xz.h5"}
BSUM = HERE.parent / "magnetics/out/summary.json"
Q_E, M_E, KB = 1.602176634e-19, 9.1093837139e-31, 1.380649e-23
M_AR = 39.948 * 1.66053906660e-27
E_IZ, E_EXC = 15.7596112, 11.5          # argon ionization threshold, excitation energy in excitation_1.dat [eV]


def parse():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    a = p.add_argument
    a("--tag", default="run", help="output folder runs/<tag>")
    a("--geometry", choices=["xz", "rz"], default="xz",
      help="xz: planar slab (x = r); rz: axisymmetric (wrong in WarpX 26.10 for r_lo > 0, see REPORT I5)")
    a("--steps", type=int, default=None, help="number of steps (overrides --t-end)")
    a("--t-end", type=float, default=10e-6, help="simulated time [s]")
    a("--dt", type=float, default=5e-11, help="time step [s]")
    a("--cell", type=float, default=0.3e-3, help="target cell size [m]")
    a("--gamma2", type=float, default=400.0, help="vacuum permittivity scaling factor")
    a("--v-anode", type=float, default=300.0, help="discharge voltage [V]")
    a("--z-cathode", type=float, default=30e-3, help="cathode plane, downstream of the exit [m]")
    a("--mdot", type=float, default=2e-6, help="argon mass flow [kg/s]")
    a("--t-gas", type=float, default=600.0, help="neutral temperature [K] (anode / gas distributor)")
    a("--n-facility", type=float, default=3e17, help="facility background neutral density [m^-3]")
    a("--n0", type=float, default=5e16, help="initial plasma density [m^-3] (uniform, as the LANDMARK benchmark; "
                                            "starting from a weak seed, cathode electrons piled up at the magnetic "
                                            "barrier and ignition would take most of 10 us)")
    a("--te0", type=float, default=10.0, help="initial electron temperature [eV]")
    a("--ppc", type=float, default=20, help="initial macroparticles per cell and species")
    a("--resample-ppc", type=int, default=80, help="thin a species (leveling thinning) when its average "
                                                    "macroparticles per cell exceed this")
    a("--i-inject", type=float, default=6.0, help="cathode electron one-way injection current [A]")
    a("--inject-depth", type=float, default=5e-3,
      help="injection plane this far upstream of the cathode plane [m]; on the boundary itself, ~97 %% of the "
           "electrons left again within their first gyro-orbit (r_L ~ 2.5 mm at 5 eV in the 20-40 G plume)")
    a("--te-inject", type=float, default=5.0, help="cathode electron temperature [eV]")
    a("--ppc-inject", type=float, default=0.15, help="injected macroparticles per cell per step (sets the "
                                                    "macroparticle weight of the cathode electrons and their "
                                                    "ionization products)")
    a("--walls", choices=["reflect", "cutoff", "absorb"], default="reflect",
      help="radial walls: reflect = both species reflected specularly (no wall losses; an axial model, like the "
           "LANDMARK benchmark's wall-free domain); cutoff = ions absorbed, electrons reflected below --wall-cut-ev "
           "of wall-normal energy; absorb = both absorbed")
    a("--wall-cut-ev", type=float, default=50.0, help="electron wall-normal energy cutoff for --walls cutoff [eV]")
    a("--alpha-in", type=float, default=1 / 160, help="anomalous collision frequency / omega_ce in the channel")
    a("--alpha-out", type=float, default=1 / 16, help="anomalous collision frequency / omega_ce in the plume")
    a("--threads", type=int, default=12, help="OpenMP threads (12 physical cores here; 24 was slower)")
    a("--diag-every", type=float, default=0.25e-6, help="field output interval [s]")
    a("--monitor-every", type=int, default=500, help="steps between current/count records")
    a("--profile-every", type=int, default=2000, help="steps between particle-profile samples")
    a("--profile-start", type=float, default=0.5, help="sample profiles after this fraction of the run")
    a("--mcc-e-every", type=int, default=2, help="electron collisions every N steps (MCC supercycling)")
    a("--mcc-i-every", type=int, default=20, help="ion collisions every N steps")
    a("--max-particles", type=float, default=6e6, help="stop early above this many macroparticles")
    a("--max-wall", type=float, default=6600.0, help="stop early after this many seconds of wall time")
    a("--chunk", type=int, default=2000, help="steps per sim.step() call (limits are checked between chunks)")
    a("--seed", type=int, default=1)
    return p.parse_args()


def neutral_density_expr(args, area):
    """Fixed neutral background: n0 in the channel (free-molecular half-Maxwellian carrying mdot), falling off
    downstream of the exit like the view factor of a slot of the channel width, plus the facility background."""
    vbar = sqrt(8 * KB * args.t_gas / (pi * M_AR))
    n_ch = 2 * (args.mdot / M_AR) / (area * vbar)              # half-Maxwellian: flux = n vbar / 2
    w = (g.R_OUT - g.R_IN) * 1e-3
    expr = f"n_ch*if(z <= 0, 1, (2/pi)*atan(w_slot/(2*z))) + n_fac"
    return expr, dict(n_ch=n_ch, w_slot=w, n_fac=args.n_facility, pi=pi), n_ch


def anomalous_expr(args):
    """nu_anom(z) = alpha(z) e |B|(z) / m_e, |B| from the fit to the magnet map at the channel mid-radius;
    alpha switches from alpha_in to alpha_out over ~1 mm at the exit plane."""
    fit = json.loads(BSUM.read_text())["bmag_fit_mid_radius"]
    bexpr = ("if(z < zpk, bfu + (bpk - bfu)*exp(-0.5*((z - zpk)/sgu)^2), "
             "bfd + (bpk - bfd)*exp(-0.5*((z - zpk)/sgd)^2))")
    alpha = "(a_in + (a_out - a_in)*0.5*(1 + tanh(z/0.001)))"
    consts = dict(bpk=fit["B_peak_T"], zpk=fit["z_peak_m"], sgu=fit["sigma_up_m"], sgd=fit["sigma_down_m"],
                  bfu=fit["B_floor_up_T"], bfd=fit["B_floor_down_T"], a_in=args.alpha_in, a_out=args.alpha_out,
                  qm=Q_E / M_E)
    nu_max = max(args.alpha_in, args.alpha_out) * Q_E / M_E * fit["B_peak_T"] * 1.05
    return f"{alpha}*qm*{bexpr}", consts, nu_max


def subst(expr, consts):
    """Write the constants into a parser expression as numbers (MCCCollisions takes no expression constants)."""
    import re
    for k in sorted(consts, key=len, reverse=True):
        expr = re.sub(rf"\b{k}\b", f"({consts[k]!r})", expr)
    return expr


def write_anomalous_xs(path, k_anom):
    """sigma(E) = k / v(E): the collision frequency n sigma v = n k is the same for every electron energy."""
    e = np.concatenate([np.geomspace(1e-3, 1.0, 400, endpoint=False), np.linspace(1.0, 5000.0, 5000)])
    v = np.sqrt(2 * e * Q_E / M_E)
    np.savetxt(path, np.column_stack([e, k_anom / v]), fmt="%.8e")


def main():
    args = parse()
    os.environ.setdefault("OMP_NUM_THREADS", str(args.threads))
    from pywarpx import callbacks, particle_containers, picmi    # after OMP_NUM_THREADS

    out = HERE / "runs" / args.tag
    out.mkdir(parents=True, exist_ok=True)
    os.chdir(out)

    # ------------------------------------------------------------------ domain and grid
    r_lo, r_hi = g.R_IN * 1e-3, g.R_OUT * 1e-3
    z_lo, z_hi = -g.DEPTH * 1e-3, args.z_cathode
    # cell counts in multiples of 8 so multigrid can coarsen (in RZ with r_lo > 0 the radial Poisson BCs are
    # unset in WarpX 26.10 and a coarsenable grid diverges; see ../REPORT.md, I5)
    nr = max(8, 8 * round((r_hi - r_lo) / args.cell / 8))
    nz = max(8, 8 * round((z_hi - z_lo) / args.cell / 8))
    planar = args.geometry == "xz"
    l_y = 2 * pi * 0.5 * (r_lo + r_hi) if planar else 1.0     # planar: per-metre-of-y quantities -> totals
    dr, dz = (r_hi - r_lo) / nr, (z_hi - z_lo) / nz
    area = pi * (r_hi ** 2 - r_lo ** 2)
    nsteps = args.steps or int(round(args.t_end / args.dt))

    Grid = picmi.Cartesian2DGrid if planar else picmi.CylindricalGrid
    grid = Grid(
        number_of_cells=[nr, nz], lower_bound=[r_lo, z_lo], upper_bound=[r_hi, z_hi],
        lower_boundary_conditions=["neumann", "dirichlet"], upper_boundary_conditions=["neumann", "dirichlet"],
        lower_boundary_conditions_particles=["reflecting" if args.walls == "reflect" else "absorbing", "absorbing"],
        upper_boundary_conditions_particles=["reflecting" if args.walls == "reflect" else "absorbing", "absorbing"],
        warpx_potential_lo_z=args.v_anode, warpx_potential_hi_z=0.0,
        warpx_blocking_factor=8, warpx_max_grid_size=32)
    solver = picmi.ElectrostaticSolver(grid=grid, method="Multigrid", required_precision=1e-6,
                                       maximum_iterations=200, warpx_self_fields_verbosity=0)

    # ------------------------------------------------------------------ species
    vth_e = sqrt(args.te0 * Q_E / M_E)
    vth_i = sqrt(KB * args.t_gas / M_AR)
    vth_inj = sqrt(args.te_inject * Q_E / M_E)
    flux_inj = args.i_inject / (Q_E * area)                    # one-way electron flux at the cathode plane
    init_e = picmi.UniformDistribution(density=args.n0, rms_velocity=[vth_e] * 3)
    init_i = picmi.UniformDistribution(density=args.n0, rms_velocity=[vth_i] * 3)
    inj_e = picmi.UniformFluxDistribution(
        flux=flux_inj, flux_normal_axis="z", surface_flux_position=z_hi - args.inject_depth, flux_direction=-1,
        rms_velocity=[vth_inj] * 3, gaussian_flux_momentum_distribution=True)
    save = {f"warpx_save_particles_at_{b}": True
            for b in (("zlo", "zhi") if args.walls == "reflect" else ("xlo", "xhi", "zlo", "zhi"))}
    # sheath placeholder on the radial walls: an electron is reflected when its wall-normal energy is below
    # wall_cut_ev and lost above it (WarpX evaluates the probability with v = wall-normal u = gamma v)
    vc2 = 2 * args.wall_cut_ev * Q_E / M_E
    wall = ({f"warpx_reflection_model_{b}": f"if(v*v < {vc2!r}, 1, 0)" for b in ("xlo", "xhi")}
            if args.walls == "cutoff" else {})
    resample = dict(warpx_do_resampling=True, warpx_resampling_algorithm="leveling_thinning",
                    warpx_resampling_trigger_max_avg_ppc=args.resample_ppc, warpx_resampling_min_ppc=10)
    electrons = picmi.Species(particle_type="electron", name="electrons", initial_distribution=[init_e, inj_e],
                              **save, **wall, **resample)
    ions = picmi.Species(particle_type="Ar", name="ions", charge="q_e", mass=M_AR, initial_distribution=init_i,
                         **save, **resample)

    # ------------------------------------------------------------------ collisions
    nn_expr, nn_consts, n_ch = neutral_density_expr(args, area)
    n_max = n_ch + args.n_facility
    coll_e = picmi.MCCCollisions(
        name="coll_e", species=electrons, background_density=subst(nn_expr, nn_consts), background_temperature=args.t_gas,
        background_mass=M_AR, max_background_density=n_max, ndt_supercycle=args.mcc_e_every,
        scattering_processes={
            "elastic": {"cross_section": str(XS / "electron_scattering.dat")},
            "excitation1": {"cross_section": str(XS / "excitation_1.dat"), "energy": E_EXC},
            "ionization": {"cross_section": str(XS / "ionization.dat"), "energy": E_IZ, "species": ions}})
    coll_i = picmi.MCCCollisions(
        name="coll_i", species=ions, background_density=subst(nn_expr, nn_consts), background_temperature=args.t_gas,
        background_mass=M_AR, max_background_density=n_max, ndt_supercycle=args.mcc_i_every,
        scattering_processes={
            "elastic": {"cross_section": str(XS / "ion_scattering.dat")},
            "elastic_back": {"cross_section": str(XS / "ion_back_scatter.dat"), "scattering_angle_model": "backward"}})
    k_anom = 1e-13                                              # [m^3/s] arbitrary; density = nu / k
    write_anomalous_xs(out / "anomalous_xs.dat", k_anom)
    nu_expr, nu_consts, nu_max = anomalous_expr(args)
    coll_anom = picmi.MCCCollisions(
        name="coll_anom", species=electrons, background_density=subst(f"({nu_expr})/k_anom", {**nu_consts, "k_anom": k_anom}),
        background_temperature=300.0, background_mass=1e-15, max_background_density=nu_max / k_anom,
        ndt_supercycle=args.mcc_e_every,
        scattering_processes={"elastic": {"cross_section": str(out / "anomalous_xs.dat")}})

    # ------------------------------------------------------------------ simulation
    n_diag = max(1, int(round(args.diag_every / args.dt)))
    diag = picmi.FieldDiagnostic(name="fields", grid=grid, period=n_diag, write_dir="diags",
                                 data_list=["phi", "Ex" if planar else "Er", "Ez", "rho_electrons", "rho_ions"],
                                 warpx_format="openpmd", warpx_openpmd_backend="h5")
    pnum = picmi.ReducedDiagnostic(diag_type="ParticleNumber", name="pnum", period=args.monitor_every, path="diags/")
    sim = picmi.Simulation(solver=solver, time_step_size=args.dt, max_steps=nsteps, verbose=1,
                           particle_shape=1, warpx_collisions=[coll_e, coll_i, coll_anom],
                           warpx_use_filter=False, warpx_random_seed=args.seed)
    sim.add_applied_field(picmi.LoadAppliedField(read_fields_from_path=str(BFILES[args.geometry]), load_E=False,
                                                 load_B=True))
    sim.add_species(electrons, layout=[picmi.PseudoRandomLayout(n_macroparticles_per_cell=args.ppc, grid=grid),
                                       picmi.PseudoRandomLayout(n_macroparticles_per_cell=args.ppc_inject, grid=grid)])
    sim.add_species(ions, layout=picmi.PseudoRandomLayout(n_macroparticles_per_cell=args.ppc, grid=grid))
    sim.add_diagnostic(diag)
    sim.add_diagnostic(pnum)

    meta = {"args": vars(args), "nr": nr, "nz": nz, "dr": dr, "dz": dz, "r_lo": r_lo, "r_hi": r_hi,
            "z_lo": z_lo, "z_hi": z_hi, "area_m2": area, "steps": nsteps, "n_neutral_channel": n_ch,
            "neutral_density_expr": nn_expr, "neutral_consts": nn_consts, "nu_anom_expr": nu_expr,
            "nu_anom_consts": nu_consts, "nu_anom_max": nu_max, "flux_inject": flux_inj,
            "mass_equivalent_current_A": args.mdot / M_AR * Q_E, "geometry": args.geometry, "l_y": l_y}
    (out / "meta.json").write_text(json.dumps(meta, indent=2, default=str))
    print(f"{args.geometry} grid {nr} x {nz}, dr {dr * 1e3:.4f} mm, dz {dz * 1e3:.4f} mm, dt {args.dt:.3g} s, {nsteps} steps; "
          f"neutrals {n_ch:.3g} m^-3 in the channel; equivalent current {meta['mass_equivalent_current_A']:.3f} A")

    # ------------------------------------------------------------------ callbacks
    inv_g2 = 1.0 / args.gamma2
    state = {"rho_scaled": 0}

    def scale_rho():
        """Permittivity scaling: the Poisson solve sees rho / gamma2, i.e. eps0 -> gamma2 * eps0."""
        mf = sim.fields.get("rho_fp", level=0)
        ng = mf.n_grow_vect
        mf.mult(inv_g2, 0, mf.n_comp, max(ng[0], ng[1]))
        state["rho_scaled"] += 1

    callbacks.installafterdeposition(scale_rho)

    import warnings
    warnings.filterwarnings("ignore", category=UserWarning, module="pywarpx")
    warnings.filterwarnings("ignore", message=".*deprecated.*")
    pcs = {s: particle_containers.ParticleContainerWrapper(s) for s in ("electrons", "ions")}
    buf = particle_containers.ParticleBoundaryBufferWrapper()
    bounds = ("x_lo", "x_hi", "z_lo", "z_hi")
    mon_path = out / "monitor.csv"
    fields_cols = ["step", "t_us", "n_e", "n_i"] + [f"I_{s}_{b}" for s in ("e", "i") for b in bounds] + [
        "I_discharge", "phi_mid_anode5", "phi_mid_exit", "phi_mid_plume15", "wall_s"]
    with open(mon_path, "w", newline="") as fh:
        csv.writer(fh).writerow(fields_cols)
    zb = np.linspace(z_lo, z_hi, nz + 1)
    prof = {k: np.zeros(nz) for k in ("we", "wi", "e_vz", "e_vr", "e_vt", "e_v2", "i_vz", "i_v2")}
    prof["samples"] = 0
    state["t0"] = state["last"] = time.time()

    def phi_line():
        from pywarpx import fields
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            phi = fields.PhiFPWrapper(0)[...]
        return phi[nr // 2, :]

    def monitor():
        step = sim.extension.warpx.getistep(0)
        if step % args.monitor_every:
            return
        dt_win = args.monitor_every * args.dt
        cur = {}
        for s, q in (("e", -Q_E), ("i", Q_E)):
            name = "electrons" if s == "e" else "ions"
            for b in bounds:
                if args.walls == "reflect" and b.startswith("x"):
                    cur[f"I_{s}_{b}"] = 0.0                     # reflecting: nothing is scraped there
                    continue
                w = buf.get_particle_boundary_buffer(name, b, "w", 0) if buf.get_particle_boundary_buffer_size(name, b) else []
                cur[f"I_{s}_{b}"] = l_y * q * float(sum(np.sum(x) for x in w)) / dt_win     # total current
        buf.clear_buffer()
        # discharge current: conventional current from the plasma into the anode, taken positive
        i_d = -(cur["I_e_z_lo"] + cur["I_i_z_lo"])
        line = phi_line()
        zn = np.linspace(z_lo, z_hi, line.size)
        pick = lambda zz: float(np.interp(zz, zn, line))                          # noqa: E731
        row = [step, step * args.dt * 1e6, pcs["electrons"].get_particle_count(), pcs["ions"].get_particle_count()]
        row += [cur[c] for c in fields_cols[4:12]] + [i_d, pick(z_lo + 5e-3), pick(0.0), pick(15e-3),
                                                       time.time() - state["t0"]]
        with open(mon_path, "a", newline="") as fh:
            csv.writer(fh).writerow(row)
        if time.time() - state["last"] > 30 or step == nsteps:
            state["last"] = time.time()
            el = time.time() - state["t0"]
            print(f"[{step}/{nsteps}] t {step * args.dt * 1e6:.3f} us, {el / max(step, 1) * 1e3:.2f} ms/step, "
                  f"N_e {row[2]}, N_i {row[3]}, I_d {i_d:.3f} A, I_beam {cur['I_i_z_hi']:.3f} A, "
                  f"phi exit {pick(0.0):.1f} V, rho scaled {state['rho_scaled']}x", flush=True)

    def profiles():
        step = sim.extension.warpx.getistep(0)
        if step % args.profile_every or step < args.profile_start * nsteps:   # late part of the run only
            return
        for s in ("electrons", "ions"):
            pc = pcs[s]
            z = np.concatenate(pc.get_particle_z(copy_to_host=True) or [np.zeros(0)])
            if z.size == 0:
                continue
            w = np.concatenate(pc.get_particle_weight(copy_to_host=True))
            ux = np.concatenate(pc.get_particle_ux(copy_to_host=True))
            uy = np.concatenate(pc.get_particle_uy(copy_to_host=True))
            uz = np.concatenate(pc.get_particle_uz(copy_to_host=True))
            if planar:                                         # x = r, y = azimuth
                vr, vt = ux, uy
            else:
                th = np.concatenate(pc.get_particle_theta(copy_to_host=True))
                vr = ux * np.cos(th) + uy * np.sin(th)
                vt = -ux * np.sin(th) + uy * np.cos(th)
            k = np.clip(np.searchsorted(zb, z) - 1, 0, nz - 1)
            bc = lambda v: np.bincount(k, weights=v, minlength=nz)                 # noqa: E731
            if s == "electrons":
                prof["we"] += bc(w)
                prof["e_vz"] += bc(w * uz)
                prof["e_vr"] += bc(w * vr)
                prof["e_vt"] += bc(w * vt)
                prof["e_v2"] += bc(w * (ux * ux + uy * uy + uz * uz))
            else:
                prof["wi"] += bc(w)
                prof["i_vz"] += bc(w * uz)
                prof["i_v2"] += bc(w * (ux * ux + uy * uy + uz * uz))
        prof["samples"] += 1
        np.savez(out / "profiles.npz", z_edges=zb, **prof)

    callbacks.installafterstep(monitor)
    callbacks.installafterstep(profiles)

    t0 = time.time()
    done, stop = 0, ""
    while done < nsteps:
        n = min(args.chunk, nsteps - done)
        sim.step(n)
        done += n
        ntot = pcs["electrons"].get_particle_count() + pcs["ions"].get_particle_count()
        if ntot > args.max_particles:
            stop = f"stopped early: {ntot} macroparticles > {args.max_particles:g}"
        elif time.time() - t0 > args.max_wall:
            stop = f"stopped early: wall time > {args.max_wall:g} s"
        if stop:
            print(stop, flush=True)
            break
    wall = time.time() - t0
    meta.update(wall_s=wall, ms_per_step=wall / done * 1e3, steps_done=done, t_end_reached=done * args.dt,
                stop_reason=stop or "completed")
    (out / "meta.json").write_text(json.dumps(meta, indent=2, default=str))
    print(f"done: {done} steps ({done * args.dt * 1e6:.3f} us) in {wall:.0f} s ({wall / done * 1e3:.2f} ms/step); "
          f"output in {out}")


if __name__ == "__main__":
    main()
