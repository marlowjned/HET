"""2-D PIC-MCC/DSMC model of Design B for breathing-mode studies (WarpX; no bplsim).

What it adds to het2d.py (the first example):
- **Kinetic neutrals.** Injected at the anode face with the propellant mass flow. Excitation and ionization are
  particle-particle collisions (DSMC), so the propellant depletes and refills, which is what breathing is.
  Elastic electron and ion collisions use a fixed neutral background (MCC) of the undepleted profile. Particle
  elastic collisions would split a new, roughly 1000x lighter neutral off every event for a negligible recoil:
  about 8000 extra neutral macroparticles per 100 steps in the argon test.
- **A quasi-neutral cathode.** Electrons are added in a thin slab upstream of the cathode plane wherever ions
  outnumber electrons there, so no virtual cathode forms.
- **A propellant table.** Krypton is the set point, once its cross sections are in xs/Kr; argon (warpx-data) is
  there for tests.

Geometry (set-point channel), magnet map (165 G), radial walls, anomalous transport and permittivity scaling are
as in het2d.py.

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
from het2d import BFILES, KB, M_E, Q_E, anomalous_expr, subst, write_anomalous_xs  # noqa: E402

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
    a("--n0", type=float, default=5e16, help="initial plasma density [m^-3] (LANDMARK start)")
    a("--te0", type=float, default=10.0)
    a("--ppc", type=float, default=20, help="initial plasma macroparticles per cell and species")
    a("--ppc-neutral", type=float, default=30, help="neutral macroparticles per cell at the channel density")
    a("--resample-ppc", type=int, default=80)
    a("--cathode-te", type=float, default=10.0, help="temperature of the cathode electrons [eV]")
    a("--inject-depth", type=float, default=5e-3, help="cathode slab centre, upstream of the cathode plane [m]")
    a("--slab", type=float, default=1e-3, help="cathode slab half-width [m]")
    a("--inject-every", type=int, default=20, help="steps between cathode corrections")
    a("--walls", choices=["reflect", "absorb"], default="reflect")
    a("--alpha-in", type=float, default=1 / 64)
    a("--alpha-out", type=float, default=1 / 16)
    a("--dsmc-e-every", type=int, default=2)
    a("--dsmc-i-every", type=int, default=20)
    a("--threads", type=int, default=12)
    a("--diag-every", type=float, default=0.5e-6)
    a("--monitor-every", type=int, default=200, help="steps between records (10 ns at the default dt)")
    a("--max-particles", type=float, default=8e6)
    a("--max-wall", type=float, default=6600.0)
    a("--chunk", type=int, default=2000)
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
    nn_expr = subst("n_ch*if(z <= 0, 1, (2/pi)*atan(w_slot/(2*z)))", dict(n_ch=n_ch, w_slot=w_slot, pi=pi))
    cell_v = dx * dz                                           # per metre of y
    w_n = n_ch * cell_v / args.ppc_neutral                     # neutral macroparticle weight (per metre of y)
    ppc_inj_n = flux_n * dx * args.dt / w_n                    # injected neutral macroparticles per cell per step
    resample = dict(warpx_do_resampling=True, warpx_resampling_algorithm="leveling_thinning",
                    warpx_resampling_trigger_max_avg_ppc=args.resample_ppc, warpx_resampling_min_ppc=10)
    save = {f"warpx_save_particles_at_{b}": True for b in ("zlo", "zhi")}
    if args.walls == "absorb":
        save.update({f"warpx_save_particles_at_{b}": True for b in ("xlo", "xhi")})
    electrons = picmi.Species(particle_type="electron", name="electrons",
                              initial_distribution=picmi.UniformDistribution(density=args.n0, rms_velocity=[vth_e] * 3),
                              **save, **resample)
    ions = picmi.Species(particle_type=args.propellant, name="ions", charge="q_e", mass=M,
                         initial_distribution=picmi.UniformDistribution(density=args.n0, rms_velocity=[vth_n] * 3),
                         **save, **resample)
    init_n = picmi.AnalyticDistribution(density_expression=nn_expr, rms_velocity=[vth_n] * 3)
    inj_n = picmi.UniformFluxDistribution(flux=flux_n, flux_normal_axis="z", surface_flux_position=z_lo + 0.05 * dz,
                                          flux_direction=+1, rms_velocity=[vth_n] * 3,
                                          gaussian_flux_momentum_distribution=True)
    neutrals = picmi.Species(particle_type=args.propellant, name="neutrals", charge=0, mass=M,
                             initial_distribution=[init_n, inj_n], warpx_do_not_deposit=True,
                             warpx_do_not_gather=True, warpx_reflection_model_zlo="1",   # the anode face re-emits
                             warpx_save_particles_at_zhi=True, warpx_do_resampling=True,
                             warpx_resampling_algorithm="leveling_thinning",
                             warpx_resampling_trigger_max_avg_ppc=int(3 * args.ppc_neutral),
                             warpx_resampling_min_ppc=10)

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
                                    background_temperature=args.t_gas, background_mass=M, max_background_density=n_ch,
                                    ndt_supercycle=args.dsmc_e_every, scattering_processes=e_mcc)
    i_proc = {name: {"cross_section": str(xs / f), "scattering_angle_model": model}
              for name, (f, model) in prop["ion"].items()}
    coll_in = picmi.MCCCollisions(name="coll_in", species=ions, background_density=nn_expr,
                                  background_temperature=args.t_gas, background_mass=M, max_background_density=n_ch,
                                  ndt_supercycle=args.dsmc_i_every, scattering_processes=i_proc)
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
    sim = picmi.Simulation(solver=solver, time_step_size=args.dt, max_steps=nsteps, verbose=1, particle_shape=1,
                           warpx_collisions=[coll_en, coll_e_el, coll_in, coll_anom], warpx_use_filter=False,
                           warpx_random_seed=args.seed)
    sim.add_applied_field(picmi.LoadAppliedField(read_fields_from_path=str(BFILES["xz"]), load_E=False, load_B=True))
    sim.add_species(electrons, layout=picmi.PseudoRandomLayout(n_macroparticles_per_cell=args.ppc, grid=grid))
    sim.add_species(ions, layout=picmi.PseudoRandomLayout(n_macroparticles_per_cell=args.ppc, grid=grid))
    sim.add_species(neutrals, layout=[picmi.PseudoRandomLayout(n_macroparticles_per_cell=args.ppc_neutral, grid=grid),
                                      picmi.PseudoRandomLayout(n_macroparticles_per_cell=ppc_inj_n, grid=grid)])
    sim.add_diagnostic(diag)

    i_mass = args.mdot / M * Q_E
    meta = {"args": vars(args), "propellant": {k: v for k, v in prop.items() if k != "dir"} | {"dir": str(prop["dir"])},
            "nr": nr, "nz": nz, "dx": dx, "dz": dz, "r_lo": r_lo, "r_hi": r_hi, "z_lo": z_lo, "z_hi": z_hi,
            "l_y": l_y, "area_m2": area, "steps": nsteps, "flux_neutral": flux_n, "n_neutral_channel": n_ch,
            "w_neutral": w_n, "ppc_inject_neutral": ppc_inj_n, "mass_equivalent_current_A": i_mass,
            "nu_anom_expr": nu_expr, "nu_anom_consts": nu_consts}
    (out / "meta.json").write_text(json.dumps(meta, indent=2, default=str))
    print(f"{args.propellant}: {nr} x {nz} cells ({dx * 1e3:.3f} x {dz * 1e3:.3f} mm), dt {args.dt:g} s, {nsteps} steps; "
          f"neutrals {n_ch:.3g} m^-3, mass-equivalent current {i_mass:.3f} A", flush=True)

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
            "I_n_out_A_equiv", "I_cathode_inj", "I_discharge", "phi_mid_exit", "wall_s"]
    with open(mon, "w", newline="") as fh:
        csv.writer(fh).writerow(cols)

    def wsum(name, b):
        if not buf.get_particle_boundary_buffer_size(name, b):
            return 0.0
        return float(sum(np.sum(x) for x in buf.get_particle_boundary_buffer(name, b, "w", 0)))

    def monitor():
        step = sim.extension.warpx.getistep(0)
        if step % args.monitor_every:
            return
        dtw = args.monitor_every * args.dt
        c = l_y * Q_E / dtw
        ie_a, ii_a = c * wsum("electrons", "z_lo"), c * wsum("ions", "z_lo")
        ie_c, ii_c = c * wsum("electrons", "z_hi"), c * wsum("ions", "z_hi")
        in_out = c * wsum("neutrals", "z_hi")
        buf.clear_buffer()
        i_inj = c * st["inj_w"]
        st["inj_w"] = 0.0
        wn = pcs["neutrals"].get_particle_weight(copy_to_host=True)
        mass_ug = (float(sum(np.sum(x) for x in wn)) if wn else 0.0) * l_y * M * 1e9
        from pywarpx import fields
        phi = fields.PhiFPWrapper(0)[...]
        zn = np.linspace(z_lo, z_hi, phi.shape[1])
        row = [step, step * args.dt * 1e6, pcs["electrons"].get_particle_count(), pcs["ions"].get_particle_count(),
               pcs["neutrals"].get_particle_count(), mass_ug, ie_a, ii_a, ie_c, ii_c, in_out, i_inj, ie_a - ii_a,
               float(np.interp(0.0, zn, phi[phi.shape[0] // 2])), time.time() - st["t0"]]
        with open(mon, "a", newline="") as fh:
            csv.writer(fh).writerow(row)
        if time.time() - st["last"] > 30 or step == nsteps:
            st["last"] = time.time()
            print(f"[{step}/{nsteps}] t {row[1]:.3f} us, {(time.time() - st['t0']) / max(step, 1) * 1e3:.2f} ms/step, "
                  f"N_e {row[2]}, N_i {row[3]}, N_n {row[4]}, neutral mass {mass_ug:.3f} ug, I_d {row[12]:.3f} A, "
                  f"I_beam {ii_c:.3f} A, cathode {i_inj:.3f} A", flush=True)

    callbacks.installafterstep(cathode)
    callbacks.installafterstep(monitor)

    t0 = time.time()
    done, stop = 0, ""
    while done < nsteps:
        n = min(args.chunk, nsteps - done)
        sim.step(n)
        done += n
        tot = sum(pcs[s].get_particle_count() for s in pcs)
        if tot > args.max_particles:
            stop = f"stopped early: {tot} macroparticles"
        elif time.time() - t0 > args.max_wall:
            stop = f"stopped early: wall time > {args.max_wall:g} s"
        if stop:
            print(stop, flush=True)
            break
    wall = time.time() - t0
    meta.update(wall_s=wall, ms_per_step=wall / max(done, 1) * 1e3, steps_done=done, stop_reason=stop or "completed")
    (out / "meta.json").write_text(json.dumps(meta, indent=2, default=str))
    print(f"done: {done} steps ({done * args.dt * 1e6:.3f} us) in {wall:.0f} s; output in {out}")


if __name__ == "__main__":
    main()
