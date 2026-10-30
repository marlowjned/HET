# Breathing-mode PIC: handoff (2026-10-30)

Read this first when resuming. It records where the breathing-mode work stands, what is verified, and the next
steps in order. Everything lives in the WSL clone `~/Projects/HET` (Ubuntu-24.04; from Windows
`\\wsl.localhost\Ubuntu-24.04\home\marlow\Projects\HET`), not in the Windows clone `C:\Users\marlo\Projects\HET`.

## Goal and ground rules

- **Goal:** a 2-D PIC (WarpX) run of BPL-1000 Design B at the set point that shows the breathing mode, compared
  with a 1-D reference (HallThruster.jl) on the same inputs, to validate the tools.
- **Set point:** krypton, 3.0 mg/s, 300 V, 3.53 A, 1060 W, 56 mN; channel r_in 20.5 mm, r_out 33.5 mm
  (width 13 mm); magnet map scaled to 165 G peak.
- **Do not put this in bplsim yet** (user, 2026-10-29): make the standalone example work first.
- **How the user wants the work done** (memory notes): plan and wait for approval before big decisions; say why each
  sub-task moves the goal and do not polish what does not; one bounded attempt at incidental breakage, then report;
  diagnose against physics before fixing; recompute every number before presenting. A reduced model will not hit the
  set point exactly: take the closest case, state the residual, keep going.

## Environment

- Drive WSL from the Windows shell: `MSYS_NO_PATHCONV=1 wsl.exe -e bash -c '...'`. Files written from Windows get
  CRLF: run `sed -i "s/\r$//" <file>` before using them in WSL. Never `pkill -f <pattern>` (it matches its own
  shell): kill by PID. Background jobs die when their `wsl.exe` exits. Use a tracked background task (up to 2 h),
  or, for longer runs, a detached Windows process (see "Long run" below). `/tmp` in WSL can be wiped.
- **CPU WarpX:** `source ~/Projects/HET/pic/activate.sh` (micromamba env `hetpic`: WarpX 26.10, numpy, scipy,
  matplotlib).
- **GPU WarpX:** `bash pic/pic2d/run_gpu.sh <het_breathing args>` (env `warpx-gpu`). It is WarpX 26.10 built from
  `~/tools/warpx_src` in `~/tools/warpx_build_gpu`:
  - CUDA 12.9, sm_120 (RTX 5080, 16 GB), dims `2;RZ`;
  - fields in double, particles in **single** precision;
  - RZ Poisson patch for r_lo > 0;
  - cupy 14.2 installed; `CUDA_PATH=<env>/targets/x86_64-linux` is set by `run_gpu.sh`, and cupy's JIT needs it;
  - py-spy installed, but ptrace_scope = 1, so profile only as the parent: `py-spy record --native -- <cmd>`.
- **Julia:** `~/.juliaup/bin/julia`, HallThruster.jl 0.23.2 in the project env `pic/hallthruster`.

## Results so far (verified)

### 1-D reference (`pic/hallthruster`, `ht_case.py`)

| run | c1 / c2 | walls | I_d [A] | thrust [mN] | breathing |
|---|---|---|---|---|---|
| `ref` (**the reference**) | 1/800, 1/16 | BN WallSheath, scale 1 | 3.97 (std 0.59) | 59.0 | 14 kHz, 15 % |
| `ref_nowall` | 1/800, 1/16 | none (`--wall-scale 0`) | 4.22 (std 2.63) | 76.8 | 15 kHz, 62 % |
| `base` (earlier scan) | 1/160, 1/16 | BN, scale 1, Kr2+ on | 5.80 | 77.8 | 8 kHz, 20 % |

- Both `ref` runs use Kr+ only (`--max-charge 1`), as the PIC does.
- The 1-D model has no plume divergence, so its thrust is axial momentum. With a typical 0.9 factor, `ref` gives
  about 53 mN against 56 mN.
- No transport setting reaches 3.53 A: the 41.4 mm channel ionizes about 100 % of the flow.
- `ref` profiles:
  - time-averaged n_e peaks at 1.15e18 m^-3 at z = -12.8 mm (exit plane at 0), with instantaneous peaks of
    6.5e18 m^-3;
  - Te peaks at 24.5 eV at z = -3.2 mm;
  - neutral density is 5.09e19 m^-3 at the anode face and 2.1e17 m^-3 at the exit plane.
- From a cold start the 1-D model needs about 400 us to settle: a 72 A ignition spike at 50-100 us, then relaxation
  swings.
- Transport scan: `scan.csv`, `results/scan.png`.

### Krypton cross sections (`pic/pic2d/xs/make_kr.py` -> `xs/Kr/`, `kr_fit.png`, `kr_fit.txt`)

They are fitted to HallThruster.jl's own krypton rate tables, so the PIC and the 1-D reference share the same rates.
The tables are indexed by eps = 3/2 Te with a Maxwellian average. Max rate error over Te 2-80 eV:

- **Ionization:** threshold law sigma = 1.006e-19 (x - 1)^1.604 / x^2.054 with x = E/14.0 eV; 4.4 %; peak 3.4e-20 m^2
  at 64 eV; between argon and xenon at every energy checked.
- **Lumped excitation** (9.99 eV): smooth non-negative curve with penalty lam = 10, the weakest that gives a single
  peak; up to 15 % at Te 2-3 eV, at most 4.2 % for Te 4-80 eV.
- **Elastic** (momentum transfer): 0.7 %.
- **Method check:** the same fit on xenon recovers the measured warpx-data ionization cross section within 10 %.
- No ion-neutral collisions, as in HallThruster.jl.

### PIC script (`pic/pic2d/het_breathing.py`, planar x-z, 40 x 240 cells, dt 5e-11 s, gamma^2 = 400)

- **Kinetic neutrals:** DSMC for ionization and excitation, so the neutrals deplete. The quasi-neutral cathode adds
  electrons in a slab 5 mm upstream of the cathode plane.
- **Anomalous transport as in the 1-D model:** nu = c(z) omega_ce, with c going linearly from 1/800 to 1/16 over
  4.1 mm (0.1 L_ch) centred on the exit plane.
- **Sheath walls** (default):
  - an electron at a radial wall is reflected if its normal energy is below a cut-off, otherwise absorbed and
    re-emitted cold (2 eV) at the same place, so energy is lost but the electron is kept, as in the 1-D wall model;
  - the cut-off, 121.1 V, loses the 1-D BN WallSheath power at Te = 20 eV (`sheath_cut`);
  - ions and neutrals reflect.
- **Start from the 1-D state** (`--init-1d ref`, default):
  - the `ref` time-averaged profiles (neutral density, n_e, Te, ion velocity) enter as piecewise-linear parser
    expressions with 1 mm knots, built as a balanced if-tree;
  - the same neutral profile is the fixed background of the elastic MCC. The old undepleted background overstated
    the classical mobility at the exit plane by about 250 times;
  - neutrals start drifting at vbar/2.
- **Checkpoints and restart:** checkpoints every 20 us (`--checkpoint-every`); continue with `--restart
  ~/Projects/HET/pic/pic2d/runs/<tag>/chk/chkNNNNNN` and otherwise the same arguments. Use an **absolute** path:
  the script changes into `runs/<tag>` before resolving it, a bug listed below. The monitor appends on restart.
- **Monitor** `runs/<tag>/monitor.csv`, one record per 200 steps: I_discharge, I_i_beam, I_e_wall, P_wall_W,
  neutral mass, phi at the exit plane.
- **Plots:** `python plot_breathing.py <tag> --skip 20` writes `results/breathing_<tag>.png`: I_d(t), neutral
  inventory, the 1-D reference, and both spectra.

### Pilot runs (`pic/pic2d/runs/`)

- `kr_pilot_cpu`, 2 us on CPU at 24 ms/step: I_d fell from about 5 A to 3.2-3.4 A, ion beam 3.2-3.9 A, wall loss
  6-36 W, neutral mass steady at 0.372 ug. There is a 68 A electron surge in the first 5 ns. The run is far shorter
  than one breathing period (about 70 us).
- `kr_gpu_smoke`: the GPU matches the CPU at 5 and 10 ns (I_d 68.5 / 4.80 A against 68.1 / 4.70 A on CPU).
- Stale runs from the GPU hang, safe to delete: `kr_gpu_noinit`, `kr_gpu_prof`, and the old `chk` in
  `kr_gpu_smoke`.

## Fixed on 2026-10-30: GPU hang at step 1

- **Symptom:** runs stopped after the first "Setting up Monte-Carlo collisions" message, with the CPU at 100 % and
  the GPU idle.
- **Cause:** WarpX `BackgroundMCCCollision::get_nu_max` scans 0-5000 eV in steps of the smallest spacing of each
  table, with `E += E_step` in particle precision (single precision on the GPU build). The anomalous-collision table
  (`het2d.write_anomalous_xs`) had a 1.7e-5 eV step. Past 512 eV that is below half a float spacing, so E never
  grows. On CPU (double precision) the scan ends.
- **Fix:** the anomalous table is now evenly spaced at 0.05 eV (`het2d.py`).
- **Rule:** every MCC/DSMC table needs its smallest energy step above about 5e-4 eV. The krypton tables use 0.01 eV.
- This edit was applied by a command that was interrupted before its GPU rerun finished; the rerun (`kr_gpu_smoke`)
  then confirmed the fix.

## Open problem: GPU speed

- The GPU run takes **about 395 ms/step**, against **24-25 ms/step on CPU** (10 threads). Measured from the
  `wall_s` column of `kr_gpu_smoke` between steps 100 and 200; the same in both windows, so it is not start-up
  loading.
- A 0.4 ms breathing run is 8 M steps: about 53 h on CPU, impossible at the current GPU speed.
- Suspects, cheapest first:
  1. **Small boxes.** `warpx_max_grid_size=32` gives 16 boxes of 16 x 16 to 24 x 32 cells. On a GPU every
     operation is a kernel launch per box, and WSL2 kernel launches are slow. Add a `--max-grid-size` argument and
     try 256 (one 40 x 240 box).
  2. **Python callbacks.** `cathode()` copies every electron and ion z and weight to the host every 20 steps. Time
     a run with the callbacks off; if they cost, sum the slab on the device (cupy, `copy_to_host=False`).
  3. **Managed memory or pinned copies under WSL2.** Run a few hundred steps to the end and read AMReX's
     TinyProfiler table at finalize, if the build prints one. Otherwise use `py-spy record --native` as the
     parent, or Nsight Systems.
- **Decision point for the user:** if the GPU cannot get well below the CPU, the choices are a multi-day CPU run
  (more threads, a larger dt if the plasma frequency allows), a coarser or shorter set-up, or more GPU work. Present
  this as a plan; do not launch a multi-day run unasked.

## Next steps (in order)

1. **GPU speed:** the suspects above, bounded to about an hour. Report the speed-up before going further.
2. **Settling check:** run about 50-100 us from `--init-1d ref` (GPU if fast, otherwise CPU) and check that the
   discharge stays near the 1-D state:
   - I_d about 4 A and ion beam about 3.4 A at most (the mass-flow equivalent is 3.454 A);
   - Te not running away;
   - wall power of the order of the 1-D figure.
   Also look at the field dumps for the mid-channel band seen in the earlier planar PIC (cause unknown).
3. **Long run**, at least 300-400 us (3-4 periods at 14 kHz after settling), started detached so it outlives the
   shell:
   ```
   Start-Process wsl.exe -WindowStyle Hidden -ArgumentList '-e','bash','-lc','bash ~/Projects/HET/pic/pic2d/run_gpu.sh --tag breathing_kr --t-end 400e-6 --diag-every 2e-6 --max-wall 200000 --max-particles 2e7 > ~/Projects/HET/pic/pic2d/runs_breathing_kr.log 2>&1'
   ```
   For CPU, swap `run_gpu.sh` for `source activate.sh; python het_breathing.py`. Checkpoints land every 20 us.
4. **Compare** with `plot_breathing.py`: breathing frequency and amplitude against `ref` (14 kHz, 15 %). The PIC's
   walls lie between `ref` and `ref_nowall`, so bracket with both.
5. Write up in `pic/REPORT.md` with plots; only then discuss moving it into bplsim.

## Known differences, PIC against the 1-D reference

- Planar x-z geometry: no 1/r. Radial walls run the full length of the domain, as the 1-D model's walls do (it keeps
  wall losses in the plume, `electron_plume_loss_scale = 1`).
- Elastic collisions use a fixed neutral background, the `ref` time average, so the breathing modulation of the
  classical mobility is missing. Ionization and excitation do see the depleting particle neutrals.
- The sheath cut-off is a single voltage matched at Te = 20 eV: it loses too little where Te is low.
- Kinetic electrons give a non-Maxwellian EEDF, so the rates differ from the Maxwellian tables by design.
- Permittivity is scaled (gamma^2 = 400); the anomalous profile uses |B| at the mid-radius (fit), while the 1-D model
  uses channel-averaged B_r (within about 5 % near the peak).

## Smaller items

- Bug: `main()` resolves `--restart` after `os.chdir(out)`, so relative restart paths break. Resolve it before the
  chdir. The restart path has not been exercised yet: test it on a short run before relying on it.
- `propellant()` in `het_breathing.py` still points to REPORT.md when the krypton tables are missing; it should say
  `python xs/make_kr.py`.
- The RZ patch in the GPU build is untested (vacuum linear-profile check). It is not on the critical path.
- Defaults set to the reference: `--alpha-in 1/800`, `--alpha-out 1/16`, `--walls sheath`, `--init-1d ref`.
