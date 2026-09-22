"""
nl_watchdog.py -- watch a nonlinear GetDP solve's Picard residuals and kill
it if a load stage stalls, instead of letting it burn hours at a residual
that isn't going anywhere (the 2-stage ramp on the 0.750in core sat at rel
~0.11 for 60 iterations / 2.4 h before this existed).

Prints one line per finished stage and one every PRINT_EVERY iterations.
Kills getdp.exe when, within a stage:
  - after >= MIN_ITERS iterations, the best rel. residual of the last WINDOW
    iterations improved by less than MIN_GAIN x on the WINDOW before it,
    while still above STALL_FLOOR (the solver tolerance); or
  - the stage passes MAX_ITERS iterations still above MAX_ITERS_FLOOR.
A healthy stage here contracts ~2.3x per iteration (reaches 1e-4 in ~10).

Usage: python nl_watchdog.py <run_dir> [ramp fractions, comma-separated]
"""
import re
import subprocess
import sys
import time
from pathlib import Path

WINDOW, MIN_GAIN, MIN_ITERS = 5, 1.4, 8
STALL_FLOOR = 1e-4   # = the solver tolerance: a crawl anywhere above it is a stall
MAX_ITERS, MAX_ITERS_FLOOR = 25, 1e-3
PRINT_EVERY = 5
LINE = re.compile(r"IFrac step (\d+), NL iter (\d+): residual abs (\S+) rel (\S+)")


def kill(reason):
    print(f"STALL -> terminating getdp: {reason}", flush=True)
    subprocess.run(["taskkill", "/IM", "getdp.exe", "/F"], capture_output=True)
    sys.exit(2)


def main():
    log = Path(sys.argv[1]) / "getdp.log"
    fracs = [float(f) for f in sys.argv[2].split(",")] if len(sys.argv) > 2 else None
    while not log.exists():
        time.sleep(2)
    stage, hist, t_stage = None, [], time.time()
    with open(log, errors="replace") as fh:
        while True:
            line = fh.readline()
            if not line:
                time.sleep(1)
                continue
            if "Stopped" in line:
                if hist:
                    print(f"stage {stage} last rel {hist[-1]:.2e} after {len(hist) - 1} iters", flush=True)
                print("getdp finished", flush=True)
                return
            m = LINE.search(line)
            if not m:
                continue
            n, k, rel = int(m.group(1)), int(m.group(2)), float(m.group(4))
            if n != stage:
                # GetDP moves on at NL_iter_max whether or not the stage converged;
                # every later stage would then start from a wrong field.
                if stage is not None and hist[-1] > STALL_FLOOR:
                    kill(f"stage {stage} ended unconverged at rel {hist[-1]:.2e} (iteration cap)")
                if stage is not None:
                    print(f"stage {stage} ({fracs[stage] if fracs else '?'}) done: "
                          f"{len(hist) - 1} iters, rel {hist[-1]:.2e}, {(time.time() - t_stage) / 60:.0f} min",
                          flush=True)
                stage, hist, t_stage = n, [], time.time()
            hist.append(rel)
            if k and k % PRINT_EVERY == 0:
                print(f"stage {n} iter {k}: rel {rel:.2e}", flush=True)
            if k >= MIN_ITERS and len(hist) >= 2 * WINDOW + 1 and rel > STALL_FLOOR:
                recent = min(hist[-WINDOW:])
                before = min(hist[-2 * WINDOW:-WINDOW])
                if before / recent < MIN_GAIN:
                    kill(f"stage {n} iter {k}: rel {rel:.2e}, best improved only "
                         f"{before / recent:.2f}x over the last {WINDOW} iterations")
            if k >= MAX_ITERS and rel > MAX_ITERS_FLOOR:
                kill(f"stage {n} passed {MAX_ITERS} iterations at rel {rel:.2e}")


if __name__ == "__main__":
    main()
