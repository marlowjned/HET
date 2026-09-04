#!/usr/bin/env python3
"""
Hall thruster sizing tool.

Implements the scaling relations from Goebel & Katz, "Fundamentals of
Electric Propulsion: Ion and Hall Thrusters" (2008), Chapter 7. Equation
numbers in comments refer to that chapter (e.g. "7.2-6" = Eq. 7.2-6).

Inputs are the usual top-level design targets (discharge power, discharge
voltage, propellant) plus a set of efficiency/scaling assumptions that you'd
normally pull from comparable existing thrusters (SPT-100, NASA-173Mv2,
etc.) or refine iteratively. Defaults are representative mid-range values.

Note: the channel-length / ionization-mean-free-path sizing (Eqs. 7.2-14 to
7.2-19) requires the Maxwellian ionization rate coefficient <sigma_i*v_e>
for the propellant at the local Te, tabulated in Appendix E of the book.
That data isn't reproduced here, so this script does not compute L from
first principles -- it only performs the Larmor-radius magnetization check
(Eqs. 7.2-1, 7.2-3) using an assumed channel length.
"""

import math
from dataclasses import dataclass

# --- Physical constants (SI) ---
E_CHARGE = 1.602176634e-19   # C
AMU = 1.66053906660e-27      # kg
G0 = 9.80665                 # m/s^2
M_ELECTRON = 9.1093837015e-31  # kg

# --- Propellant atomic mass (amu) ---
PROPELLANTS = {
    "xenon": 131.293,
    "krypton": 83.798,
    "argon": 39.948,
}


@dataclass
class ThrusterTargets:
    Pd: float                  # discharge power, W
    Vd: float                  # discharge voltage, V
    propellant: str = "xenon"
    j: float = 0.12e4          # current density, A/m^2 (0.1-0.15 A/cm^2 typical, Sec 7.2.2)
    di_do_ratio: float = 0.6   # channel inner/outer diameter ratio
    eta_b: float = 0.70        # current utilization, Ib/Id            (Eq. 7.3-10)
    eta_v: float = 0.95        # voltage utilization, Vb/Vd             (Eq. 7.3-11)
    eta_m: float = 0.90        # mass utilization, mdot_i/mdot_p        (Eq. 7.3-14)
    eta_c: float = 0.90        # cathode flow fraction, mdot_a/mdot_p   (Eq. 7.3-5)
    eta_o: float = 0.95        # electrical utilization, Pd/Pin         (Eq. 7.3-7)
    gamma: float = 0.90        # thrust correction (divergence + multi-charge)
    Te_eV: float = 25.0        # electron temp in ionization region, eV (for Larmor check)
    B_radial_T: float = 0.015  # assumed peak radial field, T (150 G is the book's example)
    L_channel_m: float = None  # magnetized plasma length; defaults to channel width w


@dataclass
class ThrusterSizing:
    M_kg: float
    Id: float
    Ae_m2: float
    Do_m: float
    Di_m: float
    w_m: float
    R_mean_m: float
    Vb: float
    Ib: float
    T_N: float
    mdot_p_kgs: float
    mdot_a_kgs: float
    Isp_s: float
    eta_T: float
    eta_anode: float
    discharge_loss_eVion: float
    r_e_Larmor_m: float
    r_i_Larmor_m: float
    L_channel_m: float
    magnetization_ok: bool


def size_thruster(t: ThrusterTargets) -> ThrusterSizing:
    M_kg = PROPELLANTS[t.propellant.lower()] * AMU

    # Discharge current
    Id = t.Pd / t.Vd

    # Channel area from (roughly constant) current density -- Eq. 7.2-20 scaling
    Ae_m2 = Id / t.j

    # Ae = (pi/4)(Do^2 - Di^2), with Di = k*Do  =>  Do = sqrt(4*Ae / (pi*(1-k^2)))
    k = t.di_do_ratio
    Do_m = math.sqrt(4 * Ae_m2 / (math.pi * (1 - k ** 2)))
    Di_m = k * Do_m
    w_m = (Do_m - Di_m) / 2
    R_mean_m = (Do_m + Di_m) / 4

    # Beam voltage / current
    Vb = t.eta_v * t.Vd            # Eq. 7.3-11
    Ib = t.eta_b * Id              # Eq. 7.3-10

    # Thrust -- Eq. 7.3-9: T = gamma * sqrt(2M/e) * Ib * sqrt(Vb)
    T_N = t.gamma * math.sqrt(2 * M_kg / E_CHARGE) * Ib * math.sqrt(Vb)

    # Mass flow -- Eqs. 7.3-13, 7.3-14, 7.3-5
    mdot_i = (M_kg / E_CHARGE) * Ib        # ion mass flow (7.3-13, using Ib = eta_b*Id)
    mdot_p_kgs = mdot_i / t.eta_m          # total propellant flow (7.3-14)
    mdot_a_kgs = t.eta_c * mdot_p_kgs      # anode flow (7.3-5)

    # Specific impulse
    Isp_s = T_N / (mdot_p_kgs * G0)

    # Efficiencies -- Eqs. 7.3-15, 7.3-16
    eta_T = t.gamma ** 2 * t.eta_b * t.eta_v * t.eta_m * t.eta_o
    eta_anode = eta_T / (t.eta_o * t.eta_c)

    # Discharge loss, eV per beam ion -- Eq. 7.3-17
    discharge_loss_eVion = t.Pd * (1 - t.eta_b * t.eta_v) / Ib

    # Magnetization check -- Eqs. 7.2-1 (electron Larmor radius) and 7.2-3 (ion)
    L_channel_m = t.L_channel_m if t.L_channel_m is not None else w_m
    Te_J = t.Te_eV * E_CHARGE
    v_th = math.sqrt(8 * Te_J / (math.pi * M_ELECTRON))
    omega_ce = E_CHARGE * t.B_radial_T / M_ELECTRON
    r_e = v_th / omega_ce                              # Eq. 7.2-1

    v_i = math.sqrt(2 * E_CHARGE * Vb / M_kg)
    omega_ci = E_CHARGE * t.B_radial_T / M_kg
    r_i = v_i / omega_ci                                # Eq. 7.2-3

    # r_e << L and r_i >> L required; using order-of-magnitude margins (not exact
    # thresholds from the text, which only states "much less/greater than")
    magnetization_ok = (r_e < 0.3 * L_channel_m) and (r_i > 3 * L_channel_m)

    return ThrusterSizing(
        M_kg=M_kg, Id=Id, Ae_m2=Ae_m2, Do_m=Do_m, Di_m=Di_m, w_m=w_m,
        R_mean_m=R_mean_m, Vb=Vb, Ib=Ib, T_N=T_N, mdot_p_kgs=mdot_p_kgs,
        mdot_a_kgs=mdot_a_kgs, Isp_s=Isp_s, eta_T=eta_T, eta_anode=eta_anode,
        discharge_loss_eVion=discharge_loss_eVion, r_e_Larmor_m=r_e,
        r_i_Larmor_m=r_i, L_channel_m=L_channel_m,
        magnetization_ok=magnetization_ok,
    )


def report(t: ThrusterTargets, s: ThrusterSizing) -> None:
    print(f"=== Hall Thruster Sizing: {t.propellant}, {t.Pd:.0f} W @ {t.Vd:.0f} V ===\n")

    print("-- Channel geometry --")
    print(f"Discharge current Id     : {s.Id:.3f} A")
    print(f"Channel area Ae          : {s.Ae_m2 * 1e4:.2f} cm^2")
    print(f"Outer diameter Do        : {s.Do_m * 1e3:.1f} mm")
    print(f"Inner diameter Di        : {s.Di_m * 1e3:.1f} mm")
    print(f"Channel width w          : {s.w_m * 1e3:.2f} mm")
    print(f"Mean channel radius R    : {s.R_mean_m * 1e3:.1f} mm\n")

    print("-- Performance --")
    print(f"Beam voltage Vb          : {s.Vb:.1f} V")
    print(f"Beam current Ib          : {s.Ib:.3f} A")
    print(f"Thrust T                 : {s.T_N * 1e3:.2f} mN")
    print(f"Total mass flow rate     : {s.mdot_p_kgs * 1e6:.3f} mg/s")
    print(f"Anode mass flow rate     : {s.mdot_a_kgs * 1e6:.3f} mg/s")
    print(f"Specific impulse Isp     : {s.Isp_s:.0f} s\n")

    print("-- Efficiency --")
    print(f"Total efficiency eta_T   : {s.eta_T:.3f}")
    print(f"Anode efficiency eta_a   : {s.eta_anode:.3f}")
    print(f"Discharge loss           : {s.discharge_loss_eVion:.1f} eV/ion\n")

    print("-- Magnetization check (Eqs. 7.2-1, 7.2-3) --")
    print(f"Assumed channel length L : {s.L_channel_m * 1e3:.2f} mm")
    print(f"Electron Larmor radius   : {s.r_e_Larmor_m * 1e3:.3f} mm  (want << L)")
    print(f"Ion Larmor radius        : {s.r_i_Larmor_m * 1e3:.1f} mm  (want >> L)")
    print(f"Criteria satisfied       : {s.magnetization_ok}")
    if not s.magnetization_ok:
        print("  -> adjust B_radial_T, channel length, or geometry")


if __name__ == "__main__":
    # Example target: SPT-100-class thruster. Edit these values (and the
    # efficiency assumptions in ThrusterTargets) for your own design point.
    targets = ThrusterTargets(
        Pd=1350,        # W
        Vd=300,         # V
        propellant="xenon",
        eta_b=0.70,
        eta_v=0.95,
    )
    sizing = size_thruster(targets)
    report(targets, sizing)


