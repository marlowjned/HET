"""
generate_regions.py

Standalone re-generation of assembly_regions_generated.pro for an arbitrary
winding excitation (turns/current/polarity per pole type), WITHOUT
re-running the expensive geometry+meshing step in build_assembly.py. Reads
winding geometry (pole type, A_cross, loc) from assembly_params.txt, which
build_assembly.py already writes once per mesh.

This is the SOLE generator of the Group/Function block now -- build_assembly.py
used to duplicate this same generation inline, but that got hard to justify
once the Iron material definition below grew to a hand-copied B-H curve
(see AISI1008_H/AISI1008_B); build_assembly.py now calls
generate_regions_pro() directly instead of keeping a second copy in sync.
magnetostatics_assembly.pro is untouched -- it just `Include`s whatever
this writes.

Pole names are "inner_coil"/"outer_coil" (the real Onshape/STEP body
names) since build_assembly.py switched to importing real coil solids --
see that script's docstring for why windings no longer carry a zdir (every
winding's axis is confirmed global Y directly from geometry, so current
sign is one constant per pole TYPE, not a per-occurrence flip).

Iron is a NONLINEAR material. The hardware is **ASTM A36**; the B-H table
below is **AISI 1008 data used as a proxy**, since no numeric A36 curve
was obtainable. It is deliberately not relabelled -- see the block above
AISI1008_H for why, which way the error runs, and roughly how big it is.
This replaced, in order: the original mur_iron=1000 linear placeholder,
then GetDP's bundled generic "SteelGeneric" dataset.

The reluctivity-interpolation derivation (nu as a function of B^2, built
via ListAlt[] + InterpolationLinear[SquNorm[$1]]{...}) is the pattern
GetDP's own Lib_Materials.pro uses, and that its examples/magnet.pro
exercises for a real nonlinear a-formulation solve -- not invented for
this project. Air/Windings stay linear (mu_r=1, no argument needed); only
Iron's nu[] definition takes an argument, which is why the Formulation
(magnetostatics_assembly.pro) splits the single Galerkin term over Domain
into a linear part (Air+Windings, nu[]) and a nonlinear part (Iron,
nu[{d a}]) -- same split GetDP's own Lib_Magnetostatics_a_phi.pro template
uses for Vol_L_Mag vs Vol_NL_Mag.
"""
import math

from params_utils import parse_params

MU0 = 4 * math.pi * 1e-7

# Kept in sync with build_assembly.py's EXCITATION -- see its comment for
# how 1.30 A was derived from the 300 G channel-exit target.
DEFAULT_EXCITATION = {
    "inner_coil": dict(turns=260, current=1.30, polarity=-1),
    "outer_coil": dict(turns=165, current=1.30, polarity=+1),
}

# DC magnetization curve. H in A/m, B in T.
#
# *** THIS IS AISI 1008 DATA USED AS A PROXY FOR ASTM A36. ***
#
# The hardware is A36 (SendCutSend hot-rolled pickled-and-oiled, the only
# grade they stock thick enough for the 0.375in plates). No numeric A36
# B-H table could be obtained -- the published A36 magnetization work is
# transformer-tank literature behind paywalls -- so the 1008 curve below
# stands in for it, deliberately NOT relabelled, because a curve wearing
# the wrong material's name is exactly the failure mode this project
# already hit once.
#
# Direction of the error: A36 carries up to 0.26% carbon against 1008's
# ~0.08%, plus up to 1.2% manganese. More carbon means more pearlite and
# more domain-wall pinning, so real A36 is LESS permeable and more
# coercive than this curve. The model is therefore optimistic about the
# iron. Offsetting that slightly, A36 here is hot-rolled rather than
# cold-rolled, which leaves a more relaxed microstructure.
#
# Size of the error: switching the generic GetDP curve for this one
# dropped mu_r at the operating point 2.7x and moved the channel field
# 3.9%. A36 is plausibly another ~1.8x less permeable again, which by the
# same empirical sensitivity suggests a further few percent -- i.e. the
# real design current is probably nearer 1.35 A than 1.30 A. That is an
# extrapolation, not a result.
#
# PROVENANCE of the 1008 points themselves: they originate in the Ansys
# Maxwell SV material library and were transcribed via a public
# engineering forum -- the primary source could not be re-fetched directly
# (403). They are NOT from a mill certificate for any real stock.
#
# What gives reasonable confidence they are genuine digitized data rather
# than something invented: every H value is an exact multiple of one
# Oersted (2, 4, 6, 8, 10, 20, 40, 60, 80, 100, 200, 400, 600, 800, 1000,
# 2000, 4000, 5000 Oe), which is how pre-SI magnetization tables were
# tabulated. The curve is monotonic in both variables, differential
# permeability rises then falls as a real magnetization curve must, and B
# at 1000 Oe is 2.165 T -- squarely in the 2.1-2.2 T band expected of
# low-carbon steel.
#
# The important caveat, which applies to A36 at least as strongly: neither
# grade is an electrical steel, so ASTM imposes no magnetic requirement on
# either. Two suppliers' stock can differ substantially and processing
# history dominates. A measured curve on the actual plate is the only way
# to do better than this.
AISI1008_H = [
    0.0, 159.2, 318.3, 477.5, 636.6, 795.8, 1591.5, 3183.1, 4774.6, 6366.2,
    7957.7, 15915.5, 31831.0, 47746.5, 63662.0, 79577.5, 159155.0, 318310.0, 397887.0,
]
AISI1008_B = [
    0.0, 0.2402, 0.8654, 1.1106, 1.2458, 1.3310, 1.5000, 1.6000, 1.6830, 1.7410,
    1.7800, 1.9050, 2.0250, 2.0850, 2.1300, 2.1650, 2.2800, 2.4850, 2.5850,
]
assert len(AISI1008_H) == len(AISI1008_B) == 19
assert all(AISI1008_H[i] < AISI1008_H[i + 1] for i in range(18)), "H must increase"
assert all(AISI1008_B[i] < AISI1008_B[i + 1] for i in range(18)), "B must increase"


def _gmsh_list(values):
    return "{" + ", ".join(repr(v) for v in values) + "}"


def generate_regions_pro(params, excitation, out_path="assembly_regions_generated.pro"):
    tag = params["TAG"]
    winding_names = [f"Winding{w['index']}" for w in params["windings"]]

    lines = []
    lines.append("// AUTO-GENERATED by generate_regions.py -- do not hand-edit.")
    lines.append("// Region tags, nonlinear iron B-H curve, and winding current-density")
    lines.append("// sources for the real BPL-700 assembly. See generate_regions.py for")
    lines.append("// derivation.")
    lines.append("")
    lines.append("Group {")
    for name, t in sorted(tag.items(), key=lambda kv: kv[1]):
        lines.append(f"  {name} = Region[{t}];")
    lines.append(f"  Windings = Region[{{{', '.join(winding_names)}}}];")
    lines.append("  DomainC  = Region[{Windings}];")
    lines.append("  DomainCC = Region[{Air, Iron}];")
    lines.append("  Domain   = Region[{DomainC, DomainCC}];")
    lines.append("}")
    lines.append("")
    lines.append("Function {")
    lines.append(f"  mu0 = {MU0!r};")
    lines.append("")
    lines.append("  // Iron B-H curve. The hardware is ASTM A36; these are AISI 1008")
    lines.append("  // points used as a PROXY -- see generate_regions.py for why, which")
    lines.append("  // way the error runs, and its rough size.")
    lines.append(f"  AISI1008_H() = {_gmsh_list(AISI1008_H)};")
    lines.append(f"  AISI1008_B() = {_gmsh_list(AISI1008_B)};")
    lines.append("  AISI1008_B2() = AISI1008_B()^2;")
    lines.append("  AISI1008_nu_list() = AISI1008_H() / AISI1008_B();")
    lines.append("  AISI1008_nu_list(0) = AISI1008_nu_list(1);  // avoid 0/0 at B=0")
    lines.append("  AISI1008_nu_b2_list() = ListAlt[AISI1008_B2(), AISI1008_nu_list()];")
    lines.append("  AISI1008_nu[] = InterpolationLinear[SquNorm[$1]]{AISI1008_nu_b2_list()};")
    lines.append("")
    lines.append("  nu[Air]  = 1 / mu0;")
    lines.append("  nu[Iron] = AISI1008_nu[$1];  // nonlinear -- see magnetostatics_assembly.pro")
    lines.append("                                   // for the Vol_L_Mag/Vol_NL_Mag Galerkin split this requires")
    lines.append("  nu[Windings] = 1 / mu0;  // copper, non-magnetic")
    lines.append("")
    for w in params["windings"]:
        i, name, A_cross, loc = w["index"], w["pole"], w["A_cross"], w["loc"]
        exc = excitation[name]
        Jmag = exc["polarity"] * exc["turns"] * exc["current"] / A_cross
        lx, lz = loc
        lines.append(f"  // Winding{i}: pole={name} turns={exc['turns']} I={exc['current']}A "
                      f"A_cross={A_cross:.6e} polarity={exc['polarity']:+d} -> Jmag={Jmag:.6e} A/m^2")
        # $IFrac: excitation load-stepping fraction, ramped 0->1 by
        # magnetostatics_assembly.pro's Resolution to keep the nonlinear
        # iron iteration stable -- see that file's header comment.
        lines.append(f"  Js[Winding{i}] = $IFrac * ({Jmag!r} / Sqrt[(X[]-({lx!r}))^2 + (Z[]-({lz!r}))^2]) * "
                      f"Vector[(Z[]-({lz!r})), 0, -(X[]-({lx!r}))];")
    lines.append("}")
    lines.append("")

    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    params = parse_params()
    generate_regions_pro(params, DEFAULT_EXCITATION)
    desc = " vs ".join(
        f"{name} {e['turns']}t/{e['current']}A" for name, e in DEFAULT_EXCITATION.items()
    )
    print(f"Wrote assembly_regions_generated.pro ({desc})")
