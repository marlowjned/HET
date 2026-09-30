"""
wire_data.py

Magnet wire dimensions and insulation properties, from NEMA MW 1000 rather
than from estimates.

Why this module exists: the winding design's turn counts are set by the
INSULATED diameter of the wire, and until 2026-09-30 that number was an
engineering estimate for a product ("single-glass-served, 200 C, AWG 20")
that turned out not to be purchasable. MATERIALS.md flagged those figures
as the softest input in the project; this replaces them with a standard's
tabulated values and records what is actually available.

Everything here is the MAXIMUM overall diameter, not the nominal. A fit
calculation wants the worst case the supplier is allowed to ship, and the
difference is not negligible.

Source for both tables: NEMA Standards Publication MW 1000, as reproduced
in the magnet wire / winding wire engineering data handbook (12/2009 Rev.1),
https://www.casa.co.nz/cables/wire/magnet/EngData-Book_66p.pdf
  - HEAVY_BUILD   -- "ROUND, HEAVY BUILD FILM-INSULATED WIRE", printed p.18
  - SERVED_GLASS  -- "ROUND, SERVED WIRE, POLYESTER GLASS & POLYIMIDE
                      TAPE", printed p.36
Transcribed by parsing the PDF text layer, not by hand, then range-checked
(see _self_test at the bottom).
"""
from __future__ import annotations

# ---------------------------------------------------------------------------
# Round heavy-build film-insulated wire (NEMA MW 1000, printed p.18).
# AWG -> (bare nominal, max overall diameter), mm.
#
# "Heavy build" is the standard two-layer film thickness and is what every
# stocked high-temperature magnet wire uses, whatever the polymer: MW 16-C
# (polyimide, 240 C), MW 35-C (polyurethane/nylon, 200 C) and MW 80-C
# (polyamideimide, 200 C) all share this dimension table. So these numbers
# are insulation-chemistry-independent -- only the temperature class and
# the price change with the polymer.
# ---------------------------------------------------------------------------
HEAVY_BUILD = {
    12: (2.052, 2.163), 13: (1.829, 1.934), 14: (1.628, 1.732),
    15: (1.450, 1.549), 16: (1.290, 1.384), 17: (1.151, 1.240),
    18: (1.024, 1.110), 19: (0.912, 0.993), 20: (0.813, 0.892),
    21: (0.724, 0.800), 22: (0.643, 0.714), 23: (0.574, 0.643),
    24: (0.511, 0.577), 25: (0.455, 0.516), 26: (0.404, 0.462),
    27: (0.361, 0.419), 28: (0.320, 0.373), 29: (0.287, 0.338),
    30: (0.254, 0.307),
}

# ---------------------------------------------------------------------------
# Round polyester-glass ("Dacron glass") SERVED wire (NEMA MW 1000, p.36).
# AWG -> (bare nominal, max OD single serving, max OD double serving), mm,
# all over BARE wire (the table's "over heavy film" columns are omitted).
#
# NOTE THE RANGE. The standard tabulates round served glass wire from 1/0
# down to AWG 12 and NO FINER. This is the fact that invalidated the old
# AWG 20 fiberglass design point: there is no such standard product.
# Glass serving is built from filaments whose thickness does not scale with
# the conductor, so the covering adds ~0.15-0.21 mm to the OD regardless of
# gauge -- which is why it stops being offered on fine wire.
#
# The AWG 12 double-serving figure is printed as 0.0986 in in the source
# table. That is a digit transposition: every other row's double serving
# adds 0.226-0.310 mm over bare, whereas 0.0986 in would add 0.452 mm.
# The intended value is 0.0896 in (0.0088 in over bare, consistent with
# AWG 11's 0.226 mm and AWG 10's 0.229 mm). Corrected here, and noted
# rather than silently fixed.
# ---------------------------------------------------------------------------
SERVED_GLASS = {
    "1/0": (8.252, 8.461, 8.562), 1: (7.348, 7.549, 7.650),
    2: (6.543, 6.736, 6.838), 3: (5.827, 6.012, 6.114),
    4: (5.189, 5.392, 5.469), 5: (4.620, 4.818, 4.895),
    6: (4.115, 4.308, 4.384), 7: (3.665, 3.853, 3.929),
    8: (3.264, 3.449, 3.526), 9: (2.906, 3.086, 3.162),
    10: (2.588, 2.741, 2.817), 11: (2.304, 2.454, 2.530),
    12: (2.052, 2.200, 2.276),   # 2.276 = corrected 0.0896 in, see note above
}


class Insulation:
    """One wire product.

    `availability` is the field that matters most here and the one a
    datasheet never states: whether a spool of it can be bought.
    """

    def __init__(self, key, label, nema, temp_C, table, od_index, availability,
                 notes=""):
        self.key, self.label, self.nema = key, label, nema
        self.temp_C, self._table, self._od_index = temp_C, table, od_index
        self.availability, self.notes = availability, notes

    @property
    def gauges(self):
        return sorted(g for g in self._table if isinstance(g, int))

    def bare_mm(self, awg):
        return self._table[awg][0]

    def od_mm(self, awg):
        """Maximum overall diameter, mm."""
        return self._table[awg][self._od_index]

    def __repr__(self):
        return f"<Insulation {self.key} {self.temp_C:.0f}C {self.availability}>"


# The two options carried through the design. POLYIMIDE is the specified
# wire; the glass entries are a spec-only fallback -- see MATERIALS.md
# "Magnet wire" for the availability findings behind that split.
POLYIMIDE = Insulation(
    "polyimide", "aromatic polyimide film, heavy build", "MW 16-C", 240.0,
    HEAVY_BUILD, 1, "stocked",
    "Off-the-shelf in 2 oz to 10 lb spools (Remington Industries, AWG 12-40). "
    "A solid film, so there is no porous covering to trap gas -- which "
    "matters in vacuum.",
)

GLASS_SINGLE = Insulation(
    "glass", "single polyester-glass serving", "MW 45", 155.0,
    SERVED_GLASS, 1, "made-to-order",
    "Class 155, NOT 200 C: the 200 C variant (Dacron Glass High Temp, MW 51) "
    "is made 'only in shaped or heavy round wire sizes'. No US retailer "
    "stocks round glass-served magnet wire in any gauge; S&W Wire and "
    "Bridgeport Insulated Wire make it to order. Finest standard round "
    "size is AWG 12.",
)

GLASS_DOUBLE = Insulation(
    "glass2", "double polyester-glass serving", "MW 45", 155.0,
    SERVED_GLASS, 2, "made-to-order", GLASS_SINGLE.notes,
)

OPTIONS = {o.key: o for o in (POLYIMIDE, GLASS_SINGLE, GLASS_DOUBLE)}

# Ground insulation between the bare steel core and the first turn: two
# half-lapped wraps of 0.05 mm polyimide tape, taken as a radial allowance.
# Not a datasheet value -- an allowance. It is broken out rather than buried
# because it costs ~3% of the inner coil's 4.76 mm build, the same order as
# the OD tolerance it sits next to.
GROUND_INSULATION_MM = 0.15


def _self_test():
    # Bare diameters must agree with the values the repo already used, which
    # came from a standard AWG table (MATERIALS.md, "Copper (windings)").
    for awg, expect in ((16, 1.291), (18, 1.024), (20, 0.812), (22, 0.644), (24, 0.511)):
        got = HEAVY_BUILD[awg][0]
        assert abs(got - expect) < 0.002, (awg, got, expect)

    # Both tables monotonic in gauge, and every covering larger than bare.
    for tbl, name in ((HEAVY_BUILD, "HEAVY_BUILD"), (SERVED_GLASS, "SERVED_GLASS")):
        ints = sorted(g for g in tbl if isinstance(g, int))
        for a in ints:
            row = tbl[a]
            assert all(row[i] > row[0] for i in range(1, len(row))), (name, a)
        for a, b in zip(ints, ints[1:]):
            assert tbl[b][0] < tbl[a][0], (name, a, b)

    # Glass serving adds a roughly gauge-independent amount, because it is
    # filaments rather than a coating. This is the check that exposed the
    # AWG 12 transposition, so keep it tight enough to catch a repeat.
    adds = [SERVED_GLASS[a][2] - SERVED_GLASS[a][0] for a in (10, 11, 12)]
    assert max(adds) - min(adds) < 0.01, adds

    # A film coating, by contrast, scales with the conductor.
    assert HEAVY_BUILD[12][1] - HEAVY_BUILD[12][0] > 2 * (HEAVY_BUILD[30][1] - HEAVY_BUILD[30][0])

    # The finest standard round served-glass size is AWG 12. If this ever
    # fails, the availability argument in MATERIALS.md needs revisiting.
    assert max(g for g in SERVED_GLASS if isinstance(g, int)) == 12


_self_test()

if __name__ == "__main__":
    print(f"{'AWG':>4}  {'bare':>6}  {'polyimide':>9} {'od/bare':>7}   {'glass 1x':>8} {'glass 2x':>8}")
    for awg in sorted(HEAVY_BUILD):
        bare, od = HEAVY_BUILD[awg]
        g = SERVED_GLASS.get(awg)
        gs = f"{g[1]:>8.3f} {g[2]:>8.3f}" if g else f"{'--':>8} {'--':>8}"
        print(f"{awg:>4}  {bare:>6.3f}  {od:>9.3f} {od / bare:>7.3f}   {gs}")
    print()
    for o in OPTIONS.values():
        print(f"{o.key:>8}  {o.nema:<8} {o.temp_C:>5.0f} C  {o.availability:<14} "
              f"AWG {min(o.gauges)}-{max(o.gauges)}")
