"""Shared winding logic for the coil part types (helper: not an entry)."""
from math import ceil, pi

from bplsim.parts import Length, Material, Param, PartType


class Winding(PartType):
    """Base class for coils: wire, turns, layers, resistance.

    The winding window is either computed from layered turns (per_layer turns per layer of insulated wire),
    or given explicitly with `build` and `length` (the legacy model's smeared spec windows)."""

    class Params:
        wire: Material = Param(needs=["resistivity", "wire_diameter", "insulated_diameter"])
        turns: int = Param(min=1)                      # per coil
        per_layer: int = Param(46, min=1)
        circuit: str = Param()                         # name of an op circuit
        build: Length | None = Param(None)             # explicit window; None = from layered turns
        length: Length | None = Param(None)

    class Derived:
        layers: int
        length: Length
        build: Length
        coil_count: int

    def winding(self, p):
        d = p.wire.insulated_diameter
        layers = ceil(p.turns / p.per_layer)
        return {"layers": layers, "length": p.length if p.length else p.per_layer * d,
                "build": p.build if p.build else layers * d}

    def coil_count(self, p):
        return 1

    def core_radius(self, p):
        raise NotImplementedError

    def mean_turn_length(self, p):
        return 2 * pi * (self.core_radius(p) + self.winding(p)["build"] / 2)

    def derived(self, p):
        return {**self.winding(p), "coil_count": self.coil_count(p)}

    def resistance(self, p, temperature):
        """Resistance of one coil at `temperature` (K)."""
        area = pi * p.wire.wire_diameter ** 2 / 4
        return p.wire.resistivity.at(temperature) * p.turns * self.mean_turn_length(p) / area
