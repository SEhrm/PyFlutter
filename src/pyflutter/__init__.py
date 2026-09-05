# Copyright (c) 2026 Simon Ehrmanntraut
# Licensed under CC BY-NC 4.0, see LICENSE for details.
"""Flutter analysis of the 3-DOF flat plate airfoil with Edwards' unsteady aerodynamics.

Performs flutter analysis of a flat plate airfoil undergoing plunging, pitching and trailing edge
flapping motion, constrained by springs and dampers. Unsteady aerodynamics are modeled by
Edwards' generalized Theodorsen's function "EDWARDS, John William. Unsteady aerodynamic modeling
and active aeroelastic control. 1977. (Thesis). https://ntrs.nasa.gov/citations/19780002074".
"""

from .core import compute_locus
from .edwards import edwards_transfer
from .edwards import generalized_theodorsen

__all__ = [
    "compute_locus",
    "edwards_transfer",
    "generalized_theodorsen",
]
