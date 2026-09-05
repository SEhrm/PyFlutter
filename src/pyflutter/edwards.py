# Copyright (c) 2026 Simon Ehrmanntraut
# Licensed under CC BY-NC 4.0, see LICENSE for details.
"""Edwards' unsteady aerodynamics."""

import numpy as np
from scipy.special import kve


def generalized_theodorsen(reduced_laplace: complex | np.ndarray) -> complex | np.ndarray:
    """Edwards' generalized Theodorsen's function.

    Args:
        reduced_laplace: Reduced laplace `sb/v`.

    Returns:
        Evaluation of the generalized Theodorsen's function.
    """
    # Use the exponentially scaled Bessel function to avoid overflow; the exponents cancel out.
    k0e = kve(0, reduced_laplace + 1e-12j)
    k1e = kve(1, reduced_laplace + 1e-12j)
    return k1e / (k1e + k0e)


def edwards_transfer(
    reduced_laplace: complex | np.ndarray,
    axis: float,
    aileron: float,
) -> np.ndarray:
    """Edwards' MIMO transfer function matrix.

    The 3x3 transfer function matrix maps the structural deflections `[h/b, α, β]` to the
    aerodynamic forces `[Pb, M_α, M_β] / (ϱb²v²)`.

    Args:
        reduced_laplace: Reduced laplace `sb/v`.
        axis: Downstream distance from mid-chord to pitch axis per semi-chord.
        aileron: Downstream distance from mid-chord to aileron hinge per semi-chord.

    Returns:
        Evaluation of Edwards's transfer function.
    """
    # Theodorsen T-expressions for the flap
    acos = np.arccos(aileron)
    sin = np.sin(acos)
    t1 = -1 / 3 * sin * (2 + aileron**2) + aileron * acos
    t3 = -(.125 + aileron**2) * acos**2 + .25 * aileron * sin * acos * (7 + 2 * aileron**2) \
         - .125 * sin**2 * (5 * aileron**2 + 4)
    t4 = -acos + aileron * sin
    t5 = -sin**2 - acos**2 + 2 * aileron * sin * acos
    t7 = -(.125 + aileron**2) * acos + .125 * aileron * sin * (7 + 2 * aileron**2)
    t8 = -1 / 3 * sin * (2 * aileron**2 + 1) + aileron * acos
    t9 = .5 * (1 / 3 * sin**3 + axis * t4)
    t10 = sin + acos
    t11 = acos * (1 - 2 * aileron) + sin * (2 - aileron)
    t12 = sin * (2 + aileron) - acos * (2 * aileron + 1)
    t13 = .5 * (-t7 - (aileron - axis) * t1)

    # Edwards' transfer function matrix
    edwards_m = np.array([
        [-np.pi, np.pi * axis, t1],
        [np.pi * axis, -np.pi * (.125 + axis**2), (t7 + (aileron - axis) * t1)],
        [t1, -2. * t13, t3 / np.pi],
    ])
    edwards_b = np.array([
        [0., -np.pi, t4],
        [0., -np.pi * (.5 - axis), -(t1 - t8 - (aileron - axis) * t4 + t11 / 2)],
        [0., -(-2 * t9 - t1 + t4 * (axis - .5)), t4 * t11 / (2 * np.pi)],
    ])
    edwards_k = np.array([
        [0., 0., 0.],
        [0., 0., -(t4 + t10)],
        [0., 0., -(t5 - t4 * t10) / np.pi],
    ])
    edwards_r = np.array([[-2. * np.pi, 2. * np.pi * (axis + .5), -t12]]).T
    edwards_s1 = np.array([[0., 1., t10 / np.pi]])
    edwards_s2 = np.array([[1., (.5 - axis), t11 / 2. / np.pi]])
    added_mass = (
        np.multiply.outer(reduced_laplace**2, edwards_m) +
        np.multiply.outer(reduced_laplace**1, edwards_b) +
        np.multiply.outer(reduced_laplace**0, edwards_k)
    )
    quasi_steady = edwards_r @ (
        np.multiply.outer(reduced_laplace**1, edwards_s2) +
        np.multiply.outer(reduced_laplace**0, edwards_s1)
    )
    return added_mass + generalized_theodorsen(reduced_laplace)[..., None, None] * quasi_steady
