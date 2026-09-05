# Copyright (c) 2026 Simon Ehrmanntraut
# Licensed under CC BY-NC 4.0, see LICENSE for details.
"""PyFlutter eigenvalue computation."""
from collections.abc import Callable

import numpy as np
from scipy.linalg import eigvals
from scipy.sparse import block_array

from .edwards import edwards_transfer


def solve_qep(
    mass_mat: np.ndarray,
    damping_mat: np.ndarray,
    stiffness_mat: np.ndarray,
) -> np.ndarray:
    """Solves the quadratic eigenvalue problem.

    Args:
        mass_mat: Structural mass matrix.
        damping_mat: Structural damping matrix.
        stiffness_mat: Structural stiffness matrix.

    Returns:
        Eigenvalues with positive imaginary part.
    """
    id_mat = np.identity(len(mass_mat))
    return eigvals(
        b=block_array([[id_mat, None], [None, mass_mat]]).toarray(),
        a=block_array([[None, id_mat], [-stiffness_mat, -damping_mat]]).toarray(),
    )


def compute_flutter_eigenvalue(
    mass_mat: np.ndarray,
    damping_mat: np.ndarray,
    stiffness_mat: np.ndarray,
    transfer: Callable[[complex], np.ndarray],
    mass_ratio: float,
    speed: float,
    guess: complex,
    max_iter: int = 20,
    atol: float = 1e-7,
    rtol: float = 1e-4,
) -> complex:
    """Solves the flutter eigenvalue problem for a given speed and initial guess.

    Iterates the nonlinear eigenproblem by successive quadratic eigenvalue problems.

    Args:
        mass_mat: Structural mass matrix.
        damping_mat: Structural damping matrix.
        stiffness_mat: Structural stiffness matrix.
        transfer: Edwards' transfer function.
        mass_ratio: Mass ratio `μ=M/(πϱb²)`.
        speed: Reduced velocity `v/(bωᵣ)`.
        guess: Eigenvalue guess `s/ωᵣ` per reference frequency.
        max_iter: Maximum number of iterations.
        atol: Absolute tolerance for convergence.
        rtol: Relative tolerance for convergence.

    Returns:
        Eigenvalue per reference frequency `s/ωᵣ`; Nan if not converged.
    """
    for _ in range(max_iter):
        aero_mat = transfer(guess / speed) * speed**2 / np.pi / mass_ratio
        qep_eigenvalues = solve_qep(
            mass_mat=mass_mat, damping_mat=damping_mat,
            stiffness_mat=stiffness_mat - aero_mat,
        )
        eigenvalue = min(qep_eigenvalues, key=lambda v: abs(v - guess))
        if eigenvalue.imag < 0:
            return np.nan
        if np.allclose(eigenvalue, guess, rtol, atol):
            return eigenvalue
        guess = eigenvalue
    return np.nan


def compute_locus(
    axis: float,
    radius_gyration: float,
    unbalance: float,
    aileron: float,
    aileron_radius_gyration: float,
    aileron_unbalance: float,
    pitching_frequency: float,
    plunging_frequency: float,
    flapping_frequency: float,
    damping_ratio: float,
    mass_ratio: float,
    speed_step: float,
    speed_max: float,
    damping_ratio_min: float,
    damping_ratio_max: float,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Compute the root-locus of the flutter modes.

    Args:
        axis: Downstream distance from mid-chord to elastic axis per semi-chord.
        radius_gyration: Airfoil radius of gyration per semi-chord.
        unbalance: Downstream distance from elastic axis to airfoil center of gravity per semi-chord.
        aileron: Downstream distance from mid-chord to aileron hinge per semi-chord.
        aileron_radius_gyration: Aileron radius of gyration per semi-chord.
        aileron_unbalance: Reduced downstream distance from aileron hinge to aileron center of
            gravity per semi-chord.
        pitching_frequency: Natural pitching frequency per reference frequency (NaN to suppresses).
        plunging_frequency: Natural plunging frequency per reference frequency (NaN to suppresses).
        flapping_frequency: Natural flapping frequency per reference frequency (NaN to suppresses).
        damping_ratio: Structural viscous damping ratio.
        mass_ratio: Mass ratio `μ=M/(πϱb²)`.
        speed_step: Reduced velocity step `Δv/(bωᵣ)` for mode tracking.
        speed_max: Stopping speed for mode tracking.
        damping_ratio_min: Stopping damping ratio for mode tracking.
        damping_ratio_max: Stopping damping ratio for mode tracking.

    Returns:
        The mode tuples, each containing the reduced velocities `v/(bωᵣ)` and the eigenvalues
            per reference frequency `s/ωᵣ`.
    """
    # Structural matrices
    sep = aileron - axis
    mass_mat = np.zeros((3, 3))
    mass_mat[0, 0] = 1.
    mass_mat[1, 1] = radius_gyration**2
    mass_mat[2, 2] = aileron_radius_gyration**2
    mass_mat[0, 1] = mass_mat[1, 0] = unbalance
    mass_mat[0, 2] = mass_mat[2, 0] = aileron_unbalance
    mass_mat[1, 2] = mass_mat[2, 1] = aileron_radius_gyration**2 + aileron_unbalance * sep
    stiffness_mat = np.zeros((3, 3))
    stiffness_mat[0, 0] = plunging_frequency**2
    stiffness_mat[1, 1] = pitching_frequency**2 * radius_gyration**2
    stiffness_mat[2, 2] = flapping_frequency**2 * aileron_radius_gyration**2
    damping_mat = np.zeros((3, 3))
    damping_mat[0, 0] = 2 * damping_ratio * plunging_frequency
    damping_mat[1, 1] = 2 * damping_ratio * pitching_frequency * radius_gyration**2
    damping_mat[2, 2] = 2 * damping_ratio * flapping_frequency * aileron_radius_gyration**2

    # Root-locus
    dof_mask = np.argwhere(np.isfinite(np.diag(stiffness_mat)))
    dof_mask = (dof_mask, dof_mask.T)
    vacuum_eigenvalues = solve_qep(
        mass_mat=mass_mat[dof_mask],
        damping_mat=damping_mat[dof_mask],
        stiffness_mat=stiffness_mat[dof_mask],
    )
    modes: list[tuple[np.ndarray, np.ndarray]] = []
    for structural_eigenvalue in vacuum_eigenvalues[vacuum_eigenvalues.imag > 0]:
        guess = structural_eigenvalue
        eigenvalues = []
        speeds = []
        for speed in np.arange(1, np.rint(speed_max / speed_step) + 1) * speed_step:
            eigenvalue = compute_flutter_eigenvalue(
                mass_mat=mass_mat[dof_mask],
                damping_mat=damping_mat[dof_mask],
                stiffness_mat=stiffness_mat[dof_mask],
                transfer=lambda v: edwards_transfer(v, axis=axis, aileron=aileron)[dof_mask],
                mass_ratio=mass_ratio, speed=speed, guess=guess,
            )
            print(f"{speed:06.3f} {eigenvalue:+.4f}")
            if np.isnan(eigenvalue):
                print("Stopped due to non-convergence; decrease speed step")
                break
            if np.isclose(eigenvalue.imag, 0.):
                print("Stopped due to non-oscillatory eigenvalue")
                break
            speeds.append(speed)
            eigenvalues.append(eigenvalue)
            if -eigenvalue.real / abs(eigenvalue) < damping_ratio_min:
                print("Stopped due to flutter")
                break
            if -eigenvalue.real / abs(eigenvalue) > damping_ratio_max:
                print("Stopped due to damping ratio")
                break
            guess = eigenvalue
        else:
            print("Stopped due to step limit")
        modes.append((np.array(speeds), np.array(eigenvalues)))
    return modes
