# Copyright (c) 2026 Simon Ehrmanntraut
# Licensed under CC BY-NC 4.0, see LICENSE for details.
"""Plotting functionalities for PyFlutter eigenvalue results."""

import numpy as np

try:
    from matplotlib import pyplot as plt
except ImportError:
    plt = None


def plot_eigenvalues(
    modes: list[tuple[np.ndarray, np.ndarray]],
) -> tuple[plt.Figure, plt.Figure]:
    """Plots the flutter diagram and root-locus diagram of the two flutter modes.

    Args:
        modes: Mode tuples, each containing the reduced velocities `v/(bωᵣ)` and the eigenvalues
            per reference frequency `s/ωᵣ`.

    Returns:
        The flutter diagram and root-locus diagram figures.

    Raises:
        RuntimeError: If matplotlib is not installed.
    """
    if plt is None:
        msg = "matplotlib is not installed"
        raise RuntimeError(msg)
    # Flutter diagram
    fig_flutter = plt.figure()
    gs = fig_flutter.add_gridspec(2, hspace=0.1)
    ax_frequency, ax_damping = gs.subplots(sharex=True)
    for speeds, eigenvalues in modes:
        ax_frequency.plot(speeds, eigenvalues.imag, ".-")
    ax_frequency.set_ylabel(r"$\mathrm{Im}(s)/{\omega_\text{r}}$")
    ax_frequency.yaxis.set_major_formatter("{x:.2f}")
    ax_frequency.margins(x=0)
    ax_frequency.set_ylim(0, None)
    ax_frequency.set_xlim(0, None)
    ax_frequency.grid()
    for speeds, eigenvalues in modes:
        ax_damping.plot(speeds, -eigenvalues.real / np.abs(eigenvalues), ".-")
    ax_damping.axhline(y=0, color="k", ls="-", lw=0.5)
    ax_damping.set_ylabel(r"$-\mathrm{Re}(s)/|s|$")
    ax_damping.set_xlabel(r"$v_\infty b/{\omega_\text{r}}$")
    ax_damping.yaxis.set_major_formatter("{x:+.2f}")
    ax_damping.margins(x=0)
    ax_damping.grid()
    # Root-locus diagram
    fig_root_locus = plt.figure()
    ax = fig_root_locus.gca()
    for _, eigenvalues in modes:
        ax.plot(eigenvalues.real, eigenvalues.imag, ".-")
    ax.axvline(x=0, color="k", ls="-", lw=0.5)
    ax.set_xlabel(r"$\mathrm{Re}(s)/{\omega_\text{r}}$")
    ax.set_ylabel(r"$\mathrm{Im}(s)/{\omega_\text{r}}$")
    ax.xaxis.set_major_formatter("{x:.2f}")
    ax.yaxis.set_major_formatter("{x:.2f}")
    ax.margins(x=0)
    ax.set_xlim(None, None)
    ax.set_ylim(0, None)
    ax.grid()
    return fig_flutter, fig_root_locus
