# Copyright (c) 2026 Simon Ehrmanntraut
# Licensed under CC BY-NC 4.0, see LICENSE for details.
"""Command-line entry for PyFlutter."""

import argparse

from . import compute_locus
from .plotting import plot_eigenvalues

try:
    from matplotlib import pyplot as plt
except ImportError:
    plt = None


class _HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def _get_default_metavar_for_optional(self, action: argparse.Action) -> str:  # noqa: PLR6301
        return action.option_strings[0].lstrip("-").upper()


def run() -> None:
    """Runs the PyFlutter command-line interface."""
    parser = argparse.ArgumentParser(
        prog="pyflutter",
        description="Compute the root-locus for the 3-DOF flat plate airfoil",
        formatter_class=_HelpFormatter,
    )
    parser.add_argument(
        "-a", dest="axis", default=-0.4, type=float,
        help="distance from mid-chord to axis per semi-chord",
    )
    parser.add_argument(
        "-r", dest="radius_gyration", default=0.5, type=float,
        help="airfoil radius of gyration per semi-chord",
    )
    parser.add_argument(
        "-x", dest="unbalance", default=0.2, type=float,
        help="distance from pitch axis to airfoil center of gravity per semi-chord",
    )
    parser.add_argument(
        "-c", dest="aileron", default=0.6, type=float,
        help="distance from mid-chord to aileron hinge per semi-chord",
    )
    parser.add_argument(
        "-rb", dest="aileron_radius_gyration", default=0.079, type=float,
        help="aileron radius of gyration per semi-chord",
    )
    parser.add_argument(
        "-xb", dest="aileron_unbalance", default=0.0125, type=float,
        help="reduced distance from aileron hinge to aileron center of gravity per semi-chord",
    )
    parser.add_argument(
        "-wa", dest="pitching", default=1.0, type=float,
        help="natural pitching frequency per reference frequency (NaN to suppress DOF)",
    )
    parser.add_argument(
        "-wh", dest="plunging", default=0.5, type=float,
        help="natural plunging frequency per reference frequency (NaN to suppress DOF)",
    )
    parser.add_argument(
        "-wb", dest="flapping", default=3.0, type=float,
        help="natural flapping frequency per reference frequency (NaN to suppress DOF)",
    )
    parser.add_argument(
        "-m", dest="mass_ratio", default=40., type=float,
        help="airfoil mass per mass of a cylinder of air of diameter equal to the chord",
    )
    parser.add_argument(
        "-d", dest="damping", default=0., type=float,
        help="structural viscous damping ratio",
    )
    parser.add_argument(
        "--v-step", dest="speed_step", default=0.1, type=float,
        help="reduced velocity step for mode tracking",
    )
    parser.add_argument(
        "--v-max", dest="speed_max", default=3.5, type=float,
        help="maximum reduced velocity for mode tracking",
    )
    parser.add_argument(
        "--d-max", dest="damping_max", default=+1., type=float,
        help="maximum damping ratio for mode tracking",
    )
    parser.add_argument(
        "--d-min", dest="damping_min", default=0., type=float,
        help="minimum damping ratio for mode tracking",
    )
    parser.add_argument("--plot", action="store_true", help="launch plots")
    args = parser.parse_args()
    modes = compute_locus(
        axis=args.axis,
        radius_gyration=args.radius_gyration,
        unbalance=args.unbalance,
        aileron=args.aileron,
        aileron_radius_gyration=args.aileron_radius_gyration,
        aileron_unbalance=args.aileron_unbalance,
        plunging_frequency=args.plunging,
        pitching_frequency=args.pitching,
        flapping_frequency=args.flapping,
        damping_ratio=args.damping,
        mass_ratio=args.mass_ratio,
        speed_step=args.speed_step,
        speed_max=args.speed_max,
        damping_ratio_min=args.damping_min,
        damping_ratio_max=args.damping_max,
    )
    if args.plot:
        plot_eigenvalues(modes)
        plt.show()


if __name__ == "__main__":
    run()
