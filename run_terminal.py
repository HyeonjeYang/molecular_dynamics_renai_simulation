#!/usr/bin/env python3
"""Terminal runner for AffinityMD (renai_simu).

Example:
    python run_terminal.py --N_male 10 --N_female 10 --T 200 --seed 7
    python run_terminal.py --help
"""
import argparse
import sys
from datetime import datetime
from pathlib import Path

from renai_simu import Simulation, DEFAULT_PARAMS, visualize, io_utils


def build_parser():
    p = argparse.ArgumentParser(
        description='AffinityMD — stochastic romance simulation (terminal)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Population
    p.add_argument('--N_male',   type=int,   default=DEFAULT_PARAMS['N_male'])
    p.add_argument('--N_female', type=int,   default=DEFAULT_PARAMS['N_female'])
    p.add_argument('--age_min',  type=float, default=DEFAULT_PARAMS['age_min'])
    p.add_argument('--age_max',  type=float, default=DEFAULT_PARAMS['age_max'])

    # Box and time
    p.add_argument('--box',      type=float, default=DEFAULT_PARAMS['box_size'],
                   dest='box_size')
    p.add_argument('--T',        type=float, default=DEFAULT_PARAMS['T_total'],
                   dest='T_total', help='Total simulation time')
    p.add_argument('--dt',       type=float, default=DEFAULT_PARAMS['dt'])

    # Langevin
    p.add_argument('--gamma',    type=float, default=DEFAULT_PARAMS['gamma'])
    p.add_argument('--kT',       type=float, default=DEFAULT_PARAMS['kT'])

    # Morse
    p.add_argument('--eps0',     type=float, default=DEFAULT_PARAMS['epsilon_0'],
                   dest='epsilon_0')
    p.add_argument('--a_morse',  type=float, default=DEFAULT_PARAMS['a_morse'])
    p.add_argument('--r0_morse', type=float, default=DEFAULT_PARAMS['r0_morse'])

    # Affinity
    p.add_argument('--chi_het',  type=float, default=DEFAULT_PARAMS['chi_het'])
    p.add_argument('--chi_hom',  type=float, default=DEFAULT_PARAMS['chi_hom'])
    p.add_argument('--sigma_age', type=float, default=DEFAULT_PARAMS['sigma_age'])

    # Relationship dynamics
    p.add_argument('--k_on',     type=float, default=DEFAULT_PARAMS['k_on'])
    p.add_argument('--k_off',    type=float, default=DEFAULT_PARAMS['k_off'])
    p.add_argument('--r_c',      type=float, default=DEFAULT_PARAMS['r_c'])
    p.add_argument('--B_thresh', type=float, default=DEFAULT_PARAMS['B_thresh'])
    p.add_argument('--B_break',  type=float, default=DEFAULT_PARAMS['B_break'])

    # Perturbation
    p.add_argument('--Q_perturb', type=float, default=DEFAULT_PARAMS['Q_perturb'])
    p.add_argument('--eps_screen', type=float, default=DEFAULT_PARAMS['eps_screen'])

    # Stochastic events
    p.add_argument('--p_event',  type=float, default=DEFAULT_PARAMS['p_event'])
    p.add_argument('--B_boost',  type=float, default=DEFAULT_PARAMS['B_boost'])
    p.add_argument('--B_damage', type=float, default=DEFAULT_PARAMS['B_damage'])

    # Reproducibility / output
    p.add_argument('--seed',     type=int,   default=DEFAULT_PARAMS['seed'])
    p.add_argument('--record_interval', type=int,
                   default=DEFAULT_PARAMS['record_interval'])
    p.add_argument('--out', type=str, default=None,
                   help=('Output directory. '
                         'Default: outputs/terminal_run_<timestamp>/'))
    p.add_argument('--quiet', action='store_true',
                   help='Suppress progress output.')
    return p


def args_to_params(args):
    """Convert argparse namespace into a params dict (only known keys)."""
    params = dict(DEFAULT_PARAMS)
    for key in list(params.keys()):
        if hasattr(args, key):
            params[key] = getattr(args, key)
    return params


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    params = args_to_params(args)

    # Unique output dir for this run (never overwrites notebook output)
    if args.out is None:
        stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        out_dir = Path('outputs') / f'terminal_run_{stamp}'
    else:
        out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Output directory: {out_dir}")

    # Run
    sim = Simulation(params)
    sim.run(verbose=not args.quiet)

    # Save text outputs and figures
    io_utils.save_all(sim, out_dir)
    visualize.plot_all(sim, out_dir)

    # Console summary
    n_c = sum(1 for i in sim.males if sim.couple[i] != -1)
    max_c = min(sim.N_m, sim.N_f)
    print()
    print(f"Final couples : {n_c}/{max_c}")
    print(f"Total events  : {len(sim.event_log)}")
    for i in sim.males:
        j = int(sim.couple[i])
        if j != -1:
            print(f"  {sim.labels[i]} <3 {sim.labels[j]}  "
                  f"(B={sim.B[i, j]:.3f})")
    print(f"\nAll outputs written to: {out_dir}")


if __name__ == '__main__':
    main()
