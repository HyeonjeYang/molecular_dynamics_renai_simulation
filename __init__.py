"""
renai_simu — AffinityMD
=======================
Stochastic particle model of romantic pair formation via overdamped Langevin dynamics.

Quick start:
    from renai_simu import Simulation, DEFAULT_PARAMS, visualize, io_utils

    sim = Simulation({'N_male': 6, 'N_female': 6, 'T_total': 100})
    sim.run()
    visualize.plot_all(sim, out_dir='outputs/my_run')
    io_utils.save_all(sim,  out_dir='outputs/my_run')
"""
from .params import DEFAULT_PARAMS
from .simulation import Simulation
from . import dynamics
from . import relationships
from . import events
from . import visualize
from . import io_utils

__all__ = [
    'Simulation',
    'DEFAULT_PARAMS',
    'dynamics',
    'relationships',
    'events',
    'visualize',
    'io_utils',
]
