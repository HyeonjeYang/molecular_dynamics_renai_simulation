"""Default parameters for AffinityMD simulation.

Users may copy this dict and modify values, or pass individual overrides
to Simulation(params_override=...).
"""

DEFAULT_PARAMS = {
    # --- Population ---
    'N_male':      8,
    'N_female':    8,
    'age_min':     20,
    'age_max':     30,

    # --- Simulation box (2D, periodic) ---
    'box_size':    50.0,

    # --- Time ---
    'T_total':     150.0,
    'dt':          0.10,

    # --- Langevin (overdamped) ---
    'gamma':       2.0,     # friction
    'kT':          1.5,     # thermal energy

    # --- Morse potential ---
    'epsilon_0':   6.0,     # base well depth
    'a_morse':     0.40,    # potential width (steepness)
    'r0_morse':    3.5,     # equilibrium distance
    'r_cutoff':    22.0,    # spatial force cutoff
    'r_repulse':   2.5,     # hard-core repulsion radius

    # --- Affinity (Flory-Huggins analogy) ---
    'chi_het':     5.0,     # M-F affinity
    'chi_hom':     0.15,    # same-sex affinity (weak)
    'sigma_age':   4.0,     # age-compatibility Gaussian width (years)

    # --- Relationship (B_ij) dynamics ---
    'k_on':        0.40,    # formation rate
    'k_off':       0.03,    # decay rate
    'r_c':         10.0,    # max distance for relationship building
    'B_thresh':    0.70,    # couple formation threshold
    'B_break':     0.35,    # couple breaking threshold (hysteresis)

    # --- Perturbation (3rd-party field distortion) ---
    'Q_perturb':   1.5,     # 3rd-party screening strength
    'eps_screen':  2.0,     # regularizer

    # --- Stochastic events (cheating / 환승 / 대쉬) ---
    'p_event':     0.0015,  # per M-F pair per step
    'B_boost':     0.25,    # B gain from event
    'B_damage':    0.35,    # B loss to existing couple bonds

    # --- Output ---
    'record_interval': 5,   # record state every N steps
    'seed':        42,
}
