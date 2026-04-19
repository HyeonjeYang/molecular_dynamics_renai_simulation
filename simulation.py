"""Main AffinityMD simulation class."""
import numpy as np

from .params import DEFAULT_PARAMS
from .dynamics import compute_forces, langevin_step
from .relationships import update_B, update_couples
from .events import apply_stochastic_events


class Simulation:
    """AffinityMD simulation.

    Usage:
        sim = Simulation(params)        # params is a dict (full or partial)
        sim.run()
        sim.history                     # recorded trajectories
        sim.coupled_pairs               # list of (i, j) final couples
    """

    def __init__(self, params=None):
        # Merge with defaults (missing keys filled from DEFAULT_PARAMS)
        self.p = dict(DEFAULT_PARAMS)
        if params is not None:
            self.p.update(params)

        self.rng = np.random.default_rng(self.p['seed'])

        # Populations
        self.N_m = int(self.p['N_male'])
        self.N_f = int(self.p['N_female'])
        self.N   = self.N_m + self.N_f

        self.types   = np.array([0] * self.N_m + [1] * self.N_f, dtype=np.int32)
        self.males   = np.where(self.types == 0)[0]
        self.females = np.where(self.types == 1)[0]

        # Ages uniform in [age_min, age_max]
        self.ages = self.rng.uniform(
            self.p['age_min'], self.p['age_max'], self.N)

        # Random initial positions inside box (keep margin from walls)
        L = self.p['box_size']
        self.pos = self.rng.uniform(1.0, L - 1.0, (self.N, 2))

        # Relationship matrix and couple pointers
        self.B      = np.zeros((self.N, self.N))
        self.couple = np.full(self.N, -1, dtype=np.int32)

        # Labels for plots and output
        self.labels = (
            [f"M{i+1}({int(round(self.ages[i]))})" for i in range(self.N_m)]
            + [f"F{i+1}({int(round(self.ages[self.N_m+i]))})"
               for i in range(self.N_f)]
        )

        # Precompute affinity matrices (constant throughout run)
        self._build_affinity_matrices()

        # Histories
        self.history = {
            'time':           [],
            'n_couples':      [],
            'avg_B_mf':       [],
            'max_B_mf':       [],
            'n_events_step':  [],
            'positions':      [],
            'B_matrices':     [],
            'couple_history': [],
        }
        self.event_log = []     # list of (time, i, j, kind)
        self.step_count = 0

    # --------------------------------------------------------------
    def _build_affinity_matrices(self):
        chi_het   = self.p['chi_het']
        chi_hom   = self.p['chi_hom']
        sigma_age = self.p['sigma_age']

        self.chi_mat = np.zeros((self.N, self.N))
        self.age_mat = np.zeros((self.N, self.N))

        for i in range(self.N):
            for j in range(self.N):
                if i == j:
                    continue
                self.chi_mat[i, j] = (chi_het if self.types[i] != self.types[j]
                                      else chi_hom)
                d_age = self.ages[i] - self.ages[j]
                self.age_mat[i, j] = np.exp(
                    -d_age ** 2 / (2.0 * sigma_age ** 2))

    # --------------------------------------------------------------
    def step(self):
        """Advance one dt step."""
        F = compute_forces(self.pos, self.B, self.couple,
                           self.chi_mat, self.age_mat, self.p)

        self.pos = langevin_step(self.pos, F, self.p, self.rng)

        # Periodic boundary (confinement by wrap)
        L = self.p['box_size']
        self.pos = np.mod(self.pos, L)

        # B update and stochastic events
        update_B(self.pos, self.B, self.chi_mat, self.age_mat,
                 self.males, self.females, self.p)

        n_ev, events = apply_stochastic_events(
            self.B, self.couple, self.males, self.females, self.p, self.rng)

        update_couples(self.B, self.couple, self.males, self.females, self.p)

        self.step_count += 1
        t = self.step_count * self.p['dt']
        if events:
            for (i, j, kind) in events:
                self.event_log.append((float(t), i, j, kind))

        return n_ev

    # --------------------------------------------------------------
    def record(self, t, n_ev):
        """Append current state to histories."""
        n_c = int(np.sum(self.couple[self.males] != -1))
        B_mf = self.B[np.ix_(self.males, self.females)]

        self.history['time'].append(float(t))
        self.history['n_couples'].append(n_c)
        self.history['avg_B_mf'].append(float(B_mf.mean()))
        self.history['max_B_mf'].append(float(B_mf.max()))
        self.history['n_events_step'].append(int(n_ev))
        self.history['positions'].append(self.pos.copy())
        self.history['B_matrices'].append(self.B.copy())
        self.history['couple_history'].append(self.couple.copy())

    # --------------------------------------------------------------
    def run(self, verbose=True):
        """Run the full simulation, returning self for chaining."""
        N_steps = int(round(self.p['T_total'] / self.p['dt']))
        rec_int = max(int(self.p['record_interval']), 1)

        if verbose:
            print(f"AffinityMD run: N_m={self.N_m}, N_f={self.N_f}, "
                  f"T={self.p['T_total']}, dt={self.p['dt']} "
                  f"({N_steps} steps)")

        # Initial state
        self.record(0.0, 0)

        milestone = max(N_steps // 10, 1)
        for step in range(N_steps):
            n_ev = self.step()
            t = (step + 1) * self.p['dt']

            if (step + 1) % rec_int == 0:
                self.record(t, n_ev)

            if verbose and (step + 1) % milestone == 0:
                pct = 100.0 * (step + 1) / N_steps
                n_c = int(np.sum(self.couple[self.males] != -1))
                max_c = min(self.N_m, self.N_f)
                print(f"  [{pct:5.1f}%] t={t:6.1f}  "
                      f"couples={n_c}/{max_c}  "
                      f"events_total={len(self.event_log)}")

        if verbose:
            print("Done.")
        return self

    # --------------------------------------------------------------
    @property
    def coupled_pairs(self):
        """List of (male_idx, female_idx) currently coupled."""
        return [(int(i), int(self.couple[i]))
                for i in self.males if self.couple[i] != -1]
