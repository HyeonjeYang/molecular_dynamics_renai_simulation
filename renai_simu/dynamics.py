"""Physics engine: Morse forces and overdamped Langevin integrator.

Overdamped Langevin equation:
    gamma * dr/dt = -grad U + eta(t),   <eta eta> = 2*gamma*kT*delta

Euler-Maruyama discretization:
    r(t+dt) = r(t) + (F/gamma)*dt + sqrt(2*kT*dt/gamma) * xi,  xi ~ N(0,1)
"""
import numpy as np


# ------------------------------------------------------------------
# Scalar force helpers
# ------------------------------------------------------------------

def morse_force_scalar(r, eps, a, r0):
    """Radial force magnitude from Morse potential.

    U(r) = eps * (1 - exp(-a(r-r0)))^2 - eps
    F(r) = -dU/dr = -2*eps*a * exp(-a(r-r0)) * (1 - exp(-a(r-r0)))

    Sign convention: positive F => repulsive along r_hat.
    """
    exp_term = np.exp(-a * (r - r0))
    return -2.0 * eps * a * exp_term * (1.0 - exp_term)


def soft_core_repulsion(r, r_rep, strength=25.0):
    """Short-range soft-core repulsion inside r_rep to prevent overlap."""
    if r < r_rep:
        return strength * (r_rep - r) / r_rep  # positive => repulsive
    return 0.0


# ------------------------------------------------------------------
# Perturbation field
# ------------------------------------------------------------------

def perturbation_P(i, j, pos, couple, Q_perturb, eps_screen):
    """Field distortion P_ij from third parties near the (i,j) midpoint.

    Analogous to dielectric screening by nearby 'charges'.
    Already-coupled third parties contribute only weakly (they're 'taken').
    Returns a scalar clamped at >= 0.05.
    """
    mid = 0.5 * (pos[i] + pos[j])
    diff = pos - mid                          # (N, 2)
    d2 = np.einsum('ij,ij->i', diff, diff)    # (N,)

    # Screening weights
    Q = np.where(couple != -1, 0.25 * Q_perturb, Q_perturb)

    # Exclude i and j themselves
    Q[i] = 0.0
    Q[j] = 0.0

    P = 1.0 - np.sum(Q / (d2 + eps_screen))
    return max(P, 0.05)


# ------------------------------------------------------------------
# Force accumulation
# ------------------------------------------------------------------

def compute_forces(pos, B, couple, chi_mat, age_mat, params):
    """Sum pairwise forces into (N, 2) array.

    Effective well depth:
        eps_eff = eps_0 * chi_ij * B_ij * f(age_ij) * P_ij(t)
    """
    N = pos.shape[0]
    F = np.zeros_like(pos)

    eps0     = params['epsilon_0']
    a        = params['a_morse']
    r0       = params['r0_morse']
    r_cut    = params['r_cutoff']
    r_rep    = params['r_repulse']
    Q        = params['Q_perturb']
    eps_scr  = params['eps_screen']

    for i in range(N):
        for j in range(i + 1, N):
            dx = pos[j, 0] - pos[i, 0]
            dy = pos[j, 1] - pos[i, 1]
            r = np.hypot(dx, dy)
            if r < 1e-8:
                # avoid singularity: random kick direction
                dx, dy = np.random.randn(2) * 1e-3
                r = np.hypot(dx, dy)
            rhat_x = dx / r
            rhat_y = dy / r

            Fmag = 0.0

            # Long-range Morse (conditional on cutoff)
            if r < r_cut:
                P = perturbation_P(i, j, pos, couple, Q, eps_scr)
                eps_eff = eps0 * chi_mat[i, j] * B[i, j] * age_mat[i, j] * P
                Fmag += morse_force_scalar(r, eps_eff, a, r0)

            # Short-range repulsion always
            Fmag += soft_core_repulsion(r, r_rep)

            Fx = Fmag * rhat_x
            Fy = Fmag * rhat_y

            # force on i is -F (since rhat points i->j and F>0 is repulsive)
            F[i, 0] -= Fx
            F[i, 1] -= Fy
            F[j, 0] += Fx
            F[j, 1] += Fy

    return F


# ------------------------------------------------------------------
# Integrator
# ------------------------------------------------------------------

def langevin_step(pos, F, params, rng):
    """Single overdamped Euler-Maruyama step (returns new position)."""
    dt    = params['dt']
    gamma = params['gamma']
    kT    = params['kT']
    noise_amp = np.sqrt(2.0 * kT * dt / gamma)
    return pos + (F / gamma) * dt + noise_amp * rng.standard_normal(pos.shape)
