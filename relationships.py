"""Relationship strength B_ij and couple bookkeeping.

dB_ij/dt = k_on * Theta(r_c - r_ij) * chi_ij * f(age_ij) * (1 - B_ij)
         - k_off * B_ij
"""
import numpy as np


def update_B(pos, B, chi_mat, age_mat, males, females, params):
    """Update B_ij for all M-F pairs in place (Euler step)."""
    dt    = params['dt']
    k_on  = params['k_on']
    k_off = params['k_off']
    r_c2  = params['r_c'] ** 2

    for i in males:
        for j in females:
            dx = pos[j, 0] - pos[i, 0]
            dy = pos[j, 1] - pos[i, 1]
            r2 = dx * dx + dy * dy

            dB = 0.0
            if r2 < r_c2:
                dB += k_on * chi_mat[i, j] * age_mat[i, j] * (1.0 - B[i, j]) * dt
            dB -= k_off * B[i, j] * dt

            new_val = B[i, j] + dB
            if new_val < 0.0:
                new_val = 0.0
            elif new_val > 1.0:
                new_val = 1.0
            B[i, j] = new_val
            B[j, i] = new_val


def update_couples(B, couple, males, females, params):
    """Break couples with B < B_break, then greedily form new couples
    from unpaired M-F pairs with B > B_thresh, highest-B first.

    `couple` is modified in place.
    """
    B_thresh = params['B_thresh']
    B_break  = params['B_break']

    # 1) Break weak couples
    for i in males:
        j = couple[i]
        if j != -1 and B[i, j] < B_break:
            couple[i] = -1
            couple[j] = -1

    # 2) Form new couples (greedy, highest B first)
    candidates = []
    for i in males:
        if couple[i] != -1:
            continue
        for j in females:
            if couple[j] != -1:
                continue
            if B[i, j] > B_thresh:
                candidates.append((B[i, j], i, int(j)))

    candidates.sort(reverse=True)
    for b, i, j in candidates:
        if couple[i] == -1 and couple[j] == -1:
            couple[i] = j
            couple[j] = i
