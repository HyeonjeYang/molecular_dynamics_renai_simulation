"""Stochastic romantic events: cheating / 환승 (switching) / 대쉬 (asking out).

At each step, every M-F pair has a small probability p_event of experiencing
a romantic 'event' that boosts their B_ij directly (not via spatial dynamics).
Existing couples of either participant take damage.

This captures rare, discrete social events that short-circuit the otherwise
smooth physical dynamics.
"""


def apply_stochastic_events(B, couple, males, females, params, rng):
    """Apply random events in place on B.

    Returns:
        n_events (int): number of events this step
        events (list of (i, j, kind)): details for logging
    """
    p      = params['p_event']
    boost  = params['B_boost']
    damage = params['B_damage']

    n_events = 0
    events = []

    for i in males:
        i_coup = couple[i] != -1
        for j in females:
            if rng.random() >= p:
                continue

            n_events += 1
            j_coup = couple[j] != -1

            # classify
            if i_coup and j_coup:
                kind = 'double_cheat'   # 양쪽 다 커플인 상태의 불꽃
            elif i_coup or j_coup:
                kind = 'cheat_or_dash'  # 한쪽만 커플 (환승 / 대쉬)
            else:
                kind = 'spark'          # 둘 다 솔로

            # Boost the new pair
            new_B = B[i, j] + boost
            if new_B > 1.0:
                new_B = 1.0
            B[i, j] = new_B
            B[j, i] = new_B

            # Damage existing couple bonds
            if i_coup:
                p_i = couple[i]
                v = B[i, p_i] - damage
                if v < 0.0:
                    v = 0.0
                B[i, p_i] = v
                B[p_i, i] = v

            if j_coup:
                p_j = couple[j]
                v = B[j, p_j] - damage
                if v < 0.0:
                    v = 0.0
                B[j, p_j] = v
                B[p_j, j] = v

            events.append((int(i), int(j), kind))

    return n_events, events
