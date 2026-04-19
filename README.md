# AffinityMD (`renai_simu`)
##分子動力学で恋してる♡

A stochastic particle simulation of romantic pair formation in a confined 2D box.

**Idea:** Hyeonje Yang.
**Implementation assistance:** Claude (Anthropic).

---

## The idea

Treat a group of single people as diffusing particles. Different-sex pairs feel
a mutual attraction whose strength evolves over time as they interact. Nearby
third parties distort the attraction like charges in a dielectric. Rare
stochastic "events" (sparks, 환승, 대쉬) jolt the relationship network.
Run it forward and see who ends up with whom.

---

## Physics

### 1. Motion — overdamped Langevin

$$\gamma \dot{\mathbf{r}}_i = -\nabla_i \sum_{j \neq i} U_{ij}(r_{ij}) + \boldsymbol{\eta}_i(t),
\qquad \langle \eta_\alpha(t)\,\eta_\beta(t') \rangle = 2\gamma k_B T \,\delta_{\alpha\beta}\,\delta(t-t')$$

Integrated by **Euler-Maruyama** in a periodic box. Friction $\gamma$ dominates
inertia, so this is the diffusion-limited regime relevant to biology/chemistry
rather than ballistic Newtonian motion.

### 2. Pair potential — Morse with dynamic well depth

$$U_{ij}(r) = \varepsilon_{ij}^{\text{eff}} \left(1 - e^{-a(r-r_0)}\right)^2 - \varepsilon_{ij}^{\text{eff}}$$

Morse is the standard soft diatomic potential — finite well, hard core, smooth
dissociation. The well depth is *not* constant; it factorises into four pieces:

$$\varepsilon_{ij}^{\text{eff}}(t) \;=\; \varepsilon_0 \cdot \chi_{ij} \cdot B_{ij}(t) \cdot f(\Delta\text{age}_{ij}) \cdot P_{ij}(t)$$

| factor | meaning | physical analogy |
|---|---|---|
| $\chi_{ij}$ | sex compatibility | **Flory-Huggins** mixing parameter — heterogeneous pairs mix favorably ($\chi_\text{het} \gg \chi_\text{hom}$), homogeneous pairs are close to inert |
| $B_{ij}$ | relationship strength (see §3) | slow reaction coordinate / bond-occupancy variable |
| $f(\Delta\text{age})$ | age compatibility | Gaussian similarity kernel: $\exp\!\left(-\Delta\text{age}^2 / 2\sigma_\text{age}^2\right)$ |
| $P_{ij}(t)$ | 3rd-party screening | **dielectric screening** — nearby "charges" weaken the pair's field |

Screening by third parties:

$$P_{ij}(t) = 1 - \sum_{k \neq i,j} \frac{Q_k}{|\mathbf{r}_k - \mathbf{r}_{ij}^{\text{mid}}|^2 + \epsilon}$$

### 3. What is $B_{ij}$?

$B_{ij} \in [0, 1]$ is a **dynamic order parameter** for the bond between $i$ and $j$.
Think of it as a bond-occupancy variable — how "fused" the pair currently is.
It follows first-order chemical kinetics:

$$\frac{dB_{ij}}{dt} \;=\; k_\text{on}\, \Theta(r_c - r_{ij})\, \chi_{ij}\, f(\Delta\text{age})\,(1 - B_{ij})
\;-\; k_\text{off}\, B_{ij}$$

- When $B_{ij}$ is high, the Morse well is deep → the pair physically attracts.
- When $B_{ij}$ is low, the potential is nearly flat → they diffuse freely past each other.
- **Position dynamics drive $B$, and $B$ deepens the well — positive feedback**, analogous to cooperative binding.

### 4. Couples — a hysteretic switch

A pair is declared a **couple** when $B_{ij} > B_\text{thresh}$, and is broken when
$B_{ij} < B_\text{break}$, with $B_\text{break} < B_\text{thresh}$.

The hysteresis gap is **physically principled**: forming and breaking a bond
both require crossing an activation barrier — exactly as in protein folding,
ligand binding, or any first-order transition. It also suppresses numerical
flickering near the threshold.

### 5. Stochastic events

With small probability per M–F pair per step, a discrete event fires: it boosts
$B_{ij}$ by $\Delta_+$ and damages any existing couple-bond of $i$ or $j$ by
$\Delta_-$. Classified as `spark` (both singles), `cheat_or_dash` (one already
coupled), or `double_cheat` (both coupled).

---

## Usage

**Terminal:**
```bash
python run_terminal.py --N_male 10 --N_female 10 --T 200 --seed 7
python run_terminal.py --help          # all knobs
```

**Notebook:** open `run_notebook.ipynb`, edit `params_override` in cell 2, run all cells.

**Output folders are separated** — terminal runs go to `outputs/terminal_run_<timestamp>/`
and notebook runs go to `outputs/notebook_run_<timestamp>/`, so they never overwrite each other.

---

## Output files

| file | description |
|---|---|
| `parmsdict.txt` | parameters + final summary (human-readable) |
| `summary.json` | same info, machine-readable |
| `event_log.tsv` | every stochastic event (time, pair, kind) |
| `time_evolution.png` | #couples, mean/max $B$, cumulative events |
| `relationship_graph.png` | final bipartite M–F attraction network |
| `B_matrix.png` | final $B$ heatmap (magma) + top-5 pair $B(t)$ |
| `snapshots.png` | particle positions at 4 time points |

---

## Layout

```
affinity_md/
├── renai_simu/           ← core package
│   ├── params.py         ← DEFAULT_PARAMS
│   ├── dynamics.py       ← forces + Langevin integrator
│   ├── relationships.py  ← B kinetics + couple bookkeeping
│   ├── events.py         ← stochastic events
│   ├── simulation.py     ← Simulation class
│   ├── visualize.py      ← plotting
│   └── io_utils.py       ← text/JSON outputs
├── run_terminal.py       ← CLI
├── run_notebook.ipynb    ← Jupyter
└── outputs/
```
