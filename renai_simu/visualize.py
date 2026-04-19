"""Visualization utilities for AffinityMD simulation."""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import networkx as nx


# ------------------------------------------------------------------
def plot_time_evolution(sim, path=None):
    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    fig.suptitle('AffinityMD — Time Evolution',
                 fontsize=14, fontweight='bold')

    t = np.asarray(sim.history['time'])

    axes[0].plot(t, sim.history['n_couples'], color='#D62728', lw=2)
    axes[0].axhline(min(sim.N_m, sim.N_f), ls='--', color='gray', alpha=0.6,
                    label='Max possible')
    axes[0].set_ylabel('# Couples')
    axes[0].legend(loc='best')
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(t, sim.history['avg_B_mf'], color='#1F77B4', lw=2,
                 label='Mean B (M-F)')
    axes[1].plot(t, sim.history['max_B_mf'], color='#1F77B4', lw=1.5,
                 ls='--', alpha=0.6, label='Max B (M-F)')
    axes[1].axhline(sim.p['B_thresh'], color='green', ls=':',
                    label=f"B_thresh={sim.p['B_thresh']}")
    axes[1].axhline(sim.p['B_break'], color='orange', ls=':',
                    label=f"B_break={sim.p['B_break']}")
    axes[1].set_ylabel('Relationship strength B')
    axes[1].set_ylim(0, 1.05)
    axes[1].legend(fontsize=8, loc='best')
    axes[1].grid(True, alpha=0.3)

    cum = np.cumsum(sim.history['n_events_step'])
    axes[2].plot(t, cum, color='purple', lw=2)
    axes[2].set_xlabel('Time')
    axes[2].set_ylabel('Cumulative stochastic\nevents')
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    if path:
        plt.savefig(path, dpi=150, bbox_inches='tight')
        plt.close(fig)
    return fig


# ------------------------------------------------------------------
def plot_relationship_graph(sim, path=None, min_B_show=0.10):
    fig, ax = plt.subplots(figsize=(9, 8))
    ax.set_title('Final Relationship Network',
                 fontsize=13, fontweight='bold')

    G = nx.Graph()
    for i in range(sim.N):
        G.add_node(i)

    node_colors = ['#4A90E2' if sim.types[i] == 0 else '#E25C87'
                   for i in range(sim.N)]

    edge_widths = []
    edge_colors = []
    edges = []
    for i in sim.males:
        for j in sim.females:
            b = sim.B[i, j]
            if b > min_B_show:
                G.add_edge(int(i), int(j), weight=float(b))
                edges.append((int(i), int(j)))
                edge_widths.append(b * 5.0)
                edge_colors.append('#D62728' if sim.couple[i] == j
                                   else '#AAAAAA')

    # Bipartite layout
    pos = {}
    for idx, i in enumerate(sim.males):
        pos[int(i)] = (-1.0, (idx - (sim.N_m - 1) / 2.0) * 1.3)
    for idx, j in enumerate(sim.females):
        pos[int(j)] = (1.0, (idx - (sim.N_f - 1) / 2.0) * 1.3)

    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=1100,
                           edgecolors='black', linewidths=1.0, ax=ax)
    nx.draw_networkx_labels(
        G, pos, labels={i: sim.labels[i] for i in range(sim.N)},
        font_size=7, ax=ax)
    if edges:
        nx.draw_networkx_edges(G, pos, edgelist=edges,
                               width=edge_widths, edge_color=edge_colors,
                               alpha=0.7, ax=ax)

    legend_elements = [
        Patch(facecolor='#4A90E2', edgecolor='black', label='Male'),
        Patch(facecolor='#E25C87', edgecolor='black', label='Female'),
        Line2D([0], [0], color='#D62728', lw=3, label='Couple'),
        Line2D([0], [0], color='#AAAAAA', lw=2,
               label=f'Attraction (B > {min_B_show})'),
    ]
    ax.legend(handles=legend_elements, loc='lower center', fontsize=9, ncol=2,
              bbox_to_anchor=(0.5, -0.06))
    ax.set_xlim(-1.8, 1.8)
    ax.axis('off')

    plt.tight_layout()
    if path:
        plt.savefig(path, dpi=150, bbox_inches='tight')
        plt.close(fig)
    return fig


# ------------------------------------------------------------------
def plot_B_heatmap(sim, path=None, top_k=5):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle('Relationship Matrix B', fontsize=14, fontweight='bold')

    # Final B heatmap (magma: dark = low B, bright = high B)
    ax = axes[0]
    im = ax.imshow(sim.B, cmap='magma', vmin=0, vmax=1, aspect='auto')
    ax.set_title('Final B matrix')
    ax.set_xticks(range(sim.N))
    ax.set_xticklabels(sim.labels, rotation=45, ha='right', fontsize=7)
    ax.set_yticks(range(sim.N))
    ax.set_yticklabels(sim.labels, fontsize=7)
    plt.colorbar(im, ax=ax, label='$B_{ij}$')

    # Cyan rectangles for couples — pops against magma's warm tones
    for i in sim.males:
        j = int(sim.couple[i])
        if j != -1:
            ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                       fill=False, edgecolor='#00E5FF',
                                       lw=2.8))
            ax.add_patch(plt.Rectangle((i - 0.5, j - 0.5), 1, 1,
                                       fill=False, edgecolor='#00E5FF',
                                       lw=2.8))

    # Top-k pair B(t) evolution
    ax2 = axes[1]
    all_pairs = [(sim.B[i, j], int(i), int(j))
                 for i in sim.males for j in sim.females]
    all_pairs.sort(reverse=True)
    t = np.asarray(sim.history['time'])

    colors = plt.cm.tab10(np.linspace(0, 1, max(top_k, 1)))
    for idx, (b, i, j) in enumerate(all_pairs[:top_k]):
        B_hist = [bm[i, j] for bm in sim.history['B_matrices']]
        ax2.plot(t, B_hist, lw=1.8, color=colors[idx],
                 label=f"{sim.labels[i]} — {sim.labels[j]}")

    ax2.axhline(sim.p['B_thresh'], color='green', ls=':', alpha=0.7,
                label='B_thresh')
    ax2.axhline(sim.p['B_break'], color='orange', ls=':', alpha=0.7,
                label='B_break')
    ax2.set_xlabel('Time')
    ax2.set_ylabel('B_ij')
    ax2.set_title(f'Top {top_k} M-F pairs (by final B)')
    ax2.legend(fontsize=7, loc='best')
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(0, 1.05)

    plt.tight_layout()
    if path:
        plt.savefig(path, dpi=150, bbox_inches='tight')
        plt.close(fig)
    return fig


# ------------------------------------------------------------------
def plot_snapshots(sim, path=None, n_snapshots=4):
    """Particle positions at equally-spaced snapshots in time."""
    rows = 2
    cols = (n_snapshots + rows - 1) // rows
    fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 5 * rows))
    axes = np.atleast_1d(axes).flatten()
    fig.suptitle('Particle Positions — Snapshots',
                 fontsize=14, fontweight='bold')

    n_rec = len(sim.history['time'])
    n_snapshots = min(n_snapshots, n_rec)
    indices = np.linspace(0, n_rec - 1, n_snapshots).astype(int)
    L = sim.p['box_size']

    last_ax_idx = 0
    for pidx, ridx in enumerate(indices):
        ax = axes[pidx]
        last_ax_idx = pidx
        t   = sim.history['time'][ridx]
        pos = sim.history['positions'][ridx]
        coup = sim.history['couple_history'][ridx]
        B_s = sim.history['B_matrices'][ridx]

        ax.set_title(f't = {t:.1f}', fontsize=11)
        ax.set_xlim(0, L); ax.set_ylim(0, L)
        ax.set_aspect('equal')
        ax.set_facecolor('#f7f7f7')
        ax.grid(True, alpha=0.25)

        # Weak attraction lines
        for i in sim.males:
            for j in sim.females:
                if B_s[i, j] > 0.25 and coup[i] != j:
                    ax.plot([pos[i, 0], pos[j, 0]],
                            [pos[i, 1], pos[j, 1]],
                            color='#4A90E2', lw=1.0,
                            alpha=min(B_s[i, j] * 0.6, 0.6),
                            zorder=1)

        # Couple lines
        for i in sim.males:
            j = int(coup[i])
            if j != -1:
                ax.plot([pos[i, 0], pos[j, 0]],
                        [pos[i, 1], pos[j, 1]],
                        color='#D62728', lw=2.2, alpha=0.9, zorder=2)

        # Particles
        for i in range(sim.N):
            c = '#4A90E2' if sim.types[i] == 0 else '#E25C87'
            m = '^' if sim.types[i] == 0 else 'o'
            coupled = coup[i] != -1
            ec = 'gold' if coupled else 'black'
            lw = 2.0 if coupled else 0.6
            ax.scatter(pos[i, 0], pos[i, 1], c=c, s=140, marker=m,
                       edgecolors=ec, linewidths=lw, zorder=3)
            ax.annotate(sim.labels[i], (pos[i, 0], pos[i, 1]),
                        fontsize=5.5, ha='center', va='center', zorder=4)

        n_c = int(np.sum(coup[sim.males] != -1))
        ax.set_xlabel(f'Couples: {n_c}/{min(sim.N_m, sim.N_f)}')

    # hide unused subplots
    for k in range(last_ax_idx + 1, len(axes)):
        axes[k].axis('off')

    legend_elements = [
        Line2D([0], [0], marker='^', color='w', markerfacecolor='#4A90E2',
               markersize=12, label='Male', markeredgecolor='black'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#E25C87',
               markersize=12, label='Female', markeredgecolor='black'),
        Line2D([0], [0], color='#D62728', lw=2.2, label='Couple bond'),
        Line2D([0], [0], color='#4A90E2', lw=1.2, alpha=0.6,
               label='Attraction (B > 0.25)'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=4,
               fontsize=9, bbox_to_anchor=(0.5, -0.01))

    plt.tight_layout(rect=[0, 0.04, 1, 1])
    if path:
        plt.savefig(path, dpi=150, bbox_inches='tight')
        plt.close(fig)
    return fig


# ------------------------------------------------------------------
def plot_all(sim, out_dir):
    """Save every figure as PNG into out_dir."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    plot_time_evolution(sim,    out_dir / 'time_evolution.png')
    plot_relationship_graph(sim, out_dir / 'relationship_graph.png')
    plot_B_heatmap(sim,         out_dir / 'B_matrix.png')
    plot_snapshots(sim,         out_dir / 'snapshots.png', n_snapshots=4)
