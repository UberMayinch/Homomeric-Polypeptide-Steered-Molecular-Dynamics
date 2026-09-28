#!/usr/bin/env python3
"""
Steered Molecular Dynamics Analysis for Homomeric Polypeptides
==============================================================
Comprehensive analysis of SMD trajectories for all 20 amino acids.
Computes PMF via Jarzynski's equality, Kuhn lengths, force-extension
curves, and generates publication-quality figures.

Usage:
    python analysis.py
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from scipy.ndimage import gaussian_filter1d
import os, re, warnings
warnings.filterwarnings('ignore')

# ============================================================
# SECTION 1: Setup & Configuration
# ============================================================

# ── All 20 amino acids ──
proteins = [
    'alanine', 'arginine', 'aspartic_acid', 'asparagine', 'cysteine',
    'glutamine', 'glutamic_acid', 'glycine', 'histidine', 'isoleucine',
    'leucine', 'lysine', 'methionine', 'phenylalanine', 'proline',
    'serine', 'threonine', 'tryptophan', 'tyrosine', 'valine'
]

# ── Physicochemical groups ──
groups = {
    'Nonpolar':           ['alanine', 'glycine', 'valine', 'leucine',
                           'isoleucine', 'proline', 'phenylalanine',
                           'methionine', 'tryptophan'],
    'Polar uncharged':    ['serine', 'threonine', 'cysteine', 'tyrosine',
                           'asparagine', 'glutamine'],
    'Positively charged': ['arginine', 'lysine', 'histidine'],
    'Negatively charged': ['aspartic_acid', 'glutamic_acid'],
}

group_colors = {
    'Nonpolar':           '#e69f00',
    'Polar uncharged':    '#56b4e9',
    'Positively charged': '#009e73',
    'Negatively charged': '#d55e00',
}

# ── Kyte-Doolittle hydrophobicity, MW, side-chain volume, charge, H-bond donors ──
aa_properties = {
    'alanine':       ( 1.8,  89.1,  67,   0,   0),
    'arginine':      (-4.5, 174.2, 148,  +1,   5),
    'aspartic_acid': (-3.5, 133.1,  91,  -1,   1),
    'asparagine':    (-3.5, 132.1,  96,   0,   2),
    'cysteine':      ( 2.5, 121.2,  86,   0,   1),
    'glutamine':     (-3.5, 146.2, 114,   0,   2),
    'glutamic_acid': (-3.5, 147.1, 109,  -1,   1),
    'glycine':       (-0.4,  75.0,  48,   0,   0),
    'histidine':     (-3.2, 155.2, 118, +0.1,  2),
    'isoleucine':    ( 4.5, 131.2, 124,   0,   0),
    'leucine':       ( 3.8, 131.2, 124,   0,   0),
    'lysine':        (-3.9, 146.2, 135,  +1,   3),
    'methionine':    ( 1.9, 149.2, 124,   0,   0),
    'phenylalanine': ( 2.8, 165.2, 135,   0,   0),
    'proline':       (-1.6, 115.1,  90,   0,   0),
    'serine':        (-0.8, 105.1,  73,   0,   1),
    'threonine':     (-0.7, 119.1,  93,   0,   1),
    'tryptophan':    (-0.9, 204.2, 163,   0,   1),
    'tyrosine':      (-1.3, 181.2, 141,   0,   2),
    'valine':        ( 4.2, 117.1, 105,   0,   0),
}

base_dir = '/home/kiranr/Documents/3-2/PSCM/PROJECT-1'
fig_dir  = os.path.join(base_dir, 'figures')
os.makedirs(fig_dir, exist_ok=True)

n_traj = 10
dt     = 0.1       # ps  (Tclfreq=50 steps × 2 fs/step)
v      = 0.01      # Å/ps  (0.00002 Å/step ÷ 0.002 ps/step)
k_B    = 1.987e-3  # kcal mol⁻¹ K⁻¹
Temp   = 300       # K
T      = k_B * Temp


def get_group(protein):
    for g, members in groups.items():
        if protein in members:
            return g
    return 'Unknown'


# ============================================================
# SECTION 2: Data Loading & Computation Functions
# ============================================================

def parse_constraints(protein):
    """Parse c1z and c2z from smd.tcl for a given protein."""
    tcl_path = f"{base_dir}/{protein}/smd.tcl"
    c1z, c2z = None, None
    with open(tcl_path, 'r') as f:
        for line in f:
            m1 = re.match(r'\s*set\s+c1z\s+([-\d.]+)', line)
            m2 = re.match(r'\s*set\s+c2z\s+([-\d.]+)', line)
            if m1: c1z = float(m1.group(1))
            if m2: c2z = float(m2.group(1))
    if c1z is None or c2z is None:
        raise ValueError(f"Could not parse c1z/c2z from {tcl_path}")
    return c1z, c2z


def load_trajectories(protein, n_traj=10):
    """Load all SMD trajectories for one protein.
    Returns arrays of shape (n_timesteps, n_traj):
      t_arr  – time (ps)
      z_arr  – z-coordinate of pulled atom (Å)
      f_arr  – z-force on pulled atom (kcal/mol/Å)
    """
    t_list, z_list, f_list = [], [], []
    for i in range(n_traj):
        fname = f"{base_dir}/{protein}/{protein}_smd_tcl_1{i}.out"
        data  = np.loadtxt(fname)
        t_list.append(data[:, 0])
        z_list.append(data[:, 1])
        f_list.append(data[:, 2])
    return np.column_stack(t_list), np.column_stack(z_list), np.column_stack(f_list)


def compute_extension_and_work(t_arr, f_arr, c1z, c2z):
    """Compute constraint-point separation c(t) and cumulative work W.
    c(t) = (c2z - c1z) + v*t  (Å)
    W    = v·dt·cumsum(-F)   (kcal/mol)
    """
    t_vec = t_arr[:, 0]
    c     = (c2z - c1z) + v * t_vec
    W     = v * dt * np.cumsum(-f_arr, axis=0)
    return c, W


def compute_pmf(W, T):
    """Three PMF estimates from Jarzynski's equality (numerically stable)."""
    exponent  = -W / T
    max_exp   = np.max(exponent, axis=1, keepdims=True)
    Fexp = -T * (np.log(np.mean(np.exp(exponent - max_exp), axis=1))
                 + max_exp.squeeze())
    F1 = np.mean(W, axis=1)
    F2 = np.mean(W, axis=1) - 1/(2*T) * (np.mean(W**2, axis=1)
                                           - np.mean(W, axis=1)**2)
    return Fexp, F1, F2


def compute_kuhn_length(c_shifted, f_arr, N_residues=10):
    """Estimate Kuhn length from the low-extension linear regime.
    In the linear regime: <z> ≈ (N * a_Kuhn^2 / (3 k_B T)) * f
    We use the first ~20% of extension range for the linear fit.
    """
    mean_f = np.mean(np.abs(f_arr), axis=1)
    n_pts  = len(c_shifted)
    n_fit  = max(10, n_pts // 5)
    # Use low-extension regime
    ext = c_shifted[:n_fit]
    frc = mean_f[:n_fit]
    # Avoid division by zero
    mask = frc > 1e-10
    if np.sum(mask) < 3:
        return np.nan
    try:
        slope, _, _, _, _ = stats.linregress(frc[mask], ext[mask])
        # slope = N * a_Kuhn^2 / (3 k_B T)
        # a_Kuhn = sqrt(3 k_B T * slope / N)
        if slope > 0:
            a_kuhn = np.sqrt(3 * T * slope / N_residues)
        else:
            a_kuhn = np.nan
    except Exception:
        a_kuhn = np.nan
    return a_kuhn


# ============================================================
# SECTION 3: Compute All Results
# ============================================================

def compute_all_results():
    """Process all amino acids and return results dictionary."""
    results = {}
    for protein in proteins:
        print(f"Processing {protein}...", end=' ')
        try:
            c1z, c2z            = parse_constraints(protein)
            t_arr, z_arr, f_arr = load_trajectories(protein, n_traj)
            c, W                = compute_extension_and_work(t_arr, f_arr, c1z, c2z)
            c_shifted           = c + c1z - c2z  # starts at 0
            Fexp, F1, F2        = compute_pmf(W, T)

            W_final   = W[-1, :]
            mean_W    = float(np.mean(W_final))
            std_W     = float(np.std(W_final))
            delta_F   = float(F2[-1] - F2[0])
            W_diss    = mean_W - delta_F

            peak_F    = np.max(np.abs(f_arr), axis=0)
            mean_rupF = float(np.mean(peak_F))
            std_rupF  = float(np.std(peak_F))

            barrier_idx = int(np.argmax(F2))
            x_barrier   = float(c_shifted[barrier_idx])

            a_kuhn = compute_kuhn_length(c_shifted, f_arr, N_residues=10)

            results[protein] = {
                'c': c_shifted, 'W': W,
                'F1': F1, 'F2': F2, 'Fexp': Fexp,
                'W_final': W_final,
                'mean_W': mean_W, 'std_W': std_W,
                'delta_F': delta_F, 'W_diss': W_diss,
                'mean_rupF': mean_rupF, 'std_rupF': std_rupF,
                'x_barrier': x_barrier,
                'z_arr': z_arr, 'f_arr': f_arr, 't_arr': t_arr,
                'a_kuhn': a_kuhn,
                'c1z': c1z, 'c2z': c2z,
            }
            print(f"c1z={c1z}, c2z={c2z}, sep={c2z-c1z:.3f} Å, "
                  f"a_Kuhn={a_kuhn:.3f} Å" if not np.isnan(a_kuhn)
                  else "a_Kuhn=N/A")
        except Exception as e:
            print(f"ERROR: {e}")

    print(f"\nDone. {len(results)}/{len(proteins)} amino acids loaded.")
    return results


# ============================================================
# SECTION 4: Individual Amino Acid Plots
# ============================================================

def plot_individual(protein, r, save=True):
    """4-panel figure for one amino acid."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'{protein.replace("_"," ").title()} — SMD Analysis',
                 fontsize=16, fontweight='bold')
    c = r['c']

    # (a) Force vs Extension
    ax = axes[0, 0]
    for j in range(r['f_arr'].shape[1]):
        ax.plot(c, r['f_arr'][:, j], alpha=0.3, lw=0.5)
    ax.plot(c, np.mean(r['f_arr'], axis=1), 'k-', lw=2, label='Mean')
    ax.set_xlabel('Extension (Å)')
    ax.set_ylabel('Force (kcal/mol/Å)')
    ax.set_title('(a) Force–Extension')
    ax.legend(fontsize=8)

    # (b) Work vs Extension
    ax = axes[0, 1]
    for j in range(r['W'].shape[1]):
        ax.plot(c, r['W'][:, j], alpha=0.3, lw=0.5)
    ax.plot(c, np.mean(r['W'], axis=1), 'k-', lw=2, label='Mean')
    ax.set_xlabel('Extension (Å)')
    ax.set_ylabel('Work (kcal/mol)')
    ax.set_title('(b) Work–Extension')
    ax.legend(fontsize=8)

    # (c) PMF estimates
    ax = axes[1, 0]
    ax.plot(c, r['F1'],   'g-',  lw=2, label='Avg Work (1st cumulant)')
    ax.plot(c, r['F2'],   'b-',  lw=2, label='2nd cumulant')
    ax.plot(c, r['Fexp'], 'k-',  lw=2, label='Jarzynski exp. avg.')
    ax.set_xlabel('Extension (Å)')
    ax.set_ylabel('PMF (kcal/mol)')
    ax.set_title('(c) PMF Estimates')
    ax.legend(fontsize=8)

    # (d) Work distribution
    ax = axes[1, 1]
    ax.hist(r['W_final'], bins='auto', edgecolor='black', alpha=0.7,
            color=group_colors.get(get_group(protein), '#888888'))
    ax.axvline(r['mean_W'], color='red', ls='--', lw=2, label=f'Mean={r["mean_W"]:.2f}')
    ax.axvline(r['delta_F'], color='blue', ls='--', lw=2, label=f'ΔF={r["delta_F"]:.2f}')
    ax.set_xlabel('Final Work (kcal/mol)')
    ax.set_ylabel('Count')
    ax.set_title('(d) Work Distribution')
    ax.legend(fontsize=8)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    if save:
        fig.savefig(os.path.join(fig_dir, f'{protein}_analysis.pdf'),
                    dpi=150, bbox_inches='tight')
    plt.close(fig)


# ============================================================
# SECTION 5: Comparative Group Analysis
# ============================================================

def plot_force_extension_by_group(results, save=True):
    """Force-extension curves grouped by physicochemical class."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Force–Extension Curves by Amino Acid Group',
                 fontsize=16, fontweight='bold')
    for idx, (gname, members) in enumerate(groups.items()):
        ax = axes[idx // 2, idx % 2]
        for p in members:
            if p in results:
                r = results[p]
                mean_f = np.mean(r['f_arr'], axis=1)
                ax.plot(r['c'], gaussian_filter1d(mean_f, sigma=5),
                        label=p.replace('_', ' ').title(), lw=1.5)
        ax.set_title(gname, fontsize=13, color=group_colors[gname])
        ax.set_xlabel('Extension (Å)')
        ax.set_ylabel('Mean Force (kcal/mol/Å)')
        ax.legend(fontsize=7, ncol=2)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    if save:
        fig.savefig(os.path.join(fig_dir, 'force_extension_by_group.pdf'),
                    dpi=150, bbox_inches='tight')
    plt.close(fig)


def plot_pmf_by_group(results, save=True):
    """PMF (2nd cumulant) grouped by physicochemical class."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('PMF (2nd Cumulant) by Amino Acid Group',
                 fontsize=16, fontweight='bold')
    for idx, (gname, members) in enumerate(groups.items()):
        ax = axes[idx // 2, idx % 2]
        for p in members:
            if p in results:
                r = results[p]
                ax.plot(r['c'], r['F2'],
                        label=p.replace('_', ' ').title(), lw=1.5)
        ax.set_title(gname, fontsize=13, color=group_colors[gname])
        ax.set_xlabel('Extension (Å)')
        ax.set_ylabel('PMF (kcal/mol)')
        ax.legend(fontsize=7, ncol=2)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    if save:
        fig.savefig(os.path.join(fig_dir, 'pmf_by_group.pdf'),
                    dpi=150, bbox_inches='tight')
    plt.close(fig)


def plot_bar_comparisons(results, save=True):
    """Bar charts comparing key metrics across all amino acids."""
    sorted_proteins = sorted(results.keys(),
                             key=lambda p: results[p]['delta_F'], reverse=True)
    colors = [group_colors[get_group(p)] for p in sorted_proteins]
    labels = [p.replace('_', ' ').title()[:8] for p in sorted_proteins]

    fig, axes = plt.subplots(2, 2, figsize=(18, 12))
    fig.suptitle('Comparative Metrics Across All 20 Amino Acids',
                 fontsize=16, fontweight='bold')

    # (a) PMF barrier height
    ax = axes[0, 0]
    vals = [results[p]['delta_F'] for p in sorted_proteins]
    ax.bar(range(len(vals)), vals, color=colors)
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(labels, rotation=60, ha='right', fontsize=7)
    ax.set_ylabel('ΔF (kcal/mol)')
    ax.set_title('(a) PMF Barrier Height (2nd cumulant)')

    # (b) Mean rupture force
    ax = axes[0, 1]
    vals = [results[p]['mean_rupF'] for p in sorted_proteins]
    errs = [results[p]['std_rupF'] for p in sorted_proteins]
    ax.bar(range(len(vals)), vals, yerr=errs, color=colors, capsize=3)
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(labels, rotation=60, ha='right', fontsize=7)
    ax.set_ylabel('Peak Force (kcal/mol/Å)')
    ax.set_title('(b) Mean Rupture Force')

    # (c) Dissipated work
    ax = axes[1, 0]
    vals = [results[p]['W_diss'] for p in sorted_proteins]
    ax.bar(range(len(vals)), vals, color=colors)
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(labels, rotation=60, ha='right', fontsize=7)
    ax.set_ylabel('W_diss (kcal/mol)')
    ax.set_title('(c) Dissipated Work')

    # (d) Kuhn length
    ax = axes[1, 1]
    vals = [results[p].get('a_kuhn', np.nan) for p in sorted_proteins]
    ax.bar(range(len(vals)), vals, color=colors)
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(labels, rotation=60, ha='right', fontsize=7)
    ax.set_ylabel('Kuhn Length (Å)')
    ax.set_title('(d) Estimated Kuhn Length')

    # Legend for groups
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=c, label=g) for g, c in group_colors.items()]
    fig.legend(handles=legend_elements, loc='lower center', ncol=4, fontsize=9)
    plt.tight_layout(rect=[0, 0.05, 1, 0.95])
    if save:
        fig.savefig(os.path.join(fig_dir, 'bar_comparisons.pdf'),
                    dpi=150, bbox_inches='tight')
    plt.close(fig)


# ============================================================
# SECTION 6: Summary Statistics Table
# ============================================================

def print_summary_table(results):
    """Print a formatted summary table and save as CSV."""
    header = (f"{'Amino Acid':<18} {'Group':<22} {'⟨W⟩':>8} {'σ_W':>8} "
              f"{'ΔF':>8} {'W_diss':>8} {'⟨F_rup⟩':>8} {'σ_F':>8} "
              f"{'a_Kuhn':>8} {'x‡':>8}")
    print("\n" + "=" * len(header))
    print("SUMMARY TABLE — All 20 Amino Acids")
    print("=" * len(header))
    print(header)
    print("-" * len(header))

    rows = []
    for p in sorted(results.keys()):
        r = results[p]
        g = get_group(p)
        row = (f"{p:<18} {g:<22} {r['mean_W']:8.2f} {r['std_W']:8.2f} "
               f"{r['delta_F']:8.2f} {r['W_diss']:8.2f} "
               f"{r['mean_rupF']:8.2f} {r['std_rupF']:8.2f} "
               f"{r.get('a_kuhn', float('nan')):8.3f} {r['x_barrier']:8.2f}")
        print(row)
        rows.append([p, g, r['mean_W'], r['std_W'], r['delta_F'],
                     r['W_diss'], r['mean_rupF'], r['std_rupF'],
                     r.get('a_kuhn', float('nan')), r['x_barrier']])
    print("=" * len(header))

    # Save CSV
    csv_path = os.path.join(fig_dir, 'summary_table.csv')
    with open(csv_path, 'w') as f:
        f.write("Amino Acid,Group,Mean Work,Std Work,Delta F,W_diss,"
                "Mean Rupture Force,Std Rupture Force,Kuhn Length,Barrier Position\n")
        for row in rows:
            f.write(','.join(str(x) for x in row) + '\n')
    print(f"\nTable saved to {csv_path}")


# ============================================================
# SECTION 7: Correlation Analysis
# ============================================================

def plot_correlations(results, save=True):
    """Scatter plots of PMF/force vs physicochemical properties."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle('Correlations with Amino Acid Properties',
                 fontsize=14, fontweight='bold')

    props = {p: aa_properties[p] for p in results if p in aa_properties}
    pnames = list(props.keys())
    hydro = [props[p][0] for p in pnames]
    mw    = [props[p][1] for p in pnames]
    scvol = [props[p][2] for p in pnames]
    colors = [group_colors[get_group(p)] for p in pnames]

    # (a) ΔF vs Hydrophobicity
    ax = axes[0]
    dF = [results[p]['delta_F'] for p in pnames]
    ax.scatter(hydro, dF, c=colors, s=60, edgecolors='k', zorder=3)
    for i, p in enumerate(pnames):
        ax.annotate(p[:3], (hydro[i], dF[i]), fontsize=6, ha='center', va='bottom')
    ax.set_xlabel('Kyte-Doolittle Hydrophobicity')
    ax.set_ylabel('ΔF (kcal/mol)')
    ax.set_title('(a) PMF Barrier vs Hydrophobicity')

    # (b) Rupture Force vs Side-chain Volume
    ax = axes[1]
    rf = [results[p]['mean_rupF'] for p in pnames]
    ax.scatter(scvol, rf, c=colors, s=60, edgecolors='k', zorder=3)
    for i, p in enumerate(pnames):
        ax.annotate(p[:3], (scvol[i], rf[i]), fontsize=6, ha='center', va='bottom')
    ax.set_xlabel('Side-chain Volume (ų)')
    ax.set_ylabel('Mean Rupture Force (kcal/mol/Å)')
    ax.set_title('(b) Rupture Force vs Side-chain Volume')

    # (c) Kuhn Length vs MW
    ax = axes[2]
    ak = [results[p].get('a_kuhn', np.nan) for p in pnames]
    valid = [i for i in range(len(ak)) if not np.isnan(ak[i])]
    if valid:
        ax.scatter([mw[i] for i in valid], [ak[i] for i in valid],
                   c=[colors[i] for i in valid], s=60, edgecolors='k', zorder=3)
        for i in valid:
            ax.annotate(pnames[i][:3], (mw[i], ak[i]),
                        fontsize=6, ha='center', va='bottom')
    ax.set_xlabel('Molecular Weight (Da)')
    ax.set_ylabel('Kuhn Length (Å)')
    ax.set_title('(c) Kuhn Length vs Molecular Weight')

    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=c, label=g) for g, c in group_colors.items()]
    fig.legend(handles=legend_elements, loc='lower center', ncol=4, fontsize=8)
    plt.tight_layout(rect=[0, 0.08, 1, 0.93])
    if save:
        fig.savefig(os.path.join(fig_dir, 'correlations.pdf'),
                    dpi=150, bbox_inches='tight')
    plt.close(fig)


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("=" * 60)
    print("SMD Analysis — All 20 Amino Acids")
    print("=" * 60)

    results = compute_all_results()

    if not results:
        print("No results computed. Check base_dir and data files.")
    else:
        print(f"\nGenerating individual plots for {len(results)} amino acids...")
        for protein in results:
            plot_individual(protein, results[protein])
            print(f"  ✓ {protein}")

        print("\nGenerating comparative plots...")
        plot_force_extension_by_group(results)
        print("  ✓ Force-extension by group")
        plot_pmf_by_group(results)
        print("  ✓ PMF by group")
        plot_bar_comparisons(results)
        print("  ✓ Bar comparisons")
        plot_correlations(results)
        print("  ✓ Correlations")

        print_summary_table(results)

        print(f"\nAll figures saved to: {fig_dir}")
        print("Done.")
