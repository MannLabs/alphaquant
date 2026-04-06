"""Data-driven ICC (intraclass correlation) correction for Stouffer z-value aggregation.

Estimates a global ICC at every aggregation level of the protein tree and
annotates each node with an ``icc_correction`` attribute so that
``aggregate_node_properties`` can apply a design-effect correction.

The tree levels processed (bottom-to-top):
  - ``frgion``: correlation among fragment ions within the same precursor
  - ``ms1_isotopes``: correlation among MS1 isotopes within the same precursor
  - ``mod_seq_charge``: correlation among ion families within the same precursor
  - ``mod_seq``: correlation among charge states within the same modified peptide
  - ``seq``: correlation among modification variants within the same sequence
  - ``gene``: correlation among peptides within the same protein

A protein contributes to the null distribution only when *both* its
protein-level p-value and the p-values of its individual group nodes
exceed the significance threshold.  This two-level filter ensures the ICC
is estimated purely from technical noise, free of biological signal at
any level.

The resulting null-median ICC is applied uniformly to all proteins.
Per-protein estimation is deliberately avoided: the ICC captures a
property of the *measurement* (shared chromatographic / detector effects),
not of the individual protein, and per-protein estimates are too noisy
to be reliable given typical group sizes.

Processing proceeds bottom-to-top: after each level's ICC is estimated
and applied, trees are re-aggregated so that higher levels see the
corrected z-values from below.
"""

import anytree
import numpy as np
import matplotlib.pyplot as plt

import alphaquant.cluster.cluster_utils as aqcluster_utils

import alphaquant.config.config as aqconfig
import alphaquant.config.variables as aqvariables
import logging
aqconfig.setup_logging()
LOGGER = logging.getLogger(__name__)

# Node types that receive ICC correction, ordered bottom-to-top.
_ICC_NODE_TYPES = ("frgion", "ms1_isotopes", "mod_seq_charge", "mod_seq", "seq", "gene")

_MIN_GROUPS = 3   # minimum number of group nodes for a reliable ICC estimate
_MIN_IONS = 6     # minimum total number of child items across those groups

# Null protein selection threshold is read at runtime from
# aqvariables.ICC_NULL_PVAL_THRESHOLD (default 0.1, configurable via
# run_pipeline icc_null_pval_threshold argument).

_N_PERMUTATIONS = 3

# Gene-level subsampled estimation: each draw picks a random subset of
# null proteins and computes ICC treating each protein as a group.
_GENE_SUBSAMPLE_SIZE = 5
_GENE_N_DRAWS = 500

_NODE_TYPE_LABELS = {
    "frgion":         "Fragment ion",
    "ms1_isotopes":   "MS1 isotope",
    "mod_seq_charge": "Ion family",
    "mod_seq":        "Precursor",
    "seq":            "Mod. peptide",
    "gene":           "Peptide",
}


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def estimate_and_apply_icc_correction(protnodes, runtime_plots=False, aggregation_mode="stouffer_icc"):
    """Estimate a global ICC at every tree level, apply uniformly, and re-aggregate.

    Processes levels bottom-to-top.  After each level's ICC is estimated and
    assigned, all trees are re-aggregated so that the next (higher) level
    sees corrected z-values from below.

    Args:
        protnodes: list of protein root nodes (anytree.Node)
        runtime_plots: if True, show ICC distribution plots for all levels
        aggregation_mode: z-value combination strategy forwarded to re-aggregation
    """
    if not protnodes:
        return

    all_level_results = {}

    for node_type in _ICC_NODE_TYPES:
        if not _has_node_type(protnodes, node_type):
            continue

        LOGGER.info(f"ICC correction: estimating for node type '{node_type}'")

        if node_type == "gene":
            null_iccs, perm_iccs, icc_median = _estimate_gene_level_icc(protnodes)
        else:
            null_iccs, perm_iccs, icc_median = _estimate_null_icc_distribution(
                protnodes, node_type
            )

        n_annotated = _assign_icc_to_all_proteins(protnodes, node_type, icc_median)

        obs_med = float(np.median(null_iccs)) if null_iccs else 0.0
        perm_med = float(np.median(perm_iccs)) if perm_iccs else 0.0
        LOGGER.info(
            f"ICC correction ({node_type}): applied={icc_median:.4f} "
            f"(obs={obs_med:.4f} - perm={perm_med:.4f}), "
            f"n_null={len(null_iccs)}, n_perm={len(perm_iccs)}, "
            f"n_annotated={n_annotated}"
        )

        all_level_results[node_type] = (null_iccs, perm_iccs, icc_median)

        _re_aggregate_trees(protnodes, aggregation_mode=aggregation_mode)

    if runtime_plots and all_level_results:
        _plot_icc_distributions_all_levels(all_level_results)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _has_node_type(protnodes, node_type):
    """Return True if at least one protein tree contains a node of *node_type*."""
    for prot in protnodes:
        matches = anytree.search.findall(
            prot, filter_=lambda n: n.type == node_type
        )
        if matches:
            return True
    return False


def _estimate_null_icc_distribution(protnodes, node_type):
    """Compute observed and permutation-based ICC for each null protein.

    A protein qualifies as null only when *both*:
      - its protein-level p_val > aqvariables.ICC_NULL_PVAL_THRESHOLD, *and*
      - only group nodes whose own p_val > threshold are used.

    Returns:
        (null_iccs, permutation_iccs, icc_median).
        If no null proteins qualify, returns ([], [], 0.0).
    """
    null_iccs = []
    null_group_zvals_list = []

    pval_threshold = aqvariables.get_icc_null_pval_threshold(node_type)
    for prot in protnodes:
        if not hasattr(prot, "p_val"):
            continue
        if prot.p_val <= pval_threshold:
            continue

        group_zvals = _collect_group_zvals(prot, node_type, node_p_val_threshold=pval_threshold)
        if not group_zvals:
            continue

        icc = _compute_icc_from_groups(group_zvals)
        if icc is not None:
            null_iccs.append(icc)
            null_group_zvals_list.append(group_zvals)

    if len(null_iccs) == 0:
        LOGGER.warning(
            f"ICC correction ({node_type}): no null proteins qualified; "
            "falling back to icc_median=0.0"
        )
        return [], [], 0.0

    permutation_iccs = _compute_permutation_null(null_group_zvals_list)

    icc_median = _null_normalized_icc(null_iccs, permutation_iccs)
    return null_iccs, permutation_iccs, icc_median


def _estimate_gene_level_icc(protnodes, seed=42):
    """Estimate ICC at the protein (gene) level using subsampled estimation.

    Since each protein has exactly one gene node, per-protein ICC is
    undefined.  Instead, we subsample groups of null proteins and compute
    ICC treating each protein as a group with its peptide z-values as items.

    The permutation null shuffles all z-values across proteins (destroying
    the protein-level grouping) and repeats the same subsampled estimation.

    Returns:
        (null_iccs, permutation_iccs, icc_median).
    """
    protein_groups = []
    gene_threshold = aqvariables.get_icc_null_pval_threshold("gene")
    for prot in protnodes:
        if not hasattr(prot, "p_val") or prot.p_val <= gene_threshold:
            continue
        children_zvals = [c.z_val for c in prot.children if hasattr(c, "z_val")]
        if len(children_zvals) >= 2:
            protein_groups.append(np.array(children_zvals))

    if len(protein_groups) < _GENE_SUBSAMPLE_SIZE:
        LOGGER.warning(
            f"ICC correction (gene): only {len(protein_groups)} null proteins "
            f"with ≥2 peptides (need {_GENE_SUBSAMPLE_SIZE}); falling back to 0.0"
        )
        return [], [], 0.0

    rng = np.random.RandomState(seed)
    null_iccs = _subsampled_icc(protein_groups, rng)

    sizes = [len(g) for g in protein_groups]
    pooled = np.concatenate(protein_groups)
    rng_perm = np.random.RandomState(seed + 1)
    rng_perm.shuffle(pooled)
    shuffled_groups = []
    offset = 0
    for s in sizes:
        shuffled_groups.append(pooled[offset:offset + s].copy())
        offset += s
    perm_iccs = _subsampled_icc(shuffled_groups, np.random.RandomState(seed))

    if not null_iccs:
        LOGGER.warning("ICC correction (gene): no valid subsamples; falling back to 0.0")
        return [], [], 0.0

    icc_median = _null_normalized_icc(null_iccs, perm_iccs)
    return null_iccs, perm_iccs, icc_median


def _subsampled_icc(protein_groups, rng,
                    subset_size=_GENE_SUBSAMPLE_SIZE,
                    n_draws=_GENE_N_DRAWS):
    """Compute ICC on random subsets of protein groups.

    Each draw selects *subset_size* proteins, treats each as a group
    (with its peptide z-values as items), and computes one ICC value.
    """
    iccs = []
    n_proteins = len(protein_groups)
    for _ in range(n_draws):
        indices = rng.choice(n_proteins, size=subset_size, replace=False)
        group_zvals = [protein_groups[i] for i in indices]
        icc = _compute_icc_from_groups(group_zvals)
        if icc is not None:
            iccs.append(icc)
    return iccs


def _compute_permutation_null(group_zvals_list, n_permutations=_N_PERMUTATIONS, seed=42):
    """Compute ICC on shuffled data to establish a permutation baseline.

    For each protein's grouped z-values, all values are pooled and randomly
    reassigned to groups of the original sizes.  This destroys any real
    intra-group correlation, so the resulting ICC distribution reflects pure
    sampling noise.
    """
    rng = np.random.RandomState(seed)
    perm_iccs = []

    for group_zvals in group_zvals_list:
        sizes = [len(g) for g in group_zvals]
        pooled = np.concatenate(group_zvals)

        for _ in range(n_permutations):
            rng.shuffle(pooled)
            shuffled_groups = []
            offset = 0
            for s in sizes:
                shuffled_groups.append(pooled[offset:offset + s].copy())
                offset += s

            icc = _compute_icc_from_groups(shuffled_groups)
            if icc is not None:
                perm_iccs.append(icc)

    return perm_iccs


def _null_normalized_icc(null_iccs, perm_iccs):
    """Return the null-normalized ICC: median(observed) - median(shuffled), clamped to ≥0."""
    obs_med = float(np.median(null_iccs))
    perm_med = float(np.median(perm_iccs)) if len(perm_iccs) > 0 else 0.0
    return max(0.0, obs_med - perm_med)


def _assign_icc_to_all_proteins(protnodes, node_type, icc_median):
    """Assign the global null-median ICC uniformly to every node of *node_type*.

    Returns:
        Number of nodes annotated.
    """
    n_annotated = 0
    for prot in protnodes:
        target_nodes = anytree.search.findall(
            prot, filter_=lambda n: n.type == node_type
        )
        for node in target_nodes:
            node.icc_correction = icc_median
            n_annotated += 1
    return n_annotated


def _collect_group_zvals(protein_node, node_type, node_p_val_threshold=None):
    """Extract grouped z-values from a protein tree.

    For each node of *node_type*, collects the z-values of its immediate
    children.  At the ion-type level (frgion, ms1_isotopes) the children
    are base ions; at higher levels they are aggregated child nodes.

    Returns:
        list[np.ndarray]: One array of child z-values per qualifying group
            node, or an empty list if insufficient data.
    """
    group_nodes = anytree.search.findall(
        protein_node, filter_=lambda n: n.type == node_type
    )
    if not group_nodes:
        return []

    if node_p_val_threshold is not None:
        group_nodes = [
            n for n in group_nodes
            if hasattr(n, "p_val") and n.p_val > node_p_val_threshold
        ]

    group_zvals = []
    for gnode in group_nodes:
        children = [
            c for c in gnode.children
            if hasattr(c, "z_val")
        ]
        if children:
            group_zvals.append(
                np.array([c.z_val for c in children])
            )

    n_groups = len(group_zvals)
    if n_groups < _MIN_GROUPS:
        return []

    n_items = sum(len(g) for g in group_zvals)
    if n_items < _MIN_IONS:
        return []

    return group_zvals


def _compute_icc_from_groups(group_zvals):
    """Compute one-way random-effects ICC from pre-collected grouped z-values.

    Returns:
        ICC (float) or None if degenerate.
    """
    n_groups = len(group_zvals)
    if n_groups < _MIN_GROUPS:
        return None

    all_vals = np.concatenate(group_zvals)
    n_items = len(all_vals)
    if n_items < _MIN_IONS:
        return None

    grand_mean = all_vals.mean()

    ss_within = 0.0
    group_sizes = np.empty(n_groups, dtype=int)
    group_means = np.empty(n_groups)
    for i, gv in enumerate(group_zvals):
        gm = gv.mean()
        group_means[i] = gm
        group_sizes[i] = len(gv)
        ss_within += np.sum((gv - gm) ** 2)

    df_within = n_items - n_groups
    if df_within < 1:
        return None
    sigma2_residual = ss_within / df_within

    expanded_group_means = np.repeat(group_means, group_sizes)
    ss_between = np.sum((expanded_group_means - grand_mean) ** 2)
    df_between = n_groups - 1
    ms_between = ss_between / df_between

    n0 = (n_items - np.sum(group_sizes ** 2) / n_items) / df_between
    if n0 < 1e-15:
        return None

    sigma2_group = (ms_between - sigma2_residual) / n0
    sigma2_group = max(sigma2_group, 0.0)

    total = sigma2_group + sigma2_residual
    if total < 1e-15:
        return None

    return sigma2_group / total


def _compute_icc_from_tree(protein_node, node_type, node_p_val_threshold=None):
    """Convenience wrapper: extract groups from tree, then compute ICC."""
    group_zvals = _collect_group_zvals(protein_node, node_type, node_p_val_threshold)
    if not group_zvals:
        return None
    return _compute_icc_from_groups(group_zvals)


def _re_aggregate_trees(protnodes, aggregation_mode="stouffer_icc"):
    """Re-aggregate all protein trees bottom-to-top after ICC annotation.

    Walks every tree from leaves upward, re-computing z-values and p-values
    at each level so that the newly annotated ``icc_correction`` attributes
    take effect.
    """
    for prot in protnodes:
        for level_nodes in aqcluster_utils.iterate_through_tree_levels_bottom_to_top(prot):
            node_types = list(set(node.type for node in level_nodes))
            if node_types == ["base"]:
                continue
            for node in level_nodes:
                if node.children:
                    aqcluster_utils.aggregate_node_properties(
                        node, only_use_mainclust=True, peptide_outlier_filtering=False,
                        aggregation_mode=aggregation_mode
                    )


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _plot_icc_distributions_all_levels(all_level_results):
    """Multi-row figure with histogram + CDF per level.

    Each row corresponds to one tree level.  Left panel shows histograms of
    observed-null and permutation-null ICC distributions; right panel shows
    cumulative distribution functions.  A vertical line marks the applied
    median ICC.

    Args:
        all_level_results: dict mapping node_type →
            (null_iccs, perm_iccs, icc_median).
    """
    ordered_types = [nt for nt in _ICC_NODE_TYPES if nt in all_level_results]
    display_order = list(reversed(ordered_types))
    n_levels = len(display_order)
    if n_levels == 0:
        return

    fig, axes = plt.subplots(n_levels, 2, figsize=(7.5, 1.8 * n_levels),
                             squeeze=False)
    fig.subplots_adjust(hspace=0.65, wspace=0.35, top=0.96, bottom=0.06,
                        left=0.12, right=0.98)

    color_obs = "#4C72B0"
    color_perm = "#DD8452"
    bins = np.linspace(0, 1, 40)

    for row, node_type in enumerate(display_order):
        null_iccs, perm_iccs, icc_median = all_level_results[node_type]
        ax_hist, ax_cdf = axes[row]
        label = _NODE_TYPE_LABELS.get(node_type, node_type)

        legend_handles = []

        if len(perm_iccs) > 0:
            perm_med = float(np.median(perm_iccs))
            ax_hist.hist(perm_iccs, bins=bins, alpha=0.45, color=color_perm,
                         edgecolor="none")
            s = np.sort(perm_iccs)
            ax_cdf.plot(s, np.arange(1, len(s) + 1) / len(s), color=color_perm)
            legend_handles.append(
                plt.matplotlib.patches.Patch(
                    facecolor=color_perm, alpha=0.6,
                    label=f"shuffled (med={perm_med:.2f})")
            )

        if len(null_iccs) > 0:
            obs_med = float(np.median(null_iccs))
            ax_hist.hist(null_iccs, bins=bins, alpha=0.45, color=color_obs,
                         edgecolor="none")
            s = np.sort(null_iccs)
            ax_cdf.plot(s, np.arange(1, len(s) + 1) / len(s), color=color_obs)
            legend_handles.append(
                plt.matplotlib.patches.Patch(
                    facecolor=color_obs, alpha=0.6,
                    label=f"observed (med={obs_med:.2f})")
            )

        ax_cdf.axvline(0.5, color="gray", linestyle=":", linewidth=0.6)

        for ax in (ax_hist, ax_cdf):
            ax.set_xlim(-0.05, 1.05)
            ax.tick_params(axis='both', which='both', length=3, pad=2)

        ax_hist.set_ylabel("count")
        ax_cdf.set_ylabel("cum. fraction")
        ax_cdf.set_ylim(-0.02, 1.02)
        ax_cdf.set_yticks([0.0, 0.5, 1.0])

        if row == n_levels - 1:
            ax_hist.set_xlabel("ICC")
            ax_cdf.set_xlabel("ICC")
        else:
            ax_hist.set_xticklabels([])
            ax_cdf.set_xticklabels([])

        ax_hist.legend(handles=legend_handles, loc='upper right',
                       fontsize=6, frameon=False, handlelength=1,
                       handletextpad=0.4, borderpad=0.2)

        pos_l = ax_hist.get_position()
        pos_r = ax_cdf.get_position()
        x_mid = (pos_l.x0 + pos_r.x1) / 2
        fig.text(x_mid, pos_l.y1 + 0.008, label,
                 ha='center', va='bottom', fontsize=10, fontweight='bold')

    fig.suptitle("ICC distributions across tree levels", fontsize=12,
                 fontweight='bold', y=0.995)
    plt.show()
