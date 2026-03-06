"""Data-driven ICC (intraclass correlation) correction for Stouffer z-value aggregation.

Estimates a global ICC from null proteins and annotates each frgion /
ms1_isotopes node with an ``icc_correction`` attribute so that
``aggregate_node_properties`` can apply a design-effect correction.

A protein contributes to the null distribution only when *both* its
gene-level p-value and the p-values of its individual fragment-type nodes
exceed the significance threshold.  This two-level filter ensures the ICC
is estimated purely from technical noise, free of biological signal at
any level.

The resulting null-median ICC is applied uniformly to all proteins.
Per-protein estimation is deliberately avoided: the ICC captures a
property of the *measurement* (shared chromatographic / detector effects),
not of the individual protein, and per-protein estimates are too noisy
to be reliable given typical group sizes.

After annotation the trees are re-aggregated bottom-to-top so that the
corrected z-values propagate upward.
"""

import anytree
import numpy as np
import matplotlib.pyplot as plt

import alphaquant.cluster.cluster_utils as aqcluster_utils

import alphaquant.config.config as aqconfig
import logging
aqconfig.setup_logging()
LOGGER = logging.getLogger(__name__)

# Node types that receive ICC correction
_ICC_NODE_TYPES = ("frgion", "ms1_isotopes")

# Minimum data requirements for a reliable ICC estimate
_MIN_GROUPS = 3   # minimum number of group nodes (e.g. precursors with base ions)
_MIN_IONS = 6     # minimum total number of base ions across those groups

# Null protein selection: proteins whose gene-level p-value exceeds this
# threshold are considered part of the technical-noise background.
_P_VAL_THRESHOLD = 0.1

# Permutation null: number of shuffles per protein for the permutation baseline
_N_PERMUTATIONS = 3


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def estimate_and_apply_icc_correction(protnodes, runtime_plots=False, aggregation_mode="stouffer_icc"):
    """Estimate a global ICC from null proteins, apply it uniformly, and re-aggregate.

    For each relevant node type (``frgion``, ``ms1_isotopes``):
      - derive a robust median ICC from null proteins (filtered at both the
        gene level and the individual node level)
      - assign this median uniformly to every node of that type

    Finally, re-aggregate all trees bottom-to-top so the corrected z-values
    propagate upward.

    Args:
        protnodes: list of protein root nodes (anytree.Node)
        runtime_plots: if True, show ICC distribution histograms
        aggregation_mode: z-value combination strategy forwarded to re-aggregation
    """
    if not protnodes:
        return

    for node_type in _ICC_NODE_TYPES:
        if not _has_node_type(protnodes, node_type):
            continue

        LOGGER.info(f"ICC correction: estimating for node type '{node_type}'")

        null_iccs, perm_iccs, icc_median = _estimate_null_icc_distribution(
            protnodes, node_type
        )

        n_annotated = _assign_icc_to_all_proteins(protnodes, node_type, icc_median)

        LOGGER.info(
            f"ICC correction ({node_type}): null median={icc_median:.4f}, "
            f"n_null={len(null_iccs)}, n_perm={len(perm_iccs)}, "
            f"n_annotated={n_annotated}"
        )

        if runtime_plots and len(null_iccs) > 0:
            _plot_icc_distributions(null_iccs, perm_iccs, icc_median, node_type)

    _re_aggregate_trees(protnodes, aggregation_mode=aggregation_mode)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _has_node_type(protnodes, node_type):
    """Return True if at least one protein tree contains a node of *node_type*.

    Args:
        protnodes (list[anytree.Node]): Protein root nodes to search.
        node_type (str): The ``type`` attribute to look for (e.g. ``"frgion"``).

    Returns:
        bool
    """
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
      - its gene-level p_val > _P_VAL_THRESHOLD, *and*
      - only fragment-type nodes whose own p_val > _P_VAL_THRESHOLD are
        used for the ICC computation.

    The permutation null shuffles z-values across groups (destroying the
    group structure) to show what ICC looks like under true independence.

    Returns:
        (null_iccs, permutation_iccs, icc_median).
        If no null proteins qualify, returns ([], [], 0.0).
    """
    null_iccs = []
    null_group_zvals_list = []

    for prot in protnodes:
        if not hasattr(prot, "p_val"):
            continue
        if prot.p_val <= _P_VAL_THRESHOLD:
            continue

        group_zvals = _collect_group_zvals(prot, node_type, node_p_val_threshold=_P_VAL_THRESHOLD)
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

    icc_median = float(np.median(null_iccs))
    return null_iccs, permutation_iccs, icc_median


def _compute_permutation_null(group_zvals_list, n_permutations=_N_PERMUTATIONS, seed=42):
    """Compute ICC on shuffled data to establish a permutation baseline.

    For each protein's grouped z-values, all values are pooled and randomly
    reassigned to groups of the original sizes.  This destroys any real
    intra-group correlation, so the resulting ICC distribution reflects pure
    sampling noise.

    Args:
        group_zvals_list: list of (list of np.ndarray) — one entry per
            qualifying null protein.
        n_permutations: number of independent shuffles per protein.
        seed: RNG seed for reproducibility.

    Returns:
        list of float ICC values from the permuted data.
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

    Returns:
        list[np.ndarray]: One array of base-ion z-values per qualifying group
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
        base_children = [
            c for c in gnode.children
            if c.type == "base" and hasattr(c, "z_val")
        ]
        if base_children:
            group_zvals.append(
                np.array([c.z_val for c in base_children])
            )

    n_groups = len(group_zvals)
    if n_groups < _MIN_GROUPS:
        return []

    n_ions = sum(len(g) for g in group_zvals)
    if n_ions < _MIN_IONS:
        return []

    return group_zvals


def _compute_icc_from_groups(group_zvals):
    """Compute one-way random-effects ICC from pre-collected grouped z-values.

    Args:
        group_zvals: list of np.ndarray, one per group.

    Returns:
        ICC (float) or None if degenerate.
    """
    n_groups = len(group_zvals)
    if n_groups < _MIN_GROUPS:
        return None

    all_vals = np.concatenate(group_zvals)
    n_ions = len(all_vals)
    if n_ions < _MIN_IONS:
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

    df_within = n_ions - n_groups
    if df_within < 1:
        return None
    sigma2_residual = ss_within / df_within

    expanded_group_means = np.repeat(group_means, group_sizes)
    ss_between = np.sum((expanded_group_means - grand_mean) ** 2)
    df_between = n_groups - 1
    ms_between = ss_between / df_between

    n0 = (n_ions - np.sum(group_sizes ** 2) / n_ions) / df_between
    if n0 < 1e-15:
        return None

    sigma2_precursor = (ms_between - sigma2_residual) / n0
    sigma2_precursor = max(sigma2_precursor, 0.0)

    total = sigma2_precursor + sigma2_residual
    if total < 1e-15:
        return None

    return sigma2_precursor / total


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

    Args:
        protnodes (list[anytree.Node]): Protein root nodes.
        aggregation_mode (str): Z-value combination strategy forwarded to
            ``aggregate_node_properties`` (default ``"stouffer_icc"``).
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


def _plot_icc_distributions(null_iccs, perm_iccs, icc_median, node_type):
    """Two-panel figure: histograms (top) and cumulative distributions (bottom).

    Compares the observed null-protein ICC distribution against the
    permutation baseline.
    """
    fig, (ax_hist, ax_cdf) = plt.subplots(2, 1, figsize=(5, 5), sharex=True)
    bins = np.linspace(0, 1, 40)

    # -- Histograms --
    if len(perm_iccs) > 0:
        ax_hist.hist(
            perm_iccs, bins=bins, alpha=0.45, color="tab:gray",
            edgecolor="white", label="permutation null", density=True,
        )
    ax_hist.hist(
        null_iccs, bins=bins, alpha=0.6, color="tab:blue",
        edgecolor="white", label="observed null", density=True,
    )
    ax_hist.axvline(
        icc_median, color="red", linestyle="--", linewidth=1.5,
        label=f"applied median = {icc_median:.3f}",
    )
    ax_hist.set_ylabel("Density")
    ax_hist.set_title(f"ICC distribution — {node_type}")
    ax_hist.legend(fontsize=8)

    # -- Cumulative distributions --
    for vals, color, label in [
        (perm_iccs, "tab:gray", "permutation null"),
        (null_iccs, "tab:blue", "observed null"),
    ]:
        if len(vals) == 0:
            continue
        sorted_vals = np.sort(vals)
        cdf = np.arange(1, len(sorted_vals) + 1) / len(sorted_vals)
        ax_cdf.step(sorted_vals, cdf, where="post", color=color, label=label)

    ax_cdf.axvline(
        icc_median, color="red", linestyle="--", linewidth=1.5,
        label=f"applied median = {icc_median:.3f}",
    )
    ax_cdf.set_xlabel("ICC")
    ax_cdf.set_ylabel("Cumulative fraction")
    ax_cdf.legend(fontsize=8)

    fig.tight_layout()
    plt.show()
