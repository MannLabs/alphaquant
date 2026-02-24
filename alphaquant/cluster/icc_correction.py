"""Data-driven ICC (intraclass correlation) correction for Stouffer z-value aggregation.

Estimates per-protein ICC from the tree's own z-values and annotates each
frgion / ms1_isotopes node with an ``icc_correction`` attribute so that
``aggregate_node_properties`` can apply a tailored design-effect correction.

The estimation follows a two-stage approach:
  1. **Null distribution** -- ICC is computed for "null" proteins (gene-level
     p > 0.1) to obtain a robust median that serves as both *cap* and
     *fallback*.
  2. **Per-protein ICC** -- Each protein gets its own ICC, capped at the
     null median.  Proteins with insufficient data receive the null median
     as a fallback.

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


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def estimate_and_apply_icc_correction(protnodes, runtime_plots=False, aggregation_mode="stouffer_icc"):
    """Estimate per-protein ICC and annotate tree nodes, then re-aggregate.

    For each relevant node type (``frgion``, ``ms1_isotopes``):
      - derive a robust median ICC from null proteins
      - compute per-protein ICC, capped at the null median
      - annotate each node with ``icc_correction`` and ``icc_is_fallback``

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

        null_iccs, icc_median = _estimate_null_icc_distribution(
            protnodes, node_type
        )

        all_iccs = []
        for prot in protnodes:
            protein_icc = _compute_and_assign_icc_for_protein(
                prot, node_type, icc_median
            )
            if protein_icc is not None:
                all_iccs.append(protein_icc)

        LOGGER.info(
            f"ICC correction ({node_type}): null median={icc_median:.4f}, "
            f"n_null={len(null_iccs)}, n_all={len(all_iccs)}"
        )

        if runtime_plots and (len(null_iccs) > 0 or len(all_iccs) > 0):
            _plot_icc_distributions(null_iccs, all_iccs, icc_median, node_type)

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
    """Compute ICC for each null protein and return the list + median.

    Null proteins are those whose gene-level p_val > _P_VAL_THRESHOLD.

    Returns:
        (null_iccs, icc_median): list of float ICC values and their median.
            If no null proteins qualify, returns ([], 0.0).
    """
    null_iccs = []

    for prot in protnodes:
        if not hasattr(prot, "p_val"):
            continue
        if prot.p_val <= _P_VAL_THRESHOLD:
            continue

        icc = _compute_icc_from_tree(prot, node_type)
        if icc is not None:
            null_iccs.append(icc)

    if len(null_iccs) == 0:
        LOGGER.warning(
            f"ICC correction ({node_type}): no null proteins qualified; "
            "falling back to icc_median=0.0"
        )
        return [], 0.0

    icc_median = float(np.median(null_iccs))
    return null_iccs, icc_median


def _compute_and_assign_icc_for_protein(protein_node, node_type, icc_median):
    """Compute ICC for one protein and assign it to the relevant nodes.

    Returns:
        The ICC value used (float), or None if no nodes of node_type exist.
    """
    target_nodes = anytree.search.findall(
        protein_node, filter_=lambda n: n.type == node_type
    )
    if not target_nodes:
        return None

    icc = _compute_icc_from_tree(protein_node, node_type)

    if icc is not None:
        final_icc = min(icc, icc_median) if icc_median > 0 else icc
        is_fallback = False
    else:
        final_icc = icc_median
        is_fallback = True

    for node in target_nodes:
        node.icc_correction = final_icc
        node.icc_is_fallback = is_fallback

    return final_icc


def _compute_icc_from_tree(protein_node, node_type):
    """Compute ICC for a single protein using variance components.

    Groups are the nodes of ``node_type`` (e.g. frgion); observations are
    the z-values of their base-ion children.

    Returns:
        ICC (float) or None if insufficient data.
    """
    group_nodes = anytree.search.findall(
        protein_node, filter_=lambda n: n.type == node_type
    )
    if not group_nodes:
        return None

    # Collect z-values grouped by parent node
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
        return None

    all_vals = np.concatenate(group_zvals)
    n_ions = len(all_vals)
    if n_ions < _MIN_IONS:
        return None

    grand_mean = all_vals.mean()

    # Within-group variance (sigma2_residual)
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

    # Between-group variance (sigma2_precursor)
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


def _plot_icc_distributions(null_iccs, all_iccs, icc_median, node_type):
    """Show overlaid histograms of null vs. all-protein ICC values."""
    fig, ax = plt.subplots(figsize=(5, 3))

    bins = np.linspace(0, 1, 40)

    if len(all_iccs) > 0:
        ax.hist(
            all_iccs, bins=bins, alpha=0.5, color="tab:orange",
            edgecolor="white", label="all proteins"
        )
    if len(null_iccs) > 0:
        ax.hist(
            null_iccs, bins=bins, alpha=0.6, color="tab:blue",
            edgecolor="white", label="null proteins"
        )

    ax.axvline(
        icc_median, color="red", linestyle="--", linewidth=1.5,
        label=f"median (cap) = {icc_median:.3f}"
    )

    ax.set_xlabel("ICC")
    ax.set_ylabel("Count")
    ax.set_title(f"ICC distribution ({node_type})")
    ax.legend(fontsize=8)
    fig.tight_layout()
    plt.show()
