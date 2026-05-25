"""Intensity summarization for hierarchical ion trees.

Sums base-ion intensities at specified tree node types before differential
analysis.  For example, specifying ``["frgion"]`` will sum all individual
fragment ions under each fragment-ion group into a single intensity per
replicate, while leaving MS1-isotope and precursor base ions untouched.

When the requested summarization level is above the ion-type level
(e.g. ``"mod_seq_charge"`` or ``"seq"``), leaves are split by their
ion-type ancestor so that fragment and MS1 intensities are never mixed.
"""

import re
from collections import defaultdict

import anytree
import numpy as np
import pandas as pd

from alphaquant.cluster.cluster_ions import REGEX_FRGIONS_ISOTOPES, LEVEL_NAMES

import alphaquant.config.config as aqconfig
import logging

aqconfig.setup_logging()
LOGGER = logging.getLogger(__name__)

ION_TYPE_NODES = {"frgion", "ms1_isotopes", "precursor"}

# Appended to an ion_type node name so the result is a valid base-ion name
# that the downstream tree builder can parse back into the hierarchy.
_NODE_TYPE_TO_COMPLETION_SUFFIX = {
    "frgion": "ION_SUM",
    "ms1_isotopes": "ISOTOPES_SUM",
    "precursor": "URSOR_SUM",
}

# When summarising above ion_type level we insert a synthetic path fragment
# between the higher-level node name and the ion-type suffix.
_LEVEL_TO_SYNTHETIC_INFIX = {
    "mod_seq_charge": "_",
    "mod_seq": "CHARGE_0_",
    "seq": "SUM_CHARGE_0_",
}

_ION_TYPE_TO_FULL_SUFFIX = {
    "frgion": "FRGION_SUM",
    "ms1_isotopes": "MS1ISOTOPES_SUM",
    "precursor": "PRECURSOR_SUM",
}


# ---------------------------------------------------------------------------
# Tree construction (lightweight, from ion-name strings only)
# ---------------------------------------------------------------------------

def build_tree_from_ion_names(protein_name, ion_names):
    """Build a hierarchical tree from base-ion name strings.

    Uses the same regex logic as
    :func:`~alphaquant.cluster.cluster_ions.create_hierarchical_ion_grouping`
    but operates on plain strings rather than ``DifferentialIon`` objects.
    """
    nodes = [
        anytree.Node(name, type="base", level="base")
        for name in ion_names
    ]

    for level_idx, level_patterns in enumerate(REGEX_FRGIONS_ISOTOPES):
        name2node = {}
        for pattern, node_type in level_patterns:
            for node in nodes:
                m = re.match(pattern, node.name)
                if m:
                    matching_name = m.group(1)
                    if matching_name not in name2node:
                        name2node[matching_name] = anytree.Node(
                            matching_name,
                            type=node_type,
                            level=LEVEL_NAMES[level_idx],
                        )
                    node.parent = name2node[matching_name]
        if name2node:
            nodes = list(name2node.values())

    root = anytree.Node(protein_name, type="gene", level="gene")
    for node in nodes:
        node.parent = root
    return root


# ---------------------------------------------------------------------------
# Naming helpers
# ---------------------------------------------------------------------------

def _make_summarized_name_for_ion_type_node(node):
    """Parseable summarized name for a node at the ion_type level."""
    return node.name + _NODE_TYPE_TO_COMPLETION_SUFFIX[node.type]


def _make_summarized_name_for_higher_node(parent_node, ion_type):
    """Parseable summarized name when summarising above ion_type level.

    Inserts a synthetic path fragment so that the downstream tree builder
    can still parse the resulting name into the full hierarchy.
    """
    infix = _LEVEL_TO_SYNTHETIC_INFIX[parent_node.level]
    suffix = _ION_TYPE_TO_FULL_SUFFIX[ion_type]
    return parent_node.name + infix + suffix


# ---------------------------------------------------------------------------
# Grouping logic
# ---------------------------------------------------------------------------

def compute_summarization_groups(pep2prot, ion_names, summarization_nodes):
    """Determine which base ions to group and what to name the summaries.

    Args:
        pep2prot: dict mapping ion name -> protein name.
        ion_names: iterable of all base-ion names present in either condition.
        summarization_nodes: list of node types to summarize
            (e.g. ``["frgion"]``).

    Returns:
        groups: list of ``(new_name, [leaf_ion_names], protein)`` tuples.
        remaining: set of ion names that stay as individual rows.
    """
    if not summarization_nodes:
        return [], set(ion_names)

    prot2ions = defaultdict(list)
    for ion in ion_names:
        prot = pep2prot.get(ion)
        if prot is not None:
            prot2ions[prot].append(ion)

    groups = []
    summarized_ions = set()

    for prot, ions in prot2ions.items():
        tree = build_tree_from_ion_names(prot, ions)

        for node_type in summarization_nodes:
            target_nodes = anytree.findall(
                tree, filter_=lambda n, nt=node_type: n.type == nt
            )

            for target_node in target_nodes:
                if node_type in ION_TYPE_NODES:
                    leaf_names = [
                        l.name for l in target_node.leaves if l.type == "base"
                    ]
                    if leaf_names:
                        new_name = _make_summarized_name_for_ion_type_node(
                            target_node
                        )
                        groups.append((new_name, leaf_names, prot))
                        summarized_ions.update(leaf_names)
                else:
                    # Above ion_type: split by ion type to avoid mixing
                    type_to_leaves = defaultdict(list)
                    for desc in anytree.PreOrderIter(target_node):
                        if desc.type in ION_TYPE_NODES:
                            for leaf in desc.leaves:
                                if leaf.type == "base":
                                    type_to_leaves[desc.type].append(leaf.name)
                    for ion_type, leaf_names in type_to_leaves.items():
                        new_name = _make_summarized_name_for_higher_node(
                            target_node, ion_type
                        )
                        groups.append((new_name, leaf_names, prot))
                        summarized_ions.update(leaf_names)

    remaining = set(ion_names) - summarized_ions
    return groups, remaining


# ---------------------------------------------------------------------------
# DataFrame summarization
# ---------------------------------------------------------------------------

def summarize_condition_df(df, groups, remaining_ions):
    """Apply summarization to a per-condition dataframe.

    Sums intensities in **linear** space for grouped ions, keeps remaining
    ions as-is.

    Args:
        df: DataFrame with log2 intensities, index = quant_id,
            columns = sample names.
        groups: list of ``(new_name, [leaf_ion_names], protein)`` tuples.
        remaining_ions: set of ion names to keep unchanged.

    Returns:
        Summarized DataFrame (same column layout, modified index).
    """
    parts = []

    remaining_in_df = df.index.intersection(remaining_ions)
    if len(remaining_in_df) > 0:
        parts.append(df.loc[remaining_in_df])

    for new_name, leaf_names, _prot in groups:
        present = [ion for ion in leaf_names if ion in df.index]
        if not present:
            continue
        subset = df.loc[present]
        linear = 2.0 ** subset
        summed = linear.sum(axis=0)
        all_nan = subset.isna().all(axis=0)
        with np.errstate(divide='ignore'):
            log2_summed = np.log2(summed)
        log2_summed[all_nan] = np.nan
        log2_summed.name = new_name
        parts.append(log2_summed.to_frame().T)

    if not parts:
        return pd.DataFrame(columns=df.columns)

    return pd.concat(parts)


# ---------------------------------------------------------------------------
# Ion quality filtering per group
# ---------------------------------------------------------------------------

def _filter_group_ions(leaf_names, df_c1, df_c2):
    """Select which leaf ions to include in a summarization group.

    Strategy:
      1. Keep only ions that have values in ALL replicates of BOTH conditions.
      2. If no ion qualifies, pick the single ion with the most non-NaN values
         across both conditions.

    Returns:
        Filtered list of ion names.
    """
    present_c1 = [ion for ion in leaf_names if ion in df_c1.index]
    present_c2 = [ion for ion in leaf_names if ion in df_c2.index]
    present_both = set(present_c1) & set(present_c2)

    if not present_both:
        all_present = set(present_c1) | set(present_c2)
        if not all_present:
            return []
        best_ion = max(all_present, key=lambda ion: (
            df_c1.loc[ion].notna().sum() if ion in df_c1.index else 0
        ) + (
            df_c2.loc[ion].notna().sum() if ion in df_c2.index else 0
        ))
        return [best_ion]

    n_cols_c1 = df_c1.shape[1]
    n_cols_c2 = df_c2.shape[1]

    complete = [
        ion for ion in present_both
        if df_c1.loc[ion].notna().sum() == n_cols_c1
        and df_c2.loc[ion].notna().sum() == n_cols_c2
    ]

    if complete:
        return complete

    best_ion = max(present_both, key=lambda ion:
        df_c1.loc[ion].notna().sum() + df_c2.loc[ion].notna().sum()
    )
    return [best_ion]


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def apply_summarization(df_c1, df_c2, pep2prot, summarization_nodes):
    """Summarize ion intensities at specified tree levels.

    This is the main entry point called from
    :func:`~alphaquant.diffquant.condpair_analysis.analyze_condpair`.

    Args:
        df_c1: Per-condition DataFrame for condition 1 (log2 intensities).
        df_c2: Per-condition DataFrame for condition 2 (log2 intensities).
        pep2prot: dict mapping ion name -> protein name.
        summarization_nodes: list of node types to summarize
            (e.g. ``["frgion"]``).

    Returns:
        ``(df_c1_new, df_c2_new, pep2prot_new)``
    """
    all_ions = set(df_c1.index) | set(df_c2.index)
    groups, remaining = compute_summarization_groups(
        pep2prot, all_ions, summarization_nodes
    )

    if not groups:
        return df_c1, df_c2, pep2prot

    filtered_groups = []
    for new_name, leaf_names, prot in groups:
        selected = _filter_group_ions(leaf_names, df_c1, df_c2)
        if selected:
            filtered_groups.append((new_name, selected, prot))

    if not filtered_groups:
        return df_c1, df_c2, pep2prot

    df_c1_new = summarize_condition_df(df_c1, filtered_groups, remaining)
    df_c2_new = summarize_condition_df(df_c2, filtered_groups, remaining)

    pep2prot_new = {ion: pep2prot[ion] for ion in remaining if ion in pep2prot}
    for new_name, _leaves, prot in filtered_groups:
        pep2prot_new[new_name] = prot

    LOGGER.info(
        "Summarization at %s: %d base ions -> %d entries "
        "(%d summarized groups, %d unchanged)",
        summarization_nodes,
        len(all_ions),
        len(remaining) + len(filtered_groups),
        len(filtered_groups),
        len(remaining),
    )

    return df_c1_new, df_c2_new, pep2prot_new
