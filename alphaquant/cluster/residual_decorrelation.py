"""Residual decorrelation: remove correlated siblings before z-score aggregation.

Overview
--------
When AlphaQuant aggregates child nodes (e.g. peptides → protein, fragments →
peptide) using Stouffer's method, inflated sibling correlations bias the
combined z-score upward.  This module identifies and prunes the most
correlated children at every level of the ion tree so that the surviving
set's pairwise correlation distribution matches a condition-shuffled null.

Algorithm
---------
1. **Residual computation** (``attach_lm_residuals``):
   For every base ion the within-condition mean intensity is subtracted,
   yielding condition-mean residuals.  Residuals for higher-level nodes are
   the row-wise mean of their children's residuals, propagated bottom-up.
   These residuals capture shared technical variation independently of the
   fold-change signal.

2. **Per-parent precomputation** (``_build_parent``, ``ParentPrecompute``):
   For each parent node at a given level the children's residual vectors are
   stacked into a matrix and a Pearson correlation matrix ``C`` is computed.
   A greedy removal order is then computed once: at each step the child with
   the highest mean pairwise correlation to the surviving set is removed.
   The maximum pairwise correlation after each removal is stored as
   ``max_r_trajectory``, making it cheap to replay any cutoff later via
   ``survivors_at``.

3. **Null distribution** (``_cross_parent_shuffle_null``):
   Rows are permuted across parents (cross-parent shuffle) to produce a
   baseline that represents what sibling correlations look like when children
   are exchanged between unrelated proteins.

4. **Level sweep** (``run_level_sweep``):
   A grid of correlation cutoffs (default 1.0 → 0.1) is scanned.  For each
   cutoff the surviving correlation values across all parents are collected
   and compared to the null via a one-sided excess-CDF distance ``D``
   (``_excess_cdf_distance``): the maximum over all ``r`` of
   ``F_null(r) − F_corrected(r)``, i.e. how much the corrected distribution
   still exceeds the null.  The lowest cutoff with ``D ≤ tolerance`` is
   chosen.  If none qualifies the tightest cutoff is used regardless.

5. **Application** (``apply_residual_decorrelation``):
   Main entry point.  Runs steps 1–4 for every ``LEVEL_PAIRS`` pair, marks
   pruned children with ``node.exclude_residual_decorrelation = True``, then
   re-aggregates node statistics bottom-up with the decorrelation-aware
   aggregation mode.

Structure
---------
Data classes
    ParentPrecompute   – precomputed correlation matrix + removal trajectory
    LevelSweepResult   – outcome of one level sweep (cutoffs, distances, traces)

Internal helpers
    _node_matches_level          – type-string → node matching
    _build_parent                – build ParentPrecompute for one parent node
    _pair_rs_from_C              – extract upper-triangle r values given survivors
    _cross_parent_shuffle_null   – permutation-based null distribution
    _aggregate_pair_rs           – collect r values across parents at a cutoff
    _excess_cdf_distance         – one-sided CDF distance metric

Public API
    attach_lm_residuals          – attach within-condition residuals to tree nodes
    run_level_sweep              – sweep cutoffs for one (parent, child) level pair
    apply_residual_decorrelation – orchestrate the full pipeline (main entry point)
    plot_level_sweep_cdfs        – CDF comparison plot for one level result
    plot_level_sweep_diagnostics – full diagnostic figure (CDF + sweep trace)
"""
from __future__ import annotations

from dataclasses import dataclass, field
import warnings

import anytree
from anytree import PreOrderIter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import alphaquant.cluster.cluster_utils as aqcluster_utils
import alphaquant.config.variables as aqvariables
import alphaquant.config.config as aqconfig
import logging

aqconfig.setup_logging()
LOGGER = logging.getLogger(__name__)


LEVEL_PAIRS = (
    ("gene", "seq"),
    ("seq", "mod_seq"),
    ("mod_seq", "mod_seq_charge"),
    ("mod_seq_charge", "ion_type"),
    ("frgion", "base"),
    ("ms1_isotopes", "base"),
)

DEFAULT_CUTOFF_GRID = tuple(round(1.0 - 0.1 * k, 2) for k in range(10))
DEFAULT_TOLERANCE = 0.10
DEFAULT_MIN_KEEP = 1


@dataclass
class ParentPrecompute:
    """Precomputed correlation structure for one parent node.

    Built once per parent by ``_build_parent``; the removal order and
    max-r trajectory are stored so that ``survivors_at`` can replay any
    cutoff in O(n) without recomputing correlations.

    Attributes
    ----------
    parent_node:      the tree node this precompute belongs to
    child_nodes:      ordered tuple of children whose residuals were used
    C:                (n × n) Pearson correlation matrix; diagonal is NaN
    remove_order:     greedy removal sequence (indices into child_nodes)
    max_r_trajectory: max pairwise r after removing k children (length n)
    """

    parent_node: anytree.Node
    child_nodes: tuple[anytree.Node, ...]
    C: np.ndarray
    remove_order: np.ndarray
    max_r_trajectory: np.ndarray

    def survivors_at(self, cutoff: float, min_keep: int) -> np.ndarray:
        """Return a boolean mask of children that survive at ``cutoff``.

        Replays the greedy removal order until ``max_r_trajectory`` drops
        below ``cutoff``, always retaining at least ``min_keep`` children.
        """
        n = self.C.shape[0]
        if n == 0:
            return np.zeros(0, dtype=bool)
        # at most n-min_keep children may be removed
        k_max = max(0, n - min_keep)
        # scan the trajectory: stop at the first step where max_r is already below cutoff
        k = 0
        while k <= k_max and self.max_r_trajectory[k] > cutoff:
            k += 1
        # clamp in case the loop overshot (shouldn't happen with monotone trajectory)
        if k > k_max:
            k = k_max
        # replay: mark the first k entries of the greedy removal order as dead
        alive = np.ones(n, dtype=bool)
        if k > 0:
            alive[self.remove_order[:k]] = False
        return alive


@dataclass
class LevelSweepResult:
    """Outcome of a full cutoff sweep for one (parent_level, child_level) pair.

    Attributes
    ----------
    level:             (parent_level, child_level) string pair
    cutoff:            chosen correlation cutoff
    d_before/d_after:  excess CDF distance before and after pruning
    n_parents:         total parents examined at this level
    parents_touched:   parents where at least one child was dropped
    children_dropped:  total children marked for exclusion
    grid_trace:        list of (cutoff, distance, dropped, touched) per grid step
    unmodified_sorted: sorted pairwise r values before pruning (for plotting)
    corrected_sorted:  sorted pairwise r values after pruning (for plotting)
    null_sorted:       sorted null distribution r values (for plotting)
    """

    level: tuple[str, str]
    cutoff: float
    d_before: float
    d_after: float
    n_parents: int
    parents_touched: int
    children_dropped: int
    grid_trace: list[tuple[float, float, int, int]] = field(default_factory=list)
    unmodified_sorted: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.float64))
    corrected_sorted: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.float64))
    null_sorted: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.float64))


def _node_matches_level(node, level: str) -> bool:
    if level == "ion_type":
        return node.type in {"frgion", "ms1_isotopes", "precursor"}
    return node.type == level


def _build_parent(parent_node, child_nodes, mat):
    """Build a ParentPrecompute by greedily computing the removal order.

    At each step the child with the highest mean pairwise correlation to the
    remaining set is removed.  The maximum pairwise r after each removal is
    recorded in ``max_r_trajectory`` and monotonically enforced (non-increasing)
    so that ``survivors_at`` can use a simple threshold scan.
    """
    if mat.shape[0] < 2:
        return None

    # compute full Pearson correlation matrix across children's residual vectors
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", r"invalid value encountered")
        C = np.corrcoef(mat)
    # constant rows produce NaN correlations; replace with 0 so they don't distort means
    if not np.all(np.isfinite(C)):
        C = np.nan_to_num(C, nan=0.0, posinf=1.0, neginf=-1.0)
    C = C.copy()
    # NaN on diagonal so nanmean excludes self-correlation when averaging rows
    np.fill_diagonal(C, np.nan)

    n = C.shape[0]
    alive = np.ones(n, dtype=bool)
    remove_order = []
    max_r = []

    # record the max pairwise r before any removal (step 0 of the trajectory)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", r"All-NaN slice encountered")
        init_max = float(np.nanmax(C)) if n >= 2 else -np.inf
    if not np.isfinite(init_max):
        init_max = -np.inf
    max_r.append(init_max)

    while alive.sum() > 1:
        # extract the submatrix of currently surviving children
        sub_idx = np.where(alive)[0]
        cc = C[np.ix_(sub_idx, sub_idx)]
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", r"Mean of empty slice")
            warnings.filterwarnings("ignore", r"All-NaN slice encountered")
            mean_r = np.nanmean(cc, axis=1)
        # the "worst" child is the one most correlated on average with its siblings;
        # replace NaN with -inf so it is never chosen as worst (avoids argmax crash)
        worst_local = int(np.nanargmax(np.where(np.isnan(mean_r), -np.inf, mean_r)))
        alive[sub_idx[worst_local]] = False
        remove_order.append(sub_idx[worst_local])
        # record the new maximum pairwise r among survivors after this removal
        if alive.sum() >= 2:
            rem = np.where(alive)[0]
            cc2 = C[np.ix_(rem, rem)]
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", r"All-NaN slice encountered")
                m = float(np.nanmax(cc2))
            if not np.isfinite(m):
                m = -np.inf
        else:
            m = -np.inf
        max_r.append(m)

    traj = np.asarray(max_r, dtype=np.float64)
    # enforce monotone non-increasing: numerical noise can cause slight upward bumps
    # which would break the threshold scan in survivors_at
    for i in range(1, traj.size):
        if traj[i] > traj[i - 1]:
            traj[i] = traj[i - 1]

    return ParentPrecompute(
        parent_node=parent_node,
        child_nodes=tuple(child_nodes),
        C=C.astype(np.float64),
        remove_order=np.asarray(remove_order, dtype=np.int64),
        max_r_trajectory=traj,
    )


def _pair_rs_from_C(C: np.ndarray, survivors: np.ndarray) -> np.ndarray:
    if survivors.sum() < 2:
        return np.empty(0, dtype=np.float64)
    idx = np.where(survivors)[0]
    sub = C[np.ix_(idx, idx)]
    # k=1 skips the diagonal, giving only unique off-diagonal pairs
    iu = np.triu_indices(sub.shape[0], k=1)
    vals = sub[iu]
    # drop NaN entries that arise from originally constant residual vectors
    vals = vals[~np.isnan(vals)]
    return vals.astype(np.float64, copy=False)


def _cross_parent_shuffle_null(mats: list[np.ndarray], rng: np.random.Generator) -> np.ndarray:
    """Build a null correlation distribution by shuffling rows across parents.

    All residual rows from every parent are pooled, randomly permuted, then
    re-assigned to groups of the original sizes.  Pairwise correlations within
    each group are computed and concatenated.  This destroys within-parent
    structure while preserving group sizes, giving the expected correlation
    distribution under the null hypothesis that children are unrelated.
    """
    if not mats:
        return np.empty(0, dtype=np.float64)
    # remember original group sizes so shuffled rows can be re-partitioned identically
    sizes = [m.shape[0] for m in mats]
    # pool all residual rows from all parents into one matrix and shuffle
    pool = np.vstack(mats)
    pool = pool[rng.permutation(pool.shape[0])]
    out = []
    idx = 0
    for size in sizes:
        # re-assign the next 'size' shuffled rows to this group
        chunk = pool[idx:idx + size]
        idx += size
        if size < 2:
            continue
        # skip constant rows: they produce NaN correlations and add no information
        keep = chunk.std(axis=1) > 0
        if keep.sum() < 2:
            continue
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", r"invalid value encountered")
            c = np.corrcoef(chunk[keep])
        iu = np.triu_indices(c.shape[0], k=1)
        vals = c[iu]
        vals = vals[~np.isnan(vals)]
        if vals.size:
            out.append(vals)
    if not out:
        return np.empty(0, dtype=np.float64)
    return np.concatenate(out).astype(np.float64, copy=False)


def _aggregate_pair_rs(parents: list[ParentPrecompute], cutoff: float, min_keep: int) -> np.ndarray:
    chunks = []
    for pp in parents:
        alive = pp.survivors_at(cutoff, min_keep)
        if alive.sum() < 2:
            continue
        chunks.append(_pair_rs_from_C(pp.C, alive))
    if not chunks:
        return np.empty(0, dtype=np.float64)
    return np.concatenate(chunks)


def _excess_cdf_distance(corrected: np.ndarray, null_sorted: np.ndarray) -> float:
    """One-sided excess CDF distance between corrected and null distributions.

    Returns max over r of ``F_null(r) − F_corrected(r)``, clipped to 0 from
    below.  A value of 0 means the corrected distribution is nowhere above the
    null; higher values indicate residual excess correlation.
    """
    if corrected.size == 0 or null_sorted.size == 0:
        return 0.0
    corr_sorted = np.sort(corrected)
    # evaluate both CDFs on the union of all observed r values
    grid = np.unique(np.concatenate([corr_sorted, null_sorted]))
    # searchsorted with side="right" gives F(r) = P(X ≤ r) at each grid point
    f_corr = np.searchsorted(corr_sorted, grid, side="right") / corr_sorted.size
    f_null = np.searchsorted(null_sorted, grid, side="right") / null_sorted.size
    # one-sided: only penalise when the corrected distribution exceeds the null
    # (F_null > F_corr means more mass above r in corrected than null → excess correlation)
    return float(np.max(np.maximum(f_null - f_corr, 0.0)))


def run_level_sweep(
    parents: list[ParentPrecompute],
    null_sorted: np.ndarray,
    *,
    cutoff_grid: tuple[float, ...] = DEFAULT_CUTOFF_GRID,
    tolerance: float = DEFAULT_TOLERANCE,
    min_keep: int = DEFAULT_MIN_KEEP,
    level: tuple[str, str] = ("", ""),
):
    """Sweep correlation cutoffs and return the mildest one within tolerance.

    For each cutoff in ``cutoff_grid`` (scanned from loose to tight), the
    surviving pairwise r values are collected and the excess CDF distance to
    the null is computed.  The first cutoff whose distance falls at or below
    ``tolerance`` is chosen.  If none qualifies the tightest cutoff is used
    unconditionally.

    Parameters
    ----------
    parents:     precomputed parent structures for this level pair
    null_sorted: sorted null pairwise r values (from cross-parent shuffle)
    cutoff_grid: correlation thresholds to scan, ordered loose → tight
    tolerance:   maximum allowed excess CDF distance D after pruning
    min_keep:    minimum children to retain per parent regardless of cutoff
    level:       (parent_level, child_level) label stored in the result

    Returns
    -------
    LevelSweepResult with the chosen cutoff and full diagnostic traces.
    """
    # measure baseline excess distance with no pruning (cutoff = 1.0 keeps everyone)
    baseline = _aggregate_pair_rs(parents, 1.0, min_keep)
    d_before = _excess_cdf_distance(baseline, null_sorted)

    chosen = None
    trace = []
    for cutoff in cutoff_grid:
        touched = 0
        dropped = 0
        chunks = []
        for pp in parents:
            alive = pp.survivors_at(cutoff, min_keep)
            n_drop = int(pp.C.shape[0] - alive.sum())
            if n_drop > 0:
                touched += 1
                dropped += n_drop
            if alive.sum() >= 2:
                chunks.append(_pair_rs_from_C(pp.C, alive))
        corrected = np.concatenate(chunks) if chunks else np.empty(0, dtype=np.float64)
        d = _excess_cdf_distance(corrected, null_sorted)
        trace.append((cutoff, d, dropped, touched))
        # take the first (loosest) cutoff that already satisfies the tolerance —
        # prefer dropping as few children as possible
        if chosen is None and d <= tolerance:
            chosen = (cutoff, d, corrected, dropped, touched)

    # if no cutoff reached tolerance, fall back to the tightest one in the grid
    if chosen is None:
        cutoff, d, dropped, touched = trace[-1]
        corrected = _aggregate_pair_rs(parents, cutoff, min_keep)
        chosen = (cutoff, d, corrected, dropped, touched)

    cutoff, d_after, corrected, dropped, touched = chosen
    return LevelSweepResult(
        level=level,
        cutoff=float(cutoff),
        d_before=float(d_before),
        d_after=float(d_after),
        n_parents=len(parents),
        parents_touched=int(touched),
        children_dropped=int(dropped),
        grid_trace=trace,
        unmodified_sorted=np.sort(baseline),
        corrected_sorted=np.sort(corrected),
        null_sorted=np.asarray(null_sorted, dtype=np.float64),
    )


def attach_lm_residuals(protnodes, df_c1_normed, df_c2_normed, min_n_per_cond=2):
    """Attach per-ion residuals from ``log2(intensity) ~ condition``.

    The AlphaQuant pipeline passes log2-normalized intensities here.  The
    saved ``*.normed.tsv`` tables are exponentiated for output, so callers that
    start from those files must convert them back to log2 before using this
    helper.
    """
    # build a single intensity matrix with all samples from both conditions
    X = pd.concat([df_c1_normed, df_c2_normed], axis=1)
    c1_cols = list(df_c1_normed.columns)
    c2_cols = list(df_c2_normed.columns)
    X = X.astype(float)

    # compute within-condition means per ion (row-wise)
    m1 = X[c1_cols].mean(axis=1, skipna=True)
    m2 = X[c2_cols].mean(axis=1, skipna=True)

    # subtract the within-condition mean: residual = intensity - condition mean
    # this removes the fold-change signal, leaving only condition-independent noise
    res = X.copy()
    res[c1_cols] = X[c1_cols].sub(m1, axis=0)
    res[c2_cols] = X[c2_cols].sub(m2, axis=0)
    # mask ions with too few valid values in either condition — their residuals are unreliable
    n1_ok = X[c1_cols].notna().sum(axis=1) >= int(min_n_per_cond)
    n2_ok = X[c2_cols].notna().sum(axis=1) >= int(min_n_per_cond)
    res.loc[~(n1_ok & n2_ok), :] = np.nan

    for protnode in protnodes:
        # initialise residuals to None on every node before filling
        for node in PreOrderIter(protnode):
            node.residuals = None

        # assign residual vectors to base (leaf) ions by matching ion name to the matrix index
        for node in PreOrderIter(protnode):
            if node.type != "base":
                continue
            if node.name in res.index:
                node.residuals = res.loc[node.name].to_numpy(dtype=float)
            else:
                node.residuals = None

        # propagate residuals bottom-up: each non-base node gets the column-wise mean
        # of its children's residual vectors, so higher-level nodes carry an averaged
        # representation of shared technical variation across their subtree
        for level_nodes in aqcluster_utils.iterate_through_tree_levels_bottom_to_top(protnode):
            for node in level_nodes:
                if node.type == "base":
                    continue
                vecs = [
                    child.residuals
                    for child in node.children
                    if isinstance(getattr(child, "residuals", None), np.ndarray)
                ]
                if not vecs:
                    node.residuals = None
                    continue
                stacked = np.vstack(vecs)
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", message="Mean of empty slice")
                    with np.errstate(all="ignore"):
                        mean_vec = np.nanmean(stacked, axis=0)
                # if all samples are NaN after averaging, treat as missing
                node.residuals = None if np.all(np.isnan(mean_vec)) else mean_vec


def _collect_level_parents(protnodes, parent_level, child_level):
    parents = []
    for protnode in protnodes:
        for node in PreOrderIter(protnode):
            if not _node_matches_level(node, parent_level):
                continue
            child_nodes = []
            vecs = []
            for child in node.children:
                if not _node_matches_level(child, child_level):
                    continue
                v = getattr(child, "residuals", None)
                # skip children without residuals or with any NaN sample —
                # NaN entries would propagate into the correlation matrix
                if v is None or not isinstance(v, np.ndarray):
                    continue
                if np.any(np.isnan(v)):
                    continue
                child_nodes.append(child)
                vecs.append(v)
            # need at least 2 children to form a correlation matrix
            if len(vecs) < 2:
                continue
            mat = np.vstack(vecs)
            pp = _build_parent(node, child_nodes, mat)
            if pp is not None:
                parents.append(pp)
    return parents


def apply_residual_decorrelation(
    protnodes,
    df_c1_normed,
    df_c2_normed,
    *,
    tolerance=DEFAULT_TOLERANCE,
    min_keep=DEFAULT_MIN_KEEP,
    cutoff_grid=DEFAULT_CUTOFF_GRID,
    aggregation_mode="stouffer_decorrelation",
    null_seed=42,
    plot_dir=None,
):
    """Main entry point: run full residual decorrelation on a list of protein nodes.

    Steps
    -----
    1. Attach within-condition mean residuals to every node (``attach_lm_residuals``).
    2. For each level pair in ``LEVEL_PAIRS``, build a cross-parent shuffle null,
       run the cutoff sweep, and mark pruned children with
       ``node.exclude_residual_decorrelation = True``.
    3. Optionally apply PTM fragment selection if ``PTM_FRAGMENT_SELECTION`` is set.
    4. Re-aggregate all node statistics bottom-up using ``aggregation_mode``.
    5. Strip residual arrays from nodes to keep the tree serializable.

    Parameters
    ----------
    protnodes:        list of root protein nodes (anytree)
    df_c1_normed:     log2-normalized intensities for condition 1 (ions × samples)
    df_c2_normed:     log2-normalized intensities for condition 2 (ions × samples)
    tolerance:        maximum excess CDF distance D allowed after pruning
    min_keep:         minimum children retained per parent at each level
    cutoff_grid:      correlation cutoffs to sweep per level pair
    aggregation_mode: z-aggregation mode used when re-aggregating after pruning
    null_seed:        random seed for the cross-parent shuffle null

    Returns
    -------
    pandas.DataFrame summarising cutoffs, distances, and drop counts per level.
    """
    # reset exclusion flags in case this function is called more than once on the same nodes
    for protnode in protnodes:
        for node in PreOrderIter(protnode):
            node.exclude_residual_decorrelation = False
            node.exclude_ptm_fragment_selection = False
            if aqvariables.RESIDUAL_DEFF_CORRECTION:
                node.icc_correction = 0.0

    # step 1: compute within-condition residuals and attach them to every node
    attach_lm_residuals(protnodes, df_c1_normed, df_c2_normed)

    rng = np.random.default_rng(null_seed)
    level_results = []

    # step 2: run the sweep for each level pair independently
    for parent_level, child_level in LEVEL_PAIRS:
        parents = _collect_level_parents(protnodes, parent_level, child_level)
        # build the null from the same residual matrices used for the real sweep
        mats = [
            np.vstack([child.residuals for child in pp.child_nodes])
            for pp in parents
        ]
        null_sorted = np.sort(_cross_parent_shuffle_null(mats, rng))
        sweep = run_level_sweep(
            parents,
            null_sorted,
            cutoff_grid=cutoff_grid,
            tolerance=tolerance,
            min_keep=min_keep,
            level=(parent_level, child_level),
        )
        level_results.append(sweep)
        msg = (
            f"residual decorrelation {parent_level}->{child_level}: "
            f"cutoff={sweep.cutoff:g} "
            f"D={sweep.d_before:.4g}->{sweep.d_after:.4g} "
            f"dropped={sweep.children_dropped:,} "
            f"parents_touched={sweep.parents_touched:,}/{sweep.n_parents:,}"
        )
        LOGGER.info(msg)
        print(msg, flush=True)

        # design-effect on the RESIDUAL correlation. The mean pairwise correlation
        # among surviving children, in EXCESS of the cross-parent shuffle null's mean,
        # is recorded on the parent so aggregation applies deff=1+(n-1)*rho. Subtracting
        # the null mean removes the finite-sample bias of correlations estimated from
        # short residual vectors (few replicates), so clean/few-rep data where survivors
        # are no more correlated than chance gets rho~0 (no-op); genuine homogeneous
        # between-child correlation that pruning cannot drop gets rho>0 (corrected).
        null_mean = float(np.mean(null_sorted)) if null_sorted.size else 0.0
        survivor_means = []
        parents_with_pairs = []
        # mark children that did not survive the chosen cutoff
        for pp in parents:
            survivors = pp.survivors_at(sweep.cutoff, min_keep)
            for keep, child in zip(survivors, pp.child_nodes):
                if not keep:
                    child.exclude_residual_decorrelation = True
            if aqvariables.RESIDUAL_DEFF_CORRECTION:
                resid_rs = _pair_rs_from_C(pp.C, survivors)
                if resid_rs.size:
                    survivor_means.append(float(np.mean(resid_rs)))
                    parents_with_pairs.append(pp)
                else:
                    pp.parent_node.icc_correction = 0.0
        # The design effect is applied ONLY at the between-peptide level (gene->seq).
        # That level carries the protein-level random effect (delta) shared across a
        # protein's peptides -- the ICC that EmpiRe's ion-variance model does NOT capture.
        # The within-peptide levels (precursor/ion-family/fragment/MS1) are already handled
        # by the EmpiRe variance model plus pruning, so applying deff there double-counts
        # and needlessly costs power on well-calibrated data.
        #
        # GATE: only correct when the peptide-level correlation actually exceeded the
        # pruning tolerance (d_before > tolerance). If peptides are no more correlated than
        # the null to begin with (clean/spike-in data, e.g. d_before ~ 0.03 < tol 0.10),
        # pruning does nothing and there is no ICC to correct -- so deff must be a no-op,
        # otherwise it needlessly costs sensitivity on well-calibrated datasets.
        gate_open = sweep.d_before > tolerance
        if aqvariables.RESIDUAL_DEFF_CORRECTION and parent_level == "gene":
            LOGGER.info("deff gate gene->seq: d_before=%.4f tolerance=%.4f -> %s",
                        sweep.d_before, tolerance, "OPEN" if gate_open else "CLOSED (deff off)")
        if aqvariables.RESIDUAL_DEFF_CORRECTION and parent_level == "gene" and gate_open:
            # Apply a single robust per-LEVEL excess correlation to every parent at this
            # level: clip ONCE on the well-estimated level mean rather than per parent.
            # Per-parent survivor means are noisy (few children/replicates); clipping each
            # at zero rectifies that noise into a large positive bias, which over-corrects
            # clean/few-replicate data. The level mean averages the noise out; subtracting
            # the shuffle-null mean removes the finite-sample floor.
            level_excess = (max(0.0, float(np.mean(survivor_means)) - null_mean)
                            if survivor_means else 0.0)
            for pp in parents_with_pairs:
                pp.parent_node.icc_correction = level_excess
            LOGGER.info(
                "deff %s->%s: mean survivor rho=%.4f  null mean=%.4f  "
                "LEVEL excess rho=%.4f (parents=%d)",
                parent_level, child_level,
                float(np.mean(survivor_means)) if survivor_means else 0.0,
                null_mean, level_excess, len(survivor_means),
            )

    # optional: save the per-level distribution diagnostics (before/after/null CDFs
    # + cutoff sweep trace) using AlphaQuant's own plotting.
    if plot_dir is not None:
        import os
        os.makedirs(plot_dir, exist_ok=True)
        for sweep in level_results:
            try:
                fig = plot_level_sweep_diagnostics(sweep)
                fig.savefig(os.path.join(plot_dir, f"decorr_{sweep.level[0]}__{sweep.level[1]}.png"),
                            dpi=120)
                plt.close(fig)
            except Exception as exc:
                LOGGER.warning("could not save decorrelation plot for %s: %s", sweep.level, exc)

    # step 3 (optional): apply PTM fragment selection on top of decorrelation exclusions
    if aqvariables.PTM_FRAGMENT_SELECTION:
        n_ptm_dropped, n_ptm_parents = aqcluster_utils.apply_ptm_fragment_selection(
            protnodes,
        )
        LOGGER.info(
            "PTM fragment low-|Z| selection after residual decorrelation: "
            "dropped %s children across %s frgion parents",
            n_ptm_dropped,
            n_ptm_parents,
        )
        print(
            "PTM fragment low-|Z| selection after residual decorrelation: "
            f"dropped {n_ptm_dropped:,} children across "
            f"{n_ptm_parents:,} frgion parents",
            flush=True,
        )

    # step 4: re-aggregate node statistics bottom-up now that exclusion flags are set
    for protnode in protnodes:
        for level_nodes in aqcluster_utils.iterate_through_tree_levels_bottom_to_top(protnode):
            for node in level_nodes:
                if node.type == "base":
                    continue
                aqcluster_utils.aggregate_node_properties(
                    node,
                    only_use_mainclust=True,
                    peptide_outlier_filtering=False,
                    aggregation_mode=aggregation_mode,
                )

    # step 5: residual vectors are only needed during the sweep; remove them before
    # downstream ML reordering / JSON export to keep the tree serializable
    for protnode in protnodes:
        for node in PreOrderIter(protnode):
            if hasattr(node, "residuals"):
                delattr(node, "residuals")

    summary = pd.DataFrame(
        [
            {
                "parent_level": result.level[0],
                "child_level": result.level[1],
                "cutoff": result.cutoff,
                "d_before": result.d_before,
                "d_after": result.d_after,
                "n_parents": result.n_parents,
                "parents_touched": result.parents_touched,
                "children_dropped": result.children_dropped,
            }
            for result in level_results
        ]
    )
    if not summary.empty:
        LOGGER.info("Residual decorrelation summary:\n%s", summary.to_string(index=False))
    return summary


def plot_level_sweep_cdfs(level_result: LevelSweepResult, *, ax=None, title: str | None = None):
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4))

    for values, label, color in (
        (level_result.null_sorted, "null", "#7f8c8d"),
        (level_result.unmodified_sorted, "before", "#d35400"),
        (level_result.corrected_sorted, "after", "#1f77b4"),
    ):
        if values.size == 0:
            continue
        y = np.arange(1, values.size + 1) / values.size
        ax.step(values, y, where="post", label=label, color=color)

    ax.set_xlabel("pairwise residual correlation r")
    ax.set_ylabel("cumulative fraction")
    if title is None:
        title = (
            f"{level_result.level[0]} -> {level_result.level[1]} | "
            f"cutoff={level_result.cutoff:.2f}, D={level_result.d_after:.3f}"
        )
    ax.set_title(title)
    ax.legend()
    ax.grid(alpha=0.25)
    return ax


def plot_level_sweep_diagnostics(level_result: LevelSweepResult):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    plot_level_sweep_cdfs(level_result, ax=axes[0])

    if level_result.grid_trace:
        cutoffs = [x[0] for x in level_result.grid_trace]
        distances = [x[1] for x in level_result.grid_trace]
        dropped = [x[2] for x in level_result.grid_trace]
        axes[1].plot(cutoffs, distances, marker="o", color="#1f77b4", label="excess CDF distance")
        axes[1].axvline(level_result.cutoff, color="#d35400", linestyle="--", label="chosen cutoff")
        axes[1].scatter([level_result.cutoff], [level_result.d_after], color="#d35400", zorder=3)
        ax2 = axes[1].twinx()
        ax2.bar(cutoffs, dropped, width=0.06, alpha=0.18, color="#2c3e50", label="children dropped")
        axes[1].set_xlabel("cutoff")
        axes[1].set_ylabel("excess CDF distance")
        ax2.set_ylabel("children dropped")
        axes[1].invert_xaxis()
        axes[1].grid(alpha=0.25)
        lines1, labels1 = axes[1].get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        axes[1].legend(lines1 + lines2, labels1 + labels2, loc="best")
        axes[1].set_title(
            f"{level_result.level[0]} -> {level_result.level[1]} sweep"
        )

    fig.tight_layout()
    return fig
