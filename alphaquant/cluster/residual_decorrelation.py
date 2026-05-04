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
    parent_node: anytree.Node
    child_nodes: tuple[anytree.Node, ...]
    C: np.ndarray
    remove_order: np.ndarray
    max_r_trajectory: np.ndarray

    def survivors_at(self, cutoff: float, min_keep: int) -> np.ndarray:
        n = self.C.shape[0]
        if n == 0:
            return np.zeros(0, dtype=bool)
        k_max = max(0, n - min_keep)
        k = 0
        while k <= k_max and self.max_r_trajectory[k] > cutoff:
            k += 1
        if k > k_max:
            k = k_max
        alive = np.ones(n, dtype=bool)
        if k > 0:
            alive[self.remove_order[:k]] = False
        return alive


@dataclass
class LevelSweepResult:
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
    if mat.shape[0] < 2:
        return None

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", r"invalid value encountered")
        C = np.corrcoef(mat)
    if not np.all(np.isfinite(C)):
        C = np.nan_to_num(C, nan=0.0, posinf=1.0, neginf=-1.0)
    C = C.copy()
    np.fill_diagonal(C, np.nan)

    n = C.shape[0]
    alive = np.ones(n, dtype=bool)
    remove_order = []
    max_r = []

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", r"All-NaN slice encountered")
        init_max = float(np.nanmax(C)) if n >= 2 else -np.inf
    if not np.isfinite(init_max):
        init_max = -np.inf
    max_r.append(init_max)

    while alive.sum() > 1:
        sub_idx = np.where(alive)[0]
        cc = C[np.ix_(sub_idx, sub_idx)]
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", r"Mean of empty slice")
            warnings.filterwarnings("ignore", r"All-NaN slice encountered")
            mean_r = np.nanmean(cc, axis=1)
        worst_local = int(np.nanargmax(np.where(np.isnan(mean_r), -np.inf, mean_r)))
        alive[sub_idx[worst_local]] = False
        remove_order.append(sub_idx[worst_local])
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
    iu = np.triu_indices(sub.shape[0], k=1)
    vals = sub[iu]
    vals = vals[~np.isnan(vals)]
    return vals.astype(np.float64, copy=False)


def _cross_parent_shuffle_null(mats: list[np.ndarray], rng: np.random.Generator) -> np.ndarray:
    if not mats:
        return np.empty(0, dtype=np.float64)
    sizes = [m.shape[0] for m in mats]
    pool = np.vstack(mats)
    pool = pool[rng.permutation(pool.shape[0])]
    out = []
    idx = 0
    for size in sizes:
        chunk = pool[idx:idx + size]
        idx += size
        if size < 2:
            continue
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
    if corrected.size == 0 or null_sorted.size == 0:
        return 0.0
    corr_sorted = np.sort(corrected)
    grid = np.unique(np.concatenate([corr_sorted, null_sorted]))
    f_corr = np.searchsorted(corr_sorted, grid, side="right") / corr_sorted.size
    f_null = np.searchsorted(null_sorted, grid, side="right") / null_sorted.size
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
        if chosen is None and d <= tolerance:
            chosen = (cutoff, d, corrected, dropped, touched)

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
    X = pd.concat([df_c1_normed, df_c2_normed], axis=1)
    c1_cols = list(df_c1_normed.columns)
    c2_cols = list(df_c2_normed.columns)
    X = X.astype(float)

    m1 = X[c1_cols].mean(axis=1, skipna=True)
    m2 = X[c2_cols].mean(axis=1, skipna=True)

    res = X.copy()
    res[c1_cols] = X[c1_cols].sub(m1, axis=0)
    res[c2_cols] = X[c2_cols].sub(m2, axis=0)
    n1_ok = X[c1_cols].notna().sum(axis=1) >= int(min_n_per_cond)
    n2_ok = X[c2_cols].notna().sum(axis=1) >= int(min_n_per_cond)
    res.loc[~(n1_ok & n2_ok), :] = np.nan

    for protnode in protnodes:
        for node in PreOrderIter(protnode):
            node.residuals = None

        for node in PreOrderIter(protnode):
            if node.type != "base":
                continue
            if node.name in res.index:
                node.residuals = res.loc[node.name].to_numpy(dtype=float)
            else:
                node.residuals = None

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
                if v is None or not isinstance(v, np.ndarray):
                    continue
                if np.any(np.isnan(v)):
                    continue
                child_nodes.append(child)
                vecs.append(v)
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
):
    for protnode in protnodes:
        for node in PreOrderIter(protnode):
            node.exclude_residual_decorrelation = False
            node.exclude_ptm_fragment_selection = False

    attach_lm_residuals(protnodes, df_c1_normed, df_c2_normed)

    rng = np.random.default_rng(null_seed)
    level_results = []

    for parent_level, child_level in LEVEL_PAIRS:
        parents = _collect_level_parents(protnodes, parent_level, child_level)
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

        for pp in parents:
            survivors = pp.survivors_at(sweep.cutoff, min_keep)
            for keep, child in zip(survivors, pp.child_nodes):
                if not keep:
                    child.exclude_residual_decorrelation = True

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

    # Residual vectors are only needed during the sweep; remove them before
    # downstream ML reordering / JSON export to keep the tree serializable.
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
