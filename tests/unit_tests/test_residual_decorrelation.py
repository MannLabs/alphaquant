from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import pytest

import anytree
import numpy as np
import pandas as pd

import alphaquant.cluster.residual_decorrelation as aq_resid


def _load_reference_module():
    path = Path("sandbox/analyses_revision_v3/paper_nbs_revision/10_alphaquant_mouse_aq/residual_correlation/auto_decorrelation.py")
    if not path.exists():
        pytest.skip("Reference sandbox file not available in this environment")
    spec = importlib.util.spec_from_file_location("auto_decorrelation_ref", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_survivors_at_drops_expected_child():
    C = np.array(
        [
            [np.nan, 0.95, 0.20],
            [0.95, np.nan, 0.80],
            [0.20, 0.80, np.nan],
        ]
    )
    pp = aq_resid.ParentPrecompute(
        parent_node=anytree.Node("parent", type="frgion"),
        child_nodes=(
            anytree.Node("c0", type="base"),
            anytree.Node("c1", type="base"),
            anytree.Node("c2", type="base"),
        ),
        C=C,
        remove_order=np.array([1, 0]),
        max_r_trajectory=np.array([0.95, 0.20, -np.inf]),
    )

    survivors = pp.survivors_at(0.5, min_keep=1)
    np.testing.assert_array_equal(survivors, np.array([True, False, True]))


def test_excess_cdf_distance_positive_when_corrected_shifted_higher():
    corrected = np.array([0.7, 0.8, 0.9])
    null = np.sort(np.array([0.1, 0.2, 0.3]))
    assert aq_resid._excess_cdf_distance(corrected, null) > 0


def test_attach_lm_residuals_subtracts_condition_means():
    root = anytree.Node("gene1", type="gene")
    frg = anytree.Node("pep1", parent=root, type="frgion")
    base = anytree.Node("ionA", parent=frg, type="base")

    vals = np.log2(np.array([10.0, 14.0, 20.0, 24.0]))
    df_c1 = pd.DataFrame([vals[:2]], index=["ionA"], columns=["c1_r1", "c1_r2"])
    df_c2 = pd.DataFrame([vals[2:]], index=["ionA"], columns=["c2_r1", "c2_r2"])

    aq_resid.attach_lm_residuals([root], df_c1, df_c2)

    expected = vals.copy()
    expected[:2] -= expected[:2].mean()
    expected[2:] -= expected[2:].mean()
    np.testing.assert_allclose(base.residuals, expected)
    np.testing.assert_allclose(frg.residuals, expected)
    np.testing.assert_allclose(root.residuals, expected)


def test_reference_port_matches_dashboard_core():
    ref = _load_reference_module()
    rng = np.random.default_rng(7)
    mats = [
        rng.normal(size=(4, 6)),
        rng.normal(size=(3, 6)),
        rng.normal(size=(5, 6)),
    ]

    our_parents = []
    ref_parents = []
    for idx, mat in enumerate(mats):
        parent = anytree.Node(f"p{idx}", type="seq")
        children = tuple(anytree.Node(f"c{idx}_{j}", parent=parent, type="mod_seq") for j in range(mat.shape[0]))
        our_pp = aq_resid._build_parent(parent, children, mat)
        ref_pp = ref._build_parent(f"g{idx}", f"p{idx}", [c.name for c in children], mat)
        assert our_pp is not None
        assert ref_pp is not None
        np.testing.assert_allclose(our_pp.C, ref_pp.C)
        np.testing.assert_array_equal(our_pp.remove_order, ref_pp.remove_order)
        np.testing.assert_allclose(our_pp.max_r_trajectory, ref_pp.max_r_trajectory)
        our_parents.append(our_pp)
        ref_parents.append(ref_pp)

    our_null = np.sort(aq_resid._cross_parent_shuffle_null(mats, np.random.default_rng(42)))
    ref_null = np.sort(ref._cross_parent_shuffle_null(mats, np.random.default_rng(42)))
    np.testing.assert_allclose(our_null, ref_null)

    our_sweep = aq_resid.run_level_sweep(
        our_parents,
        our_null,
        cutoff_grid=aq_resid.DEFAULT_CUTOFF_GRID,
        tolerance=aq_resid.DEFAULT_TOLERANCE,
        min_keep=aq_resid.DEFAULT_MIN_KEEP,
        level=("seq", "mod_seq"),
    )
    ref_sweep = ref.run_level_sweep(
        ref_parents,
        ref_null,
        cutoff_grid=ref.DEFAULT_CUTOFF_GRID,
        tolerance=ref.DEFAULT_TOLERANCE,
        min_keep=ref.DEFAULT_MIN_KEEP,
        level=("seq", "mod_seq"),
    )

    assert our_sweep.cutoff == ref_sweep.cutoff
    assert our_sweep.d_before == ref_sweep.d_before
    assert our_sweep.d_after == ref_sweep.d_after
    assert our_sweep.children_dropped == ref_sweep.children_dropped
    assert our_sweep.parents_touched == ref_sweep.parents_touched
    assert our_sweep.grid_trace == ref_sweep.grid_trace
    np.testing.assert_allclose(our_sweep.unmodified_sorted, ref_sweep.unmodified_sorted)
    np.testing.assert_allclose(our_sweep.corrected_sorted, ref_sweep.corrected_sorted)
