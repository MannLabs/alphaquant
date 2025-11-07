# For Developers: Modifying AlphaQuant

AlphaQuant is designed with modularity in mind to allow practitioners to introduce alternative numerical methods for each module. The codebase follows clear interfaces that make it straightforward to extend or replace statistical methods at different levels of the analysis pipeline.

## 1. Ion-Level Statistical Testing

**Where to modify:** `alphaquant/diffquant/diff_analysis.py`

**How it works:** Each ion (fragment, peptide, etc.) is tested independently for differential expression. The test produces three key outputs: `p_val` (p-value), `fc` (log2 fold change), and `z_val` (z-score for aggregation).

**Main class:**
- **`DifferentialIon`** - The default method that uses intensity-dependent empirical background distributions to compute p-values and z-scores. It accounts for technical variation by comparing observed fold changes against distributions derived from similarly abundant ions in the dataset. The core statistical logic is in the `_calc_diffreg_peptide()` method.

**How to extend:** We've included `DifferentialIonTTest` in the same file as example code demonstrating how to implement alternative tests. This variant uses Welch's t-test with robust variance estimation (similar to MS-EmpiRe). **To create your own method:**

1. Create a new class (e.g., `DifferentialIonMyMethod`) with the same interface:
   - `__init__()` should accept `(noNanvals_from, noNanvals_to, ...)` and any method-specific parameters
   - Set attributes: `name`, `p_val`, `fc`, `z_val`, `usable`
2. Implement your statistical test in a method (e.g., `_calc_mymethod()`)
3. Modify `alphaquant/diffquant/condpair_analysis.py` (lines 67-70) to instantiate your class
4. Optionally, add a parameter to `run_pipeline()` to select between methods

The key requirement is that your class must output `p_val`, `fc`, and `z_val` for each ion—these are used by the tree aggregation framework.

## 2. Tree-Based Ion Propagation

**Where to modify:** `alphaquant/cluster/cluster_utils.py` and `alphaquant/cluster/cluster_ions.py`

**How it works:** Statistics from child nodes (e.g., fragments) are aggregated to parent nodes (e.g., peptides → proteins) in a hierarchical tree. Z-values are combined using Stouffer's method, and fold changes are summarized using medians.

**Key functions:**
- **`aggregate_node_properties()`** - The core function that propagates statistics up the tree. It combines z-values, fold changes, and quality metrics from children to parents.
- **`sum_and_re_scale_zvalues()`** - Implements Stouffer's Z-score method: sums z-values and divides by sqrt(n), then rescales to maintain standard normal distribution.
- **`transform_znormed_to_pval()`** - Converts aggregated z-scores back to two-sided p-values.

**How to extend:** If you want to use different aggregation methods:
1. Modify `sum_and_re_scale_zvalues()` to implement your preferred meta-analysis method (e.g., Fisher's method, weighted Z-scores, etc.)
2. If your method changes the distribution, update `transform_znormed_to_pval()` accordingly
3. For fold-change aggregation, modify line 67 in `aggregate_node_properties()` where `node.fc = np.median(fcs)` is set

The tree traversal itself is in `cluster_ions.py`:
- **`cluster_along_specified_levels()`** - Iterates through tree levels bottom-to-top
- **`get_scored_clusterselected_ions()`** - Entry point for the hierarchical workflow

## 3. Multiple Testing Correction

**Where to modify:** `alphaquant/tables/diffquant_table.py` and `alphaquant/tables/proteoformtable.py`

**How it works:** FDR correction is applied separately to different result tables during output generation. The method outputs p-values in all tables, so you can always recalculate q-values from the output files.

**Key functions:**
- **Protein results** (`alphaquant/tables/diffquant_table.py`):
  - `_add_fdr_fc_based_set()` - Applies Benjamini-Hochberg to intensity-based proteins
  - `_add_fdr_counting_based_set()` - Applies adjusted Benjamini-Hochberg to proteins detected only via missing values

- **Proteoform results** (`alphaquant/tables/proteoformtable.py`):
  - `_annotate_fdr_column()` - Applies Benjamini-Hochberg to test if alternative proteoforms differ from the reference

**How to extend:**
1. Modify the relevant function to use a different method (e.g., Bonferroni, Storey's q-value, etc.)
2. Replace the `mt.multipletests(..., method='fdr_bh', ...)` call with your preferred correction
3. Alternatively, use the p-values from output tables and apply your own correction externally

## 4. Outlier Robustness

**Where to modify:** `alphaquant/diffquant/diff_analysis.py` and `alphaquant/cluster/cluster_utils.py`

**How it works:** AlphaQuant applies outlier correction at two levels to make results robust to technical variation and biological heterogeneity.

**Key functions:**
- **`calc_outlier_scaling_factor()`** (in `diff_analysis.py`) - Compares between-replicate variance to expected technical variance and inflates estimates when replicates show unusual variability
- **`remove_outlier_fragion_childs()`** (in `cluster_utils.py`) - Filters extreme fragments before aggregating to peptides (keeps the 5 most central fragments when >4 are available)

**How to extend:**
1. Modify the scaling logic in `calc_outlier_scaling_factor()` to use different robust estimators
2. Adjust `remove_outlier_fragion_childs()` to change how many fragments are retained or which criteria are used for selection
3. Set `outlier_correction=False` in `run_pipeline()` to disable this feature entirely

## 5. Main Workflow Orchestration

**Where to modify:** `alphaquant/diffquant/condpair_analysis.py`

**How it works:** The `analyze_condpair()` function coordinates the complete pipeline for comparing two conditions.

**Pipeline steps:**
1. Load and filter data for the two conditions
2. Perform normalization (within and between conditions)
3. Create empirical background distributions
4. Compute ion-level differential statistics (`DifferentialIon` or `DifferentialIonTTest`)
5. Build hierarchical trees and perform clustering to identify proteoforms
6. Apply machine learning quality scoring (if enabled)
7. Filter outlier peptides (if enabled)
8. Generate output tables with FDR correction
9. Create visualization plots

**How to extend:** This file shows how all components connect. To add custom preprocessing, normalization, or post-processing steps, modify this function or create a wrapper that calls it with modified data.

---

## Additional Resources

For general contribution guidelines, code style, and how to submit pull requests, please see [CONTRIBUTING.md](CONTRIBUTING.md).

For questions or discussions about extending AlphaQuant, please use the [GitHub Discussions](https://github.com/MannLabs/alphaquant/discussions) forum.

