import matplotlib.pyplot as plt
import sklearn.metrics
import seaborn as sns
import numpy as np


def plot_value_histogram(y, results_dir = None):
    fig, ax = plt.subplots()
    sns.histplot(y, ax = ax)
    if results_dir is not None:
        fig.savefig(f"{results_dir}/value_histogram.pdf")
    plt.show()


def scatter_ml_regression_combined(y_test, y_pred, results_dir = None):
    fig, ax = plt.subplots()

    sns.regplot(x = abs(y_test), y = y_pred, scatter_kws=dict(alpha=0.1), ax = ax)
    err = sklearn.metrics.mean_squared_error(y_test, y_pred)
    r2 = sklearn.metrics.r2_score(y_test, y_pred)

    ax.set_xlabel("true offset")
    ax.set_ylabel("predicted offset")

    ax.set_title(f"MSE: {err:.2f}, R2: {r2:.2f}")

    if results_dir is not None:
        fig.savefig(f"{results_dir}/ml_regression.pdf")
    plt.show()

def scatter_ml_regression_testsets(test_set_predictions, results_dir = None):

    fig, axes = plt.subplots(ncols=len(test_set_predictions), figsize=(2.5 * len(test_set_predictions), 4))

    for idx in range(len(test_set_predictions)):
        y_true, y_pred = test_set_predictions[idx]
        ax = axes[idx]

        sns.regplot(x = y_true, y = y_pred, scatter_kws=dict(alpha=0.1), ax = ax)
        err = sklearn.metrics.mean_squared_error(y_true, y_pred)
        r2 = sklearn.metrics.r2_score(y_true, y_pred)
        ax.set_xlabel("true offset")
        ax.set_ylabel("predicted offset")

        ax.set_title(f"MSE: {err:.2f}, R2: {r2:.2f}")

    if results_dir is not None:
        fig.savefig(f"{results_dir}/ml_regression.pdf")
    plt.show()

def plot_feature_importance_per_model(models, featurenames,top_n = np.inf, results_dir = None):
    fig, axes = plt.subplots(1, len(models), figsize=(2.5 * len(models), 4))
    for modelidx in range(len(models)):
        model = models[modelidx]
        ax = axes[modelidx]
        plot_feature_importances(model.feature_importances_, featurenames, top_n, results_dir, ax)
    fig.tight_layout()
    if results_dir is not None:
        fig.savefig(f"{results_dir}/ml_feature_importances.pdf")
        


def plot_feature_importances(coef, names, top_n = np.inf, results_dir = None, ax = None):
    importance_score,names = filter_sort_top_n(coef, names, top_n)
    if ax is None:
        fig, ax = plt.subplots()
    ax.set_title('Feature Importances')
    ax.bar(range(len(names)), importance_score, align='center')
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=90)
    if results_dir is not None:
        if 'fig' in locals():
            fig.savefig(f"{results_dir}/ml_feature_importances.pdf")
            



import sklearn.inspection
import matplotlib.pyplot as plt
import numpy as np

def compute_and_plot_feature_importances_stacked_rf(model, X_val, y_val, feature_names, top_n=np.inf, results_dir=None, n_repeats=10, random_state=42):
    # Compute permutation importance
    result = sklearn.inspection.permutation_importance(model, X_val, y_val, n_repeats=n_repeats, random_state=random_state)

    # Sort features by importance
    importance_scores = result.importances_mean
    if top_n < len(feature_names):
        sorted_idx = importance_scores.argsort()[-top_n:]
    else:
        sorted_idx = importance_scores.argsort()

    sorted_importances = importance_scores[sorted_idx]
    sorted_feature_names = np.array(feature_names)[sorted_idx]

    # Plot feature importances
    fig, ax = plt.subplots()
    ax.set_title('Feature Importances')
    ax.barh(range(len(sorted_feature_names)), sorted_importances, align='center')
    ax.set_yticks(range(len(sorted_feature_names)), sorted_feature_names)
    ax.set_yticklabels(sorted_feature_names)
    plt.tight_layout()  # Adjust layout to make room for the feature names

    if results_dir is not None:
        fig_path = f"{results_dir}/stacked_regressor_feature_importances.pdf"
        fig.savefig(fig_path)
        print(f"Feature importances plot saved to {fig_path}")
    
    plt.show()


def filter_sort_top_n(imp, names, top_n):
    tuplelist = list(zip(imp, names))
    tuplelist.sort(key = lambda x : abs(x[0]),reverse=True)
    tuplelist = tuplelist[:top_n]
    imp = [x[0] for x in tuplelist]
    names = [x[1] for x in tuplelist]
    return imp, names

