import alphaquant.config.config as aqconfig

import sklearn.ensemble
import sklearn.linear_model
import sklearn.impute
import sklearn.metrics
import sklearn.model_selection

from sklearn.model_selection import KFold, RandomizedSearchCV
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import RandomizedSearchCV, KFold
from sklearn.inspection import permutation_importance

import logging
aqconfig.setup_logging()
LOGGER = logging.getLogger(__name__)

import numpy as np

def train_random_forest_with_grid_search(X, y, shorten_features_for_speed, num_splits=5):
    y = np.abs(y)
    models = []
    test_set_predictions = []
    y_pred_cv = np.zeros_like(y)

    kf = sklearn.model_selection.KFold(n_splits=num_splits, shuffle=True, random_state=42)

    for fold_num, (train_index, test_index) in enumerate(kf.split(X), 1):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]

        # Define the parameter grid for GridSearchCV
        param_grid = {
            'n_estimators': [100, 200],
            'max_depth': [None, 5, 10],
            'min_samples_leaf': [1, 5],
            'max_features': [None, 'sqrt'],
            'max_samples': [None, 0.8]
        }

        # Adjust max_features based on the flag
        if shorten_features_for_speed:
            param_grid['max_features'] = ['sqrt']
        else:
            param_grid['max_features'] = [None, 'sqrt']

        # Initialize the Random Forest Regressor
        rf = sklearn.ensemble.RandomForestRegressor(random_state=42 + fold_num)

        # Initialize GridSearchCV
        grid_search = sklearn.model_selection.GridSearchCV(
            estimator=rf,
            param_grid=param_grid,
            cv=sklearn.model_selection.KFold(n_splits=3, shuffle=True, random_state=42),
            scoring='r2',
            n_jobs=-1,
            verbose=0,
            return_train_score=True
        )

        # Fit the GridSearchCV on the training data
        grid_search.fit(X_train, y_train)

        # Get the best estimator
        best_model = grid_search.best_estimator_

        # Print the best parameters
        LOGGER.info(f"Fold {fold_num} Best parameters found:", grid_search.best_params_)

        # Predict on the test set
        y_pred_test = best_model.predict(X_test)
        y_pred_cv[test_index] = y_pred_test

        # Collect the model
        models.append(best_model)

        # Collect test set predictions
        test_set_predictions.append((y_test, y_pred_test))

        # Evaluate performance
        fold_mse = np.mean((y_test - y_pred_test) ** 2)
        LOGGER.info(f"Fold {fold_num} MSE: {fold_mse}")

        correlation = np.corrcoef(y_test, y_pred_test)[0, 1]
        LOGGER.info(f"Overall correlation in fold {fold_num}: {correlation}")

    # Return the list of models, test set predictions, and out-of-fold predictions
    return models, test_set_predictions, y_pred_cv




def train_gradient_boosting_with_grid_search(X, y, shorten_features_for_speed, num_splits=5):
    # Do not take the absolute value of y
    y = np.abs(y)  # Ensure y retains its sign

    models = []
    test_set_predictions = []
    y_pred_cv = np.zeros_like(y)

    kf = sklearn.model_selection.KFold(n_splits=num_splits, shuffle=True, random_state=42)

    for fold_num, (train_index, test_index) in enumerate(kf.split(X), 1):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]

        # Define the parameter grid for GridSearchCV
        param_grid = {
            'n_estimators': [100, 200],
            'learning_rate': [0.01, 0.1, 0.2],
            'max_depth': [3, 5, 7],
            'min_samples_leaf': [1, 5],
            'max_features': [None, 'sqrt'],
            'subsample': [1.0, 0.8]
        }

        # Adjust max_features based on the flag
        if shorten_features_for_speed:
            param_grid['max_features'] = ['sqrt']
        else:
            param_grid['max_features'] = [None, 'sqrt']

        # Initialize the Gradient Boosting Regressor
        gbr = sklearn.ensemble.GradientBoostingRegressor(random_state=42 + fold_num)

        # Initialize GridSearchCV
        grid_search = sklearn.model_selection.GridSearchCV(
            estimator=gbr,
            param_grid=param_grid,
            cv=sklearn.model_selection.KFold(n_splits=3, shuffle=True, random_state=42),
            scoring='r2',
            n_jobs=-1,
            verbose=0,
            return_train_score=True
        )

        # Fit the GridSearchCV on the training data
        grid_search.fit(X_train, y_train)

        # Get the best estimator
        best_model = grid_search.best_estimator_

        # Print the best parameters
        LOGGER.info(f"Fold {fold_num} Best parameters found:", grid_search.best_params_)

        # Predict on the test set
        y_pred_test = best_model.predict(X_test)
        y_pred_cv[test_index] = y_pred_test

        # Collect the model
        models.append(best_model)

        # Collect test set predictions
        test_set_predictions.append((y_test, y_pred_test))

        # Evaluate performance
        fold_mse = np.mean((y_test - y_pred_test) ** 2)
        LOGGER.info(f"Fold {fold_num} MSE: {fold_mse}")

        correlation = np.corrcoef(y_test, y_pred_test)[0, 1]
        LOGGER.info(f"Overall correlation in fold {fold_num}: {correlation}")

    # Return the list of models, test set predictions, and out-of-fold predictions
    return models, test_set_predictions, y_pred_cv



def train_gradient_boosting_with_random_search(X, y, shorten_features_for_speed, num_splits=5, n_iter=10):
    # Take the absolute value of y to predict magnitudes only
    y = np.abs(y)

    models = []
    test_set_predictions = []
    y_pred_cv = np.zeros_like(y)

    kf = sklearn.model_selection.KFold(n_splits=num_splits, shuffle=True, random_state=42)

    for fold_num, (train_index, test_index) in enumerate(kf.split(X), 1):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]

        # Define the parameter distributions for RandomizedSearchCV
        param_distributions = {
            'n_estimators': [100, 200],
            'learning_rate': [0.01, 0.05, 0.1],
            'max_depth': [3, 5, 7],
            'min_samples_leaf': [1, 5],
            'max_features': [None, 'sqrt'],
            'subsample': [1.0, 0.8]
        }

        # Adjust max_features based on the flag
        if shorten_features_for_speed:
            param_distributions['max_features'] = ['sqrt']
        else:
            param_distributions['max_features'] = [None, 'sqrt']

        gbr = sklearn.ensemble.GradientBoostingRegressor(random_state=42 + fold_num)

        random_search = sklearn.model_selection.RandomizedSearchCV(
            estimator=gbr,
            param_distributions=param_distributions,
            n_iter=n_iter,
            cv=3,
            scoring='neg_mean_squared_error',
            n_jobs=-1,
            verbose=0,
            random_state=42,
            return_train_score=True
        )

        # Fit the RandomizedSearchCV on the training data
        random_search.fit(X_train, y_train)
        best_model = random_search.best_estimator_

        LOGGER.info(f"Fold {fold_num} Best parameters found:", random_search.best_params_)

        y_pred_test = best_model.predict(X_test)
        y_pred_cv[test_index] = y_pred_test

        models.append(best_model)
        test_set_predictions.append((y_test, y_pred_test))

        # Evaluate performance
        fold_mse = np.mean((y_test - y_pred_test) ** 2)
        LOGGER.info(f"Fold {fold_num} MSE: {fold_mse}")

        correlation = np.corrcoef(y_test, y_pred_test)[0, 1]
        LOGGER.info(f"Overall correlation in fold {fold_num}: {correlation}")

    return models, test_set_predictions, y_pred_cv




def train_fast_gradient_boosting(X, y, shorten_features_for_speed, num_splits=3, n_iter=10):
    LOGGER.info("Starting train_fast_gradient_boosting, no parallel processing, n=1")
    # Take the absolute value of y to predict magnitudes only
    y = np.abs(y)

    # Apply feature selection if needed
    if shorten_features_for_speed:
        # For example, select top 50% features
        from sklearn.feature_selection import SelectKBest, f_regression
        selector = SelectKBest(score_func=f_regression, k=int(X.shape[1] * 0.5))
        X = selector.fit_transform(X, y)
    else:
        selector = None

    # Define the parameter distributions for RandomizedSearchCV
    param_distributions = {
        'learning_rate': [0.01, 0.05, 0.1],
        'max_depth': [3, 5],
        'min_samples_leaf': [1, 5],
        'max_iter': [100, 200],
        'max_bins': [255],
        'max_leaf_nodes': [31],
        'l2_regularization': [0, 0.1],
        'early_stopping': [True],
        # Removed 'max_features' since it's not supported in older versions
    }

    hgb = HistGradientBoostingRegressor(random_state=42)

    random_search = RandomizedSearchCV(
        estimator=hgb,
        param_distributions=param_distributions,
        n_iter=n_iter,
        cv=3,
        scoring='neg_mean_squared_error',
        n_jobs=1,
        verbose=0,
        random_state=42,
        return_train_score=True,
        error_score='raise'  # Raise an error if a fit fails
    )
    LOGGER.info("Starting RandomizedSearchCV")
    # Fit RandomizedSearchCV on the entire dataset

    random_search.fit(X, y)
    best_params = random_search.best_params_
    LOGGER.info(f"Best parameters found: {best_params}")

    # Now do cross-validation using best_params
    models = []
    test_set_predictions = []
    y_pred_cv = np.zeros_like(y)
    kf = KFold(n_splits=num_splits, shuffle=True, random_state=42)
    LOGGER.info("Starting cross-validation")
    for fold_num, (train_index, test_index) in enumerate(kf.split(X), 1):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]

        model = HistGradientBoostingRegressor(**best_params, random_state=42 + fold_num)

        model.fit(X_train, y_train)

        # Compute permutation feature importances
        perm_importance = permutation_importance(
            model, X_test, y_test, n_repeats=5, random_state=42, scoring='neg_mean_squared_error'
        )
        feature_importances = perm_importance.importances_mean

        # Handle negative importances
        feature_importances = np.maximum(feature_importances, 0)

        # Normalize importances
        if np.sum(feature_importances) > 0:
            feature_importances = feature_importances / np.sum(feature_importances)
        else:
            feature_importances = np.zeros_like(feature_importances)
        model.feature_importances_ = feature_importances

        y_pred_test = model.predict(X_test)
        y_pred_cv[test_index] = y_pred_test

        models.append(model)
        test_set_predictions.append((y_test, y_pred_test))

        # Evaluate performance
        fold_mse = np.mean((y_test - y_pred_test) ** 2)
        LOGGER.info(f"Fold {fold_num} MSE: {fold_mse}")
        correlation = np.corrcoef(y_test, y_pred_test)[0, 1]
        LOGGER.info(f"Overall correlation in fold {fold_num}: {correlation}")

    return models, test_set_predictions, y_pred_cv





def train_random_forest_simple(X, y, shorten_features_for_speed, num_splits=5):
    kf = sklearn.model_selection.KFold(n_splits=num_splits, shuffle=True, random_state=42)
    models = []
    test_set_predictions = []
    y = np.abs(y)


    y_pred_cv = np.zeros_like(y)
    for fold_num, (train_index, test_index) in enumerate(kf.split(X), 1):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]

        model = sklearn.ensemble.RandomForestRegressor(
            n_estimators=100,
            random_state=42 + fold_num,
            n_jobs=-1,
            max_features='sqrt',
            max_depth=5,
            min_samples_leaf=5,
            bootstrap=True,
            max_samples=0.8
        )

        model.fit(X_train, y_train)
        models.append(model)

        y_pred_test = model.predict(X_test)
        y_pred_cv[test_index] = y_pred_test

        test_set_predictions.append((y_test, y_pred_test))

        # Evaluate performance
        fold_mse = np.mean((y_test - y_pred_test) ** 2)
        LOGGER.info(f"Fold {fold_num} MSE: {fold_mse}")

        correlation = np.corrcoef(y_test, y_pred_test)[0, 1]
        LOGGER.info(f"Overall correlation in fold {fold_num}: {correlation}")

    return models, test_set_predictions, y_pred_cv
