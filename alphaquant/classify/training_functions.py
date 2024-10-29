import alphaquant.config.config as aqconfig

import sklearn.ensemble
import sklearn.linear_model
import sklearn.impute
import sklearn.metrics
import sklearn.model_selection
import xgboost as xgb

from sklearn.model_selection import KFold, RandomizedSearchCV
from sklearn.feature_selection import SelectFromModel

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

def train_xgboost(X, y, shorten_features_for_speed, num_splits=5):
    # Take the absolute value of y to predict magnitudes only
    y = np.abs(y)
    models = []
    test_set_predictions = []
    y_pred_cv = np.zeros_like(y)
    
    kf = sklearn.model_selection.KFold(n_splits=num_splits, shuffle=True, random_state=42)
    
    for fold_num, (train_index, test_index) in enumerate(kf.split(X), 1):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        
        xgb_reg = xgb.XGBRegressor(
            n_estimators=1000,
            learning_rate=0.05,
            max_depth=6,
            subsample=0.8,
            colsample_bytree=0.8 if shorten_features_for_speed else 1.0,
            random_state=42 + fold_num,
            n_jobs=-1,
            eval_metric='rmse',     # Moved eval_metric here
            early_stopping_rounds=50
        )
        
        xgb_reg.fit(
            X_train, 
            y_train,
            eval_set=[(X_test, y_test)],
            verbose=False
        )
        
        y_pred_test = xgb_reg.predict(X_test)
        y_pred_cv[test_index] = y_pred_test
        models.append(xgb_reg)
        test_set_predictions.append((y_test, y_pred_test))
        
        # Evaluate performance
        fold_mse = np.mean((y_test - y_pred_test) ** 2)
        (f"Fold {fold_num} MSE: {fold_mse}")
        correlation = np.corrcoef(y_test, y_pred_test)[0, 1]
        LOGGER.info(f"Overall correlation in fold {fold_num}: {correlation}")
    
    return models, test_set_predictions, y_pred_cv



def train_xgboost_optimized(X, y, shorten_features_for_speed, num_splits=5, n_iter=20):
    # Take the absolute value of y to predict magnitudes only
    y = np.abs(y)
    
    # Feature selection
    if shorten_features_for_speed:
        model = xgb.XGBRegressor()
        model.fit(X, y)
        selector = SelectFromModel(model, prefit=True, threshold='median')
        X = selector.transform(X)
        LOGGER.info(f"Reduced features to {X.shape[1]} using feature selection.")
    
    models = []
    test_set_predictions = []
    y_pred_cv = np.zeros_like(y)
    
    kf = KFold(n_splits=num_splits, shuffle=True, random_state=42)
    
    param_distributions = {
        'n_estimators': [100, 200, 500],
        'learning_rate': [0.05, 0.1],
        'max_depth': [3, 4, 5],
        'subsample': [0.8],
        'colsample_bytree': [0.6, 0.8],
        'gamma': [0],
        'reg_alpha': [0],
        'reg_lambda': [1]
    }
    
    for fold_num, (train_index, test_index) in enumerate(kf.split(X), 1):
        LOGGER.info(f"Starting fold {fold_num}")
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        
        xgb_reg = xgb.XGBRegressor(
            objective='reg:squarederror',
            random_state=42 + fold_num,
            n_jobs=-1,
            tree_method='hist'
        )
        
        randomized_search = RandomizedSearchCV(
            estimator=xgb_reg,
            param_distributions=param_distributions,
            n_iter=n_iter,
            scoring='neg_mean_squared_error',
            cv=3,
            random_state=42,
            verbose=0,
            n_jobs=-1
        )
        
        randomized_search.fit(X_train, y_train)
        best_model = randomized_search.best_estimator_
        LOGGER.info(f"Fold {fold_num} Best parameters: {randomized_search.best_params_}")
        
        # Early stopping
        best_model.set_params(early_stopping_rounds=30, eval_metric='rmse')
        best_model.fit(
            X_train, y_train,
            eval_set=[(X_test, y_test)],
            verbose=False
        )
        
        y_pred_test = best_model.predict(X_test)
        y_pred_cv[test_index] = y_pred_test
        models.append(best_model)
        test_set_predictions.append((y_test, y_pred_test))
        
        # Evaluate performance
        fold_mse = np.mean((y_test - y_pred_test) ** 2)
        LOGGER.info(f"Fold {fold_num} MSE: {fold_mse}")
        correlation = np.corrcoef(y_test, y_pred_test)[0, 1]
        LOGGER.info(f"Fold {fold_num} Correlation: {correlation}")
        LOGGER.info(f"Fold {fold_num} Best iteration: {best_model.best_iteration}")
    
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
