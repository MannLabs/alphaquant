import numpy as np
import pytest
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor, HistGradientBoostingRegressor


from alphaquant.classify.training_functions import (
	train_random_forest_with_grid_search,
	train_gradient_boosting_with_grid_search,
	train_gradient_boosting_with_random_search,
	train_fast_gradient_boosting,
	train_random_forest_simple
)

def check_model_outputs(models, test_set_predictions, y_pred_cv, y_true, expected_model_type, num_splits):
	"""
	Helper function to perform assertions on model outputs.
	"""
	# Assert the number of models equals the number of splits
	assert len(models) == num_splits, f"Expected {num_splits} models, got {len(models)}"
	
	# Assert each model is of the expected type
	for model in models:
		assert isinstance(model, expected_model_type), f"Model is not of type {expected_model_type}"
	
	# Assert the number of test set predictions equals the number of splits
	assert len(test_set_predictions) == num_splits, f"Expected {num_splits} test set predictions, got {len(test_set_predictions)}"
	
	# Assert predictions have the same shape as the true values
	assert y_pred_cv.shape == y_true.shape, f"Predictions shape {y_pred_cv.shape} does not match true values shape {y_true.shape}"
	
	# Assert total number of test samples equals total number of samples
	total_test_samples = sum(len(y_test) for y_test, _ in test_set_predictions)
	assert total_test_samples == len(y_true), f"Total test samples {total_test_samples} does not match total samples {len(y_true)}"
	
	# Assert predictions are not all zeros
	assert not np.all(y_pred_cv == 0), "All predictions are zero, model did not make meaningful predictions"
	
	# Compute correlation between predictions and true values
	corr = np.corrcoef(y_true, y_pred_cv)[0, 1]
	
	# Assert correlation is finite
	assert np.isfinite(corr), "Correlation is not finite, possibly due to NaNs in predictions or true values"
	
	# Since data is random, correlation should be close to zero
	assert -0.15 < corr < 0.15, f"Correlation {corr} not close to zero as expected for random data"
	
	# Compute mean squared error
	mse = np.mean((y_true - y_pred_cv) ** 2)
	
	# Assert mean squared error is finite
	assert np.isfinite(mse), "Mean squared error is not finite, possibly due to NaNs in predictions or true values"
	
	# Since data is random, MSE should be close to variance of y
	y_variance = np.var(y_true)
	assert mse >= 0.7 * y_variance and mse <= 1.3 * y_variance, f"MSE {mse} should be close to variance {y_variance} for random data"


	
	# Optionally, print the MSE and correlation for debugging
	print(f"MSE: {mse}, Variance: {y_variance}, Correlation: {corr}")

@pytest.mark.parametrize("train_func, expected_model_type", [
	(train_random_forest_with_grid_search, RandomForestRegressor),
	(train_gradient_boosting_with_grid_search, GradientBoostingRegressor),
	(train_gradient_boosting_with_random_search, GradientBoostingRegressor),
	(train_fast_gradient_boosting, HistGradientBoostingRegressor),
	(train_random_forest_simple, RandomForestRegressor)
])
def test_model_training(train_func, expected_model_type):
	np.random.seed(42)
	X = np.random.rand(200, 5)
	y = np.random.rand(200)
	shorten_features_for_speed = False
	num_splits = 5

	# Some functions take the absolute value of y
	y_transformed = y.copy()
	if train_func.__name__ in [
		'train_random_forest_with_grid_search',
		'train_gradient_boosting_with_grid_search',
		'train_gradient_boosting_with_random_search',
		'train_fast_gradient_boosting',
		'train_random_forest_simple'
	]:
		y_transformed = np.abs(y)

	# Adjust parameters for functions that require n_iter
	if train_func in [train_gradient_boosting_with_random_search, train_fast_gradient_boosting]:
		n_iter = 5  # Reduced for testing purposes
		models, test_set_predictions, y_pred_cv = train_func(
			X, y, shorten_features_for_speed, num_splits=num_splits, n_iter=n_iter
		)
	else:
		models, test_set_predictions, y_pred_cv = train_func(
			X, y, shorten_features_for_speed, num_splits=num_splits
		)

	# Perform assertions
	check_model_outputs(models, test_set_predictions, y_pred_cv, y_transformed, expected_model_type, num_splits)
