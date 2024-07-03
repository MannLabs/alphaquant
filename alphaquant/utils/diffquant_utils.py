import pandas as pd
import numpy as np



def find_non_outlier_indices_ipr(data, threshold=1.5, percentile_lower = 25, percentile_upper = 75):
    
    value_lower, value_upper = np.percentile(data, [percentile_lower, percentile_upper])
    iqr = value_upper - value_lower

    # Calculate the bounds for non-outliers
    cut_off = iqr * threshold
    lowest_tolerated_value = value_lower - cut_off
    highest_tolerated_value = value_upper + cut_off

    # Identify non-outlier indices
    non_outlier_indices = np.where((data >= lowest_tolerated_value) & (data <= highest_tolerated_value))[0]

    return non_outlier_indices