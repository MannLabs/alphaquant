import pandas as pd
import numpy as np



class QualityScoreNormalizer():
    def __init__(self,  results_df, example_node):
        self.results_df = results_df
        self._example_node = example_node

        self._normalize_quality_score()
        self._invert_quality_score_if_ml()

    def _normalize_quality_score(self):
        scores = self.results_df['quality_score'].values

        # Z-Score Normalization
        mean = np.mean(scores)
        std_dev = np.std(scores)
        scores_standardized = (scores - mean) / std_dev

        # Min-Max Scaling
        min_val = np.min(scores_standardized)
        max_val = np.max(scores_standardized)
        scores_min_max_scaled = (scores_standardized - min_val) / (max_val - min_val)

        # Assigning the normalized scores back to the DataFrame
        self.results_df['quality_score'] = scores_min_max_scaled

        return self.results_df

    def _invert_quality_score_if_ml(self):
        if hasattr(self._example_node, "predscore"):
            self.results_df["quality_score"] = 1 - self.results_df["quality_score"]
        return self.results_df
