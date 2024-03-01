import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



class QualityScoreNormalizer():
    def __init__(self,  results_df, example_node):
        self.results_df = results_df
        self._example_node = example_node
        self._score_is_ml_derived = hasattr(self._example_node, "predscore")

        self._normalize_quality_score()
    

    def _normalize_quality_score(self):
        scores = self.results_df['quality_score'].values

        if self._score_is_ml_derived:
            scores = abs(scores) #the machine learning score is symmetric around 0, so we first take the absolute value


        scores_min_max_scaled = self._perform_min_max_scaling(scores) # Min Max scaling to have the scores between 0 and 1

        if self._score_is_ml_derived:
            scores_min_max_scaled = 1 - scores_min_max_scaled # The machine learning score is a "badness" score, so we want to reverse the order of the scores
        
        # Assigning the normalized scores back to the DataFrame
        self.results_df['quality_score'] = scores_min_max_scaled

        return self.results_df


    def _perform_min_max_scaling(self, scores_standardized):
        min_val = np.min(scores_standardized)
        max_val = np.max(scores_standardized)
        scores_min_max_scaled = (scores_standardized - min_val) / (max_val - min_val)
        return scores_min_max_scaled