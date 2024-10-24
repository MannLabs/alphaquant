import pandas as pd

class QualityScoreNormalizer:
    def __init__(self, results_df: pd.DataFrame):
        """
        The QualityScoreNormalizer converts the arbitrary-scaled quality scores to a normalized scale between 0 and 1, where 1
        is the best score and 0 is the worst score. 
        This normalizer assumes that higher input scores are better.
        The original scores are converted to ranks, which are then normalized to the 0-1 range.

        Args:
        results_df (pd.DataFrame): DataFrame containing the quality scores, either 'ml_score' or 'consistency_score'.
        Raises:
        ValueError: If neither 'ml_score' nor 'consistency_score' is present in the DataFrame.
        """

        self.results_df = results_df
        if "ml_score" in self.results_df.columns:
            self._normalize_quality_score('ml_score')
        elif "consistency_score" in self.results_df.columns:
            self._normalize_quality_score('consistency_score')
        else:
            raise ValueError("Quality score not recognized. Please provide a valid quality score.")

    def _normalize_quality_score(self, score_column):
        ranks = self._perform_rank_normalization(self.results_df[score_column], higher_is_better=True)
        self.results_df['quality_score'] = ranks
        self.results_df = self.results_df.drop(columns=[score_column])
    
    @staticmethod
    def _perform_rank_normalization(scores_series : pd.Series, higher_is_better: bool):
        ranks = scores_series.rank(method='average', ascending=higher_is_better).values #'average' method to handle ties
        normalized_ranks = ranks / len(ranks) # Normalize ranks to be between 0 and 1
        return normalized_ranks
