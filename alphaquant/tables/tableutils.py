import pandas as pd

class QualityScoreNormalizer:
    def __init__(self, results_df: pd.DataFrame):
        """
        The QualityScoreNormalizer converts the arbitrary-scaled quality scores to a normalized scale between 0 and 1, where 1
        is the best score and 0 is the worst score. It determines whether the quality score is
        machine learning derived ('predscore') or consistency-based ('consistency_score'), and normalizes the scores accordingly.
        For 'predscore', lower values are better. For 'consistency_score', higher values are better.
        Args:
        results_df (pd.DataFrame): DataFrame containing the quality scores, either 'predscore' or 'consistency_score'.
        Raises:
        ValueError: If neither 'predscore' nor 'consistency_score' is present in the DataFrame.
        """
        self.results_df = results_df
        if "predscore" in self.results_df.columns:
            self._normalize_quality_score_ml()
        elif "consistency_score" in self.results_df.columns:
            self._normalize_quality_score_consistency()
        else:
            raise ValueError("Quality score not recognized. Please provide a valid quality score.")

    def _normalize_quality_score_ml(self):
        self._normalize_quality_score('predscore', higher_is_better=True)  # in the ML case, lower is better

    def _normalize_quality_score_consistency(self):
        self._normalize_quality_score('consistency_score', higher_is_better=False)  # in the consistency case, higher is better

    def _normalize_quality_score(self, score_column, higher_is_better):
        ranks = self._perform_rank_normalization(self.results_df[score_column], higher_is_better=higher_is_better)
        self.results_df['quality_score'] = ranks
        self.results_df = self.results_df.drop(columns=[score_column])
    
    @staticmethod
    def _perform_rank_normalization(scores_series : pd.Series, higher_is_better: bool):
        ranks = scores_series.rank(method='average', ascending=higher_is_better).values #'average' method to handle ties
        normalized_ranks = ranks / len(ranks) # Normalize ranks to be between 0 and 1
        return normalized_ranks
