class NeedlemanWunsch:
    """
    Class to perform global sequence alignment using the Needleman-Wunsch algorithm.
    """

    def __init__(self, match_score: int = 1, mismatch_penalty: int = -1, gap_penalty: int = -2):
        """
        Initialize the Needleman-Wunsch algorithm with scoring parameters.

        Args:
            match_score (int): The score for a match (default: 1)
            mismatch_penalty (int): The penalty for a mismatch (default: -1)
            gap_penalty (int): The penalty for a gap (default: -2)
        """
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.gap_penalty = gap_penalty