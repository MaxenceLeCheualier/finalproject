class NeedlemanWunsch:
    """
    Class to perform global sequence alignment using the Needleman-Wunsch algorithm.
    """

    def __init__(self, sequences : dict, match_score: int = 1, mismatch_penalty: int = -1, gap_penalty: int = -2):
        """
        Initialize the Needleman-Wunsch algorithm with scoring parameters.

        Args:
            sequences (dict): The sequences to be aligned encoded using de fasta reader class.
            match_score (int): The score for a match (default: 1)
            mismatch_penalty (int): The penalty for a mismatch (default: -1)
            gap_penalty (int): The penalty for a gap (default: -2)
        """
        self.sequences = sequences
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.gap_penalty = gap_penalty
    
    def score_matrix(self, seq : str, var : str) :
        """
        Build the score matrix for the given sequences, using the Needleman-Wunsch algorithm.
        Args:
            seq (str): The first sequence to be aligned.
            var (str): The second sequence to be aligned.
        """

        long_seq = len(seq) + 1
        long_var = len(var) + 1

        matrix = [[0 for i in range(long_seq)] for _ in range(long_var)]

        # Initialize the first row and column of the score matrix
        for i in range(long_seq):
            matrix[0][i] = self.gap_penalty * i
        for j in range(long_var):
            matrix[j][0] = self.gap_penalty * j
        
        
    

