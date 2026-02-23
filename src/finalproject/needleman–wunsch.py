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
    
    def score_matrix(self, seq : str, var : str) -> list :
        """
        Build the score matrix for the given sequences, using the Needleman-Wunsch algorithm.
        Args:
            seq (str): The first sequence to be aligned.
            var (str): The second sequence to be aligned.
        """

        long_seq = len(seq) + 1
        long_var = len(var) + 1

        matrix = [[0 for i in range(long_seq)] for _ in range(long_var)]

        #  Firt row and first column initialized with the sum of gap penalties
        # We chose to use the first row to refer to the variant sequence and the first column to refer to the non-variant sequence

        for j in range(long_seq):
            matrix[0][j] = self.gap_penalty * j
        for i in range(long_var):
            matrix[i][0] = self.gap_penalty * i

        for i in range (1, long_var):       
            for j in range(1, long_seq):
                if seq[j-1] == var[i-1]:     #There is a match diagonal score grows by the match score 
                    diagonal_score = matrix[i-1][j-1] + self.match_score
                else :                       #There is a mismatch diagonal score grows by the mismatch penalty
                    diagonal_score = matrix[i-1][j-1] + self.mismatch_penalty
                
                #The two other scores are the up and left scores, they grow by the gap penalty
                up_score = matrix[i-1][j] + self.gap_penalty   
                left_score = matrix[i][j-1] + self.gap_penalty
                
                #The score of the cell is the maximum of the three scores calculated above
                matrix[i][j] = max(diagonal_score, up_score, left_score)
        
        return matrix
    
    def traceback(self, matrix : list, seq : str, var : str) :
        """
        Use the score matrix defined above to perform tht traceback and find an optimal global alignement between the two sequences
        Args:
            matrix (list): The score matrix calculated by the score_matrix method.
            seq (str): The first sequence to be aligned.
            var (str): The second sequence to be aligned.
        """
        aligned_seq = ""
        aligned_var = ""

        i = len(var) + 1
        j = len(seq) + 1

        while i > 0 and j > 0: 

            diagonal_score = matrix[i-1][j-1]
            up_score = matrix[i-1][j]
            left_score = matrix[i][j-1]
            max_score = max(diagonal_score, up_score, left_score)
            move_done = False

            if max_score == diagonal_score:   #We are in a match or mismatch case
                aligned_seq = seq[j-1] + aligned_seq
                aligned_var = var[i-1] + aligned_var
                i -= 1
                j -= 1
                move_done = True
            
            elif max_score == up_score and not move_done:      #We are in a gap case for the non-variant sequence
                aligned_seq = "-" + aligned_seq
                aligned_var = var[i-1] + aligned_var
                i -= 1
            
            else :                           #We are in a gap case for the variant sequence
                aligned_seq = seq[j-1] + aligned_seq
                aligned_var = "-" + aligned_var
                j -= 1
                move_done = True

            



    

