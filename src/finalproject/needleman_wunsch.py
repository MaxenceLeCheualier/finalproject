class NeedlemanWunsch:
    """
    Class to perform global sequence alignment using the Needleman-Wunsch algorithm.
    """

    def __init__(self, sequences : dict[str, list[str]], match_score: int = 1, mismatch_penalty: int = -1, gap_penalty: int = -2) -> None:
        """
        Initialize the Needleman-Wunsch algorithm with scoring parameters.

        Args:
            sequences (dict[str, list[str]]): The sequences to be aligned encoded using de fasta reader class.
            match_score (int): The score for a match (default: 1)
            mismatch_penalty (int): The penalty for a mismatch (default: -1)
            gap_penalty (int): The penalty for a gap (default: -2)
        """
        self._sequences = sequences
        self._match_score = match_score
        self._mismatch_penalty = mismatch_penalty
        self._gap_penalty = gap_penalty
    
    @property
    def sequences(self):
        return self._sequences

    @property
    def match_score(self):
        return self._match_score

    @property
    def mismatch_penalty(self):
        return self._mismatch_penalty

    @property
    def gap_penalty(self):
        return self._gap_penalty
    
    def score_matrix(self, seq : str, var : str) -> list[list[int]] :
        """
        Build the score matrix for the given sequences, using the Needleman-Wunsch algorithm.
        Args:
            seq (str): The first sequence to be aligned.
            var (str): The second sequence to be aligned.
        Returns : 
            matrix (list[list[int]]) : The score matrix
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
    
    def traceback(self, matrix : list[list[int]], seq : str, var : str) -> tuple[str, str] :
        """
        Use the score matrix defined above to perform tht traceback and find an optimal global alignement between the two sequences
        Args:
            matrix (list[list[int]]): The score matrix calculated by the score_matrix method.
            seq (str): The first sequence to be aligned.
            var (str): The second sequence to be aligned.
        Returns :
            tuple[str, str] : the tuple of the two sequences aligned
        """
        aligned_seq = ""
        aligned_var = ""

        i = len(var) 
        j = len(seq)

        while i > 0 and j > 0 : 

            current_score = matrix[i][j]
            up_score = matrix[i-1][j] + self.gap_penalty
            left_score = matrix[i][j-1] + self.gap_penalty
            diagonal_score = matrix[i-1][j-1] + self.match_score if seq[j-1] == var[i-1] else matrix[i-1][j-1] + self.mismatch_penalty
            
    
            #We give the priority to the diagonal movement, which means we need to distinguish the match and the mismatch cases.

            if current_score == diagonal_score :
                aligned_seq += seq[j-1]
                aligned_var += var[i-1]
                i -= 1
                j -= 1
            
            #The second priority is the up movement (which is an arbitrary choice)

            elif current_score == up_score :
                aligned_seq += "-"
                aligned_var += var[i-1]
                i -= 1
                
            #Here comes the left movement

            elif current_score == left_score  :
                aligned_seq += seq[j-1]
                aligned_var += "-"
                j -= 1 

        while i > 0:
            #In this cas we can only go up

            aligned_seq += "-"
            aligned_var += var[i-1]
            i -= 1

        while j > 0:
            #In this case we can only to the left 

            aligned_seq += seq[j-1]
            aligned_var += "-"
            j -= 1

        #We started from the end, then we have to reverse the string"
        return aligned_seq[::-1], aligned_var[::-1]


    
    def align(self, sequences : dict[str, list[str]]) -> dict[str, list[str]]:
        """
        Perform the global sequence alignment for all the sequences in the input dictionary and return a new dictionary with the aligned sequences.

        Args:
            sequences (dict): The sequences to be aligned encoded using de fasta reader class.
        Returns:
            aligned_sequences (dict): A dictionary where keys are sequence names and values are lists of aligned sequences (str)
        """

        aligned_sequences = {}

        for name, seqs in sequences.items():
            seq = seqs[0]
            var = seqs[1]

            matrix = self.score_matrix(seq, var)
            aligned_seq, aligned_var = self.traceback(matrix, seq, var)

            aligned_sequences[name] = [aligned_seq, aligned_var]

        return aligned_sequences
    

