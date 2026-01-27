class FastaReader:
    """""
    A class to read a fasta file and return a dictionary of sequences
    """

    def __init__(self, filepath : str) -> None:
        """""
        Initialize the FastaReader with the path to the fasta file

        Args:
            filepath (str) : The path to the fasta file 
        """ 
        self.filepath = filepath

    def read(self) -> dict(str, str):
        """
        Read the fasta file and return a dictionnary of sequences, distinguished by their names,and their variants if any.

        Args :
            None
        Returns :
            sequences (dict) : A dictionnary where keys are sequence names and values are listes of sequences (str)
        """

        sequences = {}
        name = None 
        current_sequence = []

        with open(self.filepath, 'r') as file:
            for line in file:
                line = line.strip()

                if line.startswith('>'):             # The next line is a sequence or a variant of it 
                    if name is None : 
                        name = line[1:]
                        sequence = "".join(current_sequence)