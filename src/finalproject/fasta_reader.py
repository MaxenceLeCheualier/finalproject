from finalproject.sequence import Sequence

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
        self._filepath = filepath

    @property
    def filepath(self) -> str : 
        return self._filepath

    def read(self) -> dict[str, Sequence]:
            """
            Read the fasta file and return a dictionnary of sequences, distinguished by their names,and their variants if any.
    
            Args :
                None
            Returns :
                sequences (dict[str, Sequence]) : A dictionnary where keys are sequence names and values are listes of sequences (str)
            """

            sequences = {}
            name = None 
            variant = False

            with open(self.filepath, 'r') as file:
                for line in file:
                    line = line.strip()

                    # Empty lines and comment lines are ignored 
                    if not line or line.startswith(';'): 
                        continue 

                    # Lines starting with '>' are sequence names (variant and non-variant)
                    if line.startswith('>'):
                        # Determine if the line corresponds to a variant
                        if "_var" in line:
                            variant = True
                        # Note: We don't reset 'name' here to keep the variant linked to its original sequence
                        else:
                            # We start a new original sequence
                            variant = False
                            name = line[1:].strip()
                        
                            # If this name is not yet in sequences, create a Sequence instance
                            if name not in sequences:
                                sequences[name] = Sequence(name)

                   # Lines that are sequences 
                    else : 
                        if name is not None:
                            # If name is initialized and variant is False, we are in a non-variant sequence
                            if not variant:
                                sequences[name].original = sequences[name].original + line
                            # If name is initialized and variant is True, we are in a variant sequence
                            else: 
                                sequences[name].variant = sequences[name].variant + line

            return sequences