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
            sequences (dict[str, list[str]]) : A dictionnary where keys are sequence names and values are listes of sequences (str)
        """

        sequences = {}
        name = None 

        variant = False

        with open(self.filepath, 'r') as file:
            for line in file:
                line = line.strip()

                #Empty lines and comment lines are ignored 
                if not line or line.startswith(';'): 
                    continue 

                #Lines starting with '>' are sequence names (variant and non-variant)
                if line.startswith('>'):
                    # We start a new sequence, it means we need to reinitialize the variable of the variant we have just treated
                    if name is not None and variant:
                        variant = False
                        name = None
                        
                    if name is None:
                        name = line[1:].strip()
                        sequences[name] = Sequence(name)

                        # The sequence bellow is a variant of the previous one
                    elif "_var" in line:
                        variant = True

                    #End of variant sequence case


                
                # Lines that are sequences 
                else : 
                    #If name is initialized and variant is False, we are in a non-variant sequence
                    if name is not None and not variant:
                        sequences[name].original += line
                    
                    elif name is not None and variant: 
                    #If name is initialized and variant is True, we are in a variant sequence
                        sequences[name].variant += line
        
        return sequences

                