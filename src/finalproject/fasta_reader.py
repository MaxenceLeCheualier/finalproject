class FastaReader:
    def __init__(self, filepath):
        self.filepath = filepath

    def read(self):

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