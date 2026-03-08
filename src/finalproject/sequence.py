class Sequence:
    """
    Class representing a sequence and its variant.
    """

    def __init__(self, name: str, original: str = "", variant: str = "") -> None:
        """
        Initialize a Sequence object.

        Args:
            name (str): The name of the sequence.
            original (str): The main sequence. Defaults to "".
            variant (str): The variant sequence. Defaults to "".
        """
        self._name = name
        self._original = original
        self._variant = variant

    @property
    def name(self) -> str:
        return self._name

    @property
    def original(self) -> str:
        return self._original
    
    @original.setter
    def original(self, value : str) -> None : 
        self._original = value
    

    @property
    def variant(self) -> str:
        return self._variant
    
    @variant.setter
    def variant(self, value: str) -> None:
        self._variant = value

