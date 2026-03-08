from finalproject.fasta_reader import FastaReader 

def test_fasta_reader_simple(tmp_path): 
    """""
    A test to check if the fasta reader class can read simple sequences
    Args:
        tmp_path: A pytest fixture providing a temporary directory unique to the test invocation.
    """
    
    fasta_file = tmp_path / "test.fasta"
    
    fasta_file.write_text(
        ">seq1\n"
        "GTCCAAAAATTGGGGGGAGTAGATTGACCGTTCAGGGTCTCATATTTCGTGGTGCCGACA\n"
        ">seq1_var1\n"
        "GTCCAAATATTGGGGAGAGTAGATTGATCGTTCAGGGTCTCATATTTCGGTGCCGACA\n"
        
        ">seq2\n"
        "TAGCTGCTAGTCACATTATATAACTGTTATCGCAAAAACGTGTACATTTGCACAGAGATA\n"
        ">seq2_var1\n"
        "TCCGCTTCCCATTATATAACTGTTATCGCAAAAACGTGTACACTTGCACAGAGATA\n"
        
        ">seq3\n"
        "GGGGTTGAGTTGTTGGCTCGCTCCTAGCATATTCGCACATTTGATCGGAGTGAACACAAT\n"
        ">seq3_var1\n"
        "GGGTTGAGGTGTTGGCTCGCCTAGCATATTCGCACATTAGATCGGAGTGAACACAAC\n"
        
        ">seq4\n"
        "GATGTATTGCGATTCCGGTCTTTCTTTGATGGCCCTGGCCAAGGTTACAGGTATACAGCA\n"
        ">seq4_var1\n"
        "GATGCATTGCGATTCGGTCTTTCTTTGATGTCCCTGGCCAAGGTTACAGGGTATAGCA\n"

        ">seq5\n"
        "CCCTTACTTTGAGCAGCTAGGTGGACTGTCGGATTTGTGCATGCAGCCTCCTGTATTCAG\n"
        ">seq5_var1\n"
        "CCCTAACTTTGCGCAGCTAGGTGGACTGTCGGTTTGTACATACAGCCTCCTGTATTCAG\n"
    )
        
    reader = FastaReader(fasta_file)
    sequences = reader.read()

    assert sequences["seq1"].original == "GTCCAAAAATTGGGGGGAGTAGATTGACCGTTCAGGGTCTCATATTTCGTGGTGCCGACA"
    assert sequences["seq1"].variant == "GTCCAAATATTGGGGAGAGTAGATTGATCGTTCAGGGTCTCATATTTCGGTGCCGACA"
    
    assert sequences["seq2"].original  == "TAGCTGCTAGTCACATTATATAACTGTTATCGCAAAAACGTGTACATTTGCACAGAGATA"
    assert sequences["seq2"].variant == "TCCGCTTCCCATTATATAACTGTTATCGCAAAAACGTGTACACTTGCACAGAGATA"
    
    assert sequences["seq3"].original  == "GGGGTTGAGTTGTTGGCTCGCTCCTAGCATATTCGCACATTTGATCGGAGTGAACACAAT"
    assert sequences["seq3"].variant == "GGGTTGAGGTGTTGGCTCGCCTAGCATATTCGCACATTAGATCGGAGTGAACACAAC"
    
    assert sequences["seq4"].original  == "GATGTATTGCGATTCCGGTCTTTCTTTGATGGCCCTGGCCAAGGTTACAGGTATACAGCA"
    assert sequences["seq4"].variant == "GATGCATTGCGATTCGGTCTTTCTTTGATGTCCCTGGCCAAGGTTACAGGGTATAGCA"
    
    assert sequences["seq5"].original  == "CCCTTACTTTGAGCAGCTAGGTGGACTGTCGGATTTGTGCATGCAGCCTCCTGTATTCAG"
    assert sequences["seq5"].variant == "CCCTAACTTTGCGCAGCTAGGTGGACTGTCGGTTTGTACATACAGCCTCCTGTATTCAG"

def test_fasta_reader_on_two_line(tmp_path):
    """""
    A test to check if the fasta reader class can read sequences that are writtten on two lines
    Args:
        tmp_path: A pytest fixture providing a temporary directory unique to the test invocation.
    """
    
    fasta_file = tmp_path / "test_on_two_line.fasta"
    
    fasta_file.write_text(
        ">seq1\n"
        "GTCCAAAAATTGGGGGGAGTAGATTGACCGTTCAGGG\n"
        "TCTCATATTTCGTGGTGCCGACA\n"
        ">seq1_var1\n"
        "GTCCAAATATTGGGGAGAGT\n"
        "AGATTGATCGTTCAGGGTC\n"
        "TCATATTTCGGTGCCGACA\n"
    )

    reader = FastaReader(fasta_file)
    sequences = reader.read()

    assert sequences["seq1"].original  == "GTCCAAAAATTGGGGGGAGTAGATTGACCGTTCAGGGTCTCATATTTCGTGGTGCCGACA"
    assert sequences["seq1"].variant == "GTCCAAATATTGGGGAGAGTAGATTGATCGTTCAGGGTCTCATATTTCGGTGCCGACA"

def test_fasta_reader_empty_sequence(tmp_path):
    """""
    A test to check if the fasta reader class can read an empty sequence
    Args:
        tmp_path: A pytest fixture providing a temporary directory unique to the test invocation.
    """
    
    fasta_file = tmp_path / "empty.fasta"
    fasta_file.write_text("")
    
    reader = FastaReader(str(fasta_file))
    sequences = reader.read()
    assert sequences == {}

def test_fasta_reader_commented_line(tmp_path):
    """""
    A test to check if the fasta reader class can read an line that is divided by comments 
    Args:
        tmp_path: A pytest fixture providing a temporary directory unique to the test invocation.
    """
    
    fasta_file = tmp_path / "comments.fasta"
    fasta_file.write_text(
        "; Comment\n"
        "\n"
        ">seq1\n"
        "ACGT\n"
        "; Another comment\n"
        "TGCA\n"
    )
    
    reader = FastaReader(str(fasta_file))
    sequences = reader.read()
    
    assert sequences["seq1"].original  == "ACGTTGCA"
    assert sequences["seq1"].variant == ""

def test_fasta_reader_without_variant(tmp_path):
    """""
    A test to check if the fasta reader class can read an sequence without any variant
    Args:
        tmp_path: A pytest fixture providing a temporary directory unique to the test invocation.
    """

    fasta_file = tmp_path / "without_variant.fasta"
    fasta_file.write_text(
        ">seq1\n"
        "ACGT\n"
    )

    reader = FastaReader(str(fasta_file))
    sequences = reader.read()

    assert sequences["seq1"].original  == "ACGT"
    assert sequences["seq1"].variant == ""