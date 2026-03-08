from finalproject.needleman_wunsch import NeedlemanWunsch
from finalproject.sequence import Sequence

def test_needleman_simple_match():

    sequences = {
    "seq1": Sequence(
        name="seq1",
        original="AAAA",
        variant="AAAA"
    )
    }

    nw = NeedlemanWunsch(sequences)
    aligned = nw.align(sequences)

    assert aligned["seq1"].original == "AAAA"
    assert aligned["seq1"].variant == "AAAA"

def test_needleman_one_mismatch():

    sequences = {
    "seq1": Sequence(
        name="seq1",
        original="AAAA",
        variant="TAAA"
    )
    }
    

    nw = NeedlemanWunsch(sequences)
    aligned = nw.align(sequences)

    aligned_obj = aligned["seq1"]
    seq = aligned_obj.original
    var = aligned_obj.variant 

    assert len(seq) == len(var)
    assert seq == "AAAA"
    assert var == "TAAA"

def test_needleman_one_gap():

    sequences = {
    "seq1": Sequence(
        name="seq1",
        original="ATGCTGAT",
        variant="ATGCTGT"
    )
    }

    nw = NeedlemanWunsch(sequences)
    aligned = nw.align(sequences)

    aligned_obj = aligned["seq1"]
    seq = aligned_obj.original
    var = aligned_obj.variant 

    assert seq == "ATGCTGAT"
    assert var == "ATGCTG-T"

def test_needleman_empty_variant():
    
    sequences = {
    "seq1": Sequence(name="seq1", original="ATGCTGAT", variant="")
    }

    nw = NeedlemanWunsch(sequences)
    aligned = nw.align(sequences)

    aligned_obj = aligned["seq1"]
    seq = aligned_obj.original
    var = aligned_obj.variant 

    assert seq == "ATGCTGAT"
    assert var == "--------"

def test_needleman_simple_test_multiple_seq():

    sequences = {
        "seq1": Sequence(
        name="seq1",
        original="ACGTGCTAGTACCGTATCGTAGCTAGTAC",
        variant="ACGTGCTAGTACGTATCGTAGCTAGTAC"
    ),
        "seq2": Sequence(
        name="seq2",
        original="GATCGTAGCTAGCTAGGCTATCGTAGCTA",
        variant="GATCGTAGCTAGGCTATCGTAGCTA"
    ),
        "seq3": Sequence(
        name="seq3",
        original="TTGACCGTAGCTAGCTAACGTTAGCTAGCT",
        variant="TTGACGTAGCTAGCTAACGTTAGCTAGCT"
    ),
    }
    
    nw = NeedlemanWunsch(sequences)
    aligned = nw.align(sequences)

       
    assert aligned["seq1"].original == "ACGTGCTAGTACCGTATCGTAGCTAGTAC"
    assert aligned["seq1"].variant == "ACGTGCTAGTA-CGTATCGTAGCTAGTAC"

    assert aligned["seq2"].original == "GATCGTAGCTAGCTAGGCTATCGTAGCTA"
    assert aligned["seq2"].variant == "GATC---G-TAGCTAGGCTATCGTAGCTA"

    assert aligned["seq3"].original == "TTGACCGTAGCTAGCTAACGTTAGCTAGCT"
    assert aligned["seq3"].variant == "TTGA-CGTAGCTAGCTAACGTTAGCTAGCT"
