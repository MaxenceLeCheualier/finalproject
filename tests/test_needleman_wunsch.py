import pytest
from finalproject.needleman_wunsch import NeedlemanWunsch


def test_needleman_simple_match():

    sequences = {
        "seq1": ["AAAA", "AAAA"]
    }

    nw = NeedlemanWunsch(sequences)
    aligned = nw.align(sequences)

    assert aligned["seq1"].original == "AAAA"
    assert aligned["seq1"].variant == "AAAA"

def test_needleman_one_mismatch():

    sequences = {
        "seq1": ["AAAA", "TAAA"]
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

    sequences = { "seq1" : ["ATGCTGAT", "ATGCTGT"]}

    nw = NeedlemanWunsch(sequences)
    aligned = nw.align(sequences)

    aligned_obj = aligned["seq1"]
    seq = aligned_obj.original
    var = aligned_obj.variant 

    assert seq == "ATGCTGAT"
    assert var == "ATGCTG-T"

def test_needleman_empty_variant():
    
    sequences = { "seq1" : ["ATGCTGAT", ""]}

    nw = NeedlemanWunsch(sequences)
    aligned = nw.align(sequences)

    aligned_obj = aligned["seq1"]
    seq = aligned_obj.original
    var = aligned_obj.variant 

    assert seq == "ATGCTGAT"
    assert var == "--------"

def test_needleman_simple_test_multiple_seq():

    sequences = {
        "seq1": [
            "ACGTGCTAGTACCGTATCGTAGCTAGTAC",
            "ACGTGCTAGTACGTATCGTAGCTAGTAC"
        ],
        "seq2": [
            "GATCGTAGCTAGCTAGGCTATCGTAGCTA",
            "GATCGTAGCTAGGCTATCGTAGCTA"
        ],
        "seq3": [
            "TTGACCGTAGCTAGCTAACGTTAGCTAGCT",
            "TTGACGTAGCTAGCTAACGTTAGCTAGCT"
        ],
    }
    
    nw = NeedlemanWunsch(sequences)
    aligned = nw.align(sequences)

       
    assert aligned["seq1"].original == "ACGTGCTAGTACCGTATCGTAGCTAGTAC"
    assert aligned["seq1"].variant == "ACGTGCTAGTA-CGTATCGTAGCTAGTAC"

    assert aligned["seq2"].original == "GATCGTAGCTAGCTAGGCTATCGTAGCTA"
    assert aligned["seq2"].variant == "GATC---G-TAGCTAGGCTATCGTAGCTA"

    assert aligned["seq3"].original == "TTGACCGTAGCTAGCTAACGTTAGCTAGCT"
    assert aligned["seq3"].variant == "TTGA-CGTAGCTAGCTAACGTTAGCTAGCT"
