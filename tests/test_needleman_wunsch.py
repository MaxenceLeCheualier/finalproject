import pytest
from finalproject.needleman_wunsch import NeedlemanWunsch


def test_needleman_simple_match():

    sequences = {
        "seq1": ["AAAA", "AAAA"]
    }

    nw = NeedlemanWunsch(sequences)
    aligned = nw.align(sequences)

    assert aligned["seq1"][0] == "AAAA"
    assert aligned["seq1"][1] == "AAAA"

def test_needleman_one_mismatch():

    sequences = {
        "seq1": ["AAAA", "TAAA"]
    }

    nw = NeedlemanWunsch(sequences)
    aligned = nw.align(sequences)

    seq, var = aligned["seq1"]

    assert len(seq) == len(var)
    assert seq == "AAAA"
    assert var == "TAAA"

def test_needleman_one_gap():

    sequences = { "seq1" : ["ATGCTGAT", "ATGCTGT"]}

    nw = NeedlemanWunsch(sequences)
    aligned = nw.align(sequences)

    seq, var = aligned["seq1"]

    assert seq == "ATGCTGAT"
    assert var == "ATGCTG-T"

def test_needleman_empty_variant():
    
    sequences = { "seq1" : ["ATGCTGAT", ""]}

    nw = NeedlemanWunsch(sequences)
    aligned = nw.align(sequences)

    seq, var = aligned["seq1"]

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

       
    assert aligned["seq1"][0] == "ACGTGCTAGTACCGTATCGTAGCTAGTAC"
    assert aligned["seq1"][1] == "ACGTGCTAGTA-CGTATCGTAGCTAGTAC"

    assert aligned["seq2"][0] == "GATCGTAGCTAGCTAGGCTATCGTAGCTA"
    assert aligned["seq2"][1] == "GATC---G-TAGCTAGGCTATCGTAGCTA"

    assert aligned["seq3"][0] == "TTGACCGTAGCTAGCTAACGTTAGCTAGCT"
    assert aligned["seq3"][1] == "TTGA-CGTAGCTAGCTAACGTTAGCTAGCT"
