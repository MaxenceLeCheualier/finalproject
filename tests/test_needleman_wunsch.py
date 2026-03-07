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

def test_needleman_one_gape():

    sequences = { "seq1" : ["ATGCTGAT", "ATGCTGT"]}

    nw = NeedlemanWunsch(sequences)
    aligned = nw.align(sequences)

    seq, var = aligned["seq1"]

    assert seq == "ATGCTGAT"
    assert var == "ATGCTG-T"

