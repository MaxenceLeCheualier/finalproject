import pytest
from finalproject.needleman_wunsch import NeedlemanWunsch


def test_needleman_simple_match():

    sequences = {
        "seq1": ["AAAA", "AAAA"]
    }

    nw = NeedlemanWunsch(sequences)
    aligned = nw.align()

    assert aligned["seq1"][0] == "AAAA"
    assert aligned["seq1"][1] == "AAAA"
