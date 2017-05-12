#!/usr/bin/python3
import pytest
from n_grams import Text

test_files = ["blank.txt", "repeat.txt", "short.txt"]

def test_read_text():
    assert len(Text.read_text("blank.txt")) == 0
    assert len(Text.read_text("repeat.txt")) == 7480
    assert len(Text.read_text("short.txt")) == 5
    with pytest.raises(Exception):
        # well designed code should raise an exception here as a csv file is not
        # a suitable file type for this module
        Text.read_text("unsuitable.csv")

def test_find_ngrams():
    assert len(Text("blank.txt").find_ngrams(2)) == 0

    d = Text("repeat.txt").find_ngrams(2)
    assert len(d) == 80
    assert d["a button"] == 68
    assert d["round and"] == 68*2

    d = Text("short.txt").find_ngrams(2)
    for value in d.values():
        assert value == 1

def test_word_count():
    assert Text("blank.txt").word_count() == 0
    assert Text("repeat.txt").word_count() == 7480
    assert Text("short.txt").word_count() == 5

def test_longest_words():
    assert len(Text("blank.txt").longest_words()) == 0
    assert Text("repeat.txt").longest_words()[0] == "jumblies"
    assert Text("short.txt").longest_words()[0] == "short"
