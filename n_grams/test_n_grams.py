#!/usr/bin/python3
import pytest
from n_grams import Text
import numpy

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

def test_average_word_length():
    with pytest.raises(ZeroDivisionError):
        Text("blank.txt").average_word_length()
    assert numpy.allclose(Text("repeat.txt").average_word_length(), (3.536363636363636, 3.0, 3))
    assert numpy.allclose(Text("short.txt").average_word_length(), (3.2, 4, 4))

def test_word_count():
    assert Text("blank.txt").word_count() == 0
    assert Text("repeat.txt").word_count() == 7480
    assert Text("short.txt").word_count() == 5

def test_longest_words():
    assert len(Text("blank.txt").longest_words()) == 0
    assert Text("repeat.txt").longest_words()[0] == "jumblies"
    assert Text("short.txt").longest_words()[0] == "short"

def test_common_words():
    assert len(Text("blank.txt").common_words()) == 0
    assert Text("repeat.txt").common_words()[:3] == [('sieve', 476), ('they', 408), ('sea', 340)]
    # check no word appears more than once in this file
    for word, freq in Text("short.txt").common_words():
        assert freq == 1
