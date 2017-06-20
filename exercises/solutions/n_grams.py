#!/usr/bin/python3
import string
from urllib.request import urlopen, urlretrieve
from statistics import mode, median
from collections import Counter
import sys

class Text(object):
    """
    A class for analysing texts. Includes functions for calculating various different statistics, including finding n-grams and common words.
    """

    def __init__(self, filename):
        self.text = self.read_text(filename)

    @staticmethod
    def read_text(filename):
        """
        Read the text from the file 'filename' and return lowercase.

        Parameters
        ----------
        filename : string
            Name and path of file to be analysed.

        Returns
        -------
        string list
            List of words in text, with all punctuation stripped and all text set to lowercase
        """

        if '.txt' not in filename:
            raise ValueError('Input file must be a .txt file!')

        if filename[:4] == "http": # website
            website = urlopen(filename)
            # slicing to get rid of project gutenberg preamble/license
            txt = website.read().decode('UTF-8').lower()[800:-19500]
        else: # text file
            f = open(filename)
            txt = f.read().lower()

        # strip punctuation. Have switched hyphens to a capital letter and back so that they do not get removed.
        translator = txt.maketrans('--', '  ')
        txt = txt.translate(translator)
        translator = txt.maketrans('-', 'A')
        txt = txt.translate(translator)
        translator = txt.maketrans("\n\r\t", ' '*3)
        txt = txt.translate(translator)
        translator = txt.maketrans('', '', string.punctuation + "'`’‘”“")
        txt = txt.translate(translator)
        translator = txt.maketrans('A', '-')
        txt = txt.translate(translator).split(' ')

        return [s for s in txt if s !='']

    def find_ngrams(self, n):
        """
        Find n-grams of Text object. Returns dictionary of n-grams.

        Parameters
        ----------
        n : integer
            Length of n-grams to construct

        Returns
        -------
        output : dictionary
            Dictionary of n-grams in the text. The keys are the n-grams, the values the frequency they appear in the text.
        """

        output = {}

        for i in range(len(self.text)-n+1):
            s = ' '.join(self.text[i:i+n])
            # if s is not already in dictionary, set value to 0
            output.setdefault(s, 0)
            output[s] += 1
        return output

    def average_word_length(self):
        """
        Return mean, median and mode word length. Includes only words (i.e. no numbers) in calculation.

        Returns
        -------
        float tuple
            Mean, median and mode word length
        """
        len_words_only = [len(s) if s.isalpha() else 0 for s in self.text]
        if (len_words_only == 0):
            print('Input file contains no words.')
            return 0, 0, 0
        else:
            return sum(len_words_only) / len(len_words_only), median(len_words_only), mode(len_words_only)

    def word_count(self):
        """
        Returns number of words in the text.

        Returns
        -------
        integer
            number of words in text
        """
        return len(self.text)

    def longest_words(self, n=10):
        """
        Return the n longest words in the text.

        Parameters
        ---------
        n : integer, optional
            Number of words to return. Default is 10.

        Returns
        -------
        string list
            N longest words in text (sorted with longest first)
        """
        return sorted(set(self.text), key=len, reverse=True)[:n]

    def common_words(self, n=10):
        """
        Return the n most common words in the text. Only looks for words with 3 or more letters and ignores a given set of very common words.

        Parameters
        ----------
        n : integer, option
            Number of words to return. Default is 10

        Returns
        -------
        string list
            Most common words in text (with most common first)
        """
        # remove some really common words
        ignore = ['a', 'i', 'it', 'the', 'and', 'in', 'he', 'she', 'to', 'at', 'of', 'that', 'as', 'is', 'his', 'my', 'for', 'was', 'me', 'we', 'be', 'on', 'so']
        filtered = [s for s in self.text if s not in ignore and len(s) >=3]
        dat = Counter(filtered)
        return dat.most_common(n)

    def text_report(self):
        """
        Print a report of the text, giving information about various different metrics.
        """

        word_count = self.word_count()

        print("\nThere are {} words in the text.".format(word_count))
        mean, median, mode = self.average_word_length()

        print("\nMean, median and mode word length is {}, {}, {}.".format(mean, median, mode))

        if word_count < 10:
            print("\nLongest words:")
        else:
            print("\n10 longest words:")
        for s in self.longest_words():
            print(s)

        print("\nMost common words:")
        for s in self.common_words():
            print("{} x {}".format(s[1], s[0]))

        longest_grams = []

        # find n_longest n-grams
        n_longest = 10
        # strongly doubt that there will be n-grams longer than 50
        for i in range(min(50, word_count), 1, -1):
            if len(longest_grams) >= n_longest:
                break
            grams = self.find_ngrams(i)
            grams_list = sorted(grams, key=grams.get, reverse=True)

            for g in grams_list:
                if grams[g] > 4:
                    # do not want to include n-grams which are substrings of longer n-grams
                    substring = False
                    for s in longest_grams:
                        if g in s[1]:
                            substring = True
                            break
                    if not substring:
                        longest_grams.append([grams[g], g])

        print("\nLongest n-grams:")
        for g in longest_grams:
            print("{} x {}".format(g[0], g[1]))
        print('\n')


if __name__ == "__main__":
    # If a runtime argument is given, use that as the filename. Otherwise, use the uncommented file given below
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        #filename = "http://www.gutenberg.org/files/11/11-0.txt" # alice
        #filename = "http://www.gutenberg.org/ebooks/345.txt.utf-8" # dracula
        #filename = "http://www.gutenberg.org/ebooks/1661.txt.utf-8" # sherlock
        filename = "../the_raven.txt" # poe

    txt = Text(filename)

    txt.text_report()
