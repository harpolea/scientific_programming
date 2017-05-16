#!/usr/bin/python3

import n_grams
import matplotlib.pyplot as plt
import numpy as np
plt.style.use('ggplot')

def plot_longest_words(text_dict):
    fig, ax = plt.subplots(figsize=(10, 5))

    keys = text_dict.keys()

    ax.barh(np.arange(len(keys))+0.5, [len(text_dict[k].longest_words(1)[0]) for k in keys])

    ax.set_xlabel("Length of longest word")
    ax.set_ylabel("Text")
    ax.set_yticks(np.arange(len(keys))+0.9)
    ax.set_yticklabels(keys)
    plt.savefig("longest_words.png")

if __name__ == "__main__":
    texts = {'Alice': n_grams.Text("http://www.gutenberg.org/files/11/11-0.txt"),
        'Dracula': n_grams.Text("http://www.gutenberg.org/ebooks/345.txt.utf-8"),
        'Sherlock': n_grams.Text("http://www.gutenberg.org/ebooks/1661.txt.utf-8"),
        'The Raven': n_grams.Text("the_raven.txt"),
        'Wilde': n_grams.Text("wilde.txt")}

    plot_longest_words(texts)
