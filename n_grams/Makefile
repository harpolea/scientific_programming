text_analysis.pdf: text_analysis.tex longest_words.png
	pdflatex text_analysis
	pdflatex text_analysis

longest_words.png: analyse_texts.py n_grams.py
	python analyse_texts.py

clean:
	rm -f text_analysis.aux text_analysis.log text_analysis.pdf text_analysis.synctex.gz longest_words.png
