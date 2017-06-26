CPP = g++    # Define variable (C++ compiler)

report.pdf: report.tex library.bib figures/fig.pdf
	pdflatex report
	bibtex report
	pdflatex report
	pdflatex report

figures/fig.pdf: py/plot.py data/data.csv
	cd py; python plot.py

data/data.csv: cpp/main.cpp
	${CPP} cpp/main.cpp -o main    # Uses variable defined above for compiler command
	./main
