# Scientific Programming

This repository contains a collection of notebooks and exercises for the [scientific programming course](http://ngcm.soton.ac.uk/summer-academy/sciprog.html) at the 2018 NGCM Summer Academy.

The course covers
- Testing
- Continuous integration
- Code coverage

The notebooks are best viewed on [nbviewer](https://nbviewer.jupyter.org/github/harpolea/scientific_programming).

The notebooks rely on a number of python libraries. These can be installed most easily by first creating a runtime environment:

    conda create -n ngcm_sciprog python=3 numpy scipy jupyter matplotlib pytest pytest-cov 

Activate this environment:

    source activate ngcm_sciprog

Then launch jupyter:

    jupyter notebook
