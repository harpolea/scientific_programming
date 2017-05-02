# Scientific Programming

This repository contains a collection of notebooks and exercises for the [scientific programming course](http://ngcm.soton.ac.uk/summer-academy/sciprog.html) at the 2017 NGCM Summer Academy.

The course covers
- Version control
- Testing
- Continuous integration
- Code coverage
- Documentation
- Publishing code (including make, docker and advanced plotting)

The notebooks are best viewed on [nbviewer](https://nbviewer.jupyter.org/github/harpolea/scientific_programming).

The notebooks rely on a number of python libraries. These can be installed most easily by first creating a runtime environment:

    conda create -n ngcm_sciprog python=3 numpy scipy jupyter pandas matplotlib seaborn bokeh plotly Sphinx nbsphinx silly

Activate this environment:

    source activate ngcm_sciprog

Then launch jupyter:

    jupyter notebook
