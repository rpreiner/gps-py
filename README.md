# gps-py
Python demos for the publication  
*Gaussian-Product Subdivision Surfaces*  
Website: https://cgv.tugraz.at/preiner/gps  
DOI: 10.1145/3306346.3323026

Copyright 2020 Reinhold Preiner


## Demos

This repository contains two examples from the paper, demoing the functionality of Gaussian-Product subdivision (GPS) and Gaussian inference in covariance meshes.

* **tweety-inference** (Fig. 7)  
Demoes the automatic enrichment of a given input triangle mesh with covariances based on Eq. 16, using the tweety model. The resulting covmesh leads to a GPS surface with sharper features than ordinary with loop subdivision.

* **cones** (Fig. 10)
Shows the different GPS variants of a fixed cone control mesh under varying Gaussian covariances at the apex.


## Files

The demos are provided in two versions:
* as plain python scripts (subdir 'py'). Here the subdivison results are written to ./data directory.
* as ineractive Jupyter Notebooks (subdir 'jupyter-nb'). Here the subdivision results are directly visualized using meshplot.


## Installation/Dependencies

* The scripts depend on the python bindings of libigl for mesh file i/o as well as linear subdivision (https://libigl.github.io/libigl-python-bindings/). Install via 
  
  conda install -c conda-forge igl

* To run the Jupyter scripts and perform interactive mesh visualization, you also need to install meshplot (https://skoch9.github.io/meshplot/):

  conda install -c conda-forge meshplot 

