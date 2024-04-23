# Galaxy-Assembly-Bias-in-SDSS-CMASS


Paper Files
---------------------------------------------
Clustering Measurement
xi_su_tree.c - A C program to calculate the 2-point correlation function for a lightcone distribution of galaxies. Pair counts are performed with a tree structure, only computing pairs for galaxies in the current and neighboring cells to expedite computation.
fileprep_tree.py - A Python program to build a table that will localize a galaxy in a particular cell determine by transforming the galaxy coordinates from redshift-space to Cartesian and dividing up the (X,Y,Z)-space into cells/bins dependent on the clustering scales of interest.
