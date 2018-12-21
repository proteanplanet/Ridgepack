# Ridgepack Case: JAMES\_2019\_VarRidging

This directory contains scripts used to generate figures 3 to 16, S1 and S2 in:

Roberts, A.F., E.C. Hunke, S.M. Kamal, W.H. Lipscomb, C. Horvat, W. Maslowski (2019),
Variational Method for Sea Ice Ridging in Earth System Models, *J. Adv. Model Earth Sy.*.

These scripts also test other aspects of variational ridging introduced in this paper.

## Synopsis

| Script | Description |
| -- | -- |
| ridgepack\_JAMES\_figureX.m | Produce figures 3, 4, 5, 6, 7, 10, 11, 12, 14, & 15 |
| ridgepack\_JAMES\_ridgegraph.m | Produces figures 8, 9, S1, S2 |
| ridgepack\_JAMES\_integrateg.m | Calculates redistribution in figure 16 |
| ridgepack\_grplot.m            | Compares numeric & analytic solutions of gR(h,phi) |
| ridgepack\_JAMES\_figure13.nb  | Mathematica script for generating thickness distribution for a single ridge in figure 13 |
              

## Description

Graphical output from these scripts is written to a directory called "output" adjacent to this scripts directory, except for the Mathematica script, which writes graphical output to this directory.  Each script includes a commented section explaining the code, which, in the case of MATLAB code, can be read using the help command:  

help ridgepack/cases/JAMES\_2019\_VarRidging

## Version/Cases

Ridgepack 1.0.1/JAMES\_2019\_VarRidging

## Contact

Andrew Roberts, afroberts@lanl.gov 

## History 

Author: Andrew Roberts, Naval Postgraduate School, April 2018 
Reviewer: Samy Kamal, Naval Postgraduate School, May 2018
Update: Andrew Roberts, Los Alamos National Laboratory, December 2018

## Required software

Designed to be used with MATLAB 2018b or a newer version, and Mathematica 11.3
