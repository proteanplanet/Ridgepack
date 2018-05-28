# Ridgepack Case: JAMES\_2018\_Variational\_Ridging

This directory contains scripts used to generate figures 3 to 16, S1 and S2 in:

Roberts, A.F., E.C. Hunke, S.M. Kamal, W.H. Lipscomb, C. Horvat, W. Maslowski (2018),
Variational Method for Sea Ice Ridging in Earth System Models, Part I: Theory, *submitted to J. Adv. Model Earth Sy.*.

These scripts also test other aspects of variational ridging introduced in this paper.

## Synopsis

ridgepack\_JAMES\_figureX.m    - Produce figures 3, 4, 5, 6, 7, 10, 11, 12, 14, & 15 

ridgepack\_JAMES\_ridgegraph.m - Produces figures 8, 9, S1, S2

ridgepack\_JAMES\_integrateg.m - Calculates redistribution in figure 16

ridgepack\_grplot.m            - Compares numeric & analytic solutions of gR(h,phi)

## Description

Graphical output from these scripts is written to a directory called "output" adjacent to this scripts directory.  Each script includes a commented section explaining the code, which, in the case of MATLAB code, can be read using the help command:  

help ridgepack/cases/JAMES\_2018\_Variational\_Ridging/scripts

## Authors

Andrew Roberts, Naval Postgraduate School, April 2018 (afrobert@nps.edu)

Reviewed by Samy Kamal, Naval Postgraduate School, May 2018

