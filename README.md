# Ridgepack 

Ridgepack is a MATLAB sea ice model analysis and development package designed
as part of the Regional Arctic System Model (RASM) project, funded by the 
Department of Energy, Office of Naval Research, and National Science Foundation.
It is used primarily to analyse output from RASM and the Community Earth System
Model (CESM) sea ice components, and for model developments leading to contributions
to the CICE Consortium sea ice model (https://github.com/CICE-Consortium).  

## Contents

### Library

This is a library of MATLAB functions to read netCDF data output from CICE and 
make vigorous use of the metadata within the output to aid model analysis. The library also contains development libraries for new physics being implemented in the RASM
sea ice component.  There are currently four libraries:

datastructures - ingesting, twisting and turning and outputting CICE data

graphics       - publication quality graphics functions for CICE

infrastructure - basic infrastructure for Ridgepack to operate smoothly

morphology     - library of sea ice thickness distribution functions


### Cases

The cases area of Ridgepack is designed for MATLAB scripts that make use of the 
above function library. Currently there is one set of case scripts available to 
accompany a publication under review:

JAMES\_2018\_Variational\_Ridging - Scripts used to demonstrate the functionality of the sea ice morphology library described in: Roberts, A.F., E.C. Hunke, S.M. Kamal, W.H. Lipscomb, C. Horvat, W. Maslowski (2018), Variational Method for Sea Ice Ridging in Earth System Models, Part I: Theory, *submitted to J. Adv. Model Earth Sy.*.

A second case will be released in 2018:

RASM\_Sea_Ice_Toolbox - Scripts used to analyse the dynamics and thermodynamics of the coupled RASM sea ice component.


### startup.m

This is a sample MATLAB startup file to accompany Ridgepack, and includes examples of setting several environment variables relevant to this MATLAB utility. It is anticipated that MATLAB will be run on a Unix, Linux or Mac OSX platform when using this package. 


## Authors
Andrew Roberts, Naval Postgraduate School, April 2018 (afrobert@nps.edu)

Case reviewer for JAMES\_2018\_Variational\_Ridging: Samy Kamal, Naval Postgraduate School, May 2018



