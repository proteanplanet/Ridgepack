# Ridgepack 

Ridgepack is a MATLAB sea ice model analysis and development package designed
as part of the Regional Arctic System Model (RASM) project, funded by the 
Department of Energy, Office of Naval Research, and National Science Foundation.
It has been used to analyse output from RASM, as well as the Community Earth System
Model (CESM) and Energy Exascale Earth System Model (E3SM) sea ice components, 
and for model developments leading to contributions to the CICE Consortium sea 
ice model (https://github.com/CICE-Consortium). Ridgepack libraries are readily 
applicable to other models that use CICE or its column package, called Icepack, 
including the Discrete Element Model of Sea Ice (DEMSI).

## Contents

### Libraries

These libraries form the backbone of Ridgepack to read netCDF data output 
from CICE, MPAS-Seaice and observational datasets and make vigorous use of 
the metadata within the output to aid model analysis.  The libraries include
extensive graphics and analysis function, as well as development libraries for 
new physics being implemented in Icepack.  There are four libraries:

| Library | Description |
| --- | --- |
| datastructures | Ingests, twists and turns netcdf data |
| graphics | Publication quality graphics functions |
| morphology | Library of sea ice thickness distribution functions for the paper Roberts, A.F., E.C. Hunke, S.M. Kamal, W.H. Lipscomb, C. Horvat, W. Maslowski (2019), A Variational Method for Sea Ice Ridging in Earth System Models, *J. Adv. Model Earth Sy.* |
| infrastructure | Basic functions for Ridgepack to operate smoothly |

### Cases

The cases area of Ridgepack is for MATLAB scripts that make use of the 
above function library: 

| Case | Description |
| --- | --- |
| JAMES\_2019\_VarRidging | Scripts used to demonstrate the functionality of the sea ice morphology library described in: Roberts, A.F., E.C. Hunke, S.M. Kamal, W.H. Lipscomb, C. Horvat, W. Maslowski (2019), Variational Method for Sea Ice Ridging in Earth System Models, *J. Adv. Model Earth Sy.*. |
| RASM\_Sea\_Ice\_Toolbox | Scripts used to analyse coupled RASM model components |

### Startup

startup.m is a sample MATLAB startup file to accompany Ridgepack, and includes examples of setting several environment variables relevant to this MATLAB utility. It is anticipated that MATLAB will be run on a Unix, Linux or Mac OSX platform when using this package. 

## Contact

Andrew Roberts: afroberts@lanl.gov

## History 

Author: Andrew Roberts, Naval Postgraduate School, April 2018

Reviewer: Samy Kamal, Naval Postgraduate School, May 2018 

Update: Andrew Roberts, Los Alamos National Laboratory, December 2018


## Required software

Designed to be used with MATLAB 2018b or a newer version, and for one script, Mathematica 11.3


