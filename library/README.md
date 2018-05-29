# Ridgepack Libraries

This directory includes a series of MATLAB libraries designed for development
and analysis of the sea ice components of coupled Earth System Models. The
package was developed using the CICE Consortium sea ice model (https://github.com/CICE-Consortium) as the marine cryosphere component in the Regional Arctic System Model (RASM; http://www.oc.nps.edu/NAME/RASM\_PhaseIII.html).  Ridgepack libraries are 
readily applicable to other models that use CICE, including the Community Earth
System Model (CESM; http://www.cesm.ucar.edu). 

## Contents

### datastructures

A MATLAB library for cloning netCDF data structures and using metadata in CICE 
data files to ease processing of model output, and to help understand the climate 
and evolution of the polar components of RASM and CESM.  For an introduction
to this library, use the MATLAB 'help' command, starting with understanding
'nc structures' used in this package:  'help datastructures/ridgepack\_struct'.  All
help files may be viewed using 'help datastructures' from within this library
directory.


### graphics

A MATLAB library for producing publication-quality graphics for Earth System Models. 
All code documentation may be obtained from within MATLAB by typing 'help graphics' from within this library directory.  Functions within this directory fix a number of problems with vanilla release MATLAB graphics functions when used to analyze the data on a geoid. 


### infrastructure

A repository of supporting infrastructure, currently confined to a debugger that accompanies this package. Documentation is available from within MATLAB by typing
'help infrastructure'.


### morphology

A library of sea ice ridging functions used to calculate statistics of ridges in the pack.  Functions in this library are applicable to interpreting measurements of 
sea ice topography from underneath and above sea ice.  However, the main purpose of this library is to act as a testbed for a new variational ridging method being implemented in RASM. 'help morphology' from within this library directory in MATLAB will
reveal the functions available.  A description of the science underwriting this 
morphology library is in a paper under review:

Roberts, A.F., E.C. Hunke, S.M. Kamal, W.H. Lipscomb, C. Horvat, W. Maslowski (2018),
Variational Method for Sea Ice Ridging in Earth System Models, Part I: Theory, *submitted to J. Adv. Model Earth Sy.*.

This code is published to satisfy the American Geophysical Union publication policy.  

## Author
Andrew Roberts, Naval Postgraduate School, April 2018 (afrobert@nps.edu)


## Required software

Designed to be used with MATLAB 2018a or a newer version


