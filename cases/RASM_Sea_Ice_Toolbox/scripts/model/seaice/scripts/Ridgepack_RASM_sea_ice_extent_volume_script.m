% Ridgepack_RASM_sea_ice_extent_volume_script - summary graph of RASM results
%
% This script prepares an analysis plot of RASM sea ice output equivalent
% to the example ridgepack_example_seaice_volume_graph_1980_2010_0.png.
% The script uses a wide array of observational and model data to construct the plot:
%
% 1) ICESat data from JPL
% 2) PIOMAS data from University of Washington
% 3) Ice area and extent from NSIDC NOAA CDR
% 4) RASM ice area (aice), thickness (hi) and snow depth (hs)
%
% The RASM data is assumed to be in the form of a single file for each 
% variable, organized as a timeseries using the 'rasm_cice5_series.bash'
% script on the chosen supercomputer used to integrate RASM. That script
% is located under the bash directory of the RASM Sea Ice Toolbox
% in Ridgepack. This produces files with a filename, for e.g. 'hi', as
% c6G05b.cice.h.hi.nc for the simulation 'c6G05b'.
%
% Andrew Roberts, Naval Postgraduate School, June 2018 (afrobert@nps.edu)


% set graph bounds
mintime=datenum(1980,1,1);
maxtime=datenum(2010,1,1);

% set the cases required to be analyzed, and their short notation, or 'quicknames'
rasmcases={'c6G05b','R1009RBRceap01a','R2100aRBRcaaa01a'};
quicknames={'c6G05b','RASM 1.1','RASM 2.1'};

pubdir='/Users/aroberts/science/publications/2018_RASM/Figures'

% set minimum model thickness for which concentration, extent and area are calculated
%minthick=[0.25]
minthick=[0.0]

% determine if publication-quality is required (removes some notation, and header)
%pub=true;
pub=false;

% run the function
Ridgepack_RASM_sea_ice_extent_volume(rasmcases,quicknames,minthick,mintime,maxtime,pub,pubdir)

