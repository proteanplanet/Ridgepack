function [rhoi,rhos,rhow,delrho,ghat]=ridgepack_constants

% RIDGEPACK_CONSTANTS - Gives constants of sea ice for the Ridgepack library
%
% function [rhoi,rhos,rhow,delrho,g,eincr,hincr,minthick,maxthick]=ridgepack_constants
%
% This function supplies physical constants to the Ridgepack physics library
% 
% OUTPUT:
%
% rhoi     - density of sea ice (kg/m^3)
% rhos     - density of snow (kg/m^3)
% rhow     - density of sea water (kg/m^3)
% delrho   - difference of rhow minus rhoi (kg/m^3)
% ghat     - acceleration due to gravity (m/s^2)
%
% Ridgepack Version 1.0.
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)


% density of ice (kg/m^3)
rhoi = 917.0;  

% density of snow (kg/m^3)
rhos = 330.0;  

% density of seawater (kg/m^3)
rhow = 1026.0; 

% difference density of seawater-sea ice (kg/m^3)
delrho = rhow - rhoi; 

% acceleration due to gravity (m/s^2)
ghat = 9.8;       

