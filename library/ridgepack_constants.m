function [rhoi,rhos,rhow,delrho,g,eincr,hincr,minthick,maxthick]=ridgepack_constants;

% RIDGEPACK_CONSTANT - Gives constants of sea ice for Ridgepack library
%
% This function supplies physical constants to the Ridgepack physics library
% 
% OUTPUT:
%
% rhoi     - density of sea ice (kg/m^3)
% rhos     - density of snow (kg/m^3)
% rhow     - density of sea water (kg/m^3)
% delrho   - difference of rhow minus rhoi (kg/m^3)
% g        - acceleration due to gravity (m/s^2)
% eincr    - increment of epsilon and phi on alpha-hat plane (dimensionless)
% hincr    - thickness increment on zeta-hat plane (m)
% minthick - minimum thickness on zeta-hat plane (m)
% maxthick - maximum thickness on zeta-hat plane (m)
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
g = 9.8;       

% resolution of epsilon and phi in calculating zeta hat plane (dimensionless)
eincr = 0.001;

% log thickness resolution on zeta-hat plane (m)
hincr = (log10(10)-log10(0.01))/1000;

% minimum thickness on zeta-hat plane trajectory plane (m)
minthick = 0.01;

% maximum thickness on zeta-hat plane trajectory plane (m)
maxthick = 20;

