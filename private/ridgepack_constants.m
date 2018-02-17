function [rhoi,rhos,rhow,delrho,g]=ridge_constants;

% This function declares the parameter space of variational ridging
% Written by Andrew Roberts, March 2018

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
epres = 0.001;

% maximum thickness of zeta hat trajectory plane (m)
maxthick = 50

