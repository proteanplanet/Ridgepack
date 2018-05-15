function [h]=ridgepack_astroconstants

% ridgepack_astroconstants - Provides a handle with astronomical and physical constants
%
% function [h]=ridgepack_astroconstants
%
% This function provides a handle with astronomical and physical constants 
% for use in metocean calculations. 
%
% OUTPUT:
%
% h - structure with the form:
%     h.name.descr='description of constant';
%     h.name.units='string';
%     h.name.const=(double precision value);
%
%     where 
%
%     name is the constant in the structure, and multiple names are possible
%     within the one structure.
%
% Current constants:
%
% omega  - Angular velocity of Earth's rotation
% r      - mean radius of the Earth
% ghat   - acceleration due to gravity
% rhoa   - density of air at standard temperature and pressure
% rhoi   - density of sea ice
% rhos   - density of snow on sea ice
% rhow   - density of seawater
% 
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

global debug;

if debug; disp(['Entering ',mfilename,'...']); end

h=[];

h=ridgepack_addc(h,'omega','Angular velocity of Earth''s rotation',...
                   'per second',0.00007292115);

h=ridgepack_addc(h,'r','mean radius of the Earth','m',earthRadius);

h=ridgepack_addc(h,'ghat','acceleration due to gravity','m s^{-2}',9.8);

h=ridgepack_addc(h,'rhoa','density of air at standard temperature and pressure',...
                   'kg m^{-3}',1.225);

h=ridgepack_addc(h,'rhoi','density of sea ice','kg m^{-3}',917.0);

h=ridgepack_addc(h,'rhos','density of snow on sea ice','kg m^{-3}',330.0);

h=ridgepack_addc(h,'rhow','density of seawater','kg m^{-3}',1026.0);

if debug; disp(['...Leaving ',mfilename]); end

function h=ridgepack_addc(h,name,description,units,constant)

 % This function generates the constant structure for ridgepack_astroconstants

 if isfield(h,name)
	error('reassigning constant name already allocated')
 else
	h.(name).descr=description;
	h.(name).units=units;
	h.(name).const=constant;
 end


