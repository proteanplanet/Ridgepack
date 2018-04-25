function [h]=ridgepack_astroconstants

% ridgepack_astroconstants - Provides a handle with astronomical and physical constants
%
% function [h]=ridgepack_astroconstants
%
% This function provides a handle with astronomical and physical constants 
% for use in metocean calculations. 
%
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
%     name is a constant in the structure, and multiple names are possible
%     within the one structure.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

h=[];
h=ridgepack_addc(h,'omega','Angular velocity of earth''s rotation',...
                   'per second',0.00007292115);
h=ridgepack_addc(h,'rhoa','density of air at standard temperature and pressure',...
                   'kg m^-3',1.225);
h=ridgepack_addc(h,'r','mean radius of the Earth','m',earthRadius);

if debug; disp(['...Leaving ',mfilename]); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h=ridgepack_addc(h,name,description,units,constant);

% This function generates the constant structure for ridgepack_astroconstants

if isfield(h,name)
	error('reassigning constant name already allocated')
else
	h.(name).descr=description;
	h.(name).units=units;
	h.(name).const=constant;

end

