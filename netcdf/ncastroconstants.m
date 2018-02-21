function [h]=ncastroconstants

% NCASTROCONSTANTS - Provides a handle with astronomical and physical constants
%
% function [h]=ncastroconstants
%
% This function provides a handle with astronomical and physical constants 
% for use in metocean calculations. 
%
% Output:
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
% Written by Andrew Roberts 2012
% Department of Oceanography, Naval Postgraduate School
%
% $Id$

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

h=[];
h=ncaddc(h,'omega','Angular velocity of earth''s rotation','per second',0.00007292115);
h=ncaddc(h,'rhoa','density of air at standard temperature and pressure','kg m^-3',1.225);
h=ncaddc(h,'r','mean radius of the Earth','m',earthRadius);

if debug; disp(['...Leaving ',mfilename]); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h=ncaddc(h,name,description,units,constant);

% This function generates the constant structure for ncastroconstants

if isfield(h,name)
	error('reassigning constant name already allocated')
else
	h.(name).descr=description;
	h.(name).units=units;
	h.(name).const=constant;

end

