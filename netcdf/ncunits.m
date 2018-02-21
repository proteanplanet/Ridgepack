function [units]=ncunits(nc,name)

% NCUNITS - Convert units to tex format for matlab from netcdf structure
%
% function [units]=ncunits(nc,name)
%
% Input:
% nc   - netcdf structure (see ncstruct for more details) or
%        Matlab timeseries object
% name - field of netcdf structure for which units are required
%        (if useing a netcdf structure)
%
% Ouput:
% units - character variable of tex formatted units for plotting
%
% Written by Andrew Roberts 2012
% Department of Oceanography, Naval Postgraduate School
%
% $Id$

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

interpreter=get(0,'DefaultTextInterpreter');

if isa(nc,'timeseries')
	units=nc.DataInfo.Unit;
elseif isstruct(nc)
	if isfield(nc,name)
		if isfield(nc.(name),'units') && ~strcmp(name,'time')
			units=nc.(name).units;

		else
			units='';
		end
	else
		error([name,' is not a field of nc'])
	end
else
	error('nc is not a netcdf structure or matlab timeseries')
end

% unitless
if strcmp(units,'unitless'); units=''; end
if strcmp(units,'none'); units=''; end

% longitude
if strcmp(units,'degree_east'); units='{^{\circ}}E'; end
if strcmp(units,'degrees_east'); units='{^{\circ}}E'; end
if strcmp(units,'Degrees_East'); units='{^{\circ}}E'; end
if strcmp(units,'Degrees_east'); units='{^{\circ}}E'; end

% latitude
if strcmp(units,'degree_north'); units='{^{\circ}}N'; end
if strcmp(units,'degrees_north'); units='{^{\circ}}N'; end
if strcmp(units,'Degrees_North'); units='{^{\circ}}N'; end
if strcmp(units,'Degrees_north'); units='{^{\circ}}N'; end

% temperature
if strcmp(units,'Degrees_Celsius'); units='{^{\circ}}C'; end
if strcmp(units,'degrees_Celsius'); units='{^{\circ}}C'; end
if strcmp(units,'degrees_celsius'); units='{^{\circ}}C'; end
if strcmp(units,'Degrees celsius'); units='{^{\circ}}C'; end
if strcmp(units,'Degrees Celsius'); units='{^{\circ}}C'; end
if strcmp(units,'degrees Celsius'); units='{^{\circ}}C'; end
if strcmp(units,'Degree_Celsius'); units='{^{\circ}}C'; end
if strcmp(units,'degree_Celsius'); units='{^{\circ}}C'; end
if strcmp(units,'Degree Celsius'); units='{^{\circ}}C'; end
if strcmp(units,'degree Celsius'); units='{^{\circ}}C'; end
if strcmp(units,'degrees_C'); units='{^{\circ}}C'; end
if strcmp(units,'degrees C'); units='{^{\circ}}C'; end
if strcmp(units,'degree_C'); units='{^{\circ}}C'; end
if strcmp(units,'degree C'); units='{^{\circ}}C'; end
if strcmp(units,'degC'); units='{^{\circ}}C'; end
if strcmp(units,'C'); units='{^{\circ}}C'; end
if strcmp(units,'degK'); units='K'; end

% distance
if strcmp(units,'meters'); units='m'; end

% volume
if strcmp(units,'km^3'); units='km^{3}'; end

% rate
if strcmp(units,'/s'); units='s^{-1}'; end
if strcmp(units,'/day'); units='day^{-1}'; end
if strcmpi(interpreter,'latex')
 if strcmp(units,'%/day'); units='\% day^{-1}'; end
else
 if strcmp(units,'%/day'); units='% day^{-1}'; end
end

% fraction
if strcmpi(interpreter,'latex')
 if strcmp(units,'%'); units='\%'; end
end

% speed
if strcmp(units,'m s**-1'); units=' m s^{-1}'; end
if strcmp(units,'m/s'); units='m s^{-1}'; end
if strcmp(units,'m s-1'); units='m s^{-1}'; end
if strcmp(units,'cm/s'); units='cm s^{-1}'; end
if strcmp(units,'centimeter/s'); units='cm s^{-1}'; end

% stress
if strcmp(units,'N/m'); units='N m^{-1}'; end
if strcmp(units,'N/m^2'); units='N m^{-2}'; end
if strcmp(units,'N m-2'); units='N m^{-2}'; end

% flux
if strcmp(units,'W/m^2'); units='W m^{-2}'; end
if strcmp(units,'W/m^-2'); units='W m^{-2}'; end
if strcmp(units,'W m-2'); units='W m^{-2}'; end
if strcmp(units,'kg/s'); units='kg s^{-1}'; end
if strcmp(units,'cm/day'); units='cm day^{-1}'; end
if strcmp(units,'kg/m^2/day'); units='kg m^{-2} day^{-1}'; end

% concentration or fraction
if strcmp(units,'1'); units=''; end


if debug; disp(['...Leaving ',mfilename]); end
