function [nc_type]=ridgepack_invert(mclass)

% ridgepack_invert - converts a matlab class to an equivalent netcdf type class
%
% function [nc_type]=ridgepack_invert(mclass)
%
% INPUT:
%
% class  - text string corresponding to the input:
%          int32, int8, char, int16, single, or double
%
%
% OUTPUT:
% 
% nc_type - text string of the netcdf type: 
%           NC_INT, NC_BYTE, NC_CHAR, NC_SHORT, NC_FLOAT or ...
% 	    NC_DOUBLE (respectively)
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% match the netcdf data type with the matlab
if strcmp(mclass,'int32')
        nc_type='NC_INT';
elseif strcmp(mclass,'int8')
        nc_type='NC_BYTE';
elseif strcmp(mclass,'char')
        nc_type='NC_CHAR';
elseif strcmp(mclass,'int16')
        nc_type='NC_SHORT';
elseif strcmp(mclass,'single')
        nc_type='NC_FLOAT';
elseif strcmp(mclass,'double')
        nc_type='NC_DOUBLE';
else
        error(['problem with type assignment for ',mclass]);
end

if debug; disp(['...Leaving ',mfilename]); end

