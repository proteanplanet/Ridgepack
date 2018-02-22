function [mclass]=ridgepack_onvert(nc_type)

% ridgepack_onvert - converts the netcdf type into the equivalent matlab class
%
% function [mclass]=ridgepack_onvert(nc_type)
%
% INPUT:
% 
% nc_type - text string of the netcdf type: 
%           NC_INT, NC_BYTE, NC_CHAR, NC_SHORT, NC_FLOAT or NC_DOUBLE
%
%
% OUTPUT:
%
% mclass  - text string corresponding to the input:
%           int32, int8, char, int16, single, or double (respectively)
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% match the netcdf data type with the matlab
if strcmp(nc_type,'NC_INT')
        mclass='int32';
elseif strcmp(nc_type,'NC_BYTE')
        mclass='int8';
elseif strcmp(nc_type,'NC_CHAR')
        mclass='char';
elseif strcmp(nc_type,'NC_SHORT')
        mclass='int16';
elseif strcmp(nc_type,'NC_FLOAT')
        mclass='single';
elseif strcmp(nc_type,'NC_DOUBLE')
        mclass='double';
else
        error(['problem with type assignment for ',nc_type]);
end

if debug; disp(['...Leaving ',mfilename]); end

