function ncdump(ncfile)

% NCDUMP - Lists the contents of a netcdf file.
%
% function ncdump(ncfile)
%
% Input:
% ncfile - character string of the full netcdf file name.
%
% This is a dummy function to match the name of the native
% netcdf library by Unidata that has a utility called ncdump.
%
% Written by Andrew Roberts2012
% Department of Oceanography, Naval Postgraduate School
%
% $Id$

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

ncdisp(ncfile);

if debug; disp(['...Leaving ',mfilename]); end


