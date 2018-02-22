function ridgepack_dump(ncfile)

% ridgepack_dump - Lists the contents of a netcdf file.
%
% function ridgepack_dump(ncfile)
%
% Input:
% ncfile - character string of the full netcdf file name.
%
% This is a dummy function to match the name of the native
% netcdf library by Unidata that has a utility called ridgepack_dump.
%
% Andrew Roberts, Naval Postgraduate School, March 2018  (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

ncdisp(ncfile);

if debug; disp(['...Leaving ',mfilename]); end


