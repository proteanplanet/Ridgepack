function [nc]=ridgepack_icetideoutput(nc)

% ridgepack_icetideoutput - Prepares a netcdf structure from an ice-tide output netcdf file
%
% function [nc]=ridgepack_icetideoutput(nc)
%
% This function adds a turning angle to vectors for the rotated sterographic 
% ice-tide grid.
%
% INPUT:
% 
% nc - netcdf structure prepared by ridgepack_clone
%
%
% OUTPUT:
%
% nc - netcdf strucure with certain changes applicable to ice-tide output
%
% Written Andrew Roberts 2012
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if debug; disp('Reconfiguring nc structure from ice tide model input'); end

if ~isfield(nc,'turn') & isfield(nc,'longitude')

 nc.turn.data=(90+nc.longitude.data-50);
 nc.turn.units='degrees';
 nc.turn.long_name='vector turning angle to lat-long (u,v) grid from ice-tide grid';
 nc.turn.type='NC_FLOAT';
 nc.turn.dimension={'y','x'};

end

if debug; disp(['...Leaving ',mfilename]); end
