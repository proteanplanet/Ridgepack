function [z,nc]=ridgepack_standardunits(nc,Z)

% ridgepack_standardunits - Convert units for some standard meteorological cases
%
% function [z,nc]=ridgepack_standardunits(nc,Z)
%
% INPUT:
%
% nc  - netcdf structure
% Z   - name of variable to be extracted
%
%
% OUTPUT:
%
% z   - data to be interpolated
% nc  - netcdf structure re-gridded
%
% Written Andrew Roberts 2012
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if debug; disp('Checking plotting units'); end

if isfield(nc.(Z),'units') & ischar(nc.(Z).units)

 if ~isempty(strfind(char(nc.(Z).units),'Pascals')) || ~isempty(strfind(char(nc.(Z).units),'Pa'))
	disp('Changing to hPa')
	nc.(Z).units='hPa';
	nc.(Z).data=nc.(Z).data/100;
 end

end

z=nc.(Z).data;

if debug; disp(['...Leaving ',mfilename]); end
