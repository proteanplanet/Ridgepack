function [nc]=ncpermute(nc,var,permutation)

% NCPERMUTE - Rearrange dimensions of a variable in a netcdf structure
%
% function [nc]=ncpermute(nc,var,permutation)
%
% This function creates a new permutation of an existing variable in a netcdf
% structure and write it to a new netcd structure, along with its associated
% dimension variables. 
%
% Inputs:
% nc          - netcdf structure
% var         - variable in netcdf structure to be rearranged (character)
% permutation - cell array of the dimensions in the required permutation of var
%
% Output:
% nc          - netcdf structure with rearranged var and its dimension variables.
%
% Example:
% Say you have a netcdf structure like this:
%
% netcdf structure nc {
%  attributes {global}
%     Conventions: 'CF-1.0'
%     history: '2007-11-14 22:43:55 GMT by mars2netcdf-0.92'
%  
%  4 variables, dimension=[ ], data=( )
%       
%  longitude [0 to 357.5, range 357.5]
%     long_name: 'longitude'
%         units: 'degrees_east'
%     dimension: {'longitude'}
%          data: [144x1 double]
% 	   type: 'NC_FLOAT'
%  latitude [-90 to 90, range 180]
%     long_name: 'latitude'
%         units: 'degrees_north'
%     dimension: {'latitude'}
%          data: [73x1 double]
%          type: 'NC_FLOAT'
%  time [01-Sep-1957 12:00:00 to 01-Aug-2002 12:00:00, mean step 31 days]
%     long_name: 'time'
%         units: 'hours since 1900-01-01 00:00:0.0'
%     dimension: {'time'}
%          data: [540x1 double]
%          type: 'NC_FLOAT'
%  p2t (198.5297 to 321.1956, range 122.6659, median 283.0118)
%     long_name: '2 metre temperature'
%         units: 'K'
%     dimension: {'time'  'latitude'  'longitude'}
%          data: [540x73x144 double]
%     fillvalue: -32767
%          type: 'NC_FLOAT'
% }
%
%
% If you wish to rearrange this so that the order of dimensions in p2t
% is {'latitude'  'time'  'longitude'}  instead, and thus nc.p2t.data
% would have the dimensions 73x540x144 instead of 540x73x144, then 
% one would create a new netcdf structure containing this new permutation
% with the following command:
%
% nc=ncpermute(nc,'p2t',{'latitude','time','longitude'});
%
% One reason for doing this could be to reduce the number of dimensions in p2t using
% ncreducto, for example. 
%
% Written by Andrew Roberts 2012
% Department of Oceanography, Naval Postgraduate School
%
% $Id$

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% Check data is a structure
if not(isstruct(nc));
	error([inputname(1),' is not a structure']);
end


% Check that the correct number of dimensions are specified
if length(permutation)~=length(nc.(var).dimension)
	error(['Permutation has an incorrect number of dimension for ',var])
end

% check that all dimensions listed in permutation exist in the netcdf structure
% and are dimensions of the variable var.
for i=1:length(permutation);
	if ~isfield(nc,char(permutation(i)))
		disp([char(permutation(i)),' is not in the netcdf structure provided']);
		error('Specified incorrect nc dimension');
	elseif ndims(nc.(char(permutation(i))).data)>2 | min(size(nc.(char(permutation(i))).data)>1)
		disp([char(permutation(i)),' is not a dimension in the provided netcdf structure']);
		error('Specified incorrect nc dimension');
	elseif ~any(strcmp(char(permutation(i)),nc.(var).dimension))
		disp(['The requested dimension ',char(permutation(i)),' does not exist in ',var]);
		error('Specified incorrect nc dimension');
	end
end

% Calculate the order matrix for rearranging dimensions of var 
neworder=zeros(1, length(nc.(var).dimension));
for j=1:length(permutation);
	for i=1:length(nc.(var).dimension)
		if any(strcmp(char(nc.(var).dimension(i)),char(permutation(j))))
			neworder(j)=i;
		end
	end
end

% write out the new rearranged variable
nc.(var).data=permute(nc.(var).data,neworder);
nc.(var).dimension=permutation;

% rearrange weight arrays if they exist in the structure
if isfield(nc.(var),'weight')
 nc.(var).weight=permute(nc.(var).weight,neworder);
 %disp(['weightings for ',var,' have also been permuted']);
end

if debug; disp(['...Leaving ',mfilename]); end

