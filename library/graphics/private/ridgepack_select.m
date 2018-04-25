function [ncr]=ridgepack_select(nc,X,Y,Z,dimnames,bbounds)

% ridgepack_select - For selecting dimensions and bbounds for plotting nc structured data
%
% function [ncr]=ridgepack_select(nc,X,Y,Z,dimnames,bbounds)
%
% This function is a front-end to ridgepack_reduce for selecting axes and data
% to be plotted using ncfigure.  
%
% INPUT:
%
% nc        - nc structure (see ridgepack_struct for more details)
% X         - string providing the X axis variable from nc
% Y         - string providing the Y axis or Y variable on a graph from nc
% Z         - string providing the Z Data field or Y variable on a graph in nc
%             (for the case of graphs, Z and Y should be the same string)
% dimnames  - cell array of names of dimensions over which Z (and Y for graphs) 
%             is being sliced or average. This input is explained in ridgepack_reduce.
% bbounds   - bbounds for the dimensions listed in the dimnames cell array.
%             This is a cell array and is explained in detail in ridgepack_reduce.
%
%
% OUTPUT:
%
% ncr - the reduced nc structure output from ridgepack_reduce which is checked for errors.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)


global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% Check data is a structure
if not(isstruct(nc));
 error([inputname(1),' is not a structure']);
end

% Check to see if Z has been provided, if not, assign Y
if nargin<6
	bbounds={};
elseif ~iscell(bbounds)
	error('bbounds should be a cell array')
end
if nargin<5
	dimnames={};
elseif ~iscell(dimnames)
	dimnames
	error('dimnames should be a cell array')
end
if nargin<4
	Z=Y;
	if size(Z,1)>1 && size(Z,2)>1 
	 error('Problem with Y coordinate')
        end
end

if nargin<3
	error('Input must at least include X and Y input')
end

if ~ischar(X) | ~ischar(Y) | ~ischar(Z)
	error('X, Y and Z inputs must be strings')
end

if ~isfield(nc,X);
 error(['''',X,''' is not in the nc structure'])
elseif ~isfield(nc,Y);
 error(['''',Y,''' is not in the nc structure'])
elseif ~isfield(nc,Z);
 error(['''',Z,''' is not in the nc structure'])
elseif ~isfield(nc.(X),'dimension');
 error(['''',X,''' dimension is not in the nc structure'])
elseif ~isfield(nc.(Y),'dimension');
 error(['''',Y,''' dimension is not in the nc structure'])
elseif ~isfield(nc.(Z),'dimension');
 error(['''',Z,''' dimensions are not in the nc structure'])
elseif strcmp(X,Y)
 error('X and Y are the same dimensions')
end

for i=1:length(nc.(X).dimension)
 if ~(strcmp(char(nc.(X).dimension{i}),nc.(Z).dimension))
  error([X,' dimension ',char(nc.(X).dimension{i}),' is not a dimensions of ',Z])
 end
end

for i=1:length(nc.(Y).dimension)
 if ~(strcmp(char(nc.(Y).dimension{i}),nc.(Z).dimension))
  error([Y,' dimension ',char(nc.(Y).dimension{i}),' is not a dimensions of ',Z])
 end
end


if debug
 disp('Finished dimension checking in ridgepack_select')
end

% strip type of calculation off X, Y and Z & check that X, Y and Z exist
[nc,variablenames,numbervariables]=ridgepack_sort(nc);

if ~isempty(setdiff({X},variablenames))
	error(['Unable to find X variable ',X]);
end

if ~isempty(setdiff({Y},variablenames))
	error(['Unable to find Y variable ',Y]);
end

oldZ=Z;
calctypes={'_std','_samp'};
if ~isempty(setdiff({Z},variablenames))
 for i=1:length(calctypes)
	 Z=strtok(oldZ,char(calctypes{i}));
	 unrecog=ridgepack_setdiff({Z},variablenames);
	 if isempty(unrecog)
		 break
	 elseif i==length(calctypes)
		 error(['Unable to find Z variable ',Z]);
	 end
 end
end

% Find dimensions of data
xdims=nc.(X).dimension;
ydims=nc.(Y).dimension;
zdims=nc.(Z).dimension;

% check dimensions of X and Y are in Z
if strcmp(Z,Y) & ~isempty(setdiff(xdims,ydims))
 disp(['x dimensions: ',ridgepack_cellcat(xdims)])
 disp(['y dimensions: ',ridgepack_cellcat(ydims)])
 error(['x-dimensions not found in y-dimensions :',ridgepack_cellcat(setdiff(xdims,ydims))])
elseif ~isempty(setdiff(xdims,zdims))
 disp(['x dimensions: ',ridgepack_cellcat(xdims)])
 disp(['z dimensions: ',ridgepack_cellcat(zdims)])
 error(['x-dimensions not found in z-dimensions :',ridgepack_cellcat(setdiff(xdims,zdims))])
elseif ~isempty(setdiff(ydims,zdims))
 disp(['y dimensions: ',ridgepack_cellcat(ydims)])
 disp(['z dimensions: ',ridgepack_cellcat(zdims)])
 error(['y-dimensions not found in z-dimensions :',ridgepack_cellcat(setdiff(ydims,zdims))])
end


% check that dimensions are consistent with boundaries provided
if isempty(bbounds) & ~isempty(dimnames)
 if length(xdims)>1
	 error(['No reduction dimension provided for ',X])
 elseif ~strcmp(Y,Z) & length(ydims)>1
	 error(['No reduction dimension provided for ',Y])
 end
end


% reduce dataset to a single 2D slice to be plotted
if debug
 dimnames
 bbounds
end
ncr=ridgepack_reduce(nc,dimnames,bbounds);
if debug
 size(ncr.(X).data)
 size(ncr.(Y).data)
 size(ncr.(Z).data)
end

% check that Z variable has been reduced to 2D
if length(size(squeeze(ncr.(Z).data)))>2
 error([Z,' has more than two dimensions: ', num2str(size(ncr.(Z).data))])
end

% Find dimensions of data
xdims=ncr.(X).dimension;
ydims=ncr.(Y).dimension;
zdims=ncr.(Z).dimension;

% rearrange dimensions in Z if required for 1D input dimensions
if length(zdims)>2
 
 % Ccase of 3+ dimensions with the length of the third, fourth, etc dimensions being one
 if ndims(squeeze(ncr.(Z).data))==2
  ncr=ridgepack_reduce(ncr,{char(setdiff(zdims,union(ncr.(X).dimension,ncr.(Y).dimension)))});
  return
 else % otherwise, throw an error
  disp(['Please reduce ',Z,' to only two dimensions from: ',ridgepack_cellcat(zdims)]);
  error(['Please reduce ',Z,' to only two dimensions from: ',ridgepack_cellcat(zdims)]);
 end

elseif length(xdims)==1 && length(ydims)==1 && ...
   all(strcmp({char(ydims) char(xdims)},zdims)) &&  ~strcmp(Z,Y)
 ncr=ridgepack_shuffle(ncr,{char(xdims) char(ydims)});
elseif length(xdims)==1 && ~all(strcmp({char(xdims) char(ydims)},zdims)) && ~strcmp(Z,Y)
 error('Dimensions of Z don''t mach X and Y combined')
end
 

% sort the data
%[ncr,variablenames,numbervariables]=ridgepack_sort(ncr);

if ~isempty(setdiff({X},variablenames))
 error(['Unable to find X variable ',X,' in nc output']);
end

if ~isempty(setdiff({Y},variablenames))
 error(['Unable to find Y variable ',Y,' in nc output']);
end

if ~isempty(setdiff({oldZ},variablenames))
 disp('Check that bbounds and mask allow for standard deviations and samples size arrays to be calculated');
 error(['There is no ',oldZ,' in the nc output. Check that bbounds and mask allow for std and samp.']);
end

if debug; disp(['...Leaving ',mfilename]); end

