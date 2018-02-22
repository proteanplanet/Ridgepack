function [nc]=ridgepack_shuffle(nc,lastdims)

% ridgepack_shuffle - Permute dimensions in an nc structure to a requested order
%
% function [nc]=ridgepack_shuffle(nc,lastdims)
%
% This function permutes dimensions in a netcdf structure so that the
% dimensions listed in the lastdims cell array are last for each valid
% array in the nc structure.  This is particularly useful when operating
% on certain dimensions within an array, because it can easily 
% reorganize the array.
%
% Input:
% nc       - input netcdf structure (see ridgepack_struct for more information)
%
% lastdims - cell array providing the dimensions wished to be made last
%            in arrays in the structure (see working example). The order
%            of dimensions specified in lastdims is the same order that
%            the nc structure will have.
%
% Output:
% nc       - netcdf structure with permuted dimensions.
%
% Working example:
%
% You have an existing netcdf structure like this:
%
% netcdf structure nc {
% attributes {global}
%    Conventions: 'CF-1.1'
%        history: '2007-11-14 22:43:55 GMT by mars2netcdf-0.92'
%	      
% 4 variables, dimension=[ ], data=( )
%	       
% latitude [-90 to 90, range 180, step -2.5]
%      long_name: 'latitude'
%          units: 'degrees_north'
%      dimension: {'latitude'}
%           data: [73x1 double]
%           type: 'NC_FLOAT'
%			  
% longitude [0 to 357.5, range 357.5, step 2.5]
%      long_name: 'longitude'
%          units: 'degrees_east'
%      dimension: {'longitude'}
%           data: [144x1 double]
%           type: 'NC_FLOAT'
%										     
% time [01-Sep-1957 12:00:00 to 01-Aug-2002 12:00:00, timestep 31 days]
%      long_name: 'time'
%          units: 'hours since 1900-01-01 00:00:0.0'
%      dimension: {'time'}
%           data: [540x1 double]
%           type: 'NC_FLOAT'
%
% p2t (198.5297 to 321.1956, range 122.6659, median 283.0118)
%      long_name: '2 metre temperature'
%          units: 'K'
%      dimension: {'time'  'latitude'  'longitude'}
%           data: [540x73x144 double]
%      fillvalue: -32767
%           type: 'NC_FLOAT'
%										       
% }
%
% And you want to rearrange the order of the dimensions in p2t to be 
% {'latitude'  'time'  'longitude'} which would create an array of shape 73x540x144.  
% This can be done by entering the command:
%
% nc=ridgepack_shuffle(nc,{'latitude','time','longitude'});
%
% When this is done, the structure remains much the same because time, longitude and
% latitude only have one dimension.  However p2t is changed to become:
%
% netcdf structure nc {
% attributes {global}
%    Conventions: 'CF-1.1'
%        history: '2007-11-14 22:43:55 GMT by mars2netcdf-0.92'
%	      
% 4 variables, dimension=[ ], data=( )
%	       
% ...
%
% p2t (198.5297 to 321.1956, range 122.6659, median 283.0118)
%      long_name: '2 metre temperature'
%          units: 'K'
%      dimension: {'latitude'  'time'  'longitude'}
%           data: [73x540x144 double]
%      fillvalue: -32767
%           type: 'NC_FLOAT'
%										       
% }
%
% Andrew Roberts, Naval Postgraduate School, March 2018  (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if debug; disp('Shuffling dimensions'); end

% pass through each variables of the netcdf structure to purmute the arrays
[nc,variablenames,numbervariables]=ridgepack_sort(nc);
for m = 1:numbervariables

 var=char(variablenames(m));

 % Check that name dimensions are dimensions of var and if they are,
 % build the name cell structure to create appropriate permutations.
 i=0;
 for j=1:length(lastdims)
  if any(strcmp(char(lastdims{j}),nc.(var).dimension))
   i=i+1;
   name{i}=char(lastdims{j});
  end
 end

 
 %if i>0
 %name
 %nc.(var).dimension(end-length(name)+1:end)
 %end

 if i==0 | length(nc.(var).dimension)==1

  if debug; disp([var,' does not need to be shuffled']); end

 elseif strcmp(nc.(var).dimension(end-length(name)+1:end),name)

  if debug; disp([var,' dimensions already in requested order']); end

 else

  % get old size in text
  oldsize=[' [',num2str(size(nc.(var).data)),'] '];

  % Set sizes of finaldimension and permutation cell arrays
  permutation=cell(1,length(nc.(var).dimension));
  finaldimension=cell(1,length(nc.(var).dimension)-length(name));

  % Assign dimension names for those not to be removed
  j=0;
  for i=1:length(nc.(var).dimension)
   if ~any(strcmp(char(nc.(var).dimension(i)),name))
     j=j+1;
     permutation{j}=char(nc.(var).dimension(i));
     finaldimension{j}=char(nc.(var).dimension(i));
   end
  end

  % Assign dimension names to be removed
  for i=1:length(name);
   j=j+1;
   permutation{j}=char(name(i));
  end

  if j~=length(nc.(var).dimension)
   error(['permutation does not have same length as ',var,' dimensions: j=',...
           num2str(j),' dims=',num2str(length(nc.(var).dimension))]);
  end

  % Permute var to the order specified in permutation
  nc=ridgepack_permute(nc,var,permutation);

  if debug; 
        disp([var,oldsize,'-> [',num2str(size(nc.(var).data)),...
        ']=[',ridgepack_cellcat(nc.(var).dimension),' ]']);
  end

 end

 clear name

end

if debug; disp('Finished shuffling dimensions'); end

if debug; disp(['...Leaving ',mfilename]); end


