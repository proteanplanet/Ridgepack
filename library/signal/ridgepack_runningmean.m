function [nc]=ridgepack_runningmean(nc,var,windowsize)

% ridgepack_runningmean - Calculates a running mean through time for a given variable
%
% function [nc]=ridgepack_runningmean(nc,var,windowsize)
%
% This function calculates the running mean for a given variable var in time
% assuming that the variable has is represented in the nc structure by equally
% spaced samples in time.
%
% INPUT:
%
% nc         - nc structure (see ncstruct for more details)
% var        - string providing the variable name in nc for which filtering
%              should be performed.
% windowsize - the number of samples over which each mean should be taken
% 
%
% OUTPUT:
%
% nc - nc structure with a variable var_running added to the structure which 
%      provides the running mean of the chosen variable var.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
% Naval Postgraduate School, Department of Oceanography

% Check data is a structure
if not(isstruct(nc));
    error([inputname(1),' is not a structure']);
end

% Check that var exists in ncstruct
[nc,variablenames,numbervariables]=ridgepack_sort(nc);
if ~any(strcmp(var,variablenames))
   error(['No such variable ',var,' exists in ',inputname(1)]);
end

% Check that time is a dimension of var
if ~any(strcmp('time',nc.(var).dimension))
   error(['time is not a dimension of ',var])
end

% Check that record length of time is greater than one
if length(nc.time.data)<=1 
   error('time only has 1 record or less in the structure - no filtering possible')
end

% name new variable
newvar=[var,'_running'];

% record current dimension arrangement
olddims=nc.(var).dimension;

% Arrange variable so that time is the last dimension in the structure
nc=ridgepack_shuffle(nc,{'time'});

% Reshape nc.(var).data array to be 2D
if length(nc.(var).dimension)==1
 nc.(var).data=nc.(var).data';
end
oldsize=size(nc.(var).data);
newsize=[prod(oldsize(1:end-1)), oldsize(end)];
nc.(var).data=reshape(nc.(var).data,newsize);

% populate new variable for efficient matlab memory use
nc.(newvar).data=nc.(var).data;

% Filter the data with a running mean
a=1;
b=ones(1,windowsize)/windowsize;
for i=1:newsize(1);
 %nc.(newvar).data(i,:)=filter(b,a,nc.(var).data(i,:));
 nc.(newvar).data(i,:)=filtfilt(b,a,nc.(var).data(i,:));
end

% Fill the first part of the running mean with fillvalue
nc.(newvar).fillvalue=min(-99999.9,10*min(nc.(newvar).data(:)));
for i=1:newsize(1);
 nc.(newvar).data(i,1:windowsize-1)=NaN;
 nc.(newvar).data(i,end-windowsize+1:end)=NaN;
end

% Reshape arrays back to original size
nc.(var).data=reshape(nc.(var).data,oldsize);
nc.(newvar).data=reshape(nc.(newvar).data,oldsize);
if length(nc.(var).dimension)==1
 nc.(var).data=nc.(var).data';
 nc.(newvar).data=nc.(newvar).data';
end

% Add information to new variable
nc.(newvar).long_name=[num2str(windowsize),' sample running mean ',nc.(var).long_name];
nc.(newvar).dimension=nc.(var).dimension;
if isfield(nc.(var),'units')
 nc.(newvar).units=nc.(var).units;
end
nc.(newvar).type=nc.(var).type;

% Put dimension order back to the original
nc=ridgepack_struct(ridgepack_shuffle(nc,olddims));

end % function

