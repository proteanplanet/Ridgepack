function [nc1]=ncmerge(nc1,nc2)

% NCMERGE - Merge two nc structures together, avoiding duplication of supporting data
%
% function [nc1]=ncmerge(nc1,nc2)
%
% Input:
% nc1, nc2 - two nc structures to be merged (nc2 merged into nc1)
% 
% Output:
% nc1 - new nc structure nc1 which is a merge of nc1 and nc2
%
% Written by Andrew Roberts 2012
% Department of Oceanography, Naval Postgraduate School
%
% $Id$

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% run checks on structures
[nc1,out]=ncstruct(nc1);
[nc2,out]=ncstruct(nc2);
disp(' ');

% Run a check on the global attributes
globs={'title','source','references','comment'};
firstgo=false;
for i=1:length(globs)
 globalatt=char(globs{i});
 if isfield(nc1.attributes,globalatt) & isfield(nc2.attributes,globalatt)
  if ~firstgo; disp('Attributes differ...'); firstgo=true; end
  if ~strcmpi(nc1.attributes.(globalatt),nc2.attributes.(globalatt))
   disp(['nc1 ',globalatt,': ',nc1.attributes.(globalatt)]);
   disp(['nc2 ',globalatt,': ',nc2.attributes.(globalatt)]);
   nc1.attributes.(globalatt)=input(['Enter new ',globalatt,': '],'s');
  end
 elseif isfield(nc1.attributes,globalatt)
  if ~firstgo; disp('Attributes differ...'); firstgo=true; end
  nc1=rmfield(nc1.attributes,globalatt);
  disp(['Removing ',globalatt]);
 elseif isfield(nc2.attributes,globalatt)
  if ~firstgo; disp('Attributes differ...'); firstgo=true; end
  disp(['Not merging ',globalatt]);
 end
end

% check for consistency between nc1 and nc2 and vice versa
varadd=ncconsistency(nc2,nc1);

% now add fields missing from nc1 but in nc2 into nc1.
for m=1:length(varadd)

 name=char(varadd{m});

 disp(['Merging ',name]);

 nc1.(name)=nc2.(name);

end

% Run checks on final structure
[nc1,out]=ncstruct(nc1);

% debug information
if debug; disp(['...Leaving ',mfilename]); end
