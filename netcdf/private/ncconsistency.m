function [varadd]=ncconsistency(nc1,nc2)

% NCCONSISTENCY - Check the consistency of two nc structures
%
% function [varadd]=ncconsistency(nc1,nc2)
%
% This is a sub-function of ncmerge to check that common variables between two
% nc structures are consistent in their structure and content, and to determine
% which variables need to be added to nc2 from nc1 that are not already in nc2.
%
% Input:
% nc1, nc2 - two nc structures to be merged (nc1 merged into nc2)
%
% Ouput:
% varadd - cell array of the variable names to be added from nc1 to nc2
%
% Written by Andrew Roberts 2012
% Department of Oceanography, Naval Postgraduate School
%
% $Id$

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

[variablenames1,numbervariables1]=ncname(nc1);
[variablenames2,numbervariables2]=ncname(nc2);

k=0;
for m=1:numbervariables1
 name1=char(variablenames1(m));
 if isfield(nc2,name1)
  [v1,n1]=ncname(nc1.(name1));
  for n=1:n1
   sub1=char(v1(n));
   if isfield(nc2.(name1),sub1)
     if strcmpi(sub1,'data') 
      s1=size(nc1.(name1).data);
      s2=size(nc2.(name1).data);
      if length(s1) ~= length(s2)
       error([name1,' data dimensions are inconsistent']);
      elseif s1 ~= s2
       error([name1,' data matrix sizes are different']);
      elseif nc1.(name1).data ~= nc2.(name1).data
       error([name1,' data is different']);
      end
     elseif strcmpi(sub1,'long_name')
      if nc1.(name1).long_name ~= nc2.(name1).long_name
       error([name1,' long_names are different']);
      end
     elseif strcmpi(sub1,'units')
      if nc1.(name1).units ~= nc2.(name1).units
       error([name1,' units are different']);
      end
     elseif strcmpi(sub1,'fillvalue')
      if nc1.(name1).fillvalue ~= nc2.(name1).fillvalue
       error([name1,' fillvalue is different']);
      end
     elseif strcmpi(sub1,'type')
      if nc1.(name1).type ~= nc2.(name1).type
       error([name1,' type is different']);
      end
     elseif strcmpi(sub1,'coordinates')
      if nc1.(name1).coordinates ~= nc2.(name1).coordinates
       error([name1,' coordinates are different']);
      end
     elseif strcmpi(sub1,'grid_mapping')
      if nc1.(name1).grid_mapping  ~= nc2.(name1).grid_mapping
       error([name1,' grid_mapping is different']);
      end
     elseif strcmpi(sub1,'dimension')
      if length(nc1.(name1).dimension) ~= length(nc2.(name1).dimension)
       error([name1,' dimension cell arrays are different sizes']);
      end
      for kk=1:length(nc1.(name1).dimension)
       if ~strcmpi(char(nc1.(name1).dimension{kk}),char(nc2.(name1).dimension{kk}))
        error([name1,' dimensions are inconsistent']);
       end
      end
     else
       disp(['Non-standard field ',sub1,' not changed in merged structure']);
     end
   else
     error(['Fields are not consistent for ',name1])
   end
  end
  %disp([name1,' is consistent.']); 
 else
  k=k+1;
  varadd{k}=name1;
 end
end

if k==0; varadd=[]; end

if debug; disp(['...Leaving ',mfilename]); end

