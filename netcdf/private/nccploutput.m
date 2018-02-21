function [nc]=nccploutput(nc)

% NCCPLOUTPUT - Condition nc structure from CPL7 output
%
% function [nc]=nccploutput(nc)
%
% Input:
% nc - structure input with 'CICE' in the source text
%
% Output:
% nc - structure output with a few changed attributes 
%
% Written by Andrew Roberts 2012
% Department of Oceanography, Naval Postgraduate School
%
% $Id$

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if debug; disp('Reconfiguring nc structure from CPL input'); end

if isfield(nc.attributes,'file_version') & ~isempty(strfind(char(nc.attributes.file_version),'cpl7'))

 if isfield(nc.attributes,'conventions'); nc.attributes=rmfield(nc.attributes,'conventions'); end;
 if isfield(nc.attributes,'contents'); nc.attributes=rmfield(nc.attributes,'contents'); end;
 if isfield(nc.attributes,'history'); nc.attributes=rmfield(nc.attributes,'history'); end;
 if isfield(nc.attributes,'comment1'); nc.attributes=rmfield(nc.attributes,'comment'); end;
 if isfield(nc.attributes,'comment2'); nc.attributes=rmfield(nc.attributes,'comment2'); end;
 if isfield(nc.attributes,'comment3'); nc.attributes=rmfield(nc.attributes,'comment3'); end;
 if isfield(nc.attributes,'NCO'); nc.attributes=rmfield(nc.attributes,'NCO'); end

 [variablenames,numbervariables]=ncname(nc);
 for i=1:numbervariables

   name=char(variablenames(i));

   % atmospheric variables
   if  (~isempty(strfind(name,'a2x')) | ~isempty(strfind(name,'x2a'))) 
    if isfield(nc,'doma_lat') & isfield(nc,'doma_lon')
     nc.doma_lat.dimension=ncsetdiff(nc.(name).dimension,{'time'});
     nc.doma_lon.dimension=ncsetdiff(nc.(name).dimension,{'time'});
     nc.(name).coordinates='doma_lat doma_lon';
    end
    if isfield(nc,'doma_area')
     nc.doma_area.dimension=ncsetdiff(nc.(name).dimension,{'time'});
    end
    %nc=rmfield(nc,{'doma_nx','doma_ny'});
   end % atmospheric variables

   % land variables
   if  (~isempty(strfind(name,'l2x')) | ~isempty(strfind(name,'x2l')))
    if isfield(nc,'doml_lat') & isfield(nc,'doml_lon')
     nc.doml_lat.dimension=ncsetdiff(nc.(name).dimension,{'time'});
     nc.doml_lon.dimension=ncsetdiff(nc.(name).dimension,{'time'});
     nc.(name).coordinates='doml_lat doml_lon';
    end
    if isfield(nc,'doml_area')
     nc.doml_area.dimension=ncsetdiff(nc.(name).dimension,{'time'});
    end
    %nc=rmfield(nc,{'doml_nx','doml_ny'});
   end % land variables

 end

elseif debug

 disp('This does not appear to be CPL input')

end

if debug; disp(['...Leaving ',mfilename]); end

