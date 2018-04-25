function [nc]=ridgepack_ciceoutput(nc)

% ridgepack_ciceoutput - Condition nc structure from a CICE parameters 
%
% function [nc]=ridgepack_ciceoutput(nc)
%
% INPUT:
%
% nc - structure input with 'CICE' in the source text
%
%
% OUTPUT:
%
% nc - structure output with a few changed attributes 
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if debug; disp('Reconfiguring nc structure from CICE input'); end

% set default maximum value variables may have
maxval=9.0E29;

if isfield(nc.attributes,'source') && ~isempty(strfind(char(nc.attributes.source),'CICE'))

 if isfield(nc,'time')
  nc.time.units='days since 1800-01-01 00:00:00';
 end

 if isfield(nc.attributes,'conventions'); nc.attributes=rmfield(nc.attributes,'conventions'); end;
 if isfield(nc.attributes,'contents'); nc.attributes=rmfield(nc.attributes,'contents'); end;
 if isfield(nc.attributes,'history'); nc.attributes=rmfield(nc.attributes,'history'); end;
 if isfield(nc.attributes,'comment1'); nc.attributes=rmfield(nc.attributes,'comment'); end;
 if isfield(nc.attributes,'comment2'); nc.attributes=rmfield(nc.attributes,'comment2'); end;
 if isfield(nc.attributes,'comment3'); nc.attributes=rmfield(nc.attributes,'comment3'); end;

 if isfield(nc.attributes,'NCO')
  nc.attributes=rmfield(nc.attributes,'NCO');
 end


 % rename the mask to 'mask'
 if isfield(nc,'tmask') | isfield(nc,'mask')

  if isfield(nc,'tmask')
   nc.mask=nc.tmask; 
   nc=rmfield(nc,'tmask');
  end

  nc.mask.data(nc.mask.data>maxval)=0;

  % order the structure dimensions as time,y,x
  if isfield(nc,'nc') & isfield(nc,'time') & isfield(nc,'x') & isfield(nc,'y')
   nc=ridgepack_shuffle(nc,{'time','nc','y','x'});
  elseif isfield(nc,'time') & isfield(nc,'x') & isfield(nc,'y')
   nc=ridgepack_shuffle(nc,{'time','y','x'});
  elseif isfield(nc,'nc') & isfield(nc,'time') & isfield(nc,'nj') & isfield(nc,'ni')
   nc=ridgepack_shuffle(nc,{'time','nc','nj','ni'});
  elseif isfield(nc,'time') & isfield(nc,'nj') & isfield(nc,'ni')
   nc=ridgepack_shuffle(nc,{'time','nj','ni'});
  end

  % apply the mask to each variable with the time,y & x dimensions 
  % or to supporting variables with x & y as dimensions
  mask=nc.mask.data; mask(mask==0)=NaN;

  % Mask out parts of arrays to remove extreme values
  [variablenames,numbervariables]=ridgepack_name(nc);
  for i=1:numbervariables
   name=char(variablenames(i));
   if (any(strcmp(nc.(name).dimension,'x')) & any(strcmp(nc.(name).dimension,'y'))) | ...
      (any(strcmp(nc.(name).dimension,'ni')) & any(strcmp(nc.(name).dimension,'nj')))
    if any(strcmp(nc.(name).dimension,'time')) & ndims(mask)~=ndims(nc.(name).data)
     if any(strcmp(nc.(name).dimension,'nc')) & ndims(nc.(name).data)==4
      for ncat=1:size(nc.(name).data,2)
       for ii=1:size(nc.(name).data,1)
        nc.(name).data(ii,ncat,:,:)=squeeze(nc.(name).data(ii,ncat,:,:)).*mask;
       end
      end
     else
      for ii=1:size(nc.(name).data,1)
       nc.(name).data(ii,:,:)=squeeze(nc.(name).data(ii,:,:)).*mask;
      end
     end
     nc.(name).data(nc.(name).data>maxval)=NaN;
     nc.(name).fillvalue=-99999999;
    elseif ndims(mask)==ndims(nc.(name).data) 
     nc.(name).data(nc.(name).data>maxval)=NaN;
     nc.(name).fillvalue=-99999999;
    end

   end

  end

 end

else

 error('This is not a CICE netcdf structure')

end

% Creating a turning angle if it does not exist 
[variablenames,numbervariables]=ridgepack_name(nc);
for i=1:numbervariables
 name=char(variablenames(i));
 if strcmpi(name,'angle')
    nc.turn=nc.(name);
    if strcmpi(nc.turn.units,'radians')
     nc.turn.data=rad2deg(nc.(name).data);
     nc.turn.units='degrees';
    end
 end
end


% temporary fix to grid-edge problem in RACM
%if isfield(nc,'latitude') | isfield(nc,'ULAT') | isfield(nc,'TLAT')
% disp('**TEMPORARY FIX TO GRID EDGE PROBLEM IN RACM/CICE**')
% try
%  nc.latitude.data(1,:)=nc.latitude.data(2,:);
%  nc.longitude.data(1,:)=nc.longitude.data(2,:);
% end
% try
%  nc.TLAT.data(1,:)=nc.TLAT.data(2,:);
%  nc.TLON.data(1,:)=nc.TLON.data(2,:);
% end
% try
%  nc.ULAT.data(1,:)=nc.ULAT.data(2,:);
%  nc.ULON.data(1,:)=nc.ULON.data(2,:);
% end
%end

if debug; disp(['...Leaving ',mfilename]); end


