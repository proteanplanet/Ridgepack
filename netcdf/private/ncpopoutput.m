function [nc]=ncpopoutput(nc)

% NCPOPOUTPUT - Condition an nc structure assuming it derives from a POP hist file
%
% function [nc]=ncpopoutput(nc)
%
% Input:
% nc - structure input with 'POP HIST conventions' as the convention
%
% Output:
% nc - structure output with time added from the attributes
%
% Written by Andrew Roberts 2012
% Department of Oceanography, Naval Postgraduate School
%
% $Id$

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if debug; disp('Reconfiguring nc structure from POP input'); end

if isfield(nc.attributes,'conventions') && strcmp(nc.attributes.conventions,'POP HIST conventions')

 nc.time.data=datenum(nc.attributes.iyear,nc.attributes.imonth,nc.attributes.iday); 
 nc.time.dimension={'time'};
 nc.time.units='days since 1800-01-01 00:00:00';

 nc.attributes=rmfield(nc.attributes,'iyear');
 nc.attributes=rmfield(nc.attributes,'imonth');
 nc.attributes=rmfield(nc.attributes,'iday');
 nc.attributes=rmfield(nc.attributes,'conventions');

 % - add time as a dimension if required
 % - add the turning angle for plotting vectors on a map
 [variablenames,numbervariables]=ncname(nc);
 for i=1:numbervariables

  name=char(variablenames(i));
  ii=length(nc.(name).dimension);
  if ii>1
   nc.(name).dimension{ii+1}='time';
  end

  if strcmpi(name,'angle')
   nc.turn=nc.(name)
   if strcmpi(nc.turn.units,'radians')
    nc.turn.data=rad2deg(nc.(name).data);
    nc.turn.units='degrees';
   end
  end

 end

 try
  disp('Looking to append geopositions to popcice');
  ncg=ncclone([getenv('HOME'),'/local/src/popcice/etc/popcice_scalar_grid.nc']);
  ncg.latitude.dimension={'i','j'};
  ncg.longitude.dimension=ncg.latitude.dimension;
  ncg.mask.dimension=ncg.latitude.dimension;
  ncg.rotang.dimension=ncg.latitude.dimension;
  ncg.tarea.dimension=ncg.latitude.dimension;
  ncg.latitude.data=ncg.latitude.data';
  ncg.longitude.data=ncg.longitude.data';
  ncg.mask.data=ncg.mask.data';
  ncg.rotang.data=ncg.rotang.data';
  ncg.tarea.data=ncg.tarea.data';
  nc.latitude=ncg.latitude;
  nc.longitude=ncg.longitude;
  nc.mask=ncg.mask;
  nc.rotang=ncg.rotang;
  nc.tarea=ncg.tarea;
 catch
  disp('Unable to append popcice latitude and longitude cell centers');
 end

elseif debug

 disp('This appears not to be POP model output')

end

if debug; disp(['...Leaving ',mfilename]); end

