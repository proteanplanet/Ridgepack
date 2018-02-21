function ncg=ncgeostrophic(nc,name)

% NCGEOSTROPHIC - Calculate geostropic wind from a pressure field
%
% function ncg=ncgeostrophic(nc,name)
%
% This function calculates geostrophic wind from a pressure
% field in a nc structure.
%
% Input:
% nc   - nc structure 
% name - name of pressure field in netcdf structure
%
% Output:
% ncg  - output netcdf structure containing the same
%        latitudes, longitude, x and y coordinates of the
%        original structure, and ugwind and vgwind 
%        geostrophic components for each timestep in time.
%
% Written by Andrew Roberts 2012
% Department of Oceanography, Naval Postgraduate School
%
% $Id$

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% get constants for calculations
h=ncastroconstants;

% Put into generalized coordinates if need be
[nc,calc]=ncsph2gen(nc);

% position dimensions for gradient
olddimorder=nc.(name).dimension;
ncg=ncshuffle(nc,{'y','x'});

% initialise geostrophic wind variables
ncg.ugwind=ncg.(name);
ncg.vgwind=ncg.(name);

nc

ncstruct(nc)

% perform calculations for the icepack grid
if isfield(nc,'icepack_mesh') | isfield(nc,'stereographic_mesh')

 dx=nc.x.data(2)-nc.x.data(1);
 dy=nc.y.data(2)-nc.y.data(1);

 % make sure units are in SI mks
 if strcmp('km',nc.x.units)
  dx=dx*1000.0;
  dy=dy*1000.0;
 elseif ~strcmp('m',nc.x.units)
  error('units of dx undefined for icepack mesh');
 end

 disp(['Using dx=',num2str(dx/1000.),' dy=',num2str(dy/1000.),' km']);

 if ~isempty(intersect({'time'},nc.(name).dimension))
  for i=1:length(ncg.time.data);
   disp(['Creating geostrophic wind for ',num2str(i),' ',datestr(ncg.time.data(i),0),' for ',name]);
   p=squeeze(ncg.(name).data(i,:,:));
   lat=ncg.latitude.data;
   [Px,Py]=gradient(p,dx,dy);
   rf=h.rhoa.const*2.*h.omega.const*sin(lat*pi/180);
   ncg.ugwind.data(i,:,:)=-Py./rf;
   ncg.vgwind.data(i,:,:)=+Px./rf;
  end
 else
  lat=ncg.latitude.data;
  [Px,Py]=gradient(ncg.(name).data,dx,dy);
  rf=h.rhoa.const*2.*h.omega.const*sin(lat*pi/180);
  ncg.ugwind.data(:,:)=-Py./rf;
  ncg.vgwind.data(:,:)=+Px./rf;
 end

 ncg.ugwind.units='m/s';
 ncg.ugwind.long_name='Geostrophic wind U-component from PMSL';
 ncg.ugwind.type='NC_FLOAT';

 ncg.vgwind.units='m/s';
 ncg.vgwind.long_name='Geostrophic wind V-component from PMSL';
 ncg.vgwind.type='NC_FLOAT';
 
 ncg=rmfield(ncg,name);

 try
  ncg=ncwrite('ugwind.nc',ncg);
 catch
  disp('You chose not to overwrite ugwind.nc');
  ncg=ncstruct(ncg);
 end

 ncg=ncshuffle(ncg,olddimorder);

else

 error('ncgeostrophic not set up for the supplied grid')
 
end

if debug; disp(['...Leaving ',mfilename]); end

