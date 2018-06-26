% ridgepack_extent_area_to_nc - calculate CDR sea ice extent and area timeseries
% 
% This script calculates sea ice extent and area from NOAA CDR data previously
% peeled out of NSIDC netcdf files and placed in a seperate netcdf file
% using the ridgepack_cdr_ssmi_to_nc script.
%
% Andrew Roberts, Naval Postgraduate School, June 2018 (afrobert@nps.edu) 

clear

datahome='~/data/SATELLITE/processed';

% specify hemisphere being processed
hemisphere='north'
%hemisphere='south'

% specify dataset from NSIDC (in netcdf format)
if strcmp(hemisphere,'north')
 dataset='G02202_v3_merged_conc_north_1979_2017';
else
 dataset='G02202_v3_merged_conc_south_1979_2017';
end

% add area for which the polar hole exists
areacutoff=87; % degrees north

cd(datahome)

ncgeo=ridgepack_clone(dataset,{'latitude','longitude','time'});

% read in masked areas
if strcmp(hemisphere,'north')
 ncmask=ridgepack_clone('ssmi_25km_north_area_mask_sectors');
else
 ncmask=ridgepack_clone('ssmi_25km_south_area_mask_sectors');
end
area=flipud(ncmask.area.data.*ncmask.mask.data);

% add north or south attributes
if strcmp(hemisphere,'north')
 nce.attributes.title=['Sea ice extent including pole hole for ',dataset];
else
 nce.attributes.title=['Sea ice extent for ',dataset];
end

nce.time=ncgeo.time;

% add in north/south extent long name
if strcmp(hemisphere,'north')
 nce.extent.long_name=['Sea ice extent including pole hole for ',dataset];
else
 nce.extent.long_name=['Sea ice extent for ',dataset];
end
nce.extent.units='km^2';
nce.extent.dimension=ncgeo.time.dimension;
nce.extent.data=zeros(size(ncgeo.time.data));

% add in north/south area long name
if strcmp(hemisphere,'north')
 nce.area.long_name=['Sea ice area south of ',num2str(areacutoff),'N, 1988 onwards'];
else
 nce.area.long_name=['Sea ice area'];
end
nce.area.units='km^2';
nce.area.dimension=ncgeo.time.dimension;
nce.area.data=zeros(size(ncgeo.time.data));

% add in area for where SMMR has been calculated
if strcmp(hemisphere,'north')
 nce.area84.long_name=['Sea ice area south of 84N up to 1988'];
 nce.area84.units='km^2';
 nce.area84.dimension=ncgeo.time.dimension;
 nce.area84.data=zeros(size(ncgeo.time.data));
end

numm=0;

% take pole hole into account to calculate extent
for k=1:length(ncgeo.time.data)

 nc=ridgepack_clone(dataset,{'conc','polehole'},k);

 cellextent=area(squeeze(nc.conc.data(1,:,:))>15 | squeeze(nc.polehole.data(1,:,:))==1);
 nce.extent.data(k)=sum(cellextent(:));

 if ncgeo.time.data(k)>=datenum(1988,1,1) | strcmp(hemisphere,'south')
  conc=squeeze(nc.conc.data(1,:,:))/100;
  cellarea=area.*conc;
  cellarea=cellarea(nc.latitude.data(:,:)<areacutoff & conc(:,:)>0.15);
  nce.area.data(k)=sum(cellarea(~isnan(cellarea)));
 else
  nce.area.data(k)=NaN;
 end

 if strcmp(hemisphere,'north')
  if ncgeo.time.data(k)<datenum(1988,1,1)
   conc=squeeze(nc.conc.data(1,:,:))/100;
   cellarea=area.*conc;
   cellarea=cellarea(nc.latitude.data(:,:)<84 & conc(:,:)>0.15);
   nce.area84.data(k)=sum(cellarea(~isnan(cellarea)));
  else
   nce.area84.data(k)=NaN;
  end
 end

end

ridgepack_write(nce,[dataset,'_extent_area'])


