
clear

hemisphere='north';

startyear=1979;
latestyear=2016;

r=1:361; % index of Polar Pathfinder Grid
s=1:361; % index of Polar Pathfinder Grid
r0=181;  % index of pole
s0=181;  % index of pole
C=25.067525;    % resolution in km (corresponds to 'actual' SSM/I resolution is 25km)
hem=1;   % hemisphere is 1 for north, -1 for south

%rawdir=[getenv('HOME'),'/data/SATELLITE/downloaded/nsidc0116_icemotion_vectors_v3/data/north/means/months'];
rawdir=['/Volumes/RobertsRaid3/data/SATELLITE/downloaded/nsidc0116_icemotion_vectors_v3/data/north/means'];
%processeddir=[getenv('HOME'),'/data/SATELLITE/processed'];
processeddir=['/Volumes/RobertsRaid3/data/SATELLITE/processed'];
stringtime='monthly';

nc.attributes.title=['Polar Pathfinder Gridded Sea Ice Velocity (nsidc0116) ',stringtime,' mean'];

nc.x.data=r;
nc.x.long_name='EASE-grid r index';
nc.x.dimension={'x'};
nc.x.units='';
nc.x.type='NC_FLOAT';

nc.y.data=s;
nc.y.long_name='EASE-grid s index';
nc.y.dimension={'y'};
nc.y.units='';
nc.y.type='NC_FLOAT';

[lat,lon,h,k]=ridgepack_xytoease(r,s,r0,s0,C,hem);

nc.latitude.data=lat;
nc.latitude.long_name='latitude';
nc.latitude.dimension={'x' 'y'};
nc.latitude.units='degrees_north';
nc.latitude.type='NC_FLOAT';

nc.longitude.data=lon;
nc.longitude.long_name='longitude';
nc.longitude.dimension={'x' 'y'};
nc.longitude.units='degrees_north';
nc.longitude.type='NC_FLOAT';

nc.turn.data=-nc.longitude.data;
nc.turn.units='degrees';
nc.turn.long_name='vector turning angle to lat-long (u,v) grid from EASE grid';
nc.turn.type='NC_FLOAT';
nc.turn.dimension={'x','y'};

nc.time.long_name='time';
nc.time.dimension={'time'};
nc.time.data=[];
nc.time.type='NC_DOUBLE';

nc.count.long_name='number of samples in mean';
nc.count.dimension={'x' 'y' 'time'};
nc.count.data=ones(size(lat));
nc.count.type='NC_INT';

nc.u.long_name=['Sea ice velocity u-component ',stringtime,' mean'];
nc.u.dimension={'x' 'y' 'time'};
nc.u.units='cm/s';
nc.u.type='NC_FLOAT';

nc.v.long_name=['Sea ice velocity v-component ',stringtime,' mean'];
nc.v.dimension={'x' 'y' 'time'};
nc.v.units='cm/s';
nc.v.type='NC_FLOAT';

filename=['Pathfinder_icemotion_',stringtime,'_',num2str(startyear),'_',num2str(latestyear),'_v3'];

for k=startyear:1:latestyear
for m=1:12

 disp(['______',num2str(k),'______',])

 nc.time.data=datenum(k,m,1);

 cd([rawdir,'/',num2str(k,'%4.4i')]);

 fid=fopen(['icemotion.grid.month.',num2str(k),'.',num2str(m,'%2.2i'),'.n.v3.bin']);
 [m5,count]=fread(fid,inf,'short',0,'l');
 fclose(fid);

 nc=ridgepack_struct(nc);

 nc.u.data=(reshape(m5(1:3:end-2),361,361)/10)';
 nc.v.data=(reshape(m5(2:3:end-1),361,361)/10)';
 nc.count.data=reshape(m5(3:3:end),361,361)';

 nc=ridgepack_struct(nc);

 clear m5

 cd(processeddir);

 if k==startyear
  ridgepack_write(nc,filename) 
 else
  ridgepack_write(nc,filename,{'time'},{0}) 
 end

 % mask out values with no data count
 nc.u.data(nc.count.data==0)=NaN;
 nc.v.data(nc.count.data==0)=NaN;

 % interpolate to the model grid 
 ncnew=ridgepack_regrid(nc,'u','v',7);
 ncc=ridgepack_regrid(nc,'count','',7);

 % remove all zero count areas
 ncnew.u.data(ncc.count.data<1)=NaN;
 ncnew.v.data(ncc.count.data<1)=NaN;

 % write interpolated data
 if k==startyear
  ridgepack_write(ncnew,[filename,'_RASM_CICE'])
 else
  ridgepack_write(ncnew,[filename,'_RASM_CICE'],{'time'},{0})
 end

end
end



