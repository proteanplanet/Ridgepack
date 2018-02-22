function [nc]=ridgepack_wrfoutput(nc)

% ridgepack_wrfoutput - Prepares a netcdf structure from a wrf output netcdf file
%
% function [nc]=ridgepack_wrfoutput(nc)
%
% This function reduces lat and long to two dimensions (removes the time
% dimension), and adds a turning angle for vectors for a polar stereographic
% grid.
%
% Input: 
% nc - netcdf structure prepared by ridgepack_clone
%
% Output:
% nc - netcdf strucure with certain changes applicable to the WRF output
%
% Written Andrew Roberts 2012
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if debug; disp('Reconfiguring nc structure from WRF input'); end

% deal with WRF's character arrays for time first
if isfield(nc.attributes,'SIMULATION_START_DATE') & isfield(nc,'time')

 % set WRF time units
 startdate=datenum(char(nc.attributes.SIMULATION_START_DATE));

 time=sort(nc.time.data(:));
 if startdate<time(1)
  nc.time.units=['hours since ',datestr(startdate,'yyyy-mm-dd HH:MM:SS.FFF')];
 end

end

% calculate a turning angle for WRF Polar Stereographic
if isfield(nc.attributes,'MAP_PROJ') && nc.attributes.MAP_PROJ==2 && isfield(nc,'XLONG')==2;
 disp('Adding turning angle for WRF polar stereographic grid');
 nc.turn.data=(180+nc.XLONG.data);
 nc.turn.long_name='Turning angle for vectors to lat-lon for WRF polar stereographic grid';
 nc.turn.units='degrees';
 nc.turn.dimension=nc.XLONG.dimension;
end

% strip out WRF's millions of attributes
names=fieldnames(nc.attributes);
for i=1:length(names)
 ustri=char(names{i});
 uprop=find(~isstrprop(ustri,'upper'));
 dprop=~isstrprop(ustri,'digit');
 pf=true;
 if ~isempty(uprop)
  for j=uprop
   if ~strcmp(ustri(j),'_') & dprop(j)
    pf=false;
    continue 
   end
  end
 end
 if pf
  nc.attributes=rmfield(nc.attributes,ustri);
 end
end
 
if debug; disp(['...Leaving ',mfilename]); end

