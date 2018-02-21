function [time]=nccftimeextern(time,units,CFcalendar)

% NCCFTIME - Converts CF convention time to Matlab serial time with udunits
%
% function [time]=nccftime(ncfile)
%
% Input:
% time       - time given in the CF convention units
% units      - CF convention units
% CFcalendar - type of calendar being used (e.g. gregorian, noleap etc)
%
% Output:
% time       - time given in Matlab serial time units
%
% This function invokes the system's version of ncdump for listing 
% time output in CF convention using any available calendar information.
%
% Written by Andrew Roberts 2013
% Department of Oceanography, Naval Postgraduate School
%
% $Id$

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

ncfile=['icepack.tempfile.',num2str(now)];

delim=' || ';

% Generate a dummy netcdf file with the time provided for the calendar
fid=fopen([ncfile,'.cdl'],'w+');
tstring=['netcdf timetemp \{\n dimensions\:\n time = UNLIMITED \; \n'];
fprintf(fid,'%s\n','netcdf timetemp {');
fprintf(fid,'%s\n','dimensions:');
fprintf(fid,'%s\n',' time = UNLIMITED ;');
fprintf(fid,'%s\n','variables:');
fprintf(fid,'%s\n',' double time(time) ;');
fprintf(fid,'%s\n',[' time:units = "',units,'" ;']);
fprintf(fid,'%s\n',[' time:calendar = "',CFcalendar,'" ;']);
fprintf(fid,'%s\n',' time:long_name = "Time" ;');
fprintf(fid,'%s\n','data:');

timechar=num2str(time(1));
for i=2:length(time)
 timechar=[timechar,', ',num2str(time(i))];
end

fprintf(fid,'%s\n',[' time = ',timechar,' ;']);
fprintf(fid,'%s\n','}');
fclose(fid);

% Now read that dummy netcdf file using ncdump -v time -t to interpret calendar
if isunix;

   [h,t]=unix('which ncgen');
   if h ~= 0; 
     error(['Could not find ncgen on the system: Install netcdf libraries, or check the path']);
   end 

   [h,t]=unix(['ncgen -o ',ncfile,'.nc ',ncfile,'.cdl']);
   if h ~= 0; 
     disp(['CHECK YOU ARE USING VERSION 4 OR HIGHER NETCDF LIBRARIES'])
     error(['Could not create a netcdf file from ',ncfile,'.cdl: error ',num2str(h)]);
   end
        
   [h,t]=unix('which ncdump');
   if h ~= 0; 
     error(['Could not find ncdump on the system: Install netcdf libraries, or check the path']);
   end 
        
   [h,t]=unix(['ncdump -v time -t ',ncfile,'.nc | sed -e :a -e''$!N;s/\n//;ta'' | sed -e''s/\(.*\)\(time = \"\)/\2/g'' | sed -e''s/\(\"\, *\"\)/',delim,'/g'' | sed -e''s/\(time = \"\)\(.*\)\(\" \;\}\)/\2/g''']);

   if h ~= 0; 
     disp(['CHECK YOU ARE USING VERSION 4 OR HIGHER NETCDF LIBRARIES'])
     error(['Could not extract CF convention netcdf time data: error ',num2str(h)]);
   end 

   % Cleanup
   delete([ncfile,'.cdl'],[ncfile,'.nc'])

else

   error('Can only use nccftime with Linux, UNIX or Mac OS X systems')

end

% find the occurrence of the delimiter
xpos=strfind(t,delim);

% case of single date only
if isempty(xpos); 
 timestrings{1}=t;
else  
 timestrings{1}=t(1:xpos(1)-1);
 for i=2:length(xpos)
  timestrings{i}=t(xpos(i-1)+length(delim):xpos(i)-1);
 end
 timestrings{end+1}=t(xpos(end)+length(delim):end-1);
end

% Add in colon so that datenum can interpret the date, using strcat
% to elliminate the presense of new line characters (\n)
for i=1:length(timestrings)
 if length(char(timestrings{i}))>11;
   timestrings{i}=[strcat(char(timestrings{i})),':']; 
 end
end

% Complete final stage of assigning converted time
time=datenum(timestrings);

if debug; disp(['...Leaving ',mfilename]); end

