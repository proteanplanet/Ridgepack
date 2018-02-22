function [nc]=ridgepack_timesubset(ncfile,ncvar,months,yearstart,yearend)

% ridgepack_timesubset - Generate seasonal subsets of data based on consecutive months
%
% function [nc]=ridgepack_timesubset(ncfile,ncvar,months,yearstart,yearend)
%
% This function extracts only data from specified months from a netcdf dataset
% and writes them to another netcdf file.  It then extracts the mean, sample size,
% and standard deviation from just those months specified.  The months must 
% be consecutive, such as April, May, June, or November, December, January. As in
% the latter example, the months can cross between years.  If this happens,
% The years sampled pass are staggered. In other words, if one wishes to construct
% a mean of November through January for 1991 to 1996, then November and December
% would be sampled from 1991 to 1995, and January would be sampled from 1992 to 1996
% to construct the seasonal mean.
%
% INPUT:
%
%
% ncfile    - character variable giving the name of the netcdf file to be sampled.
% ncvar     - character variable of variable for which seasonal means are required.
% months    - vector of numbers 1 through 12 for January through December indicating
%             months to be extracted.
% yearstart - first year of timeseries to be analyzed.
% yearend   - last year of timeseries to be analyzed.
%
% 
%
% OUTPUT:
%
%
% nc - netcdf structure containing the mean values for each entry in ncvar, as well
%      as the sample size, and the standard deviation.
% 
% 
% File output:
% 
% Two files are created when running this program if they do not already exist. 
% The first file is a subset of ncfile containing just records from the months
% specified.  The second file is the mean, sample size and standard deviation of 
% these files. If these files already exist in the operating director, they are 
% not re-generated, but instead they are read to give the output netcdf structure.
%
% Note on time representation:
%
% This function searches not only for a time variable, but also for a time_bounds
% variable.  If this is found, then the mid-date of the time bounds is used as the 
% indicator of the month of the data. If the time bounds spans more than 31 days,
% then the function exits with an error, since the samples already represent a mean 
% of more than a month.
%
% Example:
%
% To obtain the winter (December to Feburary) mean of all time-dependent variables in
% the netcdf dataset timeseries.nc for 1990 to 1992, one would 
% enter this command:
%
% nc=ridgepack_timesubset('timeseries.nc',{},[1 2 12],1990,1992);
%
% This would generate a mean for Dec/Jan/Feb for the 1989/1991 (Dec) and 1990/1992 (Jan/Feb)
% winters and this would be provided in the output netcdf structure. It would also create 
% two files:
% timeseries.cice.h.aice.months_1_2_12_1990_1992.nc
% and 
% timeseries.cice.h.aice.months_1_2_12_1990_1992.mean.nc
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% check input
if nargin<5
 error('Missing input information')
elseif ~ischar(ncfile);
 error([inputname(1),' is not a netcdf file name']);
elseif ~ischar(ncvar);
 error('Variable (2nd argument) is not a character variable');
elseif isempty(ncvar);
 error('Variable (2nd argument) name is empty');
elseif nnz(isstrprop(ncfile, 'alphanum'))==0;
 error('No name given for the netCDF file');
elseif ~isnumeric(months)
 error('months must be a vector')
elseif length(months)==0 | length(months)>12
 error('only specify up to 12 months in month')
elseif min(months)<1 | max(months)>12
 error('months can only be from 1 to 12')
elseif ~isnumeric(yearstart) | ~isnumeric(yearend)
 error('yearstart and yearend must be numbers')
end

% Check the variable exists in the netcdf file
nc=ridgepack_clone(ncfile,ncvar,1);

% remove ".nc" from filename if need be
ui=length(ncfile); li=max(1,ui-2);
if strcmp(ncfile(li:ui),'.nc')
 ncfile=ncfile(1:li-1);
end 

% sort months and make sure they are consecutive
months=sort(months);

% find out if months cross a year boundary
if any(abs(diff(months))>1)
 idx=find(diff(months)>1);
 % check that, for example, the user hasn't specified [11 1 2]
 % instead of [12 1 2] (the former has a one-month gap in December
 if length(idx)>1 | (months(end)-months(1))<11
  error(['months must be consecutive: ',num2str(months)])
 end 
 months1=months(1:idx);
 yearstart1=yearstart;
 yearend1=yearend;
 months2=months(idx+1:end);
 yearstart2=yearstart-1;
 yearend2=yearend-1;
else
 months1=months;
 yearstart1=yearstart;
 yearend1=yearend;
 months2=[];
 yearstart2=[];
 yearend2=[];
end

% get time data from netcdf file
try 
 nc=ridgepack_clone(ncfile,{'time_bounds'});
 times=datevec(mean(nc.time_bounds.data,1));
 timebounds=true;
catch
 disp('Using absolute time values instead of time bounds')
 nc=ridgepack_clone(ncfile,{'time'});
 times=datevec(nc.time.data);
 timebounds=false;
end

if timebounds && any(diff(nc.time_bounds.data,1)>31)
 error('Time bounds exceed one month')
end

minyear=min(times(:,1));
maxyear=max(times(:,1));

if minyear>yearstart | maxyear<yearend
 error('data does not span the requested years')
end

% obtain indices of months to be extracted
index=[];
fileout=[ncfile,'.months'];
for i=months1
 index=[index; find(times(:,1)>=yearstart1 & times(:,1)<=yearend1 & times(:,2)==i)];
 fileout=[fileout,'_',num2str(i,'%2.2i')];
end
if ~isempty(months2)
 for i=months2
  index=[index; find(times(:,1)>=yearstart2 & times(:,1)<=yearend2 & times(:,2)==i)];
  fileout=[fileout,'_',num2str(i,'%2.2i')];
 end
end
index=sort(index)';

% set fileout name, taking into account bridging across years
if ~isempty(yearstart2) & ~isempty(yearend2)
 ys=min(yearstart1,yearstart2);
 ye=max(yearend1,yearend2);
elseif isempty(yearstart2) & isempty(yearend2)
 ys=yearstart1;
 ye=yearend1;
else
 error('yearstart2 and yearend2 are inconsistent')
end
fileout=[fileout,'_',num2str(ys),'_',num2str(ye),'.',ncvar];

% check whether mean file already exists
finalfile=[fileout,'.mean'];
xdir=dir([finalfile,'.nc']);

% set variable names
stdvar=[ncvar,'_std'];
sampvar=[ncvar,'_samp'];
equivvar=[ncvar,'_equiv'];

% extract mean, standard deviation and sample size
if ~isempty(xdir)
 disp(['Reading ',finalfile])
 try
  nccheck=ridgepack_clone(finalfile);
  if ~isfield(nccheck,ncvar) | ~isfield(nccheck,stdvar) | ...
     ~isfield(nccheck,sampvar) | ~isfield(nccheck,equivvar) 
   recreatemean=true;
  else
   recreatemean=false; 
  end
 catch
  disp(['Must regenerate ',finalfile])
  recreatemean=true;
 end
else
 recreatemean=true;
end

if recreatemean

 % check whether abridged seasonal timeseries already exists
 xdir=dir([fileout,'.nc']);

 % extract data if necessary and write to netcdf file
 if isempty(index);
  error(['Unable to extract these months from ',ncfile]);
 elseif ~isempty(xdir)
  % check times are correct in the file
  disp([char(xdir.name),' already exists.'])
  nctimes=ridgepack_clone(fileout,{'time'});
  generate=false;
  if length(nctimes.time.data)~=length(index)
   disp(['Records are missing from ',fileout])
   generate=true;
  else
   generate=false;
   for t=1:length(nctimes.time.data)
    if isempty(find(nctimes.time.data(t)==nc.time.data))
     generate=true;
     datestr(nctimes.time.data(t))
     disp(['Times are not correct in ',fileout,'...regenerating file'])
     break
    end
   end
  end
 else
  generate=true;
 end

 if generate
  for i=index
   nc=ridgepack_clone(ncfile,ncvar,{'time'},{i}); % this specification sucks in time_bounds
   if i==index(1)
    nc.attributes.comment=['Reduced dataset only;',...
                           ' Months: ',num2str(months,' %2.2i'),...
                           ', Years: ',num2str(ys),'-',num2str(ye)];
    ridgepack_write(nc,fileout);
   else
    ridgepack_write(nc,fileout,{'time'},{0});
   end
  end
 end

 % extract mean, standard deviation and sample size
 disp(['Reading ',fileout])
 nc=ridgepack_reduce(ridgepack_clone(fileout),{'time'},{},true);
 nc.attributes.comment=['Seasonal mean dataset only; ',...
                        ' Months: ',num2str(months,' %2.2i'),...
                        ', Years: ',num2str(ys),'-',num2str(ye)];
 ridgepack_write(nc,finalfile);

else

 nc=nccheck;

end

