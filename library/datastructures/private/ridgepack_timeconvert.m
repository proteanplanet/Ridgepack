function [newtime,units,CFcalendar]=ridgepack_timeconvert(time,units,CFcalendar,tofrom)

% ridgepack_timeconvert - Convert CF convention time to serial time
%
% function [newtime,units,CFcalendar]=ridgepack_timeconvert(time,units,CFcalendar,tofrom)
%
% Within icepack, nc structures carry around the calendar type from which the data
% was converted as a string in nc.time.units, but the actual time values in
% nc.time.data are in Matlab serial units of "days since 0000-01-01 00:00:0.0" in 
% the Proleptic Gregorian calendar. This function carries out the latter conversion
% from the original calendar to Matlab serial UTC time, and converts back, depending 
% on certain constraints discussed here.
% 
% INPUT:
%
% time       - time given in the CF convention units
% units      - CF convention units
% CFcalendar - type of calendar being used (e.g. gregorian, noleap etc)
% tofrom     - logical for converting to (0) or from (1) serial time units 
%              to netcdf units.
%
%
% OUTPUT:
%
% newtime    - time given in Matlab serial time units
% units      - CF convention units if changed where tofrom=1
% CFcalendar - Calendar output for cases where conversion was not possible for tofrom=1
%
% The Gregorian Calendar started on October 15, 1582 as declared by Pope Gregory XIII, 
% going from the Julian Calendar, for which day 0 is Midday, January 1, 4712 B.C.. 
% Conversion to the Gregorian Calendar skipped from October 4 to October 15, 1582,
% to make up for drift of the Julian Calendar relative to the equinox. This annoying
% fact sometimes makes calendar conversions to the internal Matlab calendar, which
% is proleptic_gregorian, prone to error. The proleptic_gregorian calendar assumes
% the Gregorian calendar prior to October 15, 1582. For this reason, all conversions
% of CF calendar types 'gregorian' (also known as 'standard') and 'no_leap' (otherwise
% known as '365_day') are performed by udunits prior to 1582-10-15 via an interface
% with the ncgen/ridgepack_dump utilities that need to be installed on the computer being 
% used to run this Matlab script. More complicated conversions from 'all_leap' (otherwise
% known as '366_day'), '360_day' (12 x 30-day months) and 'julian' calendars are always
% performed using udunits. Conversion from these calendars to the Matlab
% proleptic_gregorian is subject to possible innaccuracy for dates sandwiched between
% October 4 and October 15, 1582, and across the equivalent Proleptic Gregorian dates
% February 29 0000Z to March 2 00000Z. For this reason, Icepack allows conversions
% from other calendars to the Matlab proleptic_gregorian calendar, but does not allow
% an inverse conversion where the innaccuracy may perpetuate the error. In summary:
%
% Times converted to Matlab serial time (Proleptic Gregorian) without warning:
% 'proleptic_gregorian' - all cases
% 'gregorian','standard','no_leap','365_day' - where the minimum time record >= 1582-10-15
% [These conversions are conducted internally in this function, as with their inverse]
%
% Times converted to Matlab serial time (Proleptic Gregorian) with warning:
% 'all_leap','366_day','360_day','julian' 
% 'gregorian','standard','no_leap','365_day' - where the minimum time record < 1582-10-15
% [These conversions are conducted externally in udunits via a call to ridgepack_cftime]
%
% Times not converted to Matlab serial time:
% 'none' - a preset CF calendar setting or assigned during read in of some files
% 'gregorian','standard','no_leap','365_day' where seconds, minutes, hours or days not in units
%
% Inverse time conversions from Matlab serial time not allowed:
% 'gregorian','standard','no_leap','365_day' - where the minimum time record < 1582-10-15
% 'all_leap','366_day','360_day','julian' - no inverse conversion allowed
%
% For more information please reference CF conventions for netcdf and udunits. Note that
% the ridgepack_dump and ncgen utilities are used by way of a call to a seperate function 
% ridgepack_cftime.m as part of the icepack package.  These calls are made without need for mexnc
% files, but require installation of Version 4 or later netcdf libraries on the 
% Unix, Linux and Mac OS X system being used. As already alluded to, calls to this 
% external function are not made for calendars 
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
%

global debug;
%debug=true;

if debug; disp(['Entering ',mfilename,'...']); end

if debug; disp(['Calendar is ',CFcalendar]); end

if nargin<3; error('Incorrect inputs'); end

time=double(time);

if any(strcmpi(CFcalendar,{'standard','gregorian','proleptic_gregorian','noleap','365_day'}))

   baseunits={'seconds','minutes','hours','days'};
   basefactor=[1/86400 1/1440 1/24 1];
   basetime=[];

   for i=1:length(baseunits)
    compstr=[char(baseunits{i}),' since '];
    if strcmpi(units(1:min(length(units),length(compstr))),compstr)
        basestring=units(length(compstr)+1:end);
      
        % shorten the basestring for stipulations like: "2009-10-23 0"
        if strcmp(basestring(end-1:end),' 0') & ...
           (length(basestring)<=12 & length(basestring)>=11) 
         basestring=basestring(1:end-2);
        end

	% extract correction from UTC suffix
	utcpos=strfind(basestring,' +');
	utcneg=strfind(basestring,' -');
        if length(utcpos)==1 
	 utccorrection=[basestring(utcpos+2:end),':00'];
         basestring=basestring(1:utcpos-1);
        elseif length(utcneg)==1
	 utccorrection=[basestring(utcneg+2:end),':00'];
         basestring=basestring(1:utcneg-1);
	else
	 utccorrection=[];
        end

	if debug; disp(['Base time string is ',basestring]); end

        if length(basestring)<10
         error('Base time string does not include yyyy-mm-dd')
        end

        % set basetime
        try

	 yearplus=false;

         year=str2num(basestring(1:4));

         % grab month, making allownce for this possibly being single digit
         if strcmp(basestring(7),'-')
          month=str2num(basestring(6));
          nextdigit=8;
         else
          month=str2num(basestring(6:7));
          nextdigit=9;
         end

         % grab day, making allownce for this possibly being single digit
         if strcmp(basestring(nextdigit+1),'-')
          day=str2num(basestring(nextdigit));
          nextdigit=nextdigit+1;
         else
          day=str2num(basestring(nextdigit:nextdigit+1));
          nextdigit=nextdigit+2;
         end

         % move past space(s) between date and time
         if length(basestring)>=nextdigit
          while strcmp(basestring(nextdigit),' ')
           nextdigit=nextdigit+1;
          end
         end

         % grab hours, if given
         if length(basestring)>=nextdigit
          hour=str2num(basestring(nextdigit:nextdigit+1));
          nextdigit=nextdigit+3;
         else
          hour=0;
         end

         % grab minutes, if given
         if length(basestring)>=nextdigit
          minute=str2num(basestring(nextdigit:nextdigit+1));
          nextdigit=nextdigit+3;
         else
          minute=0;
         end

         % grab seconds, if given
         if length(basestring)>=nextdigit
          second=str2num(basestring(nextdigit:end));
         else
          second=0;
         end

         basetime=datenum(year,month,day,hour,minute,second);

         if basetime==0 | strcmp(basestring,'0000-00-00 00:00:00') 

	  disp(['Base time string is ',basestring]); 
	  basetime=datenum(input('Enter the new base time (mm-dd-yyyy): ','s'));

	 elseif strcmp(basestring(1:4),'0000')

          % This section is to take into account when only the year is set to 0000 which
          % matlab has trouble dealing with, so we add one year and take it away later.
	  x=datevec(basetime);
	  x(1)=0001;
	  basetime=datenum(x);
	  basestring=['0001',basestring(5:end)];
	  units=[compstr,basestring];
	  yearplus=true;
	  if debug; disp(['Reseting time base to ',datestr(basetime)]); end

         end

        catch
         disp(['ERROR: processing ', basestring]);
         error('Unable to determine the base time')
        end

        % add correction for difference from UTC
	if length(utcpos)==1
         basetime=basetime-(datenum(utccorrection)-datenum([year(now) 1 1]));
        elseif length(utcneg)==1
         basetime=basetime+(datenum(utccorrection)-datenum([year(now) 1 1]));
        end

        disp(['Converted base time string is ',datestr(basetime)])

	% Work on specifics of calendars that have been specified.  
        % This program as the capability to convert to the gregorian calendar, 
        % also referred to as the "standard" calendar in CF convention, only if the 
        % reference date, and hence the earliest date in the time series,
        % is after October 15, 1582.  Since the proleptic_gregorian calendar is not 
        % affected by the switch over from Julian to Gregorian Calendar, 
        % any time can be used in that case.

        if (any(strcmpi(CFcalendar,{'gregorian','standard'})) && ...
            basetime>=datenum([1582 10 15])) | ...
            any(strcmpi(CFcalendar,{'proleptic_gregorian','noleap','365_day'}))

          disp(['Processing ',CFcalendar,' calendar internally in Icepack'])

          if tofrom

           if debug; 
             disp(['Converting time from serial to ',[compstr,datestr(basetime)],' UTC']); 
           end

	   if strcmpi(CFcalendar,'noleap')

            for j=length(time):-1:1
             x=datevec(time(j));
	     if x(2)==2 & x(3)==29 % check that a leap year is not in the time record
              error(['time includes ',datestr(time(j)),...
                     ' but the calendar is ',CFcalendar])
             end
             x=datevec(time(j));
             if leapyear(x(1))
              leapday=datenum(x(1),2,29);
              if j>1 & time(j)>leapday & time(j-1)<leapday
               time(j:end)=time(j:end)-1;
              end
             end
            end

	    for j=time(1)-1:-1:basetime
             x=datevec(j-1);
	     if x(2)==2 & x(3)==29 % leap year before time count (count backwords)
              time=time-1;
             end
            end

           end

           newtime=(time-basetime)/basefactor(i);

          elseif ~tofrom

           if debug; 
            disp(['Converting time to serial from ',[compstr,datestr(basetime)],' UTC']); 
           end

           if any(time>1.E36)
            error('There is a problem with the time data - out of bounds')
           end

           newtime=basetime+basefactor(i)*time;

	   % correct for the case where basetime includes year 0000 by subtracting 
           % one year from the date.
           if yearplus; 
	    if debug; disp('Correcting for year 0000 in time units'); end
	    newvec=datevec(newtime(1)); 
	    newvec(1,1)=newvec(1,1)-1; 
            newtime=newtime-(newtime(1)-datenum(newvec));
	   end 

	   % remove 29th of Feb from leap year calendar
	   if strcmpi(CFcalendar,'noleap')

	    if debug; disp('Adjusting for noleap calendar'); end

	    % Calculate number of leap years between the basetime and the initial 
            % starting time, and add these to the start time.  This must be performed 
            % recursively until no more leap years exist within this time period.
	    leapsexist=true;
            baset=basetime;
            while leapsexist
	     x=datevec(baset:1:newtime(1)-1);
             xt=datenum(x);
	     leapindex=find(x(:,2)==2 & x(:,3)==29);
	     if isempty(leapindex)
	      leapsexist=false;
             else
	      baset=newtime(1);
	      newtime=newtime+length(leapindex);
             end
            end

            % Add a day to the time after each Feb 29th is encountered within
            % the times encountered in the netcdf record.  This must also be 
            % performed recursively.
            
            % advance all times if first day is leap day
	    x=datevec(newtime);
	    if (x(1,2)==2 & x(1,3)==29)
              newtime=newtime+1;
            end

            % now search for leap days through the timeseries
            for j=2:length(newtime)
	     x1=datevec(newtime(j-1));
	     x2=datevec(newtime(j));

             % find leap years between on within two consecutive time records
             yearrange=[x1(1):1:x2(1)];
             leapflag=leapyear(yearrange);

             % data in the same leap year
             if length(yearrange(leapflag))==1 
              leapday=datenum(yearrange(leapflag),2,29);
              if newtime(j-1)<leapday & newtime(j)>=leapday
               newtime(j:end)=newtime(j:end)+1;
              end
             elseif length(yearrange(leapflag))>1 
              error('Unable to process noleap times spanning more than one leap year jump')
             end
            end  

           end % noleap

          else

	   error('Specifying the conversion direction for time')

          end

        else


         if tofrom

          disp(['WARNING: Unable to convert from Matlab proleptic_gregorian to ',CFcalendar,...
                ' for basetime ',datestr(basetime),' UTC']);
          disp(['WARNING: Changing ',CFcalendar,' to proleptic_gregorian with possible errors']);
          newtime=time-fix(min(time));
          CFcalendar='proleptic_gregorian';
          units=['days since ',datestr(min(time),'yyyy-mm-dd')];

         elseif ~tofrom

          disp(['Processing ',CFcalendar,' calendar units ',...
               [compstr,datestr(basetime)],' UTC using ridgepack_dump/ncgen/udunit'])

	  newtime=ridgepack_cftime(time,units,CFcalendar);

          if min(newtime)<datenum([1582 10 15])
           disp(['WARNING: Conversion from ',CFcalendar,' to Matlab proleptic_gregorian'])
           disp(['         may give errors for times preceding October 15 1582.'])
          end

	 end

        end
    end

   end

   if isempty(basetime)

     disp(['ERROR: Unable to find any of these time units: ',ridgepack_cellcat(baseunits)])
     error(['Unable to find any of these time units: ',ridgepack_cellcat(baseunits)])

   end

elseif any(strcmpi(CFcalendar,{'all_leap','366_day','julian','360_day'}))

   if tofrom

    disp(['WARNING: Unable to convert back from Matlab proleptic_gregorian to ',CFcalendar]);
    disp(['WARNING: Changing ',CFcalendar,' to proleptic_gregorian with possible errors']);
    newtime=time-fix(min(time));
    CFcalendar='proleptic_gregorian';
    units=['days since ',datestr(min(time),'yyyy-mm-dd')];

   elseif ~tofrom

    disp(['Processing ',CFcalendar,' calendar externally using ncgen/ridgepack_dump/udunits'])
    newtime=ridgepack_cftime(time,units,CFcalendar);
    disp(['WARNING: Conversion from ',CFcalendar,' to Matlab proleptic_gregorian may give errors '])

   end

elseif any(strcmpi(CFcalendar,{'none'}))

   disp(['WARNING: No calendar conversion']);
   newtime=time;

else

   disp(['WARNING: Unable to find an appropriate calendar specification'])
   disp(['WARNING: No calendar conversion']);
   newtime=time;

end


if debug; disp(['...Leaving ',mfilename]); end

