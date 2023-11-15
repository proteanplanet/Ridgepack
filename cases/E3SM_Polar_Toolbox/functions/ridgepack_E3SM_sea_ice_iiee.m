
clear
close all

startyear=1980;
endyear=2014;

monthsets={[1 2 3]};
%monthsets={[1 2 3],[4 5 6]};
%monthsets={[1 2 3],[4 5 6],[7 8 9],[10 11 12]};

hemisphere='nh';
%hemisphere='sh';

nruns=[1];

ensnames={'v2.NARRM.historical_0101.mpassi.daily',...
          'v2.NARRM.historical_0301.mpassi.daily',...
          'v2.LR.historical_0101.mpassi.daily',...
          'v2.LR.historical_0151.mpassi.daily'};

locations={'/Users/afroberts/data/MODEL/E3SM/v2/v2.NARRM.historical_0101/data/ice/processed',...
           '/Users/afroberts/data/MODEL/E3SM/v2/v2.NARRM.historical_0301/data/ice/processed',...
           '/Users/afroberts/data/MODEL/E3SM/v2/v2.LR.historical_0101/data/ice/processed',...
           '/Users/afroberts/data/MODEL/E3SM/v2/v2.LR.historical_0151/data/ice/processed'};

rnametags={'NARRM 101','NARRM 301','LR 101','LR 151'};

mons={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

if strcmp(hemisphere,'nh')
 hem=1;
elseif strcmp(hemisphere,'sh')
 hem=-1;
end

alphatag=['abcdefghijklmnopqrstuvwxyz'];

for ms=1:length(monthsets)
 for ns=nruns

  monthnum=['months'];
  for months=monthsets{ms};
   monthnum=[monthnum,'_',num2str(months,'%2.2i')];
  end
  yeartag=['_',num2str(startyear,'%4.4i'),'_',num2str(endyear,'%4.4i')];

  obslocation='/Users/afroberts/data/data/SATELLITE/processed/G02202_v4/';
  if strcmp(hemisphere,'nh')
   %infileobs=['G02202_v4_merged_north_r00_1979_2022.',monthnum,yeartag,'.conc'];
   infileobs=['G02202_v4_merged_north_r00_1979_2022'];
  elseif strcmp(hemisphere,'sh')
   %infileobs=['G02202_v4_merged_south_r00_1979_2022.',monthnum,yeartag,'.conc'];
   infileobs=['G02202_v4_merged_south_r00_1979_2022'];
  end

  infileobsnh=['G02202_v4_merged_north_r00_1979_2022'];
  infileobssh=['G02202_v4_merged_south_r00_1979_2022'];

  cd(obslocation)

  nctimenh=ridgepack_clone(infileobsnh,'time');
  yearsobsnh=str2num(datestr(nctimenh.time.data,'yyyy'));

  nctimesh=ridgepack_clone(infileobssh,'time');
  yearsobssh=str2num(datestr(nctimesh.time.data,'yyyy'));

Linearly interpolate across missing records

  cd(char(locations{ns}))

  %for nyear=[startyear:endyear]
%  for nyear=[startyear]

%   cd(obslocation)
%   idxonh=find(nyear==yearsobsnh);
%   ncobsnh=ridgepack_clone(infileobsnh,'extent',{'time'},{idxo(1)},{idxo(end)});

   %cd(char(locations{ns}))
   %if strcmp(hemisphere,'nh')
   % infilemod=[char(ensnames{ns}),'.nh.regrid.',monthnum,'.years',yeartag];
   %elseif strcmp(hemisphere,'sh')
   % infilemod=[char(ensnames{ns}),'.sh.regrid.',monthnum,'.years',yeartag];
   %end

   %ncmod=ridgepack_clone(infilemod,'conc',datestr(nctime.time.data(i)));

%  end

 end
end






