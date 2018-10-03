function Ridgepack_RASM_sea_ice_make_speed(rasmcases,fieldu,fieldv,hourly)

% set cases
home=getenv('HOME');
dirdata=[home,'/data'];

if nargin<4
 hourly=false
end

for j=1:length(rasmcases)

 rasmcase=char(rasmcases{j});

 if hourly
  dircase=['/Volumes/RobertsRaid3/work/processing/',rasmcase,'/ice/hourly'];
  delim='.cice.h_inst.';
 else
  dircase=['/Volumes/RobertsRaid3/work/processing/',rasmcase,'/ice/monthly'];
  delim='.cice.h.';
 end

 cd(dircase)

 fileu=[rasmcase,delim,fieldu];
 filev=[rasmcase,delim,fieldv];
 files=[rasmcase,delim,'speed'];

 nctimeu=ridgepack_clone(fileu,'time');
 nctimev=ridgepack_clone(filev,'time');

 if length(nctimeu.time.data)==length(nctimev.time.data)
  timediff=(nctimeu.time.data(:)-nctimev.time.data(:));
 else
  error(['Time length different for ',fileu,' and ',filev])
 end

 if any(timediff~=0)
  error(['Time different for ',fileu,' and ',filev])
 end

 % attempt to find speed record
 try
  nctimes=ridgepack_clone(files,'time');
  generate=false;
 catch
  generate=true;
 end

 if ~generate
  if length(nctimeu.time.data)==length(nctimes.time.data)
   timediff=(nctimeu.time.data(:)-nctimes.time.data(:));
  else
   generate=true;
  end  

  if any(timediff~=0)
   generate=true;
  end
 end

 if generate

  disp(['Generating speed file ',files])

  for i=1:length(nctimeu.time.data)

   ncu=ridgepack_clone(fileu,fieldu,{'time'},{i});
   ncv=ridgepack_clone(filev,fieldv,{'time'},{i});

   ncspeed=ncu;  
   ncspeed.speed=ncu.(fieldu);
   ncspeed=rmfield(ncspeed,fieldu);
   ncspeed.speed.data=sqrt(ncu.(fieldu).data.^2 + ncv.(fieldv).data.^2);

   if i==1
    ridgepack_write(ncspeed,files);
   else
    ridgepack_write(ncspeed,files,{'time'},{0});
   end

  end

 else

  disp([files,' already exists.'])

 end

end


