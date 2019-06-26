
clear
clf

startyear=31;
endyear=60;
cd('/Users/afroberts/SIhMSatArray/E3SM/v0highres/monthly')
%cd('/Volumes/MacBookProT3SeaIce/E3SM/highresv0/monthly')
resname='HR_v0'
setname='state';

%for month=[3,9]
for month=[1,2,3,4,5,6,7,8,9,10,11,12];

 k=0

 for year=startyear:endyear
  k=k+1;
  dncai=['aice.b1850c5_acmev0_highres.cice.h.',...
           num2str(year,'%4.4i'),'-',...
           num2str(month,'%2.2i')];
  dnchi=['hi.b1850c5_acmev0_highres.cice.h.',...
           num2str(year,'%4.4i'),'-',...
           num2str(month,'%2.2i')];
  dnchs=['hs.b1850c5_acmev0_highres.cice.h.',...
           num2str(year,'%4.4i'),'-',...
           num2str(month,'%2.2i')];
  dncuvel=['uvel.b1850c5_acmev0_highres.cice.h.',...
           num2str(year,'%4.4i'),'-',...
           num2str(month,'%2.2i')];
  dncvvel=['vvel.b1850c5_acmev0_highres.cice.h.',...
           num2str(year,'%4.4i'),'-',...
           num2str(month,'%2.2i')];

  ncai=ridgepack_clone(dncai);
  nchi=ridgepack_clone(dnchi);
  nchs=ridgepack_clone(dnchs);
  ncuvel=ridgepack_clone(dncuvel);
  ncvvel=ridgepack_clone(dncvvel);

  if k==1

   nc=ncai;
   nc.hi=nchi.hi;
   nc.hs=nchs.hs;
   nc.uvel=ncuvel.uvel;
   nc.vvel=ncvvel.vvel;

   clear ncai nchi nchs ncuvel ncvvel
 
  else

   ncnew=ncai;
   ncnew.hi=nchi.hi;
   ncnew.hs=nchs.hs;
   ncnew.uvel=ncuvel.uvel;
   ncnew.vvel=ncvvel.vvel;
 
   nc.aice.data=nc.aice.data+ncnew.aice.data;
   nc.hi.data=nc.hi.data+ncnew.hi.data;
   nc.hs.data=nc.hs.data+ncnew.hs.data;
   nc.uvel.data=nc.uvel.data+ncnew.uvel.data;
   nc.vvel.data=nc.vvel.data+ncnew.vvel.data;

   clear ncnew ncai nchi nchs ncuvel ncvvel

  end
 end

 nc.aice.data=nc.aice.data/k;
 nc.hi.data=nc.hi.data/k;
 nc.hs.data=nc.hs.data/k;
 nc.uvel.data=nc.uvel.data/k;
 nc.vvel.data=nc.vvel.data/k;

 nc=rmfield(nc,'attributes');

 nc.attributes.title=['Mean years ',num2str(startyear),...
                      ' to ',num2str(endyear)];

 outfilename=['b1850c5_acmev0_highres.',resname,'.',setname,...
              '.mean.',num2str(month,'%2.2i'),'.',...
              num2str(startyear,'%4.4i'),'_',...
              num2str(endyear,'%4.4i')];

 ridgepack_write(nc,outfilename);

end


