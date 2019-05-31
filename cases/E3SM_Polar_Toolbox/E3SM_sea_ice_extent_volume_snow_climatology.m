

clear
clf

for months=1:12;

 month=char(mons(months));

 if i==1
  ncnorth=ridgepack_clone(['G02202_v3_merged_conc_north_1979_1999_',month,'_mean.nc'])
  ncsouth=ridgepack_clone(['G02202_v3_merged_conc_south_1979_1999_',month,'_mean.nc'])
  ncnorth.month.data(i)=i;
  ncsouth.month.data(i)=i;
  concnorth(i)=ncnorth.conc.data
 else
  ncnorth=ridgepack_clone(['G02202_v3_merged_conc_north_1979_1999_',month,'_mean.nc'])
  ncsouth=ridgepack_clone(['G02202_v3_merged_conc_south_1979_1999_',month,'_mean.nc'])
  
 end


end





