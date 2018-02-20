
 % plot on figure 2
 figure(2)
 clf

 % plot VR on a log color scale
 contours=10.^[floor(log10(vr)):1:ceil(log10(vr))];
 cmap=colormap(parula(ceil(log10(vr))-floor(log10(vr))+1));
 [zindex,truecol]=ridgepack_colorindex(vr,contours,0);
 surface(epsilon,phi,-0.1*ones(size(vr)),truecol,'EdgeColor','k')
 shading flat
 hold on

 % plot the dilation field over the top as streamlines
 hstr=streamslice(epsilonsplit,phisplit,d1,d2);
 set(hstr,'Color',[1 1 1])

 % only plot up to a min strain of -0.98
 idx=find(strainp>=-0.96);
 strainp=strainp(idx);
 phip=phip(idx);
 pep=pep(idx);

...
...

