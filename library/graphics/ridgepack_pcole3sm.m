function ridgepack_pcole3sm(nc,var,ncvert,cont,mask,loglin,ref,horiz,colors,colvals)

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% check that number of inputs is sufficient
if nargin<4
 error('insufficient information for plot')
end

% check for mask
if (nargin>=5 & isempty(mask)) | nargin<5
 mask=ncvert.nCells.data;
end

if nargin<6;
 loglin='linear' ;
elseif isempty(loglin)
 loglin='';
elseif ~(ischar(loglin) | isnumeric(loglin))
 loglin
 error('loglin must be a character input')
end

if nargin<7;
 ref=0.0;
elseif ~isnumeric(ref) || length(ref)~=1
 ref
 error('ref must be a single number')
end

if nargin<8;
 horiz='vertical';
elseif ~ischar(horiz)
 horiz
 error('horiz must be a character input')
end

if nargin<9;
 colors='bluered' ;
elseif ~ischar(colors)
 colors
 error('colors must be a character input')
end

% get current axes
hmap=get(gcf,'CurrentAxes');
if ~ismap(hmap)
 error('Current axes must be a map')
end
figout=get(hmap,'OuterPosition');
fontsize=min(11,max(8,10*sqrt((figout(3).^2)+(figout(4).^2))/sqrt(2)));
set(hmap,'fontsize',fontsize);

% set colorscheme
if strcmp(loglin,'linear')
 ridgepack_colormap(cont,ref,colors);
else
 ridgepack_colormap(cont,ref,colors,true);
end

% plot patches for each contour interval
for j=1:length(cont)

 if j==1
  idx=find(nc.(var).data<cont(j+1));
 elseif j==length(cont)
  idx=find(nc.(var).data>=cont(j-1));
 else
  idx=find(nc.(var).data>=cont(j) & nc.(var).data<cont(j+1));
 end
 idxn=intersect(idx,mask);

 [zindex,truecolor]=ridgepack_colorindex(nc.(var).data(idxn),cont,ref);

 if length(idxn)>0

      lat=zeros(length(idxn),8);
      lon=zeros(length(idxn),8);

      for i=1:length(idxn)

       maxidx=ncvert.nEdgesOnCell.data(idxn(i));

       lat(i,1:maxidx)=ncvert.latitude.data(...
               ncvert.verticesOnCell.data(...
                1:maxidx,idxn(i)))*180/pi;
       lon(i,1:maxidx)=ncvert.longitude.data(...
               ncvert.verticesOnCell.data(...
                1:maxidx,idxn(i)))*180/pi;
       
       lon(i,maxidx+1:8)=lon(i,1);
       lat(i,maxidx+1:8)=lat(i,1);

      end

      [cc,dd] = mfwdtran(gcm,lat(:,:),lon(:,:));

      patch(cc',dd',truecolor(1,:),'EdgeColor','none')

 end

 clear zindex truecolor cc dd lon lat idxn idx

end

drawnow


