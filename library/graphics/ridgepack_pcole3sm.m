function ridgepack_pcole3sm(nc,var,ncvert,mask,cont,loglin,ref,horiz,colors,colvals)

% ridgepack_pcole3sm - Color fills data from E3SM/MPAS on an unstructured mesh
%
% function ridgepack_pcole3sm(nc,var,ncvert,mask,cont,loglin,ref,horiz,colors,colvals)
%
% This function generates color shaded maps of E3SM fields from MPAS components
% 
% INPUTS:
%
% nc     - A netcdf structure (see ridgepack_struct for more details). The
%          structure must contain the fields "latitude" and "longitude". 
%
% var    - Character variable giving the variable in nc to be shaded.
%
% ncvert - Netcdf grid structure from E3SM. It must include the vectors:
%          nCells, nEdgesOnCell, verticesOnCell, latitude, longitude
%
% mask   - mask of cell indices to be plotted.
%
% cont   - contour range entered as a vector [C1,C2,...,CX]. This may
%          be entered as an empty vector [] to get ncpcolor to choose
%          the contour interval for you, or it may be omitted if
%          no other proceeding arguments are required. (optional)
%
% loglin - 'linear' for linear color scaling, or for log scaling the
%          minumum order of magnitude to be plotted. When using this 
%          log option, you need to make three specifications.  Firstly
%          you need to specify the contour range differently than for
%          for a linear plots.  The numbers in the contour range correspond
%          to order of magnitude multipled by the base order of magnitude
%          which is entered in as loglin.  The sign of the contour values
%          indicates the actual sign of the contour bounds.  For example,
%          say you want to specify a range of -1.e-7 to 1.e-3 on a log10
%          scale in intervals of power of 10, with the smallest value
%          represented being +/- 10-9, then you would enter:
%          cont=[-7:1:3] and loglin=10^-9. 
%          You can adjust the reference value to emphasize color contrast
%          where it is most desired.  If, in the above example, you would
%          like to adjust the spacing to indicate minor logarithmic tick
%          marks on the color scale, you would change cont to:
%          cont=[-7:0.1:3] and loglin=10^-9. (optional)
%
% ref    - reference point about which the color is centered, and, when
%          the minimum or maximum is equal to ref, all points below or 
%          abover this value, respectively, are not filled with a color.
%          {zero is the default value} (optional)
%
% horiz  - 'horizontal' for colorbar along the bottom, 'vertical' for
%          colorbar on the right side of the plot. Enter 'none' if no
%          colorbar is required. {vertical is default} (optional)
%
% colors - This sets the color scheme required:
%          'bluered'   - highest contrast blue fading through red at ref value
%          'greenred'  - green fading through red at ref value
%          'bluegreen' - blue fading through green at ref value
%          'jet'       - standard matlab 'jet' colorscheme.  Ref has no effect.
%          'cool'      - standard matlab 'cool' color scheme.  Ref has no effect.
%          'gray'      - straight gray color scale. Ref has no effect.
%          'parula'    - standard matlab 'parula' color scheme. Ref has no effect.
%          {bluered is the default} (optional)
%
% colvals- Cell array of text for each tick mark on the colorbar if required to 
%          be different from the supplied contour values. (optional)
%          See ridgepack_colorbar if further explanation is required.
%
% Ridgepack Version 1.2
% Andrew Roberts, LANL, 2019 (afroberts@lanl.gov) 

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% check that number of inputs is sufficient
if nargin<3
 error('insufficient information for plot')
elseif ~isstruct(nc)
 error(['nc is not a structure'])
elseif ~ischar(var)
 error(['var is not a character variable'])
elseif ~isstruct(ncvert)
 error(['ncvert is not a structure'])
elseif ~isfield(nc,var)
 error([var,' is not a field in the netcdf structure nc'])
elseif ~isfield(ncvert,'nCells')
 error(['nCells is not a field in the netcdf structure ncvert'])
elseif ~isfield(ncvert,'nEdgesOnCell')
 error(['nEdgesOnCell is not a field in the netcdf structure ncvert'])
elseif ~isfield(ncvert,'verticesOnCell')
 error(['verticesOnCell is not a field in the netcdf structure ncvert'])
elseif ~isfield(ncvert,'latitude')
 error(['latitude is not a field in the netcdf structure ncvert'])
elseif ~isfield(ncvert,'longitude')
 error(['longitude is not a field in the netcdf structure ncvert'])
end

% check for mask
if (nargin>=4 & isempty(mask)) | nargin<4
 mask=ncvert.nCells.data;
end

% check for contour interval
if nargin>=5 & isempty(cont) | nargin<5
 cont=[];
end

% check whether or not a linear contour interval is used
if nargin<6;
 loglin='linear' ;
elseif isempty(loglin)
 loglin='';
elseif ~(ischar(loglin) | isnumeric(loglin))
 error('loglin must be a character input')
end

% check ref value
if nargin<7;
 ref=0.0;
elseif ~isnumeric(ref) || length(ref)~=1
 error('ref must be a single number')
end

% check orientation for colorbar
if nargin<8;
 horiz='vertical';
elseif ~ischar(horiz)
 error('horiz must be a character input')
end

% color shading
if nargin<9;
 colors='bluered' ;
elseif ~ischar(colors)
 error('colors must be a character input')
end

% get current axes
hmap=get(gcf,'CurrentAxes');
if ~ismap(hmap)
 error('Current axes must be a map')
end
figout=get(hmap,'OuterPosition');
fontsize=min(11,max(9,10*sqrt((figout(3).^2)+(figout(4).^2))/sqrt(2)));
set(hmap,'fontsize',fontsize);

% get units
if isfield(nc.(var),'units')
 units=ridgepack_units(nc,var);
else
 units=[];
end

% get contour interval if it does not exist
if isempty(cont) & strcmp(loglin,'linear')

 contda=nc.(var).data(nc.(var).data~=0);
 varstd=std(contda);
 varske=skewness(contda);
 varmed=median(contda);
 if varstd==0
  error('Data has zero variance')
 end
 mincont=min(ref,max(varmed-max(1,varske/varstd)*3*varstd,min(contda)));
 maxcont=max(ref,min(varmed+max(1,varske/varstd)*3*varstd,max(contda)));
 contint=(maxcont-mincont)/12;
 contint=str2num(num2str(contint,'%2.2g'));
 cont=[mincont:contint:maxcont];

elseif ~isempty(cont) & isnumeric(loglin)

 if loglin==0; 
   error('loglin must be greater than zero for logarithmic scaling'); 
 end

 if min(diff(cont))<1

  if any(cont>0)
   lcont=unique(round(cont(cont>=0)));
   k=0;
   for i=[min(lcont)-1 lcont]
   for j=1:9;
    k=k+1;
    pcont(k)=loglin*j*(10^(abs(i)-1));
   end
   end
   mincont=loglin*10^(abs(min(cont(cont>=0)))-1);
   if any(cont<0); mincont=floor(mincont); end
   maxcont=loglin*10^(abs(max(cont(cont>=0)))-1);
   pcont=[mincont pcont(pcont>mincont & pcont<maxcont) maxcont];
  else
   pcont=[];
  end

  if any(cont<0)
   lcont=unique(round(cont(cont<=0)));
   k=0;
   for i=[lcont min(abs(lcont))-1]
   for j=9:-1:1;
    k=k+1;
    ncont(k)=-loglin*j*(10^(abs(i)-1));
   end
   end
   mincont=-loglin*10^(abs(min(cont(cont<=0)))-1);
   maxcont=-loglin*10^(abs(max(cont(cont<=0)))-1);
   if any(cont>0); maxcont=ceil(maxcont); end
   ncont=[mincont ncont(ncont>mincont & ncont<maxcont) maxcont];
  else
   ncont=[];
  end

  cont=[ncont pcont];

 else

  cont=unique(round(cont));
  cont=sign(cont).*loglin.*10.^(abs(cont)-1);

 end

 cont=cont(cont~=0); % remove zeros from log scale
 z=nc.(var).data;
 z(abs(z)<min(abs(cont))+eps)=sign(z(abs(z)<min(abs(cont))+eps))*min(abs(cont))+eps;
 nc.(var).data=z;
 clear z

elseif isempty(cont)

 error('cont is empty')

end

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

% colorbar
if strcmp(horiz,'none')
 disp('No colorbar being plotted')
elseif nargin==11
 ridgepack_colorbar(cont,units,loglin,horiz,ref,colvals);
else
 ridgepack_colorbar(cont,units,loglin,horiz,ref);
end

% add title, removing previous title
ridgepack_title(nc,['Shading: ',char(nc.(var).long_name)],1);

% clear output if none requested
if nargout==0; clear nc; end

if debug; disp(['...Leaving ',mfilename]); end

