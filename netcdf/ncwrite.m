function ncid=ncwrite(nc,ncfile,rec,recp,recl)

% NCWRITE - Writes to a netcdf file from an nc structure
%
% function ncwrite(nc,ncfile,rec,recp)
%
% This function writes to a netcdf file.  It will clobber an existing file
% of the same name as ncfile. The function can also add a hyperslab to 
% an existing ncfile if rec and recp are specified. All netcdf attributes
% are written following the CF-1.4 convention.  
%
% Mandatory inputs:
% nc      - nc structure (defined below, but also see ncstruct).
% ncfile  - NetCDF filename to which data is to be written.
%
% Optional inputs:
% rec     - Cell array of dimensions for which grid point values are provided
%           in recp if adding a hyperslab to a netcdf file.
%           Here is an example of what rec should look like, if one were 
%           to write a slice of data for a  variable with dimensions 
%           {'time','x','y','z'} on the y and z plane, one would need to 
%           specifiy the values of time and x for which the hyperslab is defined, 
%           hence rec={'time','x'}. The exact indices of time and x that 
%           define the hyperslab must then be provided in recp. 
%           For example recp={3,10} would specify for the 3rd and 10th 
%           indices for time and x, respectively. The unlimited
%           dimension in netcdf files written is reserved only for 
%           a variable named 'time'.  No other variable can take the
%           unlimited dimension using this function.  If time is not
%           a dimension, then an unlimited dimension is not assigned.
%           If you wish to append a time-constant hyperslab to a netcdf
%           file, the you simply specify zero as the recp value for time:
%           ncwrite(nc,ncfile,{'time'},{0}).
%
% recp    - Cell array of grid points corresponding to the dimensions in rec
%           for which a hyperslab is to be slotted into the netcdf archive.
%
% recl    - Full arrays of non-time dimensions to be written when the netcdf 
%           file is first established. This is required if one intends to append
%           hyperslabs to a file that are not all for a single time slice.
%           If, for example, a variable has dimensions {'x','y','time'}, and one
%           wants to establish a netcdf file, and then write future hyperslabs
%           to that file with for hyperslabs of {'x','y'} (i.e. timeseries
%           for each (x,y) grid point), then would establish the full x and y
%           diemensions in the first write:
%
%           ncwrite(ncnew,filename,{'x','y'},{i,j},{xdata,ydata})
%
%           where ncnew only includes data for (x,y)=(i,j), but xdata and ydata
%           includes all values for those dimensions. This establishes the
%           netcdf file. Subsequent writes can be placed inside a loop too append
%           hyperslabs to the file:
%
%           for i=1:length(xdata);
%           for j=1:length(ydata)
%
%            ncwrite(ncnew,filename,{'x','y'},{i,j});
%
%           end
%           end
%
%           So in the above example, recl={xdata,ydata}
%
% Global attributes for the netcdf file appear in the nc data structure:
%
% nc.attributes.title - A succinct description of what is in the dataset.
%                       This is compulsory.
%
% nc.attributes.source - The method of production of the original data. 
%           If it is model-generated, source should name the model and its,
%           version as specifically as is useful. If it is an observational, 
%           source should characterize this (e.g., "surface observation" or 
%           "radiosonde"). This attribute should also include specifics of
%           of the instrumentation used to produce the data:
%           1) satellite/remotely sensed data: include summary 
%              of the satellite, instrument and processing algorithms.
%           2) GPS/positional data: include processing procedure and 
%              retrieval method.
%           3) directly sensed data: exact details on sensors used (model
%              numbers, sensor types)
%           4) model generated data: computer on which the data was 
%              generated.
%           5) conglomerate dataset (such as climatology or assimilated):
%              list (within reasonable bounds) all sources of data and 
%              methods/models used to compile the dataset.
%           This is compulsory.
%
% nc.attributes.references - Published or web-based references that describe 
%           the data or methods used to produce it. This need to be included 
%           if references are not applicable.  Exclusion of this will
%           will remove the entry from the netcdf file. This is optional.
%
% nc.attributes.comment - Miscellaneous information about the data or
%           methods use to produce it.  This is optional.
%
%  [Note: If attribute text is long, it is advisable to place newlines at
%  appropriate placed within the text (in netcdf files these are "\n").
%
%  This function also requires two environment variables to be supplied
%  in the UNIX/Linux shell for the global attributes when it is run: 
%
%  USER_ORG - the institution or organization in which the file is created.
%  USER_REAL - the real name of the user.  
%
%  These are supplied by environment variable rather than direct input
%  to improve portability of the function, and to remove replication input.
%  They may be set in bash using: export USER_REAL='Joe Blogs', or in csh 
%  using: setenv USER_REAL 'Joe Blogs'. These commands can be appended to 
%  .bashrc or .cshrc, files respectively.]
%
% Each variable component of the netcdf structure contains the following 
% for each variable called 'name', or other CF comliant fields (those listed
% here are most common):
% 
% nc.name.long_name - text string of the long version of the variable name
%
% nc.name.units     - text string specifying the data units.  If name='time'
%                     then the text string represents the conversion 
%                     to the CF-1.1 coordinate convention from matlab serial
%                     time when the time dimension data is written to netcdf.
%                     If no time units are specified, the default is
%                     hours since 1-1-1900 00:00:00.0 UTC.
%
% nc.name.dimension - cell arrays of texts specifying the dimensions of the 
%                     the variables. So the variables time would have
%                     nc.name.dimension={'time'}, as would a time series of, 
%                     for example, barometric pressure. Multidimensional
%                     data would have multiple dimensions, such as MSL from
%                     the ERA40 dataset: {'time','latitude','longitude'}
%                     where time, latitude and longitude are also defined
%                     variables in the netcdf structure.
%                          
% nc.name.data      - data values with the dimension allocated in line with 
%                     the defined dimensions.
%
% nc.name.fillvalue - the fillvalue of the data to be assigned during 
%                     write of the data to netcdf.  The fillvalue
%                     is assumed to be NaN in matlab, and so all NaN's are
%                     substituted for nc.name.fillvalue when writing to netcdf.
%
% nc.name.type      - type the data is to be assigned to when writing to netcdf.
%                     The data is stored in the structure nominally as double
%                     precision, but this can be reduced to single upon writing
%                     to netcdf.
%
% For further information on the structure, see documentation in ncstruct.
%
% Output:
% A netcdf file is added to or created in the directory in which you are working.
%
% ncid - Optional output giving the netcdf id of the file written to.
%
% Written by Andrew Roberts 2012
% Department of Oceanography, Naval Postgraduate School
%
% $Id$

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% check for default mode option
if nargin<2 
 error('nc and ncfile must be specified'); 
elseif not(isstruct(nc)) & ischar(nc);
 error([nc,' is not an nc structure']);
elseif not(isstruct(nc));
 error([char(inputname(1)),' is not an nc structure']);
elseif ~ischar(ncfile)
 error('File name nncorrect');
elseif nnz(isstrprop(ncfile, 'alphanum'))==0
 error('No name given for the netcdf file');
elseif nargin>=4
 if ~iscell(recp)
  error('recp should be a cell array');
 elseif isempty(recp)
  error('recp entered as an empty array - using slot in the wrong context');
 elseif ~iscell(rec)
  error('rec should be a cell array')
 elseif isempty(rec)
  error('rec entered as an empty array - using slot in the wrong context')
 elseif length(recp)~=length(rec)
  error('recp and rec are different lengths');
 end
 if nargin==5 & ~iscell(recl) 
  error('recl is not a cell array')
 elseif nargin==5 & length(recl)~=length(recp)
  error('length of recl does not equal length of recp')
 end
 for i=1:length(recp)
  if ~isnumeric(recp{i})
   error(['Element ',num2str(i),' of recp must be numeric'])
  end
  if nargin==5 & ~isnumeric(recl{i})
   error(['Element ',num2str(i),' of recl must be numeric'])
  end
  if ~ischar(rec{i})
   error(['Element ',num2str(i),' of rec must be a character string'])
  end
 end
elseif nargin~=5 & nargin~=4 & nargin~=2
 error('Incorrect number of inputs')
end

% place time in last dimension
nc=ncshuffle(nc,{'time'});

% Append .nc to file name if needed 
ui=length(ncfile); li=max(1,ui-2);
if ui < 4 || not(strcmp(ncfile(li:ui),'.nc')) ; ncfile=[ncfile,'.nc']; end

% Create or open the netcdf file
if exist(ncfile,'file') > 0 & nargin==4;

 disp(['Writing hyperslab to ',ncfile,' (netcdf v',netcdf.inqLibVers,')']);

 ncid=netcdf.open(ncfile,'NC_WRITE');

 % write to variables netcdf 
 ncvariable(ncid,nc,rec,recp);

elseif nargin==4 | nargin==5;

 if exist(ncfile,'file') & nargin==5;
  disp(['Overwriting ',ncfile,' (netcdf v',netcdf.inqLibVers,')']);
 else
  disp(['Creating hyperslab in new file ',ncfile,' (netcdf v',netcdf.inqLibVers,')']);
 end

 ncid=netcdf.create(ncfile,'NC_64BIT_OFFSET');
 netcdf.endDef(ncid);

 % write global attributes
 ncglobeatt(ncid,nc);

 % write to variables netcdf 
 if nargin==4
  ncvariable(ncid,nc,rec,recp);
 elseif nargin==5 
  ncvariable(ncid,nc,rec,recp,recl);
 end

elseif nargin==2;

 disp(['Creating or overwriting ',ncfile,' using netcdf V',netcdf.inqLibVers]);

 ncid=netcdf.create(ncfile,'NC_64BIT_OFFSET');
 netcdf.endDef(ncid);

 % write global attributes
 ncglobeatt(ncid,nc);

 % write to variables netcdf 
 ncvariable(ncid,nc);

end

% Cloase the file
disp(['Closing ',ncfile]);
netcdf.close(ncid);

% check final write
if debug; ncdump(ncfile); end

% clear output if none requested
if nargout==0; clear ncid; end

if debug; disp(['...Leaving ',mfilename]); end

