function ridgepack_fprint(form,name,fig,movie)

% ridgepack_fprint - Save a matlab figure to a file
%
% function ridgepack_fprint(form,name,fig,movie)
%
% This prints without using the zbuffer which is better for some rendering
% images with opacity built into them.
% 
% INPUT:
%
% form  - format of output (default is 'epsc'). Valid values:
%
%         ps       - PostScript for black and white printers
%         psc      - PostScript for color printers
%         ps2      - Level 2 PostScript for black and white printers
%         psc2     - Level 2 PostScript for color printers
%         eps      - Encapsulated PostScript
%         epsc     - Encapsulated Color PostScript
%         eps2     - Encapsulated Level 2 PostScript
%         epsc2    - Encapsulated Level 2 Color PostScript
%         jpeg<nn> - JPEG image, quality level of nn (figures only)
%         tiff     - TIFF with packbits (lossless run-length encoding)
%         tiffnocompression - TIFF without compression (figures only)
%         png      - Portable Network Graphic 24-bit truecolor image
%         pdf      - PDF
%         
% name  - name of the file to be output (default is 'figure1')
%
% fig   - matlab figure number (default is 1)
%
% movie - changes dots per inch to 300 for movies.  Set to 1 for
%         this option to work, otherwise set to 0 (default is 0).
%         To Set to 1200dpi, set movie to 2.
%
%
% OUTPUT:
%
% Graphics file called name.form
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)


global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% change to zbuffer
set(gcf,'renderer','zbuffer')
%set(gcf,'renderer','painters')

if nargin<2
        error('no file name')
elseif nargin<3 
	disp(['Printing 150dpi resolution ',name])
	form=[form,' -r150'];
        fig=1;
        movie=1;
elseif nargin<4 
	disp(['Printing 150dpi resolution ',name])
	form=[form,' -r150'];
        movie=1;
elseif movie==1
	disp(['Printing 300dpi resolution ',name])
	form=[form,' -r300'];
elseif movie==2
	disp(['Printing 600dpi resolution ',name])
	form=[form,' -r600'];
else
	disp(['Printing 150dpi resolution ',name])
	form=[form,' -r150'];
end

if nargin<2
	name='figure1';
end
if nargin<3
	fig=get(0,'CurrentFigure');
end

if ishandle(fig) 
 printcommand=['print -f',num2str(fig),' -d',form,' ',name];
else
 error([num2str(fig),' is not a valid figure'])
end

if debug;
 disp(printcommand)
end

eval(printcommand);

if debug; disp(['...Leaving ',mfilename]); end


