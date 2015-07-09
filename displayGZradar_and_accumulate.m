% Displaying NIMROD radar rainfall data 

clear all; close all;

TOTRAIN = zeros(2175,1725);
x1 = datenum(2007,7,20,9,0,0);


for i=1: 288,

    i
    datex1 = datevec(x1);
%nimrodfilename = sprintf('Y:\\miguel\\raddis\\displaynimrod\\metoffice-c-band-rain-radar_uk_%.4d%.2d%.2d%.2d%.2d_1km-composite.dat.gz', datex1(1), datex1(2), datex1(3), datex1(4), datex1(5) )

nimrodfilename = 'Y:\miguel\raddis\displaynimrod\metoffice-c-band-rain-radar_uk_200707200900_1km-composite.dat.gz';    % input nimrod file
mapfilename = 'Y:\miguel\raddis\displaynimrod\UKNG_m.dat';  % coastline map. Modify the directory to the correct path
%tmpfolder = '/home/public/miguel/TMP/';  % temporary directory to decompress gz files (CREATE THIS FOLDER IF NECESSARY!)
tmpfolder = 'C:\Temp\';   %   /home/public/will/tmp/

% check if file is compressed
tt=length(nimrodfilename);
if ( strcmp('.gz', nimrodfilename(tt-2:tt)) == 1)   % if the input file is compressed (i.e. with extension '.gz')
     newname = gunzip(nimrodfilename, tmpfolder); % decompress file in a temporary folder and output new file name
     nimrodfilename = newname{1,1};    % copy new name into the nimrodfilename
end

% loading NIMROD file
[MTX DATE PARAM] = readnimrod(nimrodfilename);  % reading input nimrod file (decompressed). Note that data matrix is save in MTX
if ( prod(size(MTX))>0 ) 
    TOTRAIN = TOTRAIN + MTX/12;% Accumulation 
end

x1 = x1 + 1/(24*12);    % increasing the time by 5min

end

% SAVE THE ACCUMULATION MATRIX  HERE






% 
% % loading UK's boundary 
% UKCoastLine = load(mapfilename);    % in meters
% UKCoastLine = UKCoastLine / 1000;   % in kms
% 
% % box of the data matrix
% nrows = PARAM(7);  ncols = PARAM(6);
% resolution = PARAM(4);
% easting0 = PARAM(3);
% northing0 = PARAM(2)-ncols*resolution;
% x = (0: 1: nrows-1)';  x = (x*resolution + easting0)/1000;  % in km
% y = (0: 1: ncols-1)';  y = (y*resolution + northing0)/1000;  % in km
% 
% figure(1); 
% pcolor( x,y, MTX); shading flat;
% hold on; plot( UKCoastLine(:,1), UKCoastLine(:,2), 'k' ); hold off;
% caxis([0 10]);  colorbar('vert');  
% title( sprintf('RADAR RAINFALL %.2d/%.2d/%.4d %.2d:%.2d:%.2d  ', DATE(3), DATE(2), DATE(1), DATE(4), DATE(5), DATE(6)) );
% xlabel('Easting (km)'); ylabel('Northing (km)');
% axis([-200 800 -300 1200]);   % this scale is to see only UK
% 
