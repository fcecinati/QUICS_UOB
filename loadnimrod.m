clear all, close all;

filename = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\data_north_england\01\metoffice-c-band-rain-radar_uk_200701021650_1km-composite.dat.gz';
tmpfolder = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Temp\';

% check if file is compressed
tt=length(filename);
if ( strcmp('.gz', filename(tt-2:tt)) == 1)   % if the input file is compressed (i.e. with extension '.gz')
     newname = gunzip(filename, tmpfolder); % decompress file in a temporary folder and output new file name
     filename = newname{1,1};    % copy new name into the nimrodfilename
end

[mtx date parameters] = readnimrod(filename);   % function to read nimrod Met Office rainfall data



