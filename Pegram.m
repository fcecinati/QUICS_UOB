%This script will implement the (Pegram et al, 2011) method for ensemble generation

clear all
close all
clc
rng('shuffle')

%% Parameters to define:

% Path to the folder where the radar data is 
data_path = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\data_north_england\09-2009';
% file to work on
filename = 'metoffice-c-band-rain-radar_uk_200809051500_1km-composite.dat.gz';
% How many ensembles do you want to generate
number_ens = 1;

%% Data loading and pre-processing

tt=length(filename);
if (strcmp('.gz', filename(tt-2:tt)) == 1)   % if the input file is compressed (i.e. with extension '.gz')
     newname = gunzip(filename, tmpfolder); % decompress file in a temporary folder and output new file name
     filename = newname{1,1};    % copy new name into the nimrodfilename
end

[radar, date, param] = readnimrod(filename);

radar(radar<0.5) = NaN;
Z = radar - mean(mean(radar));

% pack with zeros

%% Fourier transform


