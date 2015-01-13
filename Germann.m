%This script will implement the Germann's method for ensemble generation

clear all
close all
clc
rng('shuffle')

%% Read data
% Read rain gauges and radar data in .dat format

cd('\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\data_north_england')
G_file = load('gauge_1h_2007-2011.dat');
R_file = load('radar_1h_2007-2010.dat');

years = int16(R_file(:,1));
months = int16(R_file(:,2));
days = int16(R_file(:,3));
hours = int16(R_file(:,4));

% rain gauge data is limited to 31-12-2010-23:00 as radar data is:
Rsize = size(R_file);
G = G_file(1:Rsize(1),7:Rsize(2));
R = R_file(1:Rsize(1),7:Rsize(2));

%use of NaN for negative values and to ignore days of no rain
G(G<0) = NaN;
R(R<0) = NaN;
G(G==0) = NaN;
R(R==0) = NaN;

%% Calculation of error statistics

err = 10*(log(G./R));

mean_err = nansum(R.*err)./nansum(R);

% cov_err = nansum(R.*err)./nansum(R);
 
 