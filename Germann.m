%This script will implement the Germann's method for ensemble generation

clear all
close all
clc
rng('shuffle')
%% Parameters to define:

data_path = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\data_north_england';
rain_gauge_file = 'gauge_1h_2007-2011.dat';
corresponding_radar_file = 'radar_1h_2007-2010.dat';
scripts_and_functions_path = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\QUICS_UOB';

year_to_consider = 2007;

known_corrupted_data_stations = [8 31 78 168 178 194 212];

number_ensembles = 10;

%% Read data
% Read rain gauges and radar data in .dat format

cd(data_path)
G_file = load(rain_gauge_file);
R_file = load(corresponding_radar_file);
cd(scripts_and_functions_path)

years = int16(R_file(:,1));
months = int16(R_file(:,2));
days = int16(R_file(:,3));
hours = int16(R_file(:,4));

% Data is limited to 2007:
Rsize = size(R_file);
G = G_file((years==year_to_consider),7:Rsize(2));
R = R_file((years==year_to_consider),7:Rsize(2));

%use of NaN for negative values and to ignore days of no rain
G(G<0) = NaN;
R(R<0) = NaN;
G(G==0) = NaN;
R(R==0) = NaN;

% Delete the stations with No Data or corrupted data (change the list if
% needed)
corrupted = known_corrupted_data_stations;
G(:,corrupted) = [];
R(:,corrupted) = [];


% final size of the data matrix
sz = size(R);
t=sz(1);
x=sz(2);
    

%% Calculation of error statistics

% !!! A test to check if the errors can be assumed as gaussian must be
% introduced!!!

% Error matrix
err = 10*(log(G./R));

% Expected values of the errors
mean_e = nansum(R.*err)./nansum(R);

% Covariance matrix of the errors

diff_e_m = zeros(t,x);
for i=1:t 
    diff_e_m(i,:) = err(i,:)-mean_e;
end
cov_e = zeros(x,x);
temp1 = zeros(t,1);
temp2 = zeros(t,1);
for k = 1:x
    for l = 1:x
        for tm = 1:t
            temp1(tm) = R(tm,k)*R(tm,l)*diff_e_m(tm,k)*diff_e_m(tm,l);
            temp2(tm) = R(tm,k)*R(tm,l);
        end
        cov_e(k,l) = nansum(temp1)/nansum(temp2);
    end
end

cov_e = round(cov_e*1000)/1000;

%% SVD Decomposition

[U,D,V] = svd(cov_e);
L = U*(D^0.5);


%% Autoregressive parameters

% Autoregressive parameters r1 and r2
r1=zeros(x,1);
temp1 = zeros(t-1,1);
temp2 = zeros(t-1,1);
for k = 1:x
    for tm = 1:t-1
        temp1(tm) = R(tm,k)*R(tm+1,k)*diff_e_m(tm,k)*diff_e_m(tm+1,k);
        temp2(tm) = R(tm,k)*R(tm+1,k);
    end
    r1(k) = nansum(temp1)/nansum(temp2)*cov_e(k,k);
end
r1m = nanmean(r1);

r2=zeros(x,1);
temp1 = zeros(t-2,1);
temp2 = zeros(t-2,1);
for k = 1:x
    for tm = 1:t-2
        temp1(tm) = R(tm,k)*R(tm+2,k)*diff_e_m(tm,k)*diff_e_m(tm+2,k);
        temp2(tm) = R(tm,k)*R(tm+2,k);
    end
    r2(k) = nansum(temp1)/nansum(temp2)*cov_e(k,k);
end
r2m = nanmean(r2);

% AR(2) parameters
a1 = r1m*(1-r2m)/(1-(r1m^2));
a2 = (r2m-(r1m^2))/(1-(r1m^2));

% Scaling factor
v=((1+a2)/((1-a2)*(1-a1+a2)*(1+a1+a2)))^(-0.5);

%% Generation of the perturbed field

n = number_ensembles;
delta1 = zeros(x,t,n);
delta = zeros(x,t,n);

for ne = 1:n
    delta1(:,1,ne)= L*normrnd(0,1,x,1);
    delta1(:,2,ne)= L*normrnd(0,1,x,1);
    for tm = 3:t
        delta1(:,tm,ne) = L*normrnd(0,1,x,1) + a1*delta1(:,tm-1,ne) + a2*delta1(:,tm-2,ne);
        delta(:,tm,ne) = delta1(:,tm,ne)*v + mean_e;
    end 
end





 