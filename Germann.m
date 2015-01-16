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
cd('\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\QUICS_UOB')

years = int16(R_file(:,1));
months = int16(R_file(:,2));
days = int16(R_file(:,3));
hours = int16(R_file(:,4));

% Data is limited to 2007:
Rsize = size(R_file);
G = G_file((years==2007),7:Rsize(2));
R = R_file((years==2007),7:Rsize(2));

%use of NaN for negative values and to ignore days of no rain
G(G<0) = NaN;
R(R<0) = NaN;
G(G==0) = NaN;
R(R==0) = NaN;

% Delete the stations with No Data or corrupted data (change the list if
% needed)
corrupted = [8 31 78 168 178 194 212];
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

%% Cholesky Decomposition

% checks if the matrix is positive definite
% 
% if min(eig(cov_e))>0
%     L = chol(cov_e,'lower');    %simple Cholesky decomposition
% else
%     [R, indef, E] = cholmod(cov_e);     %modified Cholesky decomposition
%     L = R';
% end


[L,p] = chol(cov_e,'lower');

%% Autoregressive parameters

% % Autoregressive parameters r1 and r2
% r1=zeros(x,1);
% temp1 = zeros(t-1,1);
% temp2 = zeros(t-1,1);
% for k = 1:x
%     for tm = 1:t-1
%         temp1(tm) = R(tm,k)*R(tm+1,k)*diff_e_m(tm,k)*diff_e_m(tm+1,k);
%         temp2(tm) = R(tm,k)*R(tm+1,k);
%     end
%     r1(k) = nansum(temp1)/nansum(temp2)*cov_e(k,k);
% end
% 
% r2=zeros(x,1);
% temp1 = zeros(t-2,1);
% temp2 = zeros(t-2,1);
% for k = 1:x
%     for tm = 1:t-2
%         temp1(tm) = R(tm,k)*R(tm+2,k)*diff_e_m(tm,k)*diff_e_m(tm+2,k);
%         temp2(tm) = R(tm,k)*R(tm+2,k);
%     end
%     r1(k) = nansum(temp1)/nansum(temp2)*cov_e(k,k);
% end


 