%This script will implement the Germann's method for ensemble generation

clear all
close all
clc
rng('shuffle')
%% Parameters to define:

data_path = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\data_north_england';

rain_gauge_file = 'gauge_1h_2007-2011.dat';

corresponding_radar_file = 'radar_1h_2007-2010.dat';

coordinates_file = 'gaugecoordinates.dat';

scripts_and_functions_path = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\QUICS_UOB';

year_to_consider = 2007;

known_corrupted_data_stations = [8 31 78 168 178 194 212];

number_ensembles = 2;

%% Read data
% Read rain gauges and radar data in .dat format

cd(data_path)
G_file = load(rain_gauge_file);
R_file = load(corresponding_radar_file);
coord = (load(coordinates_file))';
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
coord(:,corrupted) = [];

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

% % Covariance matrix of the errors
% % Need to check the weighting method
% diff_e_m = zeros(t,x);
% for i=1:t 
%     diff_e_m(i,:) = err(i,:)-mean_e;
% end
% cov_e = zeros(x,x);
% temp1 = zeros(t,1);
% temp2 = zeros(t,1);
% for k = 1:x
%     for l = 1:x
%         for tm = 1:t
%             temp1(tm) = R(tm,k)*R(tm,l)*diff_e_m(tm,k)*diff_e_m(tm,l);
%             temp2(tm) = R(tm,k)*R(tm,l);
%         end
%         cov_e(k,l) = nansum(temp1)/nansum(temp2);
%     end
% end

% Covariance matrix of the errors without weights
diff_e_m = zeros(t,x);
for i=1:t 
    diff_e_m(i,:) = err(i,:)-mean_e;
end
cov_e = zeros(x,x);
temp1 = zeros(t,1);

for k = 1:x
    for l = 1:x
        for tm = 1:t
            temp1(tm) = diff_e_m(tm,k)*diff_e_m(tm,l);
        end
        cov_e(k,l) = nansum(temp1)/t;
    end
end

%% spatial decorrelation function and plotting
correl = zeros(x,x);
distance = zeros(x,x);

for x1=1:x
    for x2=1:x
        correl(x1,x2) = cov_e(x1,x2)/((cov_e(x1,x1)*cov_e(x2,x2))^0.5);
        distance(x1,x2) = (((coord(1,x1)-coord(1,x2))^2)+((coord(2,x1)-coord(2,x2))^2))^0.5;
    end
end

c = reshape(correl,x*x,1);
d = reshape(distance,x*x,1);
d(isnan(c)) = [];
c(isnan(c)) = [];
dmax = max(max(d));

% [fitobjectexp, gofexp] = fit(d,c,'exp1');
% coeffe =  coeffvalues(fitobjectexp);
% rsq = gofexp.rsquare;
% plx =(0:dmax/100:dmax);
% ply1 = coeffe(1)*exp(coeffe(2)*plx);
% str1 = [num2str(coeffe(1)) '*exp(' num2str(coeffe(2)) '*x)    R square = ' num2str(rsq)];

[fitobjectinv, gofinv] = fit(d,c,'rat11');
coeffi =  coeffvalues(fitobjectinv);
rsq = gofinv.rsquare;
plx =(0:dmax/100:dmax);
ply2=((coeffi(1)*plx) + coeffi(2))./(plx+coeffi(3));
str2 = ['(' num2str(coeffi(1)) '*x+' num2str(coeffi(2)) ')/(x+' num2str(coeffi(3)) ')   R square = ' num2str(rsq)];

error_corr = figure;
figure(error_corr)
scatter(d,c)
title('spatial correlation of the errors');
xlabel('meters')
hold on
% plot(plx,ply1,'r')
plot(plx,ply2,'g')
legend('error correlation coefficients',str2)
saveas(error_corr,'error_spatial_correlation','fig')
hold off
close(error_corr)

%% temporal decorrelation function
rho_t = zeros(x,2);

for j=1:x
    v1 = err(1:t-1,j);
    v2 = err(2:t,j);
    m1 = nanmean(v1);
    m2 = nanmean(v2);
    sd1 = nanstd(v1);
    sd2 = nanstd(v2);
    rho_t(j,1)=(nanmean((v1-m1).*(v2-m2)))/(sd1*sd2);
end

for j=1:x
    v1 = err(1:t-2,j);
    v2 = err(3:t,j);
    m1 = nanmean(v1);
    m2 = nanmean(v2);
    sd1 = nanstd(v1);
    sd2 = nanstd(v2);
    rho_t(j,2)=(nanmean((v1-m1).*(v2-m2)))/(sd1*sd2);
end

rho = nanmean(rho_t);
x_axis = [1 2];

rho_fig = figure;
figure(rho_fig)
plot(x_axis,rho)
title('temporal correlation of the errors');
xlabel('hours')
saveas(rho_fig,'error_temporal_correlation','fig')
close(rho_fig)

%% SVD Decomposition

[U,D,V] = svd(cov_e);
L = U*(D^0.5);


% % Test to check that L*L_transp equals the covariance matrix
% cov_test = L * L';
% c1 = reshape(cov_e, x*x, 1);
% c2 = reshape(cov_test,x*x,1);
% cov_t = figure;
% figure(cov_t)
% scatter(c1,c2)
% title('test to check that L*L^T equals the covariance matrix');
% hold on
% plot((min(c1):max(c1)),(min(c1):max(c1)),'r')
% legend('original covariance against L*L^T','covariance = L*L^T','Location','southeast')
% saveas(cov_t,'covariance test','fig')
% close(cov_t)


% %% Autoregressive parameters
% 
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
% r1m = nanmean(r1);
% 
% r2=zeros(x,1);
% temp1 = zeros(t-2,1);
% temp2 = zeros(t-2,1);
% for k = 1:x
%     for tm = 1:t-2
%         temp1(tm) = R(tm,k)*R(tm+2,k)*diff_e_m(tm,k)*diff_e_m(tm+2,k);
%         temp2(tm) = R(tm,k)*R(tm+2,k);
%     end
%     r2(k) = nansum(temp1)/nansum(temp2)*cov_e(k,k);
% end
% r2m = nanmean(r2);
% 
% Instead, taking the r1 and r2 from temporal decorrelation data
r1m = rho(1);
r2m = rho(2);

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
        delta(:,tm,ne) = delta1(:,tm,ne).*v + mean_e';
%         delta(:,tm,ne) = L*normrnd(0,1,x,1)+ mean_e';
    end 
end

% decorrelation test on generated delta
for f=1:n
    d1 = delta(:,:,f)';
    mean_d = nansum(R.*d1)./nansum(R);
    diff = zeros(t,x);
    for i=1:t 
        diff(i,:) = d1(i,:)- mean_d;
    end
    cov_d = zeros(x,x);
    temp1 = zeros(t,1);

    for k = 1:x
        for l = 1:x
            for tm = 1:t
                temp1(tm) = diff(tm,k)*diff(tm,l);
            end
            cov_d(k,l) = nansum(temp1)/t;
        end
    end

    correl_d = zeros(x,x);

    for x1=1:x
        for x2=1:x
            correl_d(x1,x2) = cov_d(x1,x2)/((cov_d(x1,x1)*cov_d(x2,x2))^0.5);
        end
    end

    c = reshape(correl_d,x*x,1);
    d = reshape(distance,x*x,1);
    d(isnan(c)) = [];
    c(isnan(c)) = [];
    dmax = max(max(d));

    [fitobjectinv_d, gofinv_d] = fit(d,c,'rat11');
    coeffi =  coeffvalues(fitobjectinv_d);
    rsq = gofinv_d.rsquare;
    plx =(0:dmax/100:dmax);
    ply3=((coeffi(1)*plx) + coeffi(2))./(plx+coeffi(3));
    str3 = ['(' num2str(coeffi(1)) '*x+' num2str(coeffi(2)) ')/(x+' num2str(coeffi(3)) ')   R square = ' num2str(rsq)];

    error_corr_d = figure;
    figure(error_corr_d)
    scatter(d,c)
    title(['spatial correlation of the generated errors' num2str(f)]);
    xlabel('meters')
    hold on
    % plot(plx,ply1,'r')
    plot(plx,ply3,'g')
    plot(plx,ply2,'r')
    legend('generated error correlation coefficients','fitting curve','fitting curve of the original error spatial correlation')
    name = ['error_correlation_d' num2str(f)];
    saveas(error_corr_d,name,'fig')
    hold off
    close(error_corr)
end



 