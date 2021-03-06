% This script generates radar rainfall ensembles with a method very similar 
% to the REAL method, described by (Germann et al. 2009*).
% The main differences are:
% 1) the weighting system adopted by Germann et al. is not used.
% 2) a normal score transform is applied to the errors
% 3) the decomposition of the covariance matrix is done with a singular 
%   value decomposition instead of a Cholewsky algorithm, eliminating the 
%   requirement for a positive definite covariance matrix.

%  Read the attached README file for information on the required data and
%  files

% Code by Francesca Cecinati, University of Bristol, Civil Eng.Dept. 
% francesca.cecinati@bristol.ac.uk

% * Germann, U., Berenguer, M., Sempere-torres, D., & Zappa, M. (2009). 
% REAL � Ensemble radar precipitation estimation for hydrology in 
% mountainous region. Q. J. R. Meteorol. Soc, 135(February), 445�456.

%% 
clear all
close all
clc
rng('shuffle')

%% INFORMATION TO PROVIDE:

% Folders:

% Path to the folder where the rain gauges, the corresponding radar and the coordinate .dat files are:
data_path = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\data_north_england\';
% Name of the rain gauges .dat file with extension:
rain_gauge_file = 'gauge_1h_2007-2011.dat';
% Name of the corresponding .dat radar file:
corresponding_radar_file = 'radar_1h_2007-2010.dat';
% Name of the coordinate files of the observation points:
coordinates_file = 'gaugecoordinates.dat';
% Path to the folder where this and the other scripts/functions are:
scripts_and_functions_path = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\QUICS_UOB\';
% Results path
results_path = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_German\Test_9\';
% Path to the complete radar data
radar_folder = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\data_north_england\09-2008\';
% Temporary folder
tmpfolder = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Temp\'; 

% Parameters:

% Which of the stations have unusable data?
known_corrupted_data_stations = [8;30;31;73;78;81;98;100;107;111;141;145;152;212;213;228;229;168;174;178;194;20;40;128;208;211];
% Which year do you want to use for the covariance matrix calculation?
year_to_consider = 2007;
% How many ensembles do you want to generate
n_ensembles = 5;
% What is the starting date of the ensembles you need?
starting = '2008-09-05 15:00:00';
% How many hours is your simulation long?
time_steps = 20;

% Optional output (Y/N):

% Do you want to check the covariance decomposition? 
cov_dec = 'Y';
% Do you want to plot the spatial correlation of the errors?
plot_spa = 'Y';
% Do you want to plot the temporal correlation of the errors?
plot_tem = 'Y';
% Do you want to check the spatial correlation of the perturbed fields?
spa_dec = 'Y';
% Do you want to check the temporal correlation of the perturbed fields?
temp_dec = 'Y';


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

% Data is limited to one year for covariance calculation purposes:
Rsize = size(R_file);
G = G_file((years==year_to_consider),7:Rsize(2));
R = R_file((years==year_to_consider),7:Rsize(2));
hours = hours(years==year_to_consider);
days = days(years==year_to_consider);
months = months(years==year_to_consider);
years = years(years==year_to_consider);

clear G_file R_file 

% Use of NaN for negative values and to ignore days of no rain
G(G<0.5) = NaN;
R(R<0.5) = NaN;
G(G==0) = NaN;
R(R==0) = NaN;

% Delete the stations with No Data or corrupted data (from the list above)
corrupted = known_corrupted_data_stations;
G(:,corrupted) = [];
R(:,corrupted) = [];
coord(:,corrupted) = [];

        % % Use of NaN for negative values and to ignore days of no rain
        % G = G(:,good);
        % R = R(:,good);
        % coord = coord(:,good);

% Final size of the data matrix
sz = size(R);
t=sz(1);
x=sz(2);
    
%% Calculation of error statistics

% Error matrix
er_nongaussian = 10*(log(G./R));

% Application of normal score transform to obtain a gaussian error matrix
R(isnan(er_nongaussian)==1) = NaN;
er_nongaussian_vector = er_nongaussian(isnan(er_nongaussian)==0);
[er_gaussian, scores] = nscore(er_nongaussian_vector);
err = er_nongaussian;
err(isnan(err)==0)=er_gaussian;

% Mean without weights
mean_e = nanmean(err);

% Covariance matrix of the errors without weights
diff = zeros(t,x);
for i=1:t 
    diff(i,:) = err(i,:)-mean_e;
end
cov_e = zeros(x,x);
temp1 = zeros(t,1);

for k = 1:x
    for l = 1:x
        for tm = 1:t
            temp1(tm) = diff(tm,k)*diff(tm,l);
        end
        cov_e(k,l) = nansum(temp1)/t;
    end
end

clear diff temp1
%% Spatial decorrelation function and plotting

%calculation of correlation coefficients and distances
correl = zeros(x,x);
distance = zeros(x,x);

for x1=1:x
    for x2=1:x
        correl(x1,x2) = cov_e(x1,x2)/((cov_e(x1,x1)*cov_e(x2,x2))^0.5);
        distance(x1,x2) = (((coord(1,x1)-coord(1,x2))^2)+((coord(2,x1)-coord(2,x2))^2))^0.5;
    end
end

% Reshaping and NaN elimination for function estimation and plotting
c = reshape(correl,x*x,1);
d = reshape(distance,x*x,1);
d(isnan(c)) = [];
c(isnan(c)) = [];
d=d/100000;         %Important! the function is calculated on hundreds of kilometers as unit

% Fit on decorrelation function
myfunction = fittype('exp(-((d/R0)^F))','independent',{'d'},'coefficients',{'R0','F'});
myoption = fitoptions('Method','NonlinearLeastSquares','StartPoint',[0.5 0.9]);
myfit = fit(d,c,myfunction,myoption);

% Plotting
if plot_spa == 'Y'
    fig = figure;
    figure(fig)
    plot(myfit,d,c)
    title('spatial correlation of the errors');
    xlabel('meters x 10^5')
    name = [results_path, 'error_spatial_correlation'];
    saveas(fig,name,'fig')
    saveas(fig,name,'jpg')
    close(fig)
end

%% Temporal decorrelation function and plotting

% Correlation calculation at 1, 2, 3 hours lag
rho_t = zeros(x,4);
rho_t(:,1) = 1;
for lag=1:3
    for j=1:x
        v1 = err(1:(t-lag),j);
        v2 = err((1+lag):t,j);
        m1 = nanmean(v1);
        m2 = nanmean(v2);
        sd1 = nanstd(v1);
        sd2 = nanstd(v2);
        rho_t(j,1+lag)=(nanmean((v1-m1).*(v2-m2)))/(sd1*sd2);
    end
end
maxrho = max(rho_t);
minrho = min(rho_t);
rho = nanmean(rho_t);

x_axis = (0:3);
clear rho_t

% Plotting
if plot_tem == 'Y'
    fig = figure;
    figure(fig)
    plot(x_axis,rho)
    title('temporal correlation of the errors');
    xlabel('hours')
    hold on
    plot(x_axis, maxrho, 'g')
    plot(x_axis, minrho, 'g')
    name = [results_path, 'error_temporal_correlation'];
    saveas(fig,name,'fig')
    saveas(fig,name,'jpg')
    close(fig)
end

%% SVD Decomposition

[U,D,V] = svd(cov_e);
L = U*(D^0.5);
clear U D V

% Test to check that L*L_transp equals the covariance matrix
if cov_dec == 'Y'
    cov_test = L * L';
    c1 = reshape(cov_e, x*x, 1);
    c2 = reshape(cov_test,x*x,1);
    fig = figure;
    figure(fig)
    scatter(c1,c2)
    title('test to check that L*L^T equals the covariance matrix');
    hold on
    plot((min(c1):max(c1)),(min(c1):max(c1)),'r')
    name = [results_path, 'covariance_test'];
    saveas(fig,name,'fig')
    saveas(fig,name,'jpg')
    close(fig)
end

%% Autoregressive parameters calculated
      
% Instead, taking the r1 and r2 from temporal decorrelation data
r1m = rho(2);
r2m = rho(3);

% AR(2) parameters
a1 = r1m*(1-r2m)/(1-(r1m^2));
a2 = (r2m-(r1m^2))/(1-(r1m^2));

% Scaling factor
v=((1+a2)/((1-a2)*(1-a1+a2)*(1+a1+a2)))^(-0.5);

%% Generation of the perturbed field
n = n_ensembles;
delta1 = zeros(x,t+2,n);
delta2 = zeros(x,t+2,n);

for ens = 1:n
    delta1(:,1,ens)= L*normrnd(0,1,x,1);
    delta1(:,2,ens)= L*normrnd(0,1,x,1);
    for tm = 3:t+2
        delta1(:,tm,ens) = L*normrnd(0,1,x,1) + a1*delta1(:,tm-1,ens) + a2*delta1(:,tm-2,ens);
        delta2(:,tm,ens) = delta1(:,tm,ens).*v + mean_e';
    end 
end

delta2 = delta2(:,3:t+2,:);
delta = inscore(delta2, scores);

clear delta1 delta2

%% Spatial Decorrelation Test on Generated Delta
if spa_dec == 'Y'
    for ens=1:n
        d1 = delta(:,:,ens)';

        % Mean calculation
        mean_d = nanmean(d1);

        % Covariance Calculation
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
                cov_d(k,l) = nansum(temp1)/(t);
            end
        end

        % Spatial correlation calculation
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
        d = d/100000;

        % Plotting
        fig = figure;
        figure(fig)
        plot(myfit, d, c)
        title(['spatial correlation of the generated errors n.' num2str(ens)]);
        xlabel('meters x 10^5')
        name = [results_path, 'spatial_correlation_delta_' num2str(ens)];
        saveas(fig,name,'fig')
        saveas(fig,name,'jpg')
        close(fig)
    end
end

%% Temporal decorrelation test
if temp_dec == 'Y'
    for ens=1:n
        d1 = delta(:,:,ens)';
        rho_t_d = zeros(x,4);
        rho_t_d(:,1) = 1;
        for lag=1:3
            for j=1:x
                v1 = d1(1:(t-lag),j);
                v2 = d1((1+lag):(t),j);
                m1 = nanmean(v1);
                m2 = nanmean(v2);
                sd1 = nanstd(v1);
                sd2 = nanstd(v2);
                rho_t_d(j,1+lag)=(nanmean((v1-m1).*(v2-m2)))/(sd1*sd2);
            end
        end
        maxrhod = max(rho_t_d);
        minrhod = min(rho_t_d);
        rho_d = nanmean(rho_t_d);

        fig = figure;
        figure(fig)
        plot(x_axis(1:3),rho_d(1:3))
        hold on
        plot(x_axis(1:3),rho(1:3),'r')
        plot(x_axis(1:3), maxrho(1:3), 'g')
        plot(x_axis(1:3), minrho(1:3), 'g')
        plot(x_axis(1:3), maxrhod(1:3), 'c')
        plot(x_axis(1:3), minrhod(1:3), 'c')
        title('temporal correlation of the deltas');
        xlabel('hours')
        name = [results_path, 'temporal_correlation_delta_' num2str(ens)];
        saveas(fig,name,'fig')
        saveas(fig,name,'jpg')
        close(fig)
    end
end

%% Generate the noise
noise = Noise(delta, coord, time_steps, results_path);

%% Apply the noise to original radar data, to generate the ensembles
[ensembles, raddata] = applyensembles(starting, time_steps, n_ensembles, radar_folder, tmpfolder, results_path);

%% Save the results
save ([results_path,'delta.mat'],'delta','coord','noise','ensembles','raddata')


