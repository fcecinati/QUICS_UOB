%This script will implement the Germann's method for ensemble generation

clear all
close all
clc
rng('shuffle')
%% Parameters to define:

% Path to the folder where the rain gauges, the corresponding radar and the coordinate .dat files are:
data_path = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\data_north_england';
% Name of the rain gauges .dat file with extension:
rain_gauge_file = 'gauge_1h_2007-2011.dat';
% Name of the corresponding .dat radar file:
corresponding_radar_file = 'radar_1h_2007-2010.dat';
% Name of the coordinate files of the observation points:
coordinates_file = 'gaugecoordinates.dat';
% Path to the folder where this and the other scripts/functions are:
scripts_and_functions_path = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\QUICS_UOB\';
% Results path
results_path = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_German\Test_7\';
% 
% Which of the stations have unusable data?
% known_corrupted_data_stations = [8;30;31;73;78;81;98;100;107;111;141;145;152;212;213;228;229; ...  % invalid data  
%                                 168;174;178;194;% no data available for these gauges in 2007
%                                 20;40;128;208;211];%other

% Good gauges:
good = [1;2;3;4;5;6;7;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29;183;206;92;126;116;181;207;144;164;103;209;137;94;186;200;62];

% Which year do you want to use for the covariance matrix calculation?
year_to_consider = 2007;
% How many ensembles do you want to generate
number_ens = 1;
% How many hours is your simulation long?
sim_hours = 10;

% Do you want to check the covariance decomposition? (Y/N)
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
G = G(:,good);
R = R(:,good);
coord = coord(:,good);

G(G<0.5) = NaN;
R(R<0.5) = NaN;
G(G==0) = NaN;
R(R==0) = NaN;

% Delete the stations with No Data or corrupted data (from the list above)
% corrupted = known_corrupted_data_stations;
% G(:,corrupted) = [];
% R(:,corrupted) = [];
% coord(:,corrupted) = [];

% Final size of the data matrix
sz = size(R);
t=sz(1);
x=sz(2);
    
%% Calculation of error statistics

% !!! A test to check if the errors can be assumed as gaussian must be
% introduced!!!

% Error matrix
er_nongaussian = 10*(log(G./R));

% application of normal score transform to obtain a gaussian error matrix
R(isnan(er_nongaussian)==1) = NaN;
er_nongaussian_vector = er_nongaussian(isnan(er_nongaussian)==0);
[er_gaussian, scores] = nscore(er_nongaussian_vector);
err = er_nongaussian;
err(isnan(err)==0)=er_gaussian;

% Mean without weights
mean_e = nansum(R.*err)./nansum(R);

% Covariance matrix of the errors without weights
diff = zeros(t,x);
for i=1:t 
    diff(i,:) = err(i,:)-mean_e;
end
cov_e = zeros(x,x);
temp1 = zeros(t,1);
temp2 = zeros(t,1);
for k = 1:x
    for l = 1:x
        for tm = 1:t
            temp1(tm) = R(tm,k)*R(tm,l)*diff(tm,k)*diff(tm,l);
            temp2(tm) = R(tm,k)*R(tm,l);
        end
        cov_e(k,l) = nansum(temp1)/nansum(temp2);
    end
end

clear temp1 temp2
%% spatial decorrelation function and plotting

%calculation of correlation coefficients and distances
correl = zeros(x,x);
distance = zeros(x,x);

for x1=1:x
    for x2=1:x
        correl(x1,x2) = cov_e(x1,x2)/((cov_e(x1,x1)*cov_e(x2,x2))^0.5);
        distance(x1,x2) = (((coord(1,x1)-coord(1,x2))^2)+((coord(2,x1)-coord(2,x2))^2))^0.5;
    end
end

%reshaping and NaN elimination for function estimation and plotting
c = reshape(correl,x*x,1);
d = reshape(distance,x*x,1);
d(isnan(c)) = [];
c(isnan(c)) = [];
d=d/100000;         %Important! the function is calculated on hundreds of kilometers as unit

% Fit on decorrelation function
myfunction = fittype('exp(-((d/R0)^F))','independent',{'d'},'coefficients',{'R0','F'});
myoption = fitoptions('Method','NonlinearLeastSquares','StartPoint',[0.5 0.9]);
myfit = fit(d,c,myfunction,myoption);

%Plotting
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

%% temporal decorrelation function and plotting

% correlation calculation at 1, 2, 3 hours lag
% Autoregressive parameters r1 and r2
r1=zeros(x,1);
temp1 = zeros(t-1,1);
temp2 = zeros(t-2,1);
for k = 1:x
    for tm = 1:t-1
        temp1(tm) = R(tm,k)*R(tm+1,k)*diff(tm,k)*diff(tm+1,k);
        temp2(tm) = R(tm,k)*R(tm+1,k);
    end
    r1(k) = nansum(temp1)/(nansum(temp2)*cov_e(k,k));
end
r1m = nanmean(r1);
r1max = max(r1);
r1min = min(r1);

r2=zeros(x,1);
temp1 = zeros(t-2,1);
temp2 = zeros(t-2,1);
for k = 1:x
    for tm = 1:t-2
        temp1(tm) = R(tm,k)*R(tm+2,k)*diff(tm,k)*diff(tm+2,k);
        temp2(tm) = R(tm,k)*R(tm+2,k);
    end
    r2(k) = nansum(temp1)/(nansum(temp2)*cov_e(k,k));
end
r2m = nanmean(r2);
r2max = max(r2);
r2min = min(r2);

r_mean = [1,r1m,r2m];
r_max = [1,r1max,r2max];
r_min = [1,r1min,r2min];
x_axis = [0, 1, 2];

%Plotting
if plot_tem == 'Y'
    fig = figure;
    figure(fig)
    plot(x_axis,r_mean)
    title('temporal correlation of the errors');
    xlabel('hours')
    hold on
    plot(x_axis, r_max, 'g')
    plot(x_axis, r_min, 'g')
    name = [results_path, 'error_temporal_correlation'];
    saveas(fig,name,'fig')
    saveas(fig,name,'jpg')
    close(fig)
end

%% SVD Decomposition

[U,D,V] = svd(cov_e);
L = U*(D^0.5);
% clear U D V

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
      
% AR(2) parameters
a1 = r1m*(1-r2m)/(1-(r1m^2));
a2 = (r2m-(r1m^2))/(1-(r1m^2));

% Scaling factor
v=((1+a2)/((1-a2)*(1-a1+a2)*(1+a1+a2)))^(-0.5);

%% Generation of the perturbed field
n = number_ens;
delta1 = zeros(x,t+2,n);
delta2 = zeros(x,t+2,n);

for ens = 1:n
    delta1(:,1,ens)= L*normrnd(0,1,x,1);
    delta1(:,2,ens)= L*normrnd(0,1,x,1);
    for tm = 3:t+2
        delta1(:,tm,ens) = L*normrnd(0,1,x,1) + a1*delta1(:,(tm-1),ens) + a2*delta1(:,(tm-2),ens);
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
        mean_d = nansum(R.*d1)./nansum(R);

        % Covariance Calculation
        diffd = zeros(t,x);
        for i=1:t 
            diffd(i,:) = d1(i,:)- mean_d;
        end
        cov_d = zeros(x,x);
        temp1 = zeros(t,1);
        temp2 = zeros(t,1);
        for k = 1:x
            for l = 1:x
                for tm = 1:t
                    temp1(tm) = R(tm,k)*R(tm,l)*diffd(tm,k)*diffd(tm,l);
                    temp2(tm) = R(tm,k)*R(tm,l);
                end
                cov_d(k,l) = nansum(temp1)/nansum(temp2);
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
        
        r1d=zeros(x,1);
        temp1 = zeros(t-1,1);
        temp2 = zeros(t-1,1);
        for k = 1:x
            for tm = 1:t-1
                temp1(tm) = R(tm,k)*R(tm+1,k)*diffd(tm,k)*diffd(tm+1,k);
                temp2(tm) = R(tm,k)*R(tm+1,k);
            end
            r1d(k) = nansum(temp1)/(nansum(temp2)*cov_d(k,k));
        end
        r1dm = nanmean(r1d);
        r1dmax = max(r1d);
        r1dmin = min(r1d);

        r2d=zeros(x,1);
        temp1 = zeros(t-2,1);
        temp2 = zeros(t-2,1);
        for k = 1:x
            for tm = 1:t-2
                temp1(tm) = R(tm,k)*R(tm+2,k)*diffd(tm,k)*diffd(tm+2,k);
                temp2(tm) = R(tm,k)*R(tm+2,k);
            end
            r2d(k) = nansum(temp1)/(nansum(temp2)*cov_d(k,k));
        end
        r2dm = nanmean(r2d);
        r2dmax = max(r2d);
        r2dmin = min(r2d);

        rd_mean = [1,r1dm,r2dm];
        rd_max = [1,r1dmax,r2dmax];
        rd_min = [1,r1dmin,r2dmin];

        fig = figure;
        figure(fig)
        plot(x_axis,rd_mean)
        hold on
        plot(x_axis,r_mean,'r')
        plot(x_axis, r_max, 'g')
        plot(x_axis, r_min, 'g')
        plot(x_axis, rd_max, 'c')
        plot(x_axis, rd_min, 'c')
        title('temporal correlation of the deltas');
        xlabel('hours')
        name = [results_path, 'temporal_correlation_delta_' num2str(ens)];
        saveas(fig,name,'fig')
        saveas(fig,name,'jpg')
        close(fig)
    end
end

%% saving the results
save ([results_path,'delta.mat'],'delta','coord')

%% Generate the ensembles
ensembles = Ensembles(delta, coord, sim_hours, results_path);
