% This script generates radar ensembles with a new method...
%% 
clear all
close all
clc
rng('shuffle')
tic
%% INFORMATION TO PROVIDE:

% Folders:

% Results path
results_path = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_NewGenerator\Test1\';
% Path to the complete radar data
radar_folder = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\data_north_england\test20\';
% Temporary folder
tmpfolder = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Temp\'; 
% Path to the folder where the rain gauges, the corresponding radar and the coordinate .dat files are:
data_path = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\data_north_england\';
% Name of the coordinate files of the observation points:
coordinates_file = 'gaugecoordinates.dat';
% Name of the rain gauges .dat file with extension:
rain_gauge_file = 'gauge_1h_2007-2011.dat';
% Name of the corresponding .dat radar file:
corresponding_radar_file = 'radar_1h_2007-2010.dat';

% Which year do you want to use for the spatial correlation calculation?
year_to_consider = 2007;
% As far as you know, which of the stations have unusable data?
known_corrupted_data_stations = [8;30;31;73;78;81;98;100;107;111;141;145;152;212;213;228;229;168;174;178;194;20;40;128;208;211];
% Description for the readme file:
description = 'First test';

% Parameters:

% How many ensembles do you want to generate
n_ensembles = 1;
% What is the starting date of the ensembles you need?
starting = '2008-05-04 00:00:00';
% How many hours is your simulation long?
time_steps = 25;
% what are the easting and nothing in meters (top left corner)?
easting = 300000;
northing = 556000;
% What is the number of pixels for the square side?
side = 256;

% Do you need to calculate the original spatial correlation (i.e. you don't know the coefficients)?
spa_corr_calc = 'Y';

% If not, what are the coefficients of the semivariogram?
R0 = 9.437809846623692;
F = 0.562698743589422;
% What are the statistical characteristics of the errors?
mu = -0.00348;
sigma = 5.6855;

% Do you need to calibrate the filter?
filter_calibr = 'Y';

% Otherwise, which are the coefficients of the filter?
a=50;
b=0.1;

%% Start script
toc

% Results folder
if ~exist(results_path, 'dir')
  mkdir(results_path);
end
cd(results_path)
Readme = fopen('README.txt', 'w');
fprintf(Readme, description);
fclose(Readme);

if spa_corr_calc == 'Y'
    %% Read data
    % Read rain gauges and radar data in .dat format
    cd(data_path)
    coord = (load(coordinates_file))';
    G_file = load(rain_gauge_file);
    R_file = load(corresponding_radar_file);
    years = int16(R_file(:,1));
    currentfullpath = mfilename('fullpath'); filenamesz = (size(mfilename)); currentdir = currentfullpath(1:(end-filenamesz(2)));
    cd(currentdir)  

    % Data is limited to one year for covariance calculation purposes:
    Rsize = size(R_file);
    G = G_file((years==year_to_consider),7:Rsize(2));
    R = R_file((years==year_to_consider),7:Rsize(2));
    clear G_file R_file 

    % Use of NaN for negative values and to ignore days of no rain
    G(G<0.5) = NaN;
    R(R<0.5) = NaN;
    G(G==0) = NaN;
    R(R==0) = NaN;

    % Delete the stations with No Data or corrupted data (from the list above)
    corrupted = known_corrupted_data_stations;
    coord(:,corrupted) = [];
    G(:,corrupted) = [];
    R(:,corrupted) = [];

    % Final size of the data matrix
    sz = size(R);
    t=sz(1);
    x=sz(2);
    toc
    %% Calculation of error statistics

    % Error matrix
    err = 10*(log(G./R));
    
    % Mean
    mean_e = nanmean(err);

    % Standard Deviation
    std_e = nanstd(err);

    % Overall variance
    var_e = nanvar(err(:));

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
    toc
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
    toc
    % Reshaping and NaN elimination for function estimation and plotting
    c = reshape(correl,x*x,1);
    d = reshape(distance,x*x,1);
    d(isnan(c)) = [];
    c(isnan(c)) = [];
    d=d/1000;         %Important! the function is calculated on kilometers as unit

    % Fit on decorrelation function
    myfunction1 = fittype('exp(-((d/R0)^F))','independent',{'d'},'coefficients',{'R0','F'});
    myoption1 = fitoptions('Method','NonlinearLeastSquares','StartPoint',[50 0.9]);
    myfit1 = fit(d,c,myfunction1,myoption1);
end
toc

R0 = myfit1.R0;
F = myfit1.F;
h=(1:200)';
corr_function = exp(-((h./R0).^F));
variog_function = var_e.*(1-corr_function);

myfunction2 = fittype('c.*(1-exp(-3.*(h.^F)./a.^F))','independent',{'h'},'coefficients',{'c','a'},'problem','F');
myoption2 = fitoptions('Method','NonlinearLeastSquares','StartPoint',[7 40]);
myfit2 = fit(h,variog_function,myfunction2,myoption2,'problem',F);

sill = myfit2.c;
range = myfit2.a;

save('\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_NewGenerator\Test1\Coefficients.m','sill','range')
toc



%% Calibration of filter
toc
cd('\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\QUICS_UOB')
h = (1:128)';
sill_matrix = zeros(10,10);
range_matrix = zeros(10,10);
i=0; j=0;
for i=1:10
    a=i*10;
    sigma1 = (sigma)/(1.09641257*(a^(-0.94003379)));
    for j=1:10
        b=j*0.002;
        variog = zeros(5,128);
        for k=1:10
            %Generate normal random fields
            d=normrnd(0,sigma1,256,256);

            % Apply a lowpass filter
            f1=fir1(a,b,'low');
            f2 = ftrans2(f1);
            dnew = filter2(f2,d);

            % Calculate the semivariogram
            [hvar, vvar] = Variogram(dnew);
            variog(k,:) = (hvar + vvar)/2;
        end
        variog = mean(variog,1)';
        
        myfunction = fittype('c.*(1-exp(-3.*(h)./(a)))','independent',{'h'},'coefficients',{'c','a'});
        myoption = fitoptions('Method','NonlinearLeastSquares','StartPoint',[7 40]);
        myfit = fit(h,variog,myfunction,myoption);

        sill_matrix(i,j) = myfit.c;
        range_matrix(i,j) = myfit.a;
        
    end
end
toc
rse_matrix = ((sill_matrix - sill).^2 + (range_matrix-range).^2).^0.5;
[abest,bbest] = find(rse_matrix == min(rse_matrix(:)));



myfit = fit(h,variog,myfunction,myoption);
        



