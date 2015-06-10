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
results_path = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_NewGenerator\TestEAWAG\';
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
%Which year do you want to consider for the error statistics?
year_to_consider = 2007;
% As far as you know, which of the stations have unusable data?
known_corrupted_data_stations = [8;30;31;73;78;81;98;100;107;111;141;145;152;212;213;228;229;168;174;178;194;20;40;128;208;211];
% Description for the readme file:
description = 'Some test on the error characteristics';

% Parameters:

% How many ensembles do you want to generate
n_ensembles = 100;
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
% If not, what are the coefficients of the semivariogram exponential function?
sill = 78.7019;
range = 57.7507;
F = 0.4715;
% What are the statistical characteristics of the errors?
mu = 0;
sigma = 5.6855;

% Do you need to calibrate the filter?
filter_calibr = 'N';
% If yes, which values you would like to try:
b_test = 30:10:150;
% Otherwise, which are the coefficients of the filter?
abest=128;
bbest=60;

%% Start script

% Results folder
if ~exist(results_path, 'dir')
  mkdir(results_path);
end
cd(results_path)
Readme = fopen('README.txt', 'w');
fprintf(Readme, description);
fclose(Readme);
toc
disp('setup')

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
    disp('data loading and preparing')
    %% Calculation of error statistics

    % Errors
    er = 10*(log(G./R));
    
    % Error Mean
    mean_e = nanmean(er(:));
    
    % Standard Deviation
    std_e = nanstd(er(:));
    
    % Distance between the stations
    distance = zeros(x,x);
    for x1=1:x
        for x2=1:x
            distance(x1,x2) = (((coord(1,x1)-coord(1,x2))^2)+((coord(2,x1)-coord(2,x2))^2))^0.5;
        end
    end
    distmax = ceil(max(distance(:))/1000);
    dist=distance/1000;
    d=1:distmax;
    
    % Semivariogram
    variog=1:distmax;
    for i=1:distmax    
        [x1,x2]=find(dist>=(i-1)&dist<i);
        mat = ((er(:,x1)-er(:,x2)).^2);
        variog(i)= nanmean(mat(:));
    end
    
    %Semivariogram fitting
    myfunction2 = fittype('c.*(1-exp(-(h./a).^F))','independent',{'h'},'coefficients',{'c','a','F'});
    myoption2 = fitoptions('Method','NonlinearLeastSquares','StartPoint',[5 50 0.5]);
    myfit2 = fit(d',variog',myfunction2,myoption2);
    sill = myfit2.c;
    range = myfit2.a;
    F = myfit2.F;
    save([results_path,'SemivariogramCoefficients.mat'],'sill','range')
end
toc
disp('original semivariogram')

%% Calibration of filter
if filter_calibr == 'Y'

    %Setup
    cd('\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\QUICS_UOB')
    h = (1:128)';
    b_size = size(b_test,2);
    v=zeros(b_size,128,20);
    j=0;
  
    
  
    %Test for every a and every b
    
    % Do 20 tests to average
    for k=1:20 
        %Generate normal random fields
        d=normrnd(mu,sigma,256,256);
        variog = zeros(b_size,128);
        for bi=b_test
            j = j + 1;
                
            % Apply a lowpass filter
            base = 1:abest;
            fun = (1/(sigma^2))*exp(-(base./bi).^F);
            f1 = [fliplr(fun),(1/(sigma^2)),fun];
            f2 = ftrans2(f1);
            dnew = filter2(f2,d);
            newstd = std(dnew(:));
            dnew = dnew*sigma/newstd;

            % Calculate the semivariogram
            [hvar, vvar] = Variogram(dnew);
            variog(j,:) = (hvar + vvar)/2;
        end
        v(:,:,k)=variog;
        j=0;
    end
    v = mean(v,3);
    j=0;
  
    % Calculate the RMSE of each obtained semivariogram
    ref_v = (sill*(1-exp(-(h./range).^F)))';
    for bi=b_test
        j = j + 1;
        rmse_m(j) = mean((v(j,:) - ref_v).^2).^0.5;
    end
    plot(b_test,rmse_m)
    xlabel('values of the parameter b')
    ylabel('RMSE')
    title('RMSE of the semivariogram tests on the filter parameter b')
    print([results_path,'RMSE min small'],'-dpng','-r300')
    
    % Identify teh lowest RMSE and the corresponding a,b values
    bbestindex = find(rmse_m == min(rmse_m(:)));
    bbest = b_test(bbestindex);
    
    fig2 = figure;
    figure(fig2)
    plot(b_test,rmse_m)
    xlabel('values of the coefficient b')
    ylabel('RMSE')
    title('values of RMSE varying the coefficient b')
    saveas(fig2,[results_path, 'RMSE average'],'fig')
    print([results_path,'RMSE min small'],'-dpng','-r300')
    
end
toc
disp('filter calibration')

% Setup to generate ensembles
currentfullpath = mfilename('fullpath'); filenamesz = (size(mfilename)); currentdir = currentfullpath(1:(end-filenamesz(2)));
cd(currentdir) 
st = datenum(starting);
raddata = zeros(side,side,time_steps);
ens = zeros(side,side,time_steps,n_ensembles);
noise = zeros(256,256,25,100);

for i = 1:time_steps
    % Find the radar data for the correct date
    date = st + (i-1)/24;
    yr = year(date);    mn = month(date);   dy = day(date);     hr = hour(date);
    datename = [num2str(yr),num2str(mn, '%02d'),num2str(dy, '%02d'),num2str(hr, '%02d')];
    filename = [radar_folder,num2str(yr),num2str(mn, '%02d'),num2str(dy, '%02d'),num2str(hr, '%02d'),'00.dat.gz'];

    % Check if file is compressed
    tt=length(filename);
    if ( strcmp('.gz', filename(tt-2:tt)) == 1)   
         newname = gunzip(filename, tmpfolder); 
         filename = newname{1,1};    
    end

    % Open the radar data
    [mtx, ~, r_param] = readnimrod(filename);   % function to read nimrod Met Office rainfall data
    mtx = flipud(mtx);
    
    % Define the domain
    x_start = (easting - r_param(3))/ 1000;
    x_end = (easting - r_param(3))/ 1000 + side - 1;
    y_start = (r_param(2) - northing)/ 1000;
    y_end = (r_param(2) - northing)/ 1000 + side - 1;

    % Trim the radar data on the domain
    raddata(:,:,i) = mtx(y_start:y_end, x_start:x_end);
    
    for en=1:n_ensembles    
        %% Generate noise fields
        %Generate normal random fields
        d=normrnd(0,sigma,256,256);
        stdd = std(d(:));

        % Apply a lowpass filter
        base = 1:abest;
        fun = (1/(stdd^2))*exp(-(base./bbest).^F);
        f1 = [fliplr(fun),(1/(stdd^2)),fun];
        f2 = ftrans2(f1);        
        dnew = filter2(f2,d);
        
        %Variance scaling
        newstd = std(dnew(:));
        dnew = dnew*stdd/newstd;
        noise(:,:,i,en)=dnew;
        
        % Combine noise and radar data
        logens = dnew + 10*(log(raddata(:,:,i)));
        ens(:,:,i,en) = exp(logens/10);
    end
end

toc
disp('ensembles generated')

save([results_path,'generatedensembles'],'raddata','ens')
save([results_path,'noise'],'noise')

toc
disp('all saved. done. goodbye!')

%% Save results in Nimrod format

% PARAM = e_param;
% 
% for en = 1:n_ensembles
%     for time = 1:time_steps
%         DATE = r_date(:,time);
%         name = [ensembles_folder, 'Ensemble_n',num2str(en),'_',sprintf('%02d',r_date(:,time)),'.dat'];
%         MTX = ens(:,:,time,en);
%         
%         save2nimrod(name,MTX,DATE,PARAM)
%         
%     end
% end




