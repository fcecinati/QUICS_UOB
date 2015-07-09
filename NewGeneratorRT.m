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
results_path = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_NewGenerator\Test25\';
% Path to the complete radar data
radar_folder = 'Z:\uk-1km\hourly\2008\';
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
% As far as you know, which of the stations have unusable data?
known_corrupted_data_stations = [8;30;31;73;78;81;98;100;107;111;141;145;152;212;213;228;229;168;174;178;194;20;40;128;208;211];
% Description for the readme file:
description = 'First attempt with the FFT-MA generator';

% Parameters:

% How many ensembles do you want to generate
n_ensembles = 100;
% What is the starting date of the ensembles you need?
starting = '2008-05-04 00:00:00';
% How many hours is your simulation long?
time_steps = 100;
% what are the easting and nothing in meters (top left corner)?
easting = 300000;
northing = 556000;
% What is the number of pixels for the square side?
side = 256;
% How many hours do you want to consider for the error characteristics?
n_hours = 12;

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
disp('setup complete')

%% Read and preprocess data
% Read rain gauges and radar data in .dat format
cd(data_path)
coord = (load(coordinates_file))';
G = load(rain_gauge_file); 
R = load(corresponding_radar_file);
years = int16(R(:,1)); months = int16(R(:,2)); days = int16(R(:,3)); hours = int16(R(:,4));
G=G(:,7:end);  R=R(:,7:end);
currentfullpath = mfilename('fullpath'); filenamesz = (size(mfilename)); currentdir = currentfullpath(1:(end-filenamesz(2)));
cd(currentdir)  

% Use of NaN for negative values and to ignore days of no rain
G(G<0.1) = NaN;
R(R<0.1) = NaN;

% Delete the stations with No Data or corrupted data (from the list above)
corrupted = known_corrupted_data_stations;
coord(:,corrupted) = [];
G(:,corrupted) = [];
R(:,corrupted) = [];

% Final size of the data matrix
sz = size(R);
t=sz(1);
x=sz(2);
G = G(1:t,:);
toc
disp('data loading and preprocessing complete')

%% Calculation of error statistics

% Errors
er = 10*(log(R./G));

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

toc
disp('statistics complete')

%% Semivariogram

% General average semivariogram
dis = min(side/2, distmax);
d=1:dis;
variog=1:dis;
for i=1:dis   
    [x1,x2]=find(dist>=(i-1)&dist<i);
    mat = ((er(:,x1)-er(:,x2)).^2);
    variog(i)= nanmean(mat(:));
end

myfunction = fittype('c.*(1-exp(-(h./a)))+b','independent',{'h'},'coefficients',{'a','b','c'});
myoption = fitoptions('Method','NonlinearLeastSquares','StartPoint',[50 1 (std_e^2)]);
myfit = fit(d',variog',myfunction,myoption);
sill = myfit.c;
range = myfit.a;
nugget = myfit.b;

fig = figure;
figure(fig)
plot(myfit,d',variog')
xlabel('distance (km)')
ylabel('semivariogram')
title('general semivariogram')
saveas(fig,[results_path, 'GeneralSemivariogram'],'fig')
close

toc
disp('general semivariogram calculated')

% Specific semivariogram
sill_spec = zeros(time_steps,1);
range_spec = zeros(time_steps,1);
nug_spec = zeros(time_steps,1);
stdev = zeros(time_steps,1);
mu = zeros(time_steps,1);
bbest = zeros(time_steps,1);
time = 0;
fit_RMSE = zeros(time_steps,1);
fit_ARS = zeros(time_steps,1);
for tt = datenum(starting):(1/24):datenum(starting)+(time_steps-1)/24
    time = time + 1;
    dtst = datestr(tt); yr = year(dtst); mt = month(dtst); dy = day(dtst); hr = hour(dtst);
    dtposition = find(years==yr & months==mt & days==dy & hours==hr);
    
    errors = er(dtposition-n_hours:dtposition-1,:);
    
    % Semivariogram
    variog=1:dis;
    for i=1:dis   
        [x1,x2]=find(dist>=(i-1)&dist<i);
        mat = ((errors(:,x1)-errors(:,x2)).^2);
        variog(i)= nanmean(mat(:));
    end
    
    % Exclude the NaN
    di = d;
    di(isnan(variog)==1)=[];
    variog(isnan(variog)==1)=[];
    di = di(di<dis/2);
    variog = variog(di<dis/2);

    %Semivariogram fitting
    try
        myfunction2 = fittype('c.*(1-exp(-(h./a)))+b','independent',{'h'},'coefficients',{'c','b','a'});
        myoption2 = fitoptions('Method','NonlinearLeastSquares','StartPoint',[sill nugget range]);
        [myfit2, gof] = fit(di',variog',myfunction2,myoption2);
        sill_spec(time) = myfit2.c;
        range_spec(time) = myfit2.a;
        nug_spec(time) = myfit2.b;
        fit_RMSE(time) = gof.rmse;
        fit_ARS(time) = gof.adjrsquare;

        if (fit_RMSE(time)/sill_spec(time)<0.4) && (size(variog,2)>20) && range_spec(time)<dis
%             fig = figure;
%             figure(fig)
%             plot(myfit2,di',variog')
%             xlabel('distance (km)')
%             ylabel('semivariogram')
%             title(['specific semivariogram',datestr(tt)])
%             saveas(fig,[results_path, 'SpecificSemivariogram',num2str(time)],'fig')
%             close
            stdev(time) = nanstd(errors(:));
            mu(time) = nanmean(errors(:));
        else
            sill_spec(time) = sill;
            range_spec(time) = range;
            nug_spec(time) = nugget;
            stdev(time) = std_e;
            mu(time) = mean_e;
        end
    catch
        sill_spec(time) = sill;
        range_spec(time) = range;  
        nug_spec(time) = nugget;
        stdev(time) = std_e;
        mu(time) = mean_e;
    end    
end

toc
disp('specific semivariogram calculated')


%% Error component generation with the FFT-MA

ec = zeros(side,side,time_steps,n_ensembles);
rot = 90; %anisotropy angle of the primary direction clockwise from North
anis = 1; %anis = fraction of range of the range in the secondary direction (perpendicular to tha primary direction) to the range of the primary direction (1 for isotropy)
x_step = 1:side;
y_step = x_step;
for i=1:time_steps
    Va=[num2str(nug_spec(i)),' Nug(0) + ',num2str(sill_spec(i)-nug_spec(i)),'  Exp(',num2str(range_spec(i)),',',num2str(rot),',',num2str(anis),')'];
    for j=1:n_ensembles
        [out,~]=fft_ma_2d(x_step,y_step,Va);
        sd_out = std(out);        out = out/sd_out*stddev(i);
        m_out = mean(out);
        out = out + mu(i) - m_out;
        ec(:,:,i,j) = out;
    end
end

save([results_path,'error_components'],'ec')

toc
disp('error component generated')

%% Setup to generate ensembles
currentfullpath = mfilename('fullpath'); filenamesz = (size(mfilename)); currentdir = currentfullpath(1:(end-filenamesz(2)));
cd(currentdir) 
st = datenum(starting);
raddata = zeros(side,side,time_steps);
ens = zeros(side,side,time_steps,n_ensembles);

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
        dnew = ec(:,:,i,en);
        
        % Combine noise and radar data
        logens = dnew + 10*(log(raddata(:,:,i)));
        ens(:,:,i,en) = exp(logens/10);
    end
end

toc
disp('ensembles generated')

save([results_path,'generatedensembles'],'raddata','ens')

for ts = 1:time_steps
    RadarPlot(raddata,ens, results_path,ts)
end

toc
disp('all saved. done. goodbye!')

% % Save results in Nimrod format
% 
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




