function [] = NewGeneratorRT(starting)

% This script generates radar ensembles with a new method...
%% 
rng('shuffle')
tic
%% INFORMATION TO PROVIDE:

% Folders:

% Results path
% results_path = ['F:\Francesca\HBV',num2str(year(starting)), num2str(month(starting),'%02d'),'\'];
results_path = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_NewGenerator\ForPics';
% Path to the complete radar data
radar_folder = ['Z:\uk-1km\hourly\',num2str(year(starting)),'\'];
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
description = ['Data generated for the test in Nanding model. Time: ',num2str(month(starting),'%02d'),'-',num2str(year(starting)),'. Space: 320-500 km east - 370-500 km north. 1km scale. 12 hours window for error estimation. 100 ensembles.'];

% Parameters:

% How many ensembles do you want to generate
n_ensembles = 100;
% What is the starting date of the ensembles you need?
%starting = '2007-12-26 00:00:00';
% How many hours is your simulation long?
time_steps = eomday(year(starting), month(starting))*24;
% what are the easting and nothing in meters (top left corner)?
easting = 320000;
northing = 550000;
% What is the number of pixels for the square side?
side = 180;
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

myfunction = fittype('c.*(1-exp(-(3.*h./a)))','independent',{'h'},'coefficients',{'a','c'});
myoption = fitoptions('Method','NonlinearLeastSquares','StartPoint',[50 (std_e^2)]);
myfit = fit(d',variog',myfunction,myoption);
sill = myfit.c;
range = myfit.a;

% f = figure;
% figure(f)
% set(f,'Position',[100,100,500,300])
% plot(myfit,d',variog')
% xlabel(['distance (km) - range = ',num2str(range)])
% ylabel(['semivariogram - sill = ',num2str(sill)])
% title('general semivariogram')
% saveas(f,[results_path, 'GeneralSemivariogram'],'fig')


toc
disp('general semivariogram calculated')

% Specific semivariogram
sill_spec = zeros(time_steps,1);
range_spec = zeros(time_steps,1);
stdev = zeros(time_steps,1);
mu = zeros(time_steps,1);
fit_RMSE = zeros(time_steps,1);
time = 0;
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
        variog(i)= nanmean(mat(:))./2;
    end
    
    % Exclude the NaN
    di = d;
    di(isnan(variog)==1)=[];
    variog(isnan(variog)==1)=[];
    di = di(di<dis/2);
    variog = variog(di<dis/2);

    %Semivariogram fitting
    try
        myfunction2 = fittype('c.*(1-exp(-(3.*h./a)))','independent',{'h'},'coefficients',{'c','a'});
        myoption2 = fitoptions('Method','NonlinearLeastSquares','StartPoint',[sill range]);
        [myfit2, gof] = fit(di',variog',myfunction2,myoption2);
        sill_spec(time) = myfit2.c;
        range_spec(time) = myfit2.a;
        fit_RMSE(time) = gof.rmse;
            f = figure;
            figure(f)
            set(f,'Position',[100,100,500,300])
            plot(myfit2,di',variog')
            xlabel(['distance (km) - range = ',num2str(range_spec(time))])
            ylabel(['semivariogram - sill = ',num2str(sill_spec(time))])
            title(['specific semivariogram ',datestr(tt)])
            saveas(f,[results_path, 'Semivariogram',num2str(time)],'fig')
%             close

        if (fit_RMSE(time)/sill_spec(time)<0.4) && (size(variog,2)>20) && range_spec(time)<dis
            stdev(time) = nanstd(errors(:));
            mu(time) = nanmean(errors(:));

        else
            sill_spec(time) = sill;
            range_spec(time) = range;
            stdev(time) = std_e;
            mu(time) = mean_e;
        end
    catch
        sill_spec(time) = sill;
        range_spec(time) = range;  
        stdev(time) = std_e;
        mu(time) = mean_e;
    end    
end

save([results_path,'statistics'],'sill_spec','range_spec','mu','stdev', 'fit_RMSE','sill','range')

toc
disp('specific semivariogram calculated')

clear G R coordinates_file corresponding_radar_file corrupted currentdir currentfulpath d description di dis distance dtposition dtst er errors filenamesz hours known_corrupted_data_stations mat months rain_gauge_file sz tt variog x x1 x2 x_step y_step years yr

%% Generate ensembles

% Variable setup
rot = 90; %anisotropy angle of the primary direction clockwise from North
anis = 1; %anis = fraction of range of the range in the secondary direction (perpendicular to tha primary direction) to the range of the primary direction (1 for isotropy)
x_step = 1:side;
y_step = x_step;
st = datenum(starting);
raddata = zeros(side,side);
ens = zeros(side,side,n_ensembles);
err = zeros(side,side,n_ensembles);
missing = zeros(time_steps,1);

for i = 1:time_steps
    try
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
        raddata = mtx(y_start:y_end, x_start:x_end);
        
        % Set up the corresponding semivariogram
        Va=[num2str(sill_spec(i)),'  Exp(',num2str(range_spec(i)),',',num2str(rot),',',num2str(anis),')'];
        for en=1:n_ensembles    
            
            % Generates the error components with the fft_ma
            [out,~]=fft_ma_2d(x_step,y_step,Va);
            
            % Adjusts mean and variance
            sd_out = std(out(:));        out = out/sd_out*stdev(i);
            m_out = mean(out(:));        out = out + mu(i) - m_out;
            err(:,:,en) = out;
            
            % Combine noise and radar data
            logens = out + 10*(log(raddata));
            ens(:,:,en) = exp(logens/10);                        
        end
        save([results_path,'ensemble-time',num2str(i)],'ens','err','raddata')
    catch
        % Counts the missing files/steps
        missing(i)=1;        
    end
end
save('missing','missing')    
toc
disp('ensembles generated')

%% Variance and Mean Inflation Correction
for i = 1:time_steps
    try
        load([results_path,'ensemble-time',num2str(i)])
        newens = VarInflation_adjustment(raddata,ens);
        save([results_path,'ensemble-time',num2str(i)],'ens','err','raddata','newens')
    catch
    end
end

toc
disp('ensembles corrected')

%% Empty Temp
delete('\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Temp\*.dat')

toc
disp('all done. goodbye!')