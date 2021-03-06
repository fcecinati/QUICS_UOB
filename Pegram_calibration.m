%This script will implement calibration for the the (Pegram et al, 2011) method for ensemble generation

clear all
close all
clc
rng('shuffle')

%% Parameters to define:

%Paths:

% Path to the folder where the radar data is:
data_path = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\data_north_england\05-2008';
% Temporary folder:
tmpfolder = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Temp\'; 
% Results folder:
results_folder = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_Pegram\Test_18\';
% Brief description for this test:
description = 'extensive calibration in reduced number of points';

% Parameters and information:

% Starting date and time:
starting = '2008-05-01 04:00:00';
% Ending date and time:
ending = '2008-05-01 16:00:00';
% Time step in minutes:
tstep = 5;
% You need to work on a square area with a multiple of 2 number of pixels.
% what are the easting and nothing in meters (top left corner)?
easting = 300000;
northing = 556000;
% What is the number of pixels for the square side?
side = 256;
% Which wavelenghts you want to try?
testWL = [3 6 8 10 12];
%   how many time steps you want to test the autocorrelation on (it must be
%   smallet than the total timesteps tsteps)?
calibsteps = 12;

% Results folder
if ~exist(results_folder, 'dir')
  mkdir(results_folder);
end
cd(results_folder)
Readme = fopen('README.txt', 'w');
fprintf(Readme, description);
fclose(Readme);

ntests = size(testWL,2);
autocorrel = zeros(calibsteps,ntests);
tic
for test=1:ntests

    threshold = testWL(test);
    
    % Dates:
    s_date = datenum(starting);
    e_date = datenum(ending);
    tst = (1/1440)*tstep;     %transforms the tstep unit in days

    toc 
    t=0;
    for time = s_date:tst:e_date

        t=t+1;

        % File to work on
        datename = [sprintf('%4i',year(time)),sprintf('%02i',month(time)),sprintf('%02i',day(time)),sprintf('%02i',hour(time)),sprintf('%02i',minute(time))];
        filename = ['metoffice-c-band-rain-radar_uk_',datename,'_1km-composite.dat.gz'];
        cd(data_path)
        if exist(filename, 'file')
            tt = length(filename);
            if (strcmp('.gz', filename(tt-2:tt)) == 1)   % if the input file is compressed (i.e. with extension '.gz')
                 newname = gunzip(filename, tmpfolder); % decompress file in a temporary folder and output new file name
                 filename = newname{1,1};    % copy new name into the nimrodfilename
            end
            currentfullpath = mfilename('fullpath'); filenamesz = (size(mfilename)); currentdir = currentfullpath(1:(end-filenamesz(2)));
            cd(currentdir)        
            [radar, date, param] = readnimrod(filename);

            % Data preprocessing
            radar(radar<0.05) = NaN;
            lograd = log(radar);
            m = nanmean(nanmean(lograd));
            Z = lograd - m;

            % Trim at a multiple of 2 size
            x = (param(3):param(5):param(3)+(param(7)*param(5)-param(5)));
            y = (param(2):-param(4):(param(2)-(param(6)*param(4))+param(4)));
            res = param(4)/1000;    % in km
            eastingmax = easting + side*param(5);
            northingmin = northing - side*param(4);
            indx = (x>easting&x<eastingmax);
            indy = (y>northingmin&y<northing);
            R = Z(indy,indx);
            sz = side;

            % Creates a rain/norain mask
            mask = zeros(size(R));
            mask(isnan(R)==0)=1;

            % Substitute NaNs with zeros for the fft
            R(isnan(R)==1) = 0;

            %% Pegram method

            % Set frequencies and wavelength
            fx = (1/res).*(-(sz/2):(sz/2)); fx = (fx(2:end)-0.5)./sz; %km^-1
            fy = fx;
            [FX, FY] = meshgrid (fx, fy);
            F = (FX.^2 + FY.^2).^0.5;
            WL = 1./F;                               % km
            WLlog = log(WL);

            % Fast Fourier Transform
            RF = fft2(R);     RF = fftshift(RF);
            P = (abs(RF)).^2;
            Plog = log(P);

            % Line fitting and Ns identification
            Pcenter = P(WL>(threshold-3)&WL<(threshold+3));
            WLcenter = WL(WL>(threshold-3)&WL<(threshold+3));
            Pc = log(Pcenter);
            WLc = log(WLcenter);
            linearCoef = polyfit(WLc,Pc,1);
            Pline = WLlog.*linearCoef(1) + linearCoef(2);
            Plineexp = exp(Pline);
            ns = log(threshold)*linearCoef(1) + linearCoef(2); ns = exp(ns);

            % Separation of Signal and Noise
            Ps = NaN(sz,sz);
            Pn = NaN(sz,sz);
            Ps(P>ns) = P(P>ns);
            Pn(P<ns) = P(P<ns);
            [Pavg, rad] = nanradialavg(P, (sz/2));
            [Psavg, ~] = nanradialavg(Ps, (sz/2));
            [Pnavg, ~] = nanradialavg(Pn, (sz/2));
            Pavg(1) = Pavg(2); Psavg(1) = Psavg(2); Pnavg(1) = Pnavg(2);

            % Create signal component in the frequency domain
            position = round(F.*sz);
            SF = zeros(sz,sz);
            for x=1:sz
                for y=1:sz
                    if position(x,y)<=(sz/2)
                        SF(x,y) = RF(x,y) * (Psavg(position(x,y))/Pavg(position(x,y)));
                    end
                end
            end

            NF(:,:,t) = RF - SF;
            NF(:,:,t) = fftshift(NF(:,:,t));
            N(:,:,t,test) = ifft2(NF(:,:,t), 'symmetric');   
            N(:,:,t,test) = (exp(N(:,:,t,test)));

            toc
        end
    end
    subN = zeros(100, 145);
    for i=1:100
        subX = randi(sz);    subY = randi(sz);
        subN(i,:) = squeeze(N(subX,subY,:,test));
    end
    subN = mean(subN, 1);
    [acf(:,test),~,~] = autocorr(subN, 12);
end

toc
save ([results_folder,'calibrationtest'])

Nm = squeeze(mean(mean(N,1),2));
for test=1:5
    [acf2(:,test),~,~] = autocorr(Nm(:,test), 12);
end



leg = cell(ntests,1);
timeline = (1:size(acf2,1))*tstep;
linecolours = colormap('jet');
fig = figure;
figure(fig)
hold on
for test = 1:ntests
    plot(timeline',acf2(:,test),'color',linecolours(floor((64/ntests)*test),:))
    leg{test} = num2str(testWL(test));
end
legend(leg)
saveas(fig,[results_folder,'autocorrelation.jpg'])

toc

