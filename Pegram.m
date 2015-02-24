%This script will implement the (Pegram et al, 2011) method for ensemble generation

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
results_folder = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_Pegram\Test_11\';
% Brief description for this test:
description = 'This test is for the calibration part. it should return a plot of autocorrelation to decide which noise threshold to use.';

% Parameters and information:

% Starting date and time:
starting = '2008-05-05 15:00:00';
% Ending date and time:
ending = '2008-05-05 16:30:00';
% Time step in minutes:
tstep = 5;
% How many ensembles do you want to generate
number_ens = 5;
% You need to work on a square area with a multiple of 2 number of pixels.
% what are the easting and nothing in meters (top left corner)?
easting = 300000;
northing = 556000;
% What is the number of pixels for the square side?
side = 256;

% Calibration mode (Y/N):
cm = 'Y';
% if not in calibration mode, what is the wavelenght threshold between
% signal and noise in km?
threshold = 8;
% if in calibration mode, which wavelenghts you want to try?
testWL = [4 6 8 10 12];
%   how many time steps you want to test the autocorrelation on (it must be
%   smallet than the total timesteps tsteps)?
calibsteps = 12;

%% Script

tic

% Results folder
if ~exist(results_folder, 'dir')
  mkdir(results_folder);
end
cd(results_folder)
Readme = fopen('README.txt', 'w');
fprintf(Readme, description);
fclose(Readme);

if cm~='Y'
    % Dates:
    s_date = datenum(starting);
    e_date = datenum(ending);
    tst = (1/1440)*tstep;     %transforms the tstep unit in days
    t=0;
    
    for time = s_date:tst:e_date
        
        t=t+1;
        
        % File to work on
        datename = [sprintf('%4i',year(time)),sprintf('%02i',month(time)),sprintf('%02i',day(time)),sprintf('%02i',hour(time)),sprintf('%02i',minute(time))];
        filename = ['metoffice-c-band-rain-radar_uk_',datename,'_1km-composite.dat.gz'];
        cd(data_path)
        tt = length(filename);
        if (strcmp('.gz', filename(tt-2:tt)) == 1)   % if the input file is compressed (i.e. with extension '.gz')
             newname = gunzip(filename, tmpfolder); % decompress file in a temporary folder and output new file name
             filename = newname{1,1};    % copy new name into the nimrodfilename
        end
        currentfullpath = mfilename('fullpath'); filenamesz = (size(mfilename)); currentdir = currentfullpath(1:(end-filenamesz(2)));
        cd(currentdir)        
        [radar, date, param] = readnimrod(filename);
        
        % Data preprocessing
        radar(radar<0.5) = NaN;
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

        % Create noise component in the frequency domain
        NF = zeros(sz,sz,number_ens);
        for en=1:number_ens
            WF = normrnd(0,1,sz,sz);
            for x=1:sz
                for y=1:sz
                    if position(x,y)<=(sz/2)
                        if x>y
                            NF(x,y,en) = WF(x,y) * (Pnavg(position(x,y)).^0.5);   
                        else
                            NF(x,y,en) = WF(y,x) * (Pnavg(position(x,y)).^0.5);
                        end
                    else
                        if x>y
                            NF(x,y,en) = WF(x,y) * (Pnavg(sz/2).^0.5);   
                        else
                            NF(x,y,en) = WF(y,x) * (Pnavg(sz/2).^0.5);
                        end
                    end
                end
            end
        end
        
        % Generate the ensemble
        EF = zeros(sz,sz,number_ens);
        for en=1:number_ens
            EF(:,:,en) = SF(:,:) + NF(:,:,en);
        end

        % Inverse the Fourier transform
        E = zeros(sz,sz,number_ens);
        for en=1:number_ens
            EF(:,:,en) = fftshift(EF(:,:,en));
            E(:,:,en) = ifft2(EF(:,:,en), 'symmetric');   
            E(:,:,en) = exp(E(:,:,en) + m).* mask;
        end
        E(E<0) = 0;
        
        toc
        
        % Print out the ensembles and the radar
        Rain =exp(R + m).*mask;

        cmap = colormap('jet');
        n = size(cmap, 1);
        zerovalue = [1 1 1];
        minimum = 0.0001;
        maximum = max(max(log(Rain(:))), max(log(E(:))));
        dmap = (maximum-minimum)/n;
        for en = 1:number_ens
            fig=figure;
            figure(fig)
            imagesc(log(E(:,:,en)), [0 maximum])
            set(fig, 'Position', [200 200 800 800])
            colormap([zerovalue;cmap]);
            caxis([minimum-dmap maximum]);
            colorbar
            saveas (fig, [results_folder, 'ensemble',num2str(en),'_',datename,'.jpg'])
            close(fig)
        end
        fig=figure;
        figure(fig)
        imagesc(log(Rain), [0 maximum])
        set(fig, 'Position', [200 200 800 800])
        colormap([zerovalue;cmap]);
        caxis([minimum-dmap maximum]);
        colorbar
        saveas (fig, [results_folder, 'radar','_',datename,'.jpg'])
        close(fig)
        
        toc
        
        save([results_folder,'Ensemble_',datename,'.mat'])
        
        toc
    
    end
    
%     for  time = s_date:tstep:e_date
%         plotfilename = [results_folder,'Ensemble_',datename,'.mat'];
%         PlotPower(plotfilename)
%     end
    
    toc
else
    
    ntests = size(testWL,2);
    autocorrel = zeros(calibsteps,ntests);
    
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
            tt = length(filename);
            if (strcmp('.gz', filename(tt-2:tt)) == 1)   % if the input file is compressed (i.e. with extension '.gz')
                 newname = gunzip(filename, tmpfolder); % decompress file in a temporary folder and output new file name
                 filename = newname{1,1};    % copy new name into the nimrodfilename
            end
            currentfullpath = mfilename('fullpath'); filenamesz = (size(mfilename)); currentdir = currentfullpath(1:(end-filenamesz(2)));
            cd(currentdir)        
            [radar, date, param] = readnimrod(filename);

            % Data preprocessing
            radar(radar<0.5) = NaN;
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
            masktemp = zeros(size(R));
            masktemp(isnan(R)==0)=1;
            mask(:,:,t) = masktemp;

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
            
            N = zeros(sz,sz,t);
            for tm=1:t
                NF(:,:,tm) = fftshift(NF(:,:,tm));
                N(:,:,tm) = ifft2(NF(:,:,tm), 'symmetric');   
                N(:,:,tm) = (exp(N(:,:,tm))).*mask(:,:,tm);
            end
                      
            toc
        end

        toc
        
        timesteps = size(N,3);
        mu = mean(N,3);
        sigma = std(N,0,3);

        diff = zeros(sz,sz,timesteps);
        
        for tm=1:timesteps
            diff(:,:,tm) = N(:,:,tm) - mu;
        end

        for lag=1:calibsteps
            diff_t1 = diff(:,:,1:(end-lag));
            diff_t2 = diff(:,:,(1+lag):end);
            temp = zeros(sz,sz,t-lag);
            for tm=1:timesteps-lag
                temp(:,:,tm) = diff_t1(:,:,tm).*diff_t2(:,:,tm);
            end
            ac = mean(temp, 3)./ sigma;
            autocorrel(lag, test) = nanmean(nanmean(ac));
        end
        
        toc

    end
    
    toc
    
    save ([results_folder,'calibrationtest'])
    
    timeline = (1:12).*tstep;
    linecolours = colormap('jet');
    fig = figure;
    figure(fig)
    hold on
    for test = 1:ntests
        plot(timeline',autocorrel(:,test),'color',linecolours(floor((64/ntests)*test),:))
    end
    saveas(fig,[results_folder,'autocorrelation.jpg'])
    
    toc
    
end

    
        
        
    

