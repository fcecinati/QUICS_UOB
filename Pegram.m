%This script will implement the (Pegram et al, 2011) method for ensemble generation

clear all
close all
clc
rng('shuffle')

%% Parameters to define:

%Paths:

% Path to the folder where the radar data is:
data_path = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\data_north_england\test20\';
% Temporary folder:
tmpfolder = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Temp\'; 
% Results folder:
results_folder = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_Pegram\Test_20\';
% Brief description for this test:
description = 'Test on hourly data to be compared with Germanns test 20 and rain gauge data';

% Parameters and information:

% Starting date and time:
starting = '2008-05-04 00:00:00';
% Ending date and time:
ending = '2008-05-05 00:00:00';
% Time step in minutes:
tstep = 60;
% How many ensembles do you want to generate
number_ens = 100;
% You need to work on a square area with a multiple of 2 number of pixels.
% what are the easting and nothing in meters (top left corner)?
easting = 300000;
northing = 556000;
% What is the number of pixels for the square side?
side = 256;
% if not in calibration mode, what is the wavelenght threshold between
% signal and noise in km?
threshold = 8;


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


% Dates:
s_date = datenum(starting);
e_date = datenum(ending);
tst = (1/1440)*tstep;     %transforms the tstep unit in days
t=0;

for time = s_date:tst:e_date

    t=t+1;

    % File to work on
%     datename = [sprintf('%4i',year(time)),sprintf('%02i',month(time)),sprintf('%02i',day(time)),sprintf('%02i',hour(time)),sprintf('%02i',minute(time))];
%     filename = ['metoffice-c-band-rain-radar_uk_',datename,'_1km-composite.dat.gz'];
    datename = [sprintf('%4i',year(time)),sprintf('%02i',month(time)),sprintf('%02i',day(time)),sprintf('%02i',hour(time)),sprintf('%02i',minute(time))];
    filename = [datename,'.dat.gz'];

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
    radar(radar<0.05) = NaN;
    lograd = log(radar);
    m = nanmean(nanmean(lograd));
    Z = lograd - m;

    % Trim at a multiple of 2 size
    x = ((param(3)+param(5)/2):param(5):(param(3)+(param(7)*param(5)-param(5)/2)));
    y = ((param(2)-param(4)/2):-param(4):(param(2)-(param(6)*param(4))+param(4)/2));   y = fliplr(y);
    res = param(4)/1000;    % in km
    eastingmax = easting + side*param(5);
    northingmin = northing - side*param(4);
    indx = (x>easting&x<eastingmax);
    indy = (y>northingmin&y<northing);
    R = Z(indy,indx);   R = flipud(R);
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
    Ps = NaN(sz,sz);    WLs = NaN(sz,sz);
    Pn = NaN(sz,sz);    WLn = NaN(sz,sz);
    Ps(P>ns) = P(P>ns); WLs(P>ns) = WL(P>ns);
    Pn(P<ns) = P(P<ns); WLn(P<ns) = WL(P<ns);
    [Pavg, rad] = nanradialavg(P, (sz/2));
    [Psavg, ~] = nanradialavg(Ps, (sz/2));
    [Pnavg, ~] = nanradialavg(Pn, (sz/2));
    Pavg(1) = Pavg(2); Psavg(1) = Psavg(2); Pnavg(1) = Pnavg(2);
    
%     % Plotting
%     xh = [1 max(WL(:))]; yh = [ns ns];
%     xv = [threshold threshold]; yv = [0.00001 max(P(:))];
%     fig=figure;
%     figure(fig)
%     scatter(WLs(:),Ps(:),'b','.')
%     set(fig, 'Position', [200 200 1400 700])
%     set(gca, 'XScale','log')
%     set(gca, 'YScale','log')
%     hold on
%     scatter(WLn(:),Pn(:),'r','.')
%     plot(WL(:),Plineexp(:),'k')
%     plot(xh,yh,'--k')
%     plot(xv,yv,'--k')
%     saveas (fig, [results_folder, 'logP.jpg'])
%     close(fig)

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
    
    NF = RF - SF;
    NF = fftshift(NF);
    N = ifft2(NF, 'symmetric');   
    N = exp(N).* mask;
    
    % Create noise component in the frequency domain
    NFk = zeros(sz,sz,number_ens);
    for en=1:number_ens
        WF = normrnd(0,1,sz,sz);
        for x=1:sz
            for y=1:sz
                if position(x,y)<=(sz/2)
                    if x>y
                        NFk(x,y,en) = WF(x,y) * (Pnavg(position(x,y)).^0.5);   
                    else
                        NFk(x,y,en) = WF(y,x) * (Pnavg(position(x,y)).^0.5);
                    end
                else
                    if x>y
                        NFk(x,y,en) = WF(x,y) * (Pnavg(sz/2).^0.5);   
                    else
                        NFk(x,y,en) = WF(y,x) * (Pnavg(sz/2).^0.5);
                    end
                end
            end
        end
    end

    % Generate the ensemble
    EF = zeros(sz,sz,number_ens);
    for en=1:number_ens
        EF(:,:,en) = SF(:,:) + NFk(:,:,en);
    end

    % Inverse the Fourier transform
    E = zeros(sz,sz,number_ens);
    for en=1:number_ens
        EF(:,:,en) = fftshift(EF(:,:,en));
        E(:,:,en) = ifft2(EF(:,:,en), 'symmetric');   
        E(:,:,en) = exp(E(:,:,en) + m).* mask;
    end
    E(E<0) = 0;
    
    
    % Print out the ensembles and the radar
    Rain =exp(R + m).*mask;
    save([results_folder,'Ensemble_',datename,'.mat'], 'E')
    
    
%     cmap = colormap('jet');
%     n = size(cmap, 1);
%     zerovalue = [1 1 1];
%     minimum = -1.5;
%     maximum = 2.2;
%     dmap = (maximum-minimum)/n;
%     for en = 1:number_ens
%         Eprint = E(:,:,en);
%         fig=figure;
%         figure(fig)
%         imagesc(log(Eprint), [0 maximum])
%         set(fig, 'Position', [200 200 800 800])
%         colormap([zerovalue;cmap]);
%         caxis([minimum-dmap maximum]);
%         colorbar
%         title([starting,' - Ensemble ',num2str(en)])
%         saveas (fig, [results_folder, 'ensemble',num2str(en),'_',datename,'.jpg'])
%         close(fig)
%     end
%     fig=figure;
%     figure(fig)
%     imagesc(log(Rain), [0 maximum])
%     set(fig, 'Position', [200 200 800 800])
%     colormap([zerovalue;cmap]);
%     caxis([minimum-dmap maximum]);
%     colorbar
%     title([starting, ' - Radar'])
%     saveas (fig, [results_folder, 'radar','_',datename,'.jpg'])
%     close(fig)


end



    
        
        
    

