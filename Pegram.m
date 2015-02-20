%This script will implement the (Pegram et al, 2011) method for ensemble generation

clear all
close all
clc
rng('shuffle')

%% Parameters to define:

% Path to the folder where the radar data is 
data_path = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\data_north_england\09-2008';
% file to work on
filename = 'metoffice-c-band-rain-radar_uk_200809051700_1km-composite.dat.gz';
% Temporary folder
tmpfolder = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Temp\'; 
% Script folder
script_folder = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\QUICS_UOB';
% Results folder
results_folder = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_Pegram\Test_6\';
% How many ensembles do you want to generate
number_ens = 5;
% Resolution in km
res = 1;

%% Data loading and pre-processing

if ~exist(results_folder, 'dir')
  mkdir(results_folder);
end

cd(data_path)
tt=length(filename);
if (strcmp('.gz', filename(tt-2:tt)) == 1)   % if the input file is compressed (i.e. with extension '.gz')
     newname = gunzip(filename, tmpfolder); % decompress file in a temporary folder and output new file name
     filename = newname{1,1};    % copy new name into the nimrodfilename
end
cd(script_folder)
[radar, date, param] = readnimrod(filename);

radar(radar<0.5) = NaN;
m = nanmean(nanmean(radar));
lograd = log(radar);
Z = lograd - m;

% trim at a multiple of 2 size
x = (param(3):1000:param(3)+(param(7)*1000 -1000));
y = (param(2):-1000:(param(2)-(param(6)*1000)+1000));
indx = (x>350000&x<478000);
indy = (y>768000&y<896000);
R = Z(indy,indx);

sz = size(R); sz = sz(1);

% Creates a rain/norain mask
mask = zeros(size(R));
mask(isnan(R)==0)=1;

R(isnan(R)==1) = 0;

%% Pegram method

threshold = 8;      %km wavelenght

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
phase = angle(RF);
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

% Plotting
xh = [1 max(WL(:))]; yh = [ns ns];
xv = [threshold threshold]; yv = [0.1 max(P(:))];
fig=figure;
figure(fig)
scatter(WL(:),P(:),'k','.')
set(fig, 'Position', [200 200 1400 700])
set(gca, 'XScale','log')
set(gca, 'YScale','log')
hold on
scatter(WLcenter(:),Pcenter(:),'r','.')
plot(WL(:),Plineexp(:))
plot(xh,yh,'--k')
plot(xv,yv,'--k')
saveas (fig, [results_folder, 'logP.jpg'])
close(fig)

% Separation of Signal and Noise

Ps = NaN(sz,sz);
Pn = NaN(sz,sz);
Ps(P>ns) = P(P>ns);
Pn(P<ns) = P(P<ns);
[Pavg, rad] = nanradialavg(P, (sz/2));
[Psavg, ~] = nanradialavg(Ps, (sz/2));
[Pnavg, ~] = nanradialavg(Pn, (sz/2));
Pavg(1) = Pavg(2); Psavg(1) = Psavg(2); Pnavg(1) = Pnavg(2);

fr = rad./2;
wl = 1./fr;
fig=figure;
figure(fig)
plot(wl,Pavg,'k')
set(fig, 'Position', [200 200 1400 700])
set(gca, 'XScale','log')
set(gca, 'YScale','log')
hold on
plot(wl,Psavg,'b')
plot(wl,Pnavg,'r')
saveas (fig, [results_folder, 'avgP.jpg'])
close(fig)

position = round(F.*sz);
SF = zeros(sz,sz);

for x=1:sz
    for y=1:sz
        if position(x,y)<=(sz/2)
            SF(x,y) = RF(x,y) * (Psavg(position(x,y))/Pavg(position(x,y)));
        end
    end
end

% Generate the random noise for the ensembles
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

% Print out the ensembles and the radar
Rain =exp(R + m).*mask;

maximum = max(max(Rain(:)), max(E(:)));
for en = 1:number_ens
    fig=figure;
    figure(fig)
    imagesc(E(:,:,en), [0 maximum])
    set(fig, 'Position', [200 200 800 800])
    colorbar
    saveas (fig, [results_folder, 'ensemble',num2str(en),'.jpg'])
    close(fig)
end


fig=figure;
figure(fig)
imagesc(Rain, [0 maximum])
set(fig, 'Position', [200 200 800 800])
colorbar
saveas (fig, [results_folder, 'radar.jpg'])
close(fig)
    

%% EMD and Fourier transform

% IMS = bemd(R);
% N = IMS(:,:,1);
% S = IMS(:,:,4)+IMS(:,:,2)+IMS(:,:,3);
% 
% N = N.*mask;
% S = S.*mask;
% 
% NF = fft2(N);       NF = fftshift(NF);
% SF = fft2(S);       SF = fftshift(SF);
% 
% PN = (abs(NF)).^2;
% PS = (abs(SF)).^2;
% 
% [radavg, ~] = radialavg(P, (sz/2));
% [noiseavg, ~] = radialavg(PN, (sz/2));
% [signalavg, ~] = radialavg(PS, (sz/2));
% 
% f = (1:(sz/2))./(sz/2); % km^-1
% wl = 2*pi./f;           % km
% 
% fig=figure;
% figure(fig)
% loglog(wl, radavg,'k')
% hold on
% loglog(wl, noiseavg, 'r')
% loglog(wl, signalavg, 'b')
% set(fig, 'Position', [0 0 1400 1000])
% legend('total field','noise','signal')
% saveas(fig,[results_folder,'components.jpg'])
% 
 save([results_folder, 'Test16'])
