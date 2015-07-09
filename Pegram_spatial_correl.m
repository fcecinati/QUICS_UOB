%% spatial correlation analysis

cd('\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_Pegram\Test_19')
dates = '2008-05-04 02:00:00';
daten = datenum(dates);
subnum = 50;
ensnum = 100;

% File to work on
datename = [sprintf('%4i',year(time)),sprintf('%02i',month(time)),sprintf('%02i',day(time)),sprintf('%02i',hour(time)),sprintf('%02i',minute(time))];
ensemblefilename = ['Ensemble_', datename, '.mat'];
load(ensemblefilename)
sz = 256;

for i=1:subnum
    subX(i) = randi(sz);    subY(i) = randi(sz);
end

for ens = 1:ensnum
    Nk(:,:,ens) = ifft2(fftshift(NFk(:,:,ens)), 'symmetric');   
    Nk(:,:,ens) = exp(Nk(:,:,ens));
    for i=1:subnum
        subN(i,ens) = squeeze(Nk(subX(i),subY(i),ens));        
    end
end

% Mean without weights
mean_e = nanmean(subN,2);
sigma = nanstd(subN,1,2);
diff=zeros(subnum,ensnum);
for ens=1:ensnum
    diff(:,ens) = subN(:,ens) - mean_e;
end

%calculation of correlation coefficients and distances
correl = zeros(subnum,subnum);
distance = zeros(subnum,subnum);

for x1=1:subnum
    for x2=1:subnum
        correl(x1,x2) = nanmean(diff(x1,:).*diff(x2,:))/(sigma(x1)*sigma(x2));
        distance(x1,x2) = sqrt((subX(x1)-subX(x2))^2 + (subY(x1)-subY(x2))^2);
    end
end

% Reshaping and NaN elimination for function estimation and plotting
d = distance(:);
c = correl(:);

        % % Fit on decorrelation function
        % myfunction = fittype('exp(-((d/R0)^F))','independent',{'d'},'coefficients',{'R0','F'});
        % myoption = fitoptions('Method','NonlinearLeastSquares','StartPoint',[0.5 0.9]);
        % myfit = fit(d,c,myfunction,myoption);

% Plotting
fig = figure;
figure(fig)
set(fig, 'Position', [200 200 800 800])
scatter(d,c)
title(['spatial correlation of the noise field in the ensembles ',datename]);
xlabel('meters x 10^5')
name = ['noise_spatial_corr_',datename];
saveas(fig,name,'jpg')
close(fig)


