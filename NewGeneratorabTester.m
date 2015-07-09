% Tests how the semivariogram changes with N/b/sigma variation
res_folder = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_NewGenerator\Test15\';
close all
test = '1';
%setup
mu = 0;
F = 0.4715;
sill = 78.7019;
range = 57.7507;
h = (1:128)';
ref_v = sill.*(1-exp(-(h./range).^F));

%testing intervals
Nvec = 30:10:150;
bvec = 10:10:100;
sigmavec = 4:0.4:6;

% % test N
% sigma = 5.6;
% b = 70;
% 
% fig0 = figure;
% figure(fig0)
% linecolours = colormap('jet');
% set(fig0,'Position', [500 700 700 250])
% title ('semivariogram changing N')
% hold on
% R = normrnd(mu,sigma,256,256);
% for N=Nvec
%     base = 1:N;
%     fun = (1/(sigma^2)).*exp(-(base./b).^F);
%     f1 = [fliplr(fun),1/(sigma^2),fun];
%     f2 = ftrans2(f1);
%     Rnew = filter2(f2,R);
%     newvar = var(Rnew(:));
%     Rnew = Rnew*sigma/(newvar^0.5);
%     [hvar, vvar] = Variogram(Rnew);
%     variog = (hvar + vvar)/2;
%     h = (1:128)';
%     plot(h,variog,'color',linecolours(ceil((N-min(Nvec))*63/(max(Nvec)-min(Nvec)))+1,:))   
%     
% end
% plot(ref_v,'k')
% set(gcf,'PaperPositionMode','auto')
% xlabel('distance(km)')
% ylabel('semivariogram')
% colorbar
% caxis([min(Nvec) max(Nvec)])
% hold off
% saveas(fig0,[res_folder,'semivariogram_N.fig'])
% print(fig0,[res_folder,'sv_N'],'-dpng','-r300')

% test b
variogvec = zeros(10,128);
sigma = 5.6;
N = 128;

fig0 = figure;
figure(fig0)
set(fig0,'Position', [500 400 700 250])
title ('semivariogram changing b')
hold on
R = normrnd(mu,sigma,256,256);
for b=bvec
    base = 1:N;
    fun = (1/(sigma^2)).*exp(-(base./b).^F);
    f1 = [fliplr(fun),1/(sigma^2),fun];
    f2 = ftrans2(f1);
    Rnew = filter2(f2,R);
    newvar = var(Rnew(:));
    Rnew = Rnew*sigma/(newvar^0.5);
    [hvar, vvar] = Variogram(Rnew);
    variog = (hvar + vvar)/2;
    h = (1:128)';
    plot(h,variog,'color',linecolours(ceil((b-min(bvec))*63/(max(bvec)-min(bvec)))+1,:))        
end
plot(ref_v,'k')
set(gcf,'PaperPositionMode','auto')
colorbar
caxis([min(bvec) max(bvec)])
xlabel('distance(km)')
ylabel('semivariogram')
hold off
saveas(fig0,[res_folder,'semivariogram_b.fig'])
print(fig0,[res_folder,'sv_b'],'-dpng','-r300')

%test sigma
N = 128;
b = 90;

fig0 = figure;
figure(fig0)
set(fig0,'Position', [500 80 700 250])
title ('semivariogram changing sigma')
hold on
for sigma=sigmavec
    base = 1:N;
    fun = (1/(sigma^2)).*exp(-(base./b).^F);
    f1 = [fliplr(fun),1/(sigma^2),fun];
    f2 = ftrans2(f1);
    for k=1:10
        R = normrnd(mu,sigma,256,256);
        Rnew = filter2(f2,R);
        newvar = var(Rnew(:));
        Rnew = Rnew*sigma/(newvar^0.5);
        [hvar, vvar] = Variogram(Rnew);
        variogvec(k,:) = (hvar + vvar)/2;
    end
    variog = mean(variogvec);
    h = (1:128)';
    plot(h,variog,'color',linecolours(ceil((sigma-min(sigmavec))*63/(max(sigmavec)-min(sigmavec)))+1,:))        
end
plot(ref_v,'k')
set(gcf,'PaperPositionMode','auto')
colorbar
caxis([min(sigmavec) max(sigmavec)])
xlabel('distance(km)')
ylabel('semivariogram')
hold off
saveas(fig0,[res_folder,'semivariogram_sigma.fig'])
print(fig0,[res_folder,'sv_sigma'],'-dpng','-r300')

