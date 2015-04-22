function [] = NewGeneratorTester(b) 
% a is the filter order, b is the cut-off frequency.
res_folder = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_NewGenerator\Test3\';
close all
tic
mu = 0;
sigma = 5.6855;

R = normrnd(mu,sigma,256,256);
F = 0.566129633372927;
toc
a = 80;
c = sigma^2;

% f1 = fir1(a,b,'low');
base = 1:a;
fun = exp(-3.*(base./b).^F);
f1 = [fliplr(fun),1,fun];
f2 = ftrans2(f1);
Rnew = filter2(f2,R);
toc

newstd = std(Rnew(:));
stdratio = newstd/sigma;

Rnew = Rnew/stdratio;
newstd = std(Rnew(:));
stdratio = newstd/sigma;
m = mean(Rnew(:));

fig0 = figure;
figure(fig0)
set(fig0,'Position',[100 80 700 300])
plot(f1)
title('new filter')
saveas(fig0,[res_folder,'new_1D_', num2str(a),'_',num2str(b),'.jpg'])

fig1=figure;
figure(fig1)
set(fig1,'Position', [10 500 600 500])
freqz(f1)
title('Filter')
saveas(fig1,[res_folder,'new_mag-phase_', num2str(a),'_',num2str(b),'.jpg'])

fig2=figure;
figure(fig2)
set(fig2,'Position', [650 500 600 500])
freqz2(f2)
title('Filter')
saveas(fig2,[res_folder,'new_2D_', num2str(a),'_',num2str(b),'.jpg'])

fig3=figure;
figure(fig3)
set(fig3,'Position', [1300 500 600 500])
imagesc(Rnew)
title(['Generated Noise; stdratio = ',num2str(stdratio, 4), '; mean = ',num2str(m,4)])
colorbar
saveas(fig3,[res_folder,'new_noise_', num2str(a),'_',num2str(b),'.jpg'])

[hvar, vvar] = Variogram(Rnew);

variog = (hvar + vvar)/2;
h = (1:128)';

myfunction = fittype('c.*(1-exp(-3.*(h.^F)./a.^F))','independent',{'h'},'coefficients',{'c','a'},'problem','F');
myoption = fitoptions('Method','NonlinearLeastSquares','StartPoint',[0.0007 40]);
myfit = fit(h,variog,myfunction,myoption,'problem',F);

fig4 = figure;
figure(fig4)
hold on
set(fig4,'Position', [1000 80 700 300])
plot(myfit,h,variog)
title ('semivariogram')
hold off
saveas(fig4,[res_folder,'new_semivariogram_', num2str(a),'_',num2str(b),'.jpg'])



