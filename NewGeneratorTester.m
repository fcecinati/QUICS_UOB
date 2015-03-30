function [] = NewGeneratorTester(a,b) 
% a is the filter order, b is the cut-off frequency.
close all

mu = 0;
sigma = 5.6855;
sigma1 = (sigma)/(1.09641257*(a^(-0.94003379)));
k = (mu^2)/(sigma1^2);
th = (sigma1^2)/mu;

R = random('gam',k,th,256,256);

% R=normrnd(mu,sigma1,256,256);
F = 0.566129633372927;

f1=fir1(a,b,'low');
f2 = ftrans2(f1);
Rnew = filter2(f2,R);

fig1=figure;
figure(fig1)
set(fig1,'Position', [10 500 600 500])
freqz(f1)
title('Filter')

fig2=figure;
figure(fig2)
set(fig2,'Position', [650 500 600 500])
freqz2(f2)
title('Filter')

fig3=figure;
figure(fig3)
set(fig3,'Position', [1300 500 600 500])
imagesc(Rnew)
title('Generated Noise')
colorbar

[hvar, vvar] = Variogram(Rnew);

variog = (hvar + vvar)/2;
h = (1:128)';

myfunction = fittype('c.*(1-exp(-3.*(h.^F)./a.^F))','independent',{'h'},'coefficients',{'c','a'},'problem','F');
myoption = fitoptions('Method','NonlinearLeastSquares','StartPoint',[0.0007 40]);
myfit = fit(h,variog,myfunction,myoption,'problem',F);

fig4 = figure;
figure(fig4)
hold on
set(fig4,'Position', [250 80 700 300])
plot(myfit,h,variog)
title ('semivariogram')
hold off



