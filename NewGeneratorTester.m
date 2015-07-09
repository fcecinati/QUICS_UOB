% function [] = NewGeneratorTester(b) 
% a is the filter order, b is the cut-off frequency.
res_folder = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_NewGenerator\Test15\';
close all

mu = 0;
sigma = 5.6855;
a=6;
b = 3000;

R = normrnd(mu,sigma,256,256);

base = 1:a;
fun = exp(-(base./b));
f1 = [fliplr(fun),(1),fun];
H0 = sum(f1);
f1=f1/H0;
f2 = ftrans2(f1);
Rnew = filter2(f2,R);
newstd = std(Rnew(:));
Rnew = Rnew*sigma/newstd;

fig0 = figure;
figure(fig0)
set(fig0,'Position',[10 80 600 300])
plot(f1)
title('new filter')


fig1=figure;
figure(fig1)
set(fig1,'Position', [10 500 1200 500])
freqz(f1)
title('Filter')
% saveas(fig1,[res_folder,'fir_mag-phase_', num2str(N),'_',num2str(b),'.jpg'])
% 
% fig2=figure;
% figure(fig2)
% set(fig2,'Position', [650 500 600 500])
% freqz2(f2)
% title('Filter')
% saveas(fig2,[res_folder,'fir_2D_', num2str(N),'_',num2str(b),'.jpg'])


fig6=figure;
figure(fig6)
set(fig6,'Position', [1300 700 400 320])
imagesc(R)
title('Random Noise')
set(gcf,'PaperPositionMode','auto')
colorbar
% saveas(fig6,[res_folder,'randoimnoise', test, '.fig'])
% print(fig6,[res_folder,'randomnoise', test, ''],'-dpng','-r300')
% hold off
% 

% 
[variog,x,y] = Variogram(Rnew);
variog=variog(1:256);
h = (1:size(variog,2));

myfunction = fittype('c.*(1-exp(-(h./a)))+b','independent',{'h'},'coefficients',{'c','b','a'});
myoption = fitoptions('Method','NonlinearLeastSquares');
myfit = fit(h',variog',myfunction,myoption);

fig3=figure;
figure(fig3)
set(fig3,'Position', [1300 200 400 320])
imagesc(Rnew)
title('Generated Noise')
set(gcf,'PaperPositionMode','auto')
colorbar
hold on
scatter(x,y,'p')

fig4 = figure;
figure(fig4)
hold on
set(fig4,'Position', [650 80 600 320])
plot(myfit,h,variog)
title ('semivariogram')
% % plot(ref_v,'k')
% set(gcf,'PaperPositionMode','auto')
% legend('calculated semivariogram','fitting exponential curve','target exponential curve','Location','southeast')
% xlabel('distance (km)')
% ylabel('semivariogram')
% saveas(fig4,[res_folder,'semiv', test, '.fig'])
% print(fig4,[res_folder,'semiv', test, ''],'-dpng','-r300')
% hold off




