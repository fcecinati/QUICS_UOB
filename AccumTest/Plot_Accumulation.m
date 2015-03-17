cd('\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Comparison_Test_20\Accumulation_tests')
radius = 2;
sqkm = (radius*2+1)^2;
load('Gauge.mat')
Gername = ['Ger_acc_',num2str(sqkm),'km2.mat'];
Pegname = ['Peg_acc_',num2str(sqkm),'km2.mat'];
load(Gername)
load(Pegname)

germax1 = max(accger1');
germin1 = min(accger1');
pegmax1 = max(accpeg1');
pegmin1 = min(accpeg1');
germax2 = max(accger2');
germin2 = min(accger2');
pegmax2 = max(accpeg2');
pegmin2 = min(accpeg2');

fig3=figure;
figure(fig3)
set(fig3, 'Position', [100 100 1200 800])
hold on
a1 = area(germax1);
set(a1,'LineStyle','none'); set(a1,'FaceColor',[0.9 0.9 0.5]);
a2 = area(germin1);
set(a2,'LineStyle','none'); set(a2,'FaceColor',[1 1 1]);
plot(accrad1, 'r','Linewidth',3)
plot(accG1, 'k','Linewidth',3)
legend('Germann','','Radar','Gauge','Location','northwest');
title('Germann ensemble and rain gauge1 accumulation 24h')
saveas(fig3,['Ger_Accum1_',num2str(sqkm),'km2.jpg'])
close(fig3)

fig4=figure;
figure(fig4)
set(fig4, 'Position', [100 100 1200 800])
hold on
a1 = area(pegmax1);
set(a1,'LineStyle','none'); set(a1,'FaceColor',[0.5 0.9 0.5]);
a2 = area(pegmin1);
set(a2,'LineStyle','none'); set(a2,'FaceColor',[1 1 1]);
plot(accrad1, 'r','Linewidth',3)
plot(accG1, 'k','Linewidth',3)
legend('Pegram','','Radar','Gauge','Location','northwest');
title('Comparison of ensemble and rain gauge1 accumulation 24h')
saveas(fig4,['Peg_Accum1_',num2str(sqkm),'km2.jpg'])
close(fig4)

fig3=figure;
figure(fig3)
set(fig3, 'Position', [100 100 1200 800])
hold on
a1 = area(germax2);
set(a1,'LineStyle','none'); set(a1,'FaceColor',[0.9 0.9 0.5]);
a2 = area(germin2);
set(a2,'LineStyle','none'); set(a2,'FaceColor',[1 1 1]);
plot(accrad2, 'r','Linewidth',3)
plot(accG2, 'k','Linewidth',3)
legend('Germann','','Radar','Gauge','Location','northwest');
title('Germann ensemble and rain gauge2 accumulation 24h')
saveas(fig3,['Ger_Accum2_',num2str(sqkm),'km2.jpg'])
close(fig3)

fig4=figure;
figure(fig4)
set(fig4, 'Position', [100 100 1200 800])
hold on
a1 = area(pegmax2);
set(a1,'LineStyle','none'); set(a1,'FaceColor',[0.5 0.9 0.5]);
a2 = area(pegmin2);
set(a2,'LineStyle','none'); set(a2,'FaceColor',[1 1 1]);
plot(accrad2, 'r','Linewidth',3)
plot(accG2, 'k','Linewidth',3)
legend('Pegram','','Radar','Gauge','Location','northwest');
title('Comparison of ensemble and rain gauge2 accumulation 24h')
saveas(fig4,['Peg_Accum2_',num2str(sqkm),'km2.jpg'])
close(fig4)