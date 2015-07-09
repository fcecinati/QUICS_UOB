rad = 1:25;
en = 1:25;
res_folder = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_NewGenerator\Test16\';

for t=1:25
    R = squeeze(raddata(:,:,t));
    E = squeeze(ens(:,:,t,:));
    E = mean(E,3);
    rad(t) = mean(R(:));
    en(t) = mean(E(:));
end

fig=figure;
figure(fig)
plot(rad,'b')
hold on
plot(en,'r')
set(fig,'Position', [1000 80 600 300])
title ('Mean inflation for the new method')
set(gcf,'PaperPositionMode','auto')
legend('radar','ensemble mean')
xlabel('time (h)')
ylabel('mean')
saveas(fig,[res_folder,'meaninflation.fig'])
print(fig,[res_folder,'meaninflation'],'-dpng','-r300')