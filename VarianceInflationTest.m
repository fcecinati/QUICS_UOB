function VarianceInflationTest(raddata, ens, newens, res_folder)
sz = size(ens);
t_steps = sz(3);
n_ens = sz(4);

for t=1:t_steps
    R = squeeze(raddata(:,:,t));
    rad(t) = var(R(:));
    for ensem = 1:n_ens
        E = squeeze(ens(:,:,t,ensem));
        en(t,ensem) = var(E(:));
        EN = squeeze(newens(:,:,t,ensem));
        nen(t,ensem) = var(EN(:));
    end
end
en = mean(en,2); en = squeeze(en);
nen = mean(nen,2); nen = squeeze(nen);

fig=figure;
figure(fig)
plot(rad,'b')
hold on
plot(en,'r')
plot(nen,'k')
set(fig,'Position', [1000 80 600 300])
title ('variance inflation test')
set(gcf,'PaperPositionMode','auto')
legend('radar','ensemble mean','corrected ensemble mean')
xlabel('time (h)')
ylabel('variance')
saveas(fig,[res_folder,'varianceinflation.fig'])
print(fig,[res_folder,'varianceinflation'],'-dpng','-r300')


for t=1:t_steps
    R = squeeze(raddata(:,:,t));
    rad(t) = mean(R(:));
    for ensem = 1:n_ens
    E = squeeze(ens(:,:,t,ensem));
    en(t,ensem) = mean(E(:));
    EN = squeeze(newens(:,:,t,ensem));
    nen(t,ensem) = mean(EN(:));
    end
end
en = mean(en,2); en = squeeze(en);
nen = mean(nen,2); nen = squeeze(nen);

fig=figure;
figure(fig)
plot(rad,'b')
hold on
plot(en,'r')
plot(nen,'k')
set(fig,'Position', [1000 80 600 300])
title ('mean inflation test')
set(gcf,'PaperPositionMode','auto')
legend('radar','ensemble mean','corrected ensemble mean')
xlabel('time (h)')
ylabel('mean')
saveas(fig,[res_folder,'meaninflation.fig'])
print(fig,[res_folder,'meaninflation'],'-dpng','-r300')