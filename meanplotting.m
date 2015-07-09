me = mean(ens,1);   me = mean(me,2); me = squeeze(me);
mr = mean(raddata,1);   mr = mean(mr,2); mr = squeeze(mr);

f1 = figure;
figure(f1)
set(f1,'Position',[100,100,1500,800])
hold on
for i=1:100
    plot(me(:,i),'r')
end
plot(mr,'b')
title('Comparison of the spatial mean of the 100 ensembles and the original radar')
xlabel('Time (h)')
ylabel('Mean (mm/h)')
saveas(f1,'\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Stat_comparison\MeanNew3.jpg')


