function []=RadarPlot(raddata,ens, res_folder,ts)
close all

Eprint = raddata(:,:,ts);
name = 'Rad';
fig=figure;
figure(fig)
set(gca,'Position',[0 0 1 1])
set(fig, 'Position', [100 100 1000 900])
subplot(3,3,1)


cmap = colormap('jet');
n = size(cmap, 1);
zerovalue = [1 1 1];
minimum = -1.5;
maximum = 2.2;
imagesc(log(Eprint), [0 maximum])
dmap = (maximum-minimum)/n;
colormap([zerovalue;cmap]);
caxis([minimum-dmap maximum]);
% colorbar
title(['Radar - mean=',num2str(mean(Eprint(:))),' - std=',num2str(std(Eprint(:)))])


for i=2:9

    Eprint = ens(:,:,ts,i+10);
    name = ['Ens',num2str(i),'- mean=',num2str(mean(Eprint(:))),' - std=',num2str(std(Eprint(:)))];
    subplot(3,3,i)
    cmap = colormap('jet');
    n = size(cmap, 1);
    zerovalue = [1 1 1];
    minimum = -1.5;
    maximum = 2.2;
    dmap = (maximum-minimum)/n;
    imagesc(log(Eprint), [0 maximum])
    colormap([zerovalue;cmap]);
    caxis([minimum-dmap maximum]);
%     colorbar
    title(name)
    
end

    saveas(fig,[res_folder, 'ensembles', '.fig'])
    print(fig,[res_folder, 'ensembles'],'-dpng','-r300')