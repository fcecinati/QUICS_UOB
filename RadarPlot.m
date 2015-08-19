function []=RadarPlot(raddata,ens, name)

    close all

    Eprint = raddata(:,:);
    fig=figure;
    figure(fig)
    set(gca,'Position',[0 0 1 1])
    set(fig, 'Position', [100 100 500 450])
    colorbar
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

        Eprint = ens(:,:,i);
        nm = ['Ens',num2str(i-1),'- mean=',num2str(mean(Eprint(:))),' - std=',num2str(std(Eprint(:)))];
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
        title(nm)

    end
    hp4 = get(subplot(3,3,9),'Position');
    colorbar('Position',[hp4(1)+hp4(3)+0.02  hp4(2)  0.01  hp4(2)+hp4(3)*3.2]);

       

    saveas(fig,[name, '.fig'])
    print(fig,name,'-dpng','-r300')
