function []=PlotPower(plotfilename)

load(plotfilename)
datename = [sprintf('%4i',year(time)),sprintf('%02i',month(time)),sprintf('%02i',day(time)),sprintf('%02i',hour(time)),sprintf('%02i',minute(time))];
% Plotting puntual power
        xh = [1 max(WL(:))]; yh = [ns ns];
        xv = [threshold threshold]; yv = [0.1 max(P(:))];
        fig=figure;
        figure(fig)
        scatter(WL(:),P(:),'k','.')
        set(fig, 'Position', [200 200 1400 700])
        set(gca, 'XScale','log')
        set(gca, 'YScale','log')
        hold on
        scatter(WLcenter(:),Pcenter(:),'r','.')
        plot(WL(:),Plineexp(:))
        plot(xh,yh,'--k')
        plot(xv,yv,'--k')
        saveas (fig, [results_folder, 'PuntualP_',datename,'.jpg'])
        close(fig)
% Plotting radially averaged power
        fr = rad./2;
        wl = 1./fr;
        fig=figure;
        figure(fig)
        plot(wl,Pavg,'k')
        set(fig, 'Position', [200 200 1400 700])
        set(gca, 'XScale','log')
        set(gca, 'YScale','log')
        hold on
        plot(wl,Psavg,'b')
        plot(wl,Pnavg,'r')
        saveas (fig, [results_folder, 'AvgP_',datename,'.jpg'])
        close(fig)
        