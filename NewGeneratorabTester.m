% Tests how the semivariogram changes with a/b variation
res_folder = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_NewGenerator\Test4\';
close all
tic
mu = 0;
sigma = 5.6855;
R = normrnd(mu,sigma,256,256);
F = 0.566129633372927;
c=sigma^2;
oldvar = var(R(:));
toc

% test a for 3 different b values
b3 = [30,50,80];

for b=b3
    linecolours = colormap('jet');
    fig0 = figure;
    figure(fig0)
    set(fig0,'Position', [1000 80 700 300])
    title ('semivariogram')
    hold on
    for a=20:20:200
        base = 1:a;
        fun = c -(c.*(1-exp(-3.*(base./b).^F)));
        f1 = [fliplr(fun),c,fun];
        f2 = ftrans2(f1);
        Rnew = filter2(f2,R)./a;
        newvar = var(Rnew(:));
        varratio = newvar/oldvar;
        Rnew = Rnew/(varratio^0.5);
        [hvar, vvar] = Variogram(Rnew);
        variog = (hvar + vvar)/2;
        h = (1:128)';
        plot(h,variog,'color',linecolours(a/4,:))        
    end
    hold off
    saveas(fig0,[res_folder,'semivariogram_b_', num2str(b),'_various_a.jpg'])
end


toc

a3 = [50,100,150];

for a=a3
    linecolours = colormap('jet');
    fig0 = figure;
    figure(fig0)
    set(fig0,'Position', [1000 80 700 300])
    title ('semivariogram')
    hold on
    for b=10:10:100
        base = 1:a;
        fun = c -(c.*(1-exp(-3.*(base./b).^F)));
        f1 = [fliplr(fun),c,fun];
        f2 = ftrans2(f1);
        Rnew = filter2(f2,R)./a;
        newvar = var(Rnew(:));
        varratio = newvar/oldvar;
        Rnew = Rnew/(varratio^0.5);
        [hvar, vvar] = Variogram(Rnew);
        variog = (hvar + vvar)/2;
        h = (1:128)';
        plot(h,variog,'color',linecolours(b/2,:))        
    end
    hold off
    saveas(fig0,[res_folder,'semivariogram_a_', num2str(a),'_various_b.jpg'])
end