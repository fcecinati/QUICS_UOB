clear all
close all

load('\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_NewGenerator\Test9\noisenew.mat')
load('F:\Francesca\Test_German\Test_20\noise.mat')

for i=1:25
    ng = noise(:,:,i,:);
    nn = noisenew(:,:,i,:);
    
    mg(i) = mean(ng(:));
    mn(i) = mean(nn(:));
    
    stg(i) = std(ng(:));
    stn(i) = std(nn(:));
    
    mxg(i) = max(ng(:));
    mxn(i) = max(nn(:));
end


fig1 = figure;
figure(fig1);
set(fig1,'Position',[100,50,1500,250])
hold on
plot(mg,'b');
plot(mn,'r');
title('mean comparison - Germann blue, New red')
xlabel('time')
ylabel('mean')    
saveas(fig1,'\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Stat_comparison\NoiseComparison\mean.jpg')

fig2 = figure;
figure(fig2);
set(fig2,'Position',[100,400,1500,250])
hold on
plot(stg,'b');
plot(stn,'r');
title('std comparison - Germann blue, New red')
xlabel('time')
ylabel('std')    
saveas(fig2,'\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Stat_comparison\NoiseComparison\std.jpg')

fig3 = figure;
figure(fig3);
set(fig3,'Position',[100,750,1500,250])
hold on
plot(mxg,'b');
plot(mxn,'r');
title('max comparison - Germann blue, New red')
xlabel('time')
ylabel('max')    
saveas(fig3,'\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Stat_comparison\NoiseComparison\max.jpg')

