
res_folder = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_NewGenerator\Test16\';
load('\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_NewGenerator\Test16\noise.mat')

% Set frequencies and wavelength
fx = -128:128; fx = (fx(2:end)-0.5)/128; 
fy = fx;
[FX, FY] = meshgrid (fx, fy);
F = (FX.^2 + FY.^2).^0.5; fm = max(F(:)); F=F./fm;

R = raddata(:,:,10);
RF = fft2(R);     RF = fftshift(RF);
P = (abs(RF)).^2;
Plog = log(P);
[Pm,~] = radialavg(Plog,128);
x=(1:128)/128;

E = zeros(256,256,3); E2 = E;
T = [4,5,6];

fig1=figure;
figure(fig1)
set(fig1, 'Position', [50 50 1500 400])
for i=1:3
E(:,:,i) = noise(:,:,10,T(i));
EF = fft2(E(:,:,i));     EF = fftshift(EF);
Pe = (abs(EF)).^2;
Pelog = log(Pe);
[Em,~] = radialavg(Pelog,128);
% subplot(3,4,i)
% scatter(F(:),Plog(:),'.','b')
% hold on
% scatter(F(:),Pelog(:),'.','r')
% title(['power radar-ens',num2str(T(i))])
% xlabel('Normalised Frequency')
% ylabel('log(power)')

% subplot(3,4,i+4)
% plot(x,Pm,'b')
% hold on
% plot(x,Em,'r')
% title(['avg power radar-ens',num2str(T(i))])
% xlabel('Normalised Frequency')
% ylabel('log(power)')

ratio = Em./Pm; ratio(isnan(ratio)==1)=1;
fitcoeff = polyfit(x,ratio,1);
% subplot(3,4,i+8)
% plot(x,ratio,'k')
% hold on
% plot(x,fit,'r')
% title(['ratio of avg power radar-ens',num2str(T(i))])
% xlabel('Normalised Frequency')
% ylabel('log(power)')
subplot(1,3,i)
EF2 = EF./(exp((F.*fitcoeff(1))+fitcoeff(2)));
Pe2 = (abs(EF2)).^2;
Pelog2 = log(Pe2);
scatter(F(:),Pelog(:),'.','m')
hold on
scatter(F(:),Pelog2(:),'.','r')
scatter(F(:),Plog(:),'.','b')
legend('Ens power not corrected','Ens power corrected','Radar power')

EF2 = fftshift(EF2);
E2(:,:,i) = ifft2(EF2, 'symmetric');
end

fig2=figure;
figure(fig2)
set(fig2, 'Position', [50 50 1500 600])
cmap = colormap('jet');
n = size(cmap, 1);
zerovalue = [1 1 1];
minimum = -1.5;
maximum = 2.2;
for j=1:3
subplot(2,3,j)
E2(E2<0.0001)=0.0001;
imagesc(log(E2(:,:,j)), [0 maximum])
dmap = (maximum-minimum)/n;
colormap([zerovalue;cmap]);
caxis([minimum-dmap maximum]);
colorbar
title(['corrected ensemble',num2str(j)])
subplot(2,3,j+3)
E(E==0)=0.0001;
imagesc(log(E(:,:,j)), [0 maximum])
dmap = (maximum-minimum)/n;
colormap([zerovalue;cmap]);
caxis([minimum-dmap maximum]);
colorbar
title(['original ensemble',num2str(j)])
end




% saveas(fig,[res_folder,'powertest1.fig'])
% print(fig,[res_folder,'powertest8'],'-dpng','-r300')
