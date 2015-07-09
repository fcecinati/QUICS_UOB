close all
n=12;

r=raddata(:,:,n);
e1 = ens(:,:,n,1);
e2 = ens(:,:,n,2);
e3 = ens(:,:,n,3);
mx = max([max(r(:)),max(e1(:)),max(e2(:)),max(e3(:))]);

mr = mean(r(:));
vr = var(r(:));
mmr = max(r(:));
f1 = figure;
figure(f1)
set(f1,'Position',[100 600 700 400])
imagesc(r)
caxis([0 mx])
colorbar
title(['radar - mean = ',num2str(mr,4),'; var = ', num2str(vr,4),'; max = ',num2str(mmr,4)])
%saveas(f1,['\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_NewGenerator\Test6\Radar_t',num2str(n),'.jpg'])

me1 = mean(e1(:));
ve1 = var(e1(:));
mme1 = max(e1(:));
f2 = figure;
figure(f2)
set(f2,'Position',[1000 600 700 400])
imagesc(e1)
caxis([0 mx])
colorbar
title(['ens1 - mean = ',num2str(me1,4),'; var = ', num2str(ve1,4),'; max = ',num2str(mme1,4)])
%saveas(f2,['\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_NewGenerator\Test6\Ens1_t',num2str(n),'.jpg'])

me2 = mean(e2(:));
ve2 = var(e2(:));
mme2 = max(e2(:));
f3 = figure;
figure(f3)
set(f3,'Position',[1000 80 700 400])
imagesc(e2)
caxis([0 mx])
colorbar
title(['ens2 - mean = ',num2str(me2,4),'; var = ', num2str(ve2,4),'; max = ',num2str(mme2,4)])
%saveas(f3,['\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_NewGenerator\Test6\Ens2_t',num2str(n),'.jpg'])

me3 = mean(e3(:));
ve3 = var(e3(:));
mme3 = max(e3(:));
f4 = figure;
figure(f4)
set(f4,'Position',[100 80 700 400])
imagesc(e3)
caxis([0 mx])
colorbar
title(['ens3 - mean = ',num2str(me3,4),'; var = ', num2str(ve3,4),'; max = ',num2str(mme3,4)])
%saveas(f4,['\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_NewGenerator\Test6\Ens3_t',num2str(n),'.jpg'])