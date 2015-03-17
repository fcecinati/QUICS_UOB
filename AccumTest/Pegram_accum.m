cd('\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Comparison_Test_20\Accumulation_tests')
load('Pegram_Results.mat')

radius = 1;
coord1 = [94, 131];
coord2 = [174 212];

peg1 = ensembles((coord1(2)-radius):(coord1(2)+radius),(coord1(1)-radius):(coord1(1)+radius),:,:);
peg2 = ensembles((coord2(2)-radius):(coord2(2)+radius),(coord2(1)-radius):(coord2(1)+radius),:,:);

peg1m = squeeze(mean(mean(peg1)));
peg2m = squeeze(mean(mean(peg2)));

accpeg1 = zeros(25,100);
accpeg2 = zeros(25,100);

accpeg1(1,:) = peg1m(1,:);
accpeg2(1,:) = peg2m(1,:);

for t=2:25
    accpeg1(t,:)=accpeg1((t-1),:)+peg1m(t,:);
    accpeg2(t,:)=accpeg2((t-1),:)+peg2m(t,:);    
end

sqkm = (radius*2+1)^2;
save(['Peg_acc_',num2str(sqkm),'km2.mat'],'coord1','coord2','radius','accpeg1','accpeg2')

