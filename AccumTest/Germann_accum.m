cd('\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Comparison_Test_20\Accumulation_tests')
load('Germann_Results.mat')
radius = 1;
coord1 = [94, 131];
coord2 = [174 212];

ger1 = ensembles((coord1(2)-radius):(coord1(2)+radius),(coord1(1)-radius):(coord1(1)+radius),:,:);
ger2 = ensembles((coord2(2)-radius):(coord2(2)+radius),(coord2(1)-radius):(coord2(1)+radius),:,:);

rad1 = raddata((coord1(2)-radius):(coord1(2)+radius),(coord1(1)-radius):(coord1(1)+radius),:,:);
rad2 = raddata((coord2(2)-radius):(coord2(2)+radius),(coord2(1)-radius):(coord2(1)+radius),:,:);

ger1m = squeeze(mean(mean(ger1)));
ger2m = squeeze(mean(mean(ger2)));

rad1m = squeeze(mean(mean(rad1)));
rad2m = squeeze(mean(mean(rad2)));

accger1 = zeros(25,100);
accger2 = zeros(25,100);

accrad1 = zeros(25,1);
accrad2 = zeros(25,1);

accger1(1,:) = ger1m(1,:);
accger2(1,:) = ger2m(1,:);

accrad1(1) = rad1m(1);
accrad2(1) = rad2m(1);

for t=2:25
    accger1(t,:)=accger1((t-1),:)+ger1m(t,:);
    accger2(t,:)=accger2((t-1),:)+ger2m(t,:);
    
    accrad1(t)=accrad1((t-1))+rad1m(t);
    accrad2(t)=accrad2((t-1))+rad2m(t);
    
end

sqkm = (radius*2+1)^2;
save(['Ger_acc_',num2str(sqkm),'km2.mat'],'coord1','coord2','radius','accger1','accger2','accrad1','accrad2')

