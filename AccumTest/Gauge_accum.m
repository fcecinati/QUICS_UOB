cd('\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\data_north_england\test20');
load('goodgauges.mat')

G1 = goodgauges(:,76);
G2 = goodgauges(:,161);

accG1 = zeros(25,1);
accG2 = zeros(25,1);
accG1(1) = G1(1);
accG2(1) = G2(1);

for t=2:25
    accG1(t)=accG1((t-1))+G1(t);
    accG2(t)=accG2((t-1))+G2(t);
end

save('\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Comparison_Test_20\Accumulation_tests\Gauge.mat','G1','G2','accG1','accG2')

