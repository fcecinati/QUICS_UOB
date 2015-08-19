
% Load the masks
load('\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\8 - Other data\Nanding\catchment boundary\catchment_masks_ensembles.mat')

%set up variables
lun_series = zeros(100,8784);
swa_series = zeros(100,8784);
rib_series = zeros(100,8784);
missing = zeros(8784,1);
start = '2007-10-01 00:00:00';
finish = '2008-09-30 23:00:00';
i=0;

%loop for each time step: opens the file, and for each ensemble takes the average of each
%basin using the masks
for t = datenum(start):(1/24):datenum(finish)
    i=i+1;
    file = ['F:\Francesca\HBV',num2str(year(t),'%04d'),num2str(month(t),'%02d'),'\ensemble-time',num2str(day(t)*24 - (23 - hour(t))),'.mat'];
    try
        load(file)
        for en=1:100
            mat = newens(:,:,en);
            lun_series(en,i) = mean(mat(lun));
            swa_series(en,i) = mean(mat(swa));
            rib_series(en,i) = mean(mat(rib));
        end
    catch
        missing(i)=1;
    end    
end

% Interpolate missing values
Xall = 1:8784;
Xknown = Xall(missing ==0);
Xquery = Xall(missing ==1);

for en=1:100
    V_lun = squeeze(lun_series(en,:)); Vknown_lun = V_lun(missing == 0);
    V_swa = squeeze(swa_series(en,:)); Vknown_swa = V_swa(missing == 0);
    V_rib = squeeze(rib_series(en,:)); Vknown_rib = V_rib(missing == 0);
    
    Vquery_lun = interp1(Xknown,Vknown_lun,Xquery,'cubic');
    Vquery_swa = interp1(Xknown,Vknown_swa,Xquery,'cubic');
    Vquery_rib = interp1(Xknown,Vknown_rib,Xquery,'cubic');
    
    lun_series(missing == 1) = Vquery_lun;
    swa_series(missing == 1) = Vquery_swa;
    rib_series(missing == 1) = Vquery_rib;
end

% save files
for en=1:100
    dlun = cell(8784,2);dswa = cell(8784,2);drib = cell(8784,2);
    i=0;
    for t = datenum(start):(1/24):datenum(finish)
        i=i+1;
        dlun(i,1) = cellstr([num2str(day(t),'%02d'),'.',num2str(month(t),'%02d'),'.',num2str(year(t),'%04d'),' ',num2str(hour(t),'%02d'),':',num2str(minute(t),'%02d'),]);
        dlun(i,2) = cellstr(num2str(lun_series(en,i),'%3.3f'));
        dswa(i,1) = cellstr([num2str(day(t),'%02d'),'.',num2str(month(t),'%02d'),'.',num2str(year(t),'%04d'),' ',num2str(hour(t),'%02d'),':',num2str(minute(t),'%02d'),]);
        dswa(i,2) = cellstr(num2str(swa_series(en,i),'%3.3f'));
        drib(i,1) = cellstr([num2str(day(t),'%02d'),'.',num2str(month(t),'%02d'),'.',num2str(year(t),'%04d'),' ',num2str(hour(t),'%02d'),':',num2str(minute(t),'%02d'),]);
        drib(i,2) = cellstr(num2str(rib_series(en,i),'%3.3f'));
    end
    nameswa = ['F:\Francesca\EnsemblesHBV\swa_ens',num2str(en),'.txt'];
    namelun = ['F:\Francesca\EnsemblesHBV\lun_ens',num2str(en),'.txt'];
    namerib = ['F:\Francesca\EnsemblesHBV\rib_ens',num2str(en),'.txt'];
    dlmcell(nameswa,dswa)
    dlmcell(namerib,drib)
    dlmcell(namelun,dlun)
end    