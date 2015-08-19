t = size(gauges,1);
x = size(gauges,2);

lun_avg = zeros(t,1);
rib_avg = zeros(t,1);
swa_avg = zeros(t,1);
for i=1:t
    area = interp_gauges(:,:,i); area = area(lun_mask);
    lun_avg(i) = mean(area);
    area = interp_gauges(:,:,i); area = area(rib_mask);
    rib_avg(i) = mean(area);
    area = interp_gauges(:,:,i); area = area(swa_mask);
    swa_avg(i) = mean(area);    
end
save('\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\8 - Other data\Nanding\CalibrationR\lun_avg','lun_avg')
save('\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\8 - Other data\Nanding\CalibrationR\rib_avg','rib_avg')
save('\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\8 - Other data\Nanding\CalibrationR\swa_avg','swa_avg')