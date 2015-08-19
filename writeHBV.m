start = '2007-07-01 00:00:00';
until = '2008-06-30 23:00:00';

s = datenum(start);
u = datenum(until);
mat = zeros(8784,7);
d = cell(8784,2);
i=0;
for t = s:(1/24):u
    i=i+1;
    d(i,1) = cellstr([num2str(day(t),'%02d'),'.',num2str(month(t),'%02d'),'.',num2str(year(t),'%04d'),' ',num2str(hour(t),'%02d'),':',num2str(minute(t),'%02d'),]);
    d(i,2) = cellstr(num2str(value(i),'%3.3f'));
end

dlmcell('\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\8 - Other data\Nanding\CalibrationR\LUNflow.txt',d)



    
    