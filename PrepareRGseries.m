start = '2008-10-01 00:00:00';
finish = '2010-09-30 23:00:00';

sy = year(start); sm = month(start); sd = day(start); sh = hour(start);
fy = year(finish); fm = month(finish); fd = day(finish); fh = hour(finish);

G = load('\\ads.bris.ac.uk\filestore\myfiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\data_north_england\gauge_1h_2007-2011.dat');
C = load('\\ads.bris.ac.uk\filestore\myfiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\data_north_england\gaugecoordinates.dat');

si = find(G(:,1)==sy & G(:,2)==sm & G(:,3)==sd & G(:,4)== sh);
fi = find(G(:,1)==fy & G(:,2)==fm & G(:,3)==fd & G(:,4)== fh);

G = G(si:fi,7:end); G(G<0) = NaN;
t_steps = size(G,1); n_stations = size(G,2);

Lshape = shaperead('\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\8 - Other data\Nanding\catchment boundary\Lune\Lune.shp');
Sshape = shaperead('\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\8 - Other data\Nanding\catchment boundary\Swale\Swale1.shp');
Rshape = shaperead('\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\8 - Other data\Nanding\catchment boundary\Ribble\Ribble1.shp');
LX = (Lshape.X(1:end-1)-320000)/1000; LY = (Lshape.Y(1:end-1)-370000)/1000;
SX = (Sshape.X(1:end-1)-320000)/1000; SY = (Sshape.Y(1:end-1)-370000)/1000;
RX = (Rshape.X(1:end-1)-320000)/1000; RY = (Rshape.Y(1:end-1)-370000)/1000;
L = poly2mask(LX,LY,180,180);
S = poly2mask(SX,SY,180,180);
R = poly2mask(RX,RY,180,180);

xq = 320500:1000:499500;
yq = 370500:1000:549500;
[Xq,Yq] = meshgrid(xq,yq);

Lun = zeros(t_steps,1); Swa = zeros(t_steps,1); Rib = zeros(t_steps,1);
dlun = cell(t_steps,2);dswa = cell(t_steps,2);drib = cell(t_steps,2);
for i=1:t_steps
    v = G(i,:);
    cx = C(:,1); cy = C(:,2);
    cx(isnan(v)==1)=[]; cy(isnan(v)==1)=[];
    v(isnan(v)==1)=[];
    vq = griddata(cx,cy,v,Xq,Yq,'v4'); 
    vq(vq<0)=0;
    Lun(i) = mean(vq(L));
    Swa(i) = mean(vq(S));
    Rib(i) = mean(vq(R));
    
    t = datenum(start) + (i-1)/24; 
    
    dlun(i,1) = cellstr([num2str(day(t),'%02d'),'.',num2str(month(t),'%02d'),'.',num2str(year(t),'%04d'),' ',num2str(hour(t),'%02d'),':',num2str(minute(t),'%02d'),]);
    dlun(i,2) = cellstr(num2str(Lun(i),'%3.3f'));
    dswa(i,1) = cellstr([num2str(day(t),'%02d'),'.',num2str(month(t),'%02d'),'.',num2str(year(t),'%04d'),' ',num2str(hour(t),'%02d'),':',num2str(minute(t),'%02d'),]);
    dswa(i,2) = cellstr(num2str(Swa(i),'%3.3f'));
    drib(i,1) = cellstr([num2str(day(t),'%02d'),'.',num2str(month(t),'%02d'),'.',num2str(year(t),'%04d'),' ',num2str(hour(t),'%02d'),':',num2str(minute(t),'%02d'),]);
    drib(i,2) = cellstr(num2str(Rib(i),'%3.3f'));
end

save('\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\8 - Other data\Nanding\CalibrationR\ValTimeSeries','Lun','Swa','Rib')
dlmcell('\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\8 - Other data\Nanding\CalibrationR\SWAVal.txt',dswa)
dlmcell('\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\8 - Other data\Nanding\CalibrationR\RIBVal.txt',drib)
dlmcell('\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\8 - Other data\Nanding\CalibrationR\LUNVal.txt',dlun)