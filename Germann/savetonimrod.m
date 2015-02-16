    
nrows=300; ncols=200;
MTX = randn(nrows,ncols);
MTX = imfilter(MTX ,fspecial('gaussian',[25 25],5));
MTX = 10*(MTX + 0.2);
figure(1), pcolor(MTX), shading flat;
colorbar('vert');

ytopleft=500000;    % in meters
xtopleft=150000;    % in meters
grid=1000;    % grid resolution in meters

date1 = datenum(2012,1,1,3,15,0);
[year month day hr min sec]=datevec(date1);
fname = sprintf('File_%.4d%.2d%.2d-%.2d%.2d.dat',year,month,day,hr,min);
DATE=[year;month;day;hr;min;sec];
PARAM=[0; ytopleft; xtopleft ; grid; grid; nrows; ncols];                       
save2nimrod(fname, MTX, DATE, PARAM);
