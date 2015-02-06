

clear all, close all;

filename = 'metoffice-c-band-rain-radar_uk_200707200900_1km-composite.dat';
[mtx date parameters] = readnimrod(filename);   % function to read nimrod Met Office rainfall data

if ( prod(size(mtx))>0 )    % if the data matrix is NOT empty
    figure(1), pcolor( mtx ); shading flat; colorbar('vert'); caxis([0 10]);
    title(sprintf('Radar rainfall %.4d/%.2d/%.2d %.2d:%.2d:%.2d', date(1),date(2),date(3),date(4),date(5),date(6) ));
    xlabel('Range [km]'); ylabel('Range [km]'); 
else
    sprintf('FILE NOT FOUND!!!\n')
end


