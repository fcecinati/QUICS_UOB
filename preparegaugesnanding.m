t = size(gauges,1);
x = size(gauges,2);
intx = 340000:1000:470000; xs = size(intx,2);
inty = 430000:1000:550000; ys = size(inty,2);
[X,Y] = meshgrid(intx,inty); 
interp_gauges = zeros(ys,xs,t);
for i=1:t
    g = gauges(i,:)';
    c = coord; c(isnan(g)==1,:) = []; g(isnan(g)==1,:) = [];    
    interp_gauges(:,:,i) = griddata(c(:,1),c(:,2),g,X,Y,'v4');
end
save('\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\8 - Other data\Nanding\CalibrationR\interp_gauges','interp_gauges')