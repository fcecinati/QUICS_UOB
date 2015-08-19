function tseries = tseries(ens,boundaryfile, easting, northing)

sz = size(ens);
side = sz(1); tsteps = sz(3); nens = sz(3);
S = shaperead(boundaryfile);
xq = easting:easting+(side*1000);
yq = northing:-1000:northing-(side*1000);
xv = S.X;
yv = S.Y;
in = inpolygon(xq,yq,xv,yv);

tseries = zeros(tsteps,nens);
for t=1:tsteps
    for en=1:nens
        tseries(t,en) = mean(R(in));
    end
end
