function ensembles = Ensembles(delta, coord, sim_hours) % Start and end date to be specified with time too
% Generates the ensembles given the deltas from the german method and the
% start and end date of the event to model

sz = size(delta);
x = sz(1);
t = sz(2);
n = sz(3);

% Define the grid to interpolate
coord = coord/1000;
minx = min(coord(1,:));   minx = floor(minx)-0.5;
maxx = max(coord(1,:));   maxx = ceil(maxx)+0.5;
miny = min(coord(2,:));   miny = floor(miny)-0.5;
maxy = max(coord(2,:));   maxy = ceil(maxy)+0.5;
x_interp = (minx:maxx);
y_interp = (miny:maxy);
[x_int, y_int] = meshgrid(x_interp, y_interp);
X = size(x_interp);
Y = size(y_interp);

ensembles = zeros(X,Y,sim_hours,n);

for ens=1:n
    % Random selection of the starting point
    rand_start = randi(t - sim_hours);
    deltas = delta(:,rand_start:(rand_start + sim_hour - 1),ens);

    % interpolate the deltas
    for i=1:sim_hours
        d = deltas(:,i,n);
        perturbation = griddata(coord(1,:),coord(2,:),d,x_int,y_int, 'v4');
        ensembles(:,:,i,n) = perturbation;
    end
end

% Write the output


    
        

  