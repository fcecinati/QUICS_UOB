function noise = Noise(delta, coord, sim_hours, results_path) 
% Generates the noise component for the ensembles given the deltas from the Germann method

% Size of the deltas
sz = size(delta);
x = sz(1);
t = sz(2);
if numel(sz)>2
    n = sz(3);
else
    n = 1;
end

% Define the grid to interpolate
coord = coord/1000;
minx = min(coord(1,:));   minx = floor(minx)-0.5;
maxx = max(coord(1,:));   maxx = ceil(maxx)+0.5;
miny = min(coord(2,:));   miny = floor(miny)-0.5;
maxy = max(coord(2,:));   maxy = ceil(maxy)+0.5;
x_interp = (minx:maxx);
y_interp = (miny:maxy);
[x_int, y_int] = meshgrid(x_interp, y_interp);
X = size(x_interp);     X = X(2);
Y = size(y_interp);     Y = Y(2);

noise = zeros(Y,X,sim_hours,n);

for ens=1:n
    % Random selection of the starting point
    rand_start = randi(t - sim_hours);
    deltas = delta(:,rand_start:(rand_start + sim_hours - 1),ens);

    % interpolate the deltas
    for i=1:sim_hours
        d = deltas(:,i);
        perturbation = griddata(coord(1,:),coord(2,:),d,x_int,y_int, 'v4');
        noise(:,:,i,ens) = perturbation;
    end
end

% Write the output
name = [results_path, 'Ensembles'];
save(name,'ensembles')

PARAM = [0; maxy*1000; minx*1000; 1000; 1000; Y; X];
for ens = 1:n
    for time = 1:sim_hours
        DATE = [2000;1;1;time;0;0];
        name = [results_path, 'Ensemble_n',num2str(ens),'_step',num2str(time),'.dat'];
        MTX = noise(:,:,time,ens);
        
        save2nimrod(name,MTX,DATE,PARAM)
        
    end
end


    
        

  