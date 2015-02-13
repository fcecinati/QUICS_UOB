clear all, close all;

%% INFORMATION TO PROVIDE: 

starting = '2008-09-05 15:00:00';
time_steps = 20;                    % in hours
test = 10;
radar_folder = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\data_north_england\09-2008\';
ensembles_folder = ['\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_German\Test_',num2str(test),'\'];
tmpfolder = '\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Temp\';
mkdir('\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_German\Test_10\','Ensembles');
results_folder = ['\\ads.bris.ac.uk\filestore\MyFiles\Staff3\fc14509\Documents\QUICS\MATLAB_rep\Results\Test_German\Test_',num2str(test),'\Ensembles\'];
n_ensembles = 5;

%% Load Radar
st = datenum(starting);

for i = 1:time_steps
    date = st + (i-1)/24;
    yr = year(date);    mn = month(date);   dy = day(date);     hr = hour(date);
    filename = [radar_folder,'metoffice-c-band-rain-radar_uk_',num2str(yr),num2str(mn, '%02d'),num2str(dy, '%02d'),num2str(hr, '%02d'),'00_1km-composite.dat.gz'];

    % check if file is compressed
    tt=length(filename);
    if ( strcmp('.gz', filename(tt-2:tt)) == 1)   % if the input file is compressed (i.e. with extension '.gz')
         newname = gunzip(filename, tmpfolder); % decompress file in a temporary folder and output new file name
         filename = newname{1,1};    % copy new name into the nimrodfilename
    end
    [mtx(:,:,i), r_date, r_param] = readnimrod(filename);   % function to read nimrod Met Office rainfall data
        
    for en = 1:n_ensembles
        ensname =  [ensembles_folder, 'Ensemble_n',num2str(en),'_step' num2str(i),'.dat'];
        [noise, e_date, e_param] = readnimrod(ensname); 

        x_start = (e_param(3) - r_param(3))/ 1000;
        x_end = ((e_param(3) + e_param(7)*1000) - r_param(3))/ 1000 -1 ;
        y_start = ((e_param(2) - e_param(6)*1000) - (r_param(2) - r_param(6)*1000))/ 1000 +1;
        y_end = ((e_param(2)) - (r_param(2) - r_param(6)*1000))/ 1000;

        raddata(:,:,i) = mtx(y_start:y_end, x_start:x_end, i);

        logens = noise + 10*(log(raddata(:,:,i)));
        ens(:,:,i,en) = exp(logens/10);

        fig1 = figure;
        fig2 = figure;

        figure(fig1)
        imagesc(raddata(:,:,i))
        colorbar
        name = [results_folder,'Radar_data_step',num2str(i),'.jpg'];
        saveas(fig1, name)
        close(fig1)

        figure(fig2)
        imagesc(ens(:,:,i))
        colorbar
        name = [results_folder,'Ensemble',num2str(en),'_step',num2str(i),'.jpg'];
        saveas(fig2, name)  
        close(fig2)
    end
end

save ([results_folder,'Test_',num2str(test)], 'ens', 'raddata')
