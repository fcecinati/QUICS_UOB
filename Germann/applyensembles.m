function [ens, raddata] = applyensembles(starting, time_steps, n_ensembles, radar_folder, tmpfolder, results_folder)

%% Creates output folder
mkdir(results_folder,'Ensembles');
ensembles_folder = [results_folder,'\Ensembles\'];

%% Ensemble generation

st = datenum(starting);

for i = 1:time_steps
    % Find the radar data for the correct date
    date = st + (i-1)/24;
    yr = year(date);    mn = month(date);   dy = day(date);     hr = hour(date);
    filename = [radar_folder,'metoffice-c-band-rain-radar_uk_',num2str(yr),num2str(mn, '%02d'),num2str(dy, '%02d'),num2str(hr, '%02d'),'00_1km-composite.dat.gz'];

    % Check if file is compressed
    tt=length(filename);
    if ( strcmp('.gz', filename(tt-2:tt)) == 1)   
         newname = gunzip(filename, tmpfolder); 
         filename = newname{1,1};    
    end
    
    % Open the radar data
    [mtx(:,:,i), r_date(:,i), r_param] = readnimrod(filename);   % function to read nimrod Met Office rainfall data
    
    % Generates Ensembles
    for en = 1:n_ensembles
        % Find the noise files
        noisename =  [results_folder, 'Noise_n',num2str(en),'_step' num2str(i),'.dat'];
        [noise, e_date, e_param] = readnimrod(noisename); 
        
        % Define the domain
        x_start = (e_param(3) - r_param(3))/ 1000;
        x_end = ((e_param(3) + e_param(7)*1000) - r_param(3))/ 1000 -1 ;
        y_start = ((e_param(2) - e_param(6)*1000) - (r_param(2) - r_param(6)*1000))/ 1000 +1;
        y_end = ((e_param(2)) - (r_param(2) - r_param(6)*1000))/ 1000;
        
        % Trim the radar data on the domain
        raddata(:,:,i) = mtx(y_start:y_end, x_start:x_end, i);
        
        % Combine noise and radar data
        logens = noise + 10*(log(raddata(:,:,i)));
        ens(:,:,i,en) = exp(logens/10);
        
        % Generates figures
        fig1 = figure;
        fig2 = figure;

        figure(fig1)
        imagesc(raddata(:,:,i))
        colorbar
        name = [ensembles_folder,'Radar_data_',sprintf('%02d',r_date(:,i)),'.jpg'];
        saveas(fig1, name)
        close(fig1)

        figure(fig2)
        imagesc(ens(:,:,i))
        colorbar
        name = [ensembles_folder,'Ensemble',num2str(en),'_',sprintf('%02d',r_date(:,i)),'.jpg'];
        saveas(fig2, name)  
        close(fig2)
    end
end

%% Save results in Nimrod format

PARAM = e_param;

for en = 1:n_ensembles
    for time = 1:time_steps
        DATE = r_date(:,time);
        name = [ensembles_folder, 'Ensemble_n',num2str(en),'_',sprintf('%02d',r_date(:,time)),'.dat'];
        MTX = ens(:,:,time,en);
        
        save2nimrod(name,MTX,DATE,PARAM)
        
    end
end
