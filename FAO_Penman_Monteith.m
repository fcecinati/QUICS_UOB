%-- FAO Penman_monteith equation

% plot  air moisture, temperature, pressure
%
% -- inputs: Hyrex-> modelling->
% brue_all_rain_flow_soil_weather_MOSES(1993-200_hourly.xls

% -- put them to Matlab as
%  air_97           Solar radiation, Net radiation,  Wet temp, dry temp,  wind speed
%                   Wind direction, rainfall, pressure

%--- load SMD values of MORECS and MOSES
%save 'C:\Han\Matlab_Prog\Soil Moisture\data_soil_94_00_all' rain94 flow94 smd94 soil94 rain95 flow95 smd95 soil95 rain96 flow96 smd96 soil96 rain97 flow97 smd97 soil97 rain98 flow98 smd98 soil98 rain99 flow99 smd99 soil99 rain00 flow00 smd00 soil00
%load 'C:\Han\Matlab_Prog\Soil Moisture\data_soil_94_00_all';

%--- distribute data into colums
air=air_95;

solar=air(:,1); %   Solar radiation, W/m^2
Nsolar=air(:,2); %   Net Solar radiation, W/m^2
Wet_temp=air(:,3); %   Wet temperature, C
Dry_temp=air(:,4); %   Dry temperature, C
Wind_speed=air(:,5); % wind speed m/s
Wind_direction=air(:,6); % wind direction, degree
air_rain=air(:,7); % rainfall at the station, mm
Pressure=air(:,8); %  air pressure, mb

%-- month spans
year=2000;
dm1=[0 31 60 91 121 152 182 213 244 274 305 335, 366]; % leap year
dm2=[0 31 59 90 120 151 181 212 243 273 304 334, 365]; % normal year
month_str={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','All'};
istarthour=ones(13,1);iendhour=31*ones(13,1); %-- ini array
if year==floor(year/4)*4   % leap year
    for i=1:12,istarthour(i)=dm1(i)*24+1; iendhour(i)=dm1(i+1)*24-1;end
else % normal year
    for i=1:12,istarthour(i)=dm2(i)*24+1; iendhour(i)=dm2(i+1)*24-1; end
end
%-- month 13 is for the whole year
istarthour(13)=1; iendhour(13)=length(Pressure);


%----------------------------------------------------------
%--- in daily unit
%-- temperature
N_day=365;
if year==floor(year/4)*4   % leap year
    N_day=366;
end

day_wet_temp=reshape(Wet_temp,24,N_day); %-- each column is a day
day_dry_temp=reshape(Dry_temp,24,N_day); %-- each column is a day
day_max_dry_temp=max(day_dry_temp);
day_min_dry_temp=min(day_dry_temp);
day_mean_temp=mean(day_dry_temp); %  mean daily temperature

day_max_wet_temp=max(day_wet_temp);
day_min_wet_temp=min(day_wet_temp);

%--- air pressure ---
day_pressure=reshape(Pressure,24,N_day)/10; %-- from mb -> KPa
day_mean_pressure=mean(day_pressure);



% 1) Mean saturation vapour pressure (es)
es=(weather_vapourpressure(day_max_dry_temp)+weather_vapourpressure(day_min_dry_temp))/2;

% 2)Slope of saturation vapour pressure curve (delta)
Tm=(day_max_dry_temp+day_min_dry_temp)/2;
delta=4098*0.6108*exp(17.27*Tm./(Tm+237.3))./(Tm+237.3).^2;

% 3) Actual vapour pressure (ea) derived from psychrometric data
% actual vapour pressure is almost constant during the day, so just
% calculate the one at the max temperature
%-- if actual air pressure is not measured, use elevation to estimate
%z=20; % elevation in m
%P=101.3*((293-0.0065*z)/293)^5.26; % mean air pressure in KPa
% gama: psychrometric constant of the instrument 
gama = 0.0008*day_mean_pressure;
ea=weather_vapourpressure(day_max_wet_temp)-gama.*(day_max_dry_temp-day_max_wet_temp)

% 4) Determination of vapour pressure deficit es-ea
ed=es-ea;

% 5) daily Net solar radiation
%-- Net radiation
Rn=mean(reshape(Nsolar,24,N_day))*0.0864; %-- W/m^2 -> MJ m-2 day-1

% 6) Soil heat flux (G) 
%  ignored for daily data

%7) wind
U2=mean(reshape(Wind_speed,24,N_day));
%-- U2 should be no less than 0.5m
U2(U2<0.5)=0.5;

%8) Penman-Monteith euqation

et0=(0.408*delta.*Rn+gama*900.*U2.*ed./(day_mean_temp+273))./(delta+gama.*(1+0.34*U2));
et0=et0'; %-- into a column

%-- distribute into 24 hours
% convert daily Eo data to hourly
 nday=length(et0); %-- 
 et0_hour=zeros(nday*24,1);%  initialise the array
 %-- daily to hourly 
 k=0;
%-- 1st day, use the SMD at end of the day
     for j=1:24
         k=k+1;
         et0_hour(k)=et0(1);
     end     

 %-- from 2nd day till the end, interpolate
 for i=2:nday
     a1=et0(i-1);
     a2=et0(i);
     dd=(a2-a1)/24;
     for j=1:24
         k=k+1;
         et0_hour(k)=a1+dd*j;
     end     
 end
 et0_hour=et0_hour/24;
figure(1)
 plot(et0_hour);
 title('Reference Evapotranspiration in 2000 converted from daily to hour')
xlabel('Time (h)')
ylabel('Eo(mm/h)')
set(gcf,'color','w')

