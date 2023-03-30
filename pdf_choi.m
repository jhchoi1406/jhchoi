%% calculating IABP drift speed %%
clear;
data0 = load('buoys_1979-2021.txt');

data = sortrows(data0,5);

lat0 = data(:,6);
lon0 = data(:,7);

% error1
idx_error = find(lat0 > 90 | lat0 < 68.5 | lon0 > 360 | lon0 < -180);
data(idx_error,:) = [];

year0 = data(:,1);
month0 = data(:,2);
day0 = data(:,3);
hour0 = data(:,4);
id0 = data(:,5);
lat0 = data(:,6);
lon0 = data(:,7);

idx = find(hour0 == 0);

year = year0(idx);
month = month0(idx);
day = day0(idx);
hour = hour0(idx);
id = id0(idx);
lat = lat0(idx);
lon1 = lon0(idx);

%% calculating total DOY & month data %%

% doy

y365 = [0 31 28 31 30 31 30 31 31 30 31 30 31];

sum365 = cumsum(y365);

doy = zeros(length(id),1);
for i = 1:length(id);
    doy(i) = sum365(month(i)) + day(i);
end

yy_day = 365*ones(2021-1979+1,1);
yy_day(1:1:end) = 365;
sumday0 = cumsum(yy_day);
sumday = [0; sumday0];

cum_doy = sumday(year-1979+1)+doy;

%% calculating the daily velocity %% 
% velocity %

idx_error1 = find( diff(cum_doy) ~= 1);

lon = lon1;
idxlon = find(lon1 < 0);
lon(idxlon) = lon(idxlon)+360;

dt = 24*60*60; 

lon_b = lon(1:end-1); 
lon_f = lon(2:end); 

lat_b = lat(1:end-1); 
lat_f = lat(2:end); 

R0 = 6378*1000; 
r0 = 2*pi*R0; 
r = 2*pi*R0*cosd(lat_b); 

% U & V velocity (east-west) % 

Ui = (r/360).*(lon_f - lon_b)/dt;
Vi = (r0/360)*(lat_f - lat_b)/dt;

Ui(idx_error1) = NaN;
Vi(idx_error1) = NaN;

idxa = find(abs(Ui) > 1);
Ui(idxa) = NaN;
 
idxb = find(abs(Vi) > 1);
Vi(idxb) = NaN;

idx_zero = find(Ui==0 & Vi==0);
Ui(idx_zero) = NaN;
Vi(idx_zero) = NaN;

%% importing wind data - u10, v10 %%

ncfile_u10='u10_1979-2020_ERA5.nc';

u10=ncread(ncfile_u10,'u10');
u10lat = ncread(ncfile_u10, 'lat');
u10lon = ncread(ncfile_u10, 'lon');

ncfile_v10='v10_1979-2020_ERA5.nc';

v10=ncread(ncfile_v10,'v10');
v10lat = ncread(ncfile_v10, 'lat');
v10lon = ncread(ncfile_v10, 'lon');

%% finding out grid location - u10,v10 %%

idx_extra = find(cum_doy > 15330); 
cum_doy(idx_extra) = [];
Vi(idx_extra) = [];
Ui(idx_extra) = [];
year(idx_extra) = [];
month(idx_extra) = [];
lat(idx_extra) = [];
lon(idx_extra) = [];

lon1 = lon;
idxlon = find(lon < 0);
lon1(idxlon) = lon(idxlon)+360;

glat = NaN*zeros(size(lat));
glon = NaN*zeros(size(lon1));

u10_up = NaN*zeros(size(Ui));

for i = 1:length(Ui)
    time = cum_doy(i);
    flat = lat(i);
    flon = lon1(i);
    
    if (isnan(flat)==1 | isnan(flon) ==1)
        glat(i) = NaN;
        glon(i) = NaN;
    else 
        glat(i) = max( find( abs(u10lat-flat) == min(abs(u10lat-flat)) ) );
        glon(i) = max( find( abs(u10lon-flon) == min(abs(u10lon-flon)) ) );
        
        if (isnan(time)~=1)
            u10_up(i) = u10(glon(i),glat(i),time);
        end
    end
    
end


v10_up = NaN*zeros(size(Vi));

for i = 1:length(Vi)
    time = cum_doy(i);
    flat = lat(i);
    flon = lon1(i);
    
    if (isnan(flat)==1 | isnan(flon) ==1)
        glat(i) = NaN;
        glon(i) = NaN;
    else 
        glat(i) = max( find( abs(v10lat-flat) == min(abs(v10lat-flat)) ) );
        glon(i) = max( find( abs(v10lon-flon) == min(abs(v10lon-flon)) ) );
        
        if (isnan(time)~=1)
            v10_up(i) = v10(glon(i),glat(i),time);
        end
    end
    
end

clear u10;
clear v10;

%% calculating angle %%

% buoys
 
angle_buoy = atan2d((Ui.*v10_up-Vi.*u10_up),(Ui.*u10_up+Vi.*v10_up));

% 2007년 이전

angle_buoy_bf = angle_buoy(find(year <= 2007));

histogram(angle_buoy_bf,'Normalization','probability')

% 2007년 이후

angle_buoy_aft = angle_buoy(find(year()>2007));



