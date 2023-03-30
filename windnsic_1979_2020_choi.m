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

%% total DOY %%
y365 = [0 31 28 31 30 31 30 31 31 30 31 30 31];
y366 = [0 31 29 31 30 31 30 31 31 30 31 30 31];

sum365 = cumsum(y365);
sum366 = cumsum(y366);

doy = zeros(length(id),1);
for i = 1:length(id);
    if rem(year(i),4) == 0;
        doy(i) = sum366(month(i)) + day(i);
    else
        doy(i) = sum365(month(i)) + day(i);
    end
end

yy_day = 365*ones(2021-1979+1,1);
yy_day(2:4:end) = 366;
sumday0 = cumsum(yy_day);
sumday = [0; sumday0];

cum_doy = sumday(year-1979+1)+doy;

idx_error1 = find( diff(cum_doy) ~= 1);

lon = lon1;
idxlon = find(lon1 < 0);
lon(idxlon) = lon(idxlon)+360;

%% importing wind & sic data %%

ncfile='u10_1979-2020_ERA5.nc';

u10=ncread(ncfile,'u10');
ulat = ncread(ncfile, 'lat');
ulon = ncread(ncfile, 'lon');

ncfile='v10_1979-2020_ERA5.nc';

v10=ncread(ncfile,'v10');
vlat = ncread(ncfile, 'lat');
vlon = ncread(ncfile, 'lon');

ncfile='sic_1979-2020_NSIDC.nc';

sic = ncread(ncfile,'cdr_seaice_conc');
clat = ncread(ncfile,'lat');
clon = ncread(ncfile,'lon');
time = ncread(ncfile, 'time');

%% calculating the daily velocity %% 
% velocity %

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

%% finding out grid location %%

idx_extra = find(cum_doy > 15330); 

cum_doy(idx_extra) = [];
Vi(idx_extra) = [];
Ui(idx_extra) = [];
year(idx_extra) = [];
lat(idx_extra) = [];
lon(idx_extra) = [];

glat = NaN*zeros(size(lat));
glon = NaN*zeros(size(lon));

u10_up = NaN*zeros(size(Ui));

for i = 1:length(Ui)
    time = cum_doy(i);
    flat = lat(i);
    flon = lon(i);
    
    if (isnan(flat)==1 | isnan(flon) ==1)
        glat(i) = NaN;
        glon(i) = NaN;
    else 
        glat(i) = max( find( abs(ulat-flat) == min(abs(ulat-flat)) ) );
        glon(i) = max( find( abs(ulon-flon) == min(abs(ulon-flon)) ) );
        
        if (isnan(time)~=1)
            u10_up(i) = u10(glon(i),glat(i),time);
        end
    end
    
end


v10_up = NaN*zeros(size(Vi));

for i = 1:length(Vi)
    time = cum_doy(i);
    flat = lat(i);
    flon = lon(i);
    
    if (isnan(flat)==1 | isnan(flon) ==1)
        glat(i) = NaN;
        glon(i) = NaN;
    else 
        glat(i) = max( find( abs(vlat-flat) == min(abs(vlat-flat)) ) );
        glon(i) = max( find( abs(vlon-flon) == min(abs(vlon-flon)) ) );
        
        if (isnan(time)~=1)
            v10_up(i) = v10(glon(i),glat(i),time);
        end
    end
    
end

sic_up = NaN*zeros(size(Vi));

for i = 1:length(sic_up)
    time = cum_doy(i);
    flat = lat(i);
    flon = lon(i);
    
    if (isnan(flat)==1 | isnan(flon) ==1)
        glat(i) = NaN;
        glon(i) = NaN;
    else 
        glat(i) = max( find( abs(clat-flat) == min(abs(clat-flat)) ) );
        glon(i) = max( find( abs(clon-flon) == min(abs(clon-flon)) ) );
        
        if (isnan(time)~=1)
           sic_up(i) = sic(glon(i),glat(i),time);
        end
    end
    
end

idxnan = find(isnan(Ui) == 1);
sic_up(idxnan) = NaN;

%% wind speed & ice speed %%

icespeed = sqrt(Ui.^2 + Vi.^2);
wdspeed = sqrt(u10_up.^2+v10_up.^2);

idxw = find(wdspeed >= 1.0);
% idxs = find(sic_up > 0.9);

wdspeed_idxw = wdspeed(idxw);
icespeed_idxw = icespeed(idxw);

yearx = year(idxw);
yeary = 1979:1:2020;

wns = [];
for i = 1979:1:2020
     idxy = find(yearx == i);
     frac_speed = 100*(icespeed_idxw./wdspeed_idxw);
     frac_speed0 = frac_speed(idxy);
     mean_frac_speed = nanmean(frac_speed0);
     wns=[wns; mean_frac_speed];
end

figure(1)
plot(yeary, wns, 'LineWidth', 2); grid on;
xlim([1978 2021])
ylim([1.2 2.8])
title('Annual-mean fractional sea ice speed (m/s)', 'FontSize', 13)
