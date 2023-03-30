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

idx_error1 = find( diff(cum_doy) ~= 1);

lon = lon1;
idxlon = find(lon1 < 0);
lon(idxlon) = lon(idxlon)+360;

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

%% importing wind data - u10, v10 %%

ncfile='u10_1979-2020_ERA5.nc';

u10=ncread(ncfile,'u10');
u10lat = ncread(ncfile, 'lat');
u10lon = ncread(ncfile, 'lon');

ncfile='v10_1979-2020_ERA5.nc';

v10=ncread(ncfile,'v10');
v10lat = ncread(ncfile, 'lat');
v10lon = ncread(ncfile, 'lon');


%% importing ice drift speed data - 인공위성(nsidc) / 단위 (cm/s) %%

ncfile='uvel_1979-2020_NSIDC.nc';

un=ncread(ncfile,'u');
unlat = ncread(ncfile, 'lat');
unlon = ncread(ncfile, 'lon');

ncfile='vvel_1979-2020_NSIDC.nc';

vn=ncread(ncfile,'v');
vnlat = ncread(ncfile, 'lat');
vnlon = ncread(ncfile, 'lon');

%% importing ice drift speed data - 수치모델(cesm) %%

ncfile='uvel_1979-2020_CESM.nc';

uc=ncread(ncfile,'uvel');
uclat = ncread(ncfile, 'lat');
uclon = ncread(ncfile, 'lon');

ncfile='vvel_1979-2020_CESM.nc';

vc=ncread(ncfile,'vvel');
vclat = ncread(ncfile, 'lat');
vclon = ncread(ncfile, 'lon');

%% finding out grid location - u10,v10 %%

idx_extra = find(cum_doy > 15330); 
cum_doy(idx_extra) = [];
Vi(idx_extra) = [];
Ui(idx_extra) = [];
year(idx_extra) = [];
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

%% finding out grid location - nsidc %%

idx_extra = find(cum_doy > 15330); 
cum_doy(idx_extra) = [];
Vi(idx_extra) = [];
Ui(idx_extra) = [];
year(idx_extra) = [];
lat(idx_extra) = [];
lon(idx_extra) = [];

lon1 = lon;
idxlon = find(lon < 0);
lon1(idxlon) = lon(idxlon)+360;

glat = NaN*zeros(size(lat));
glon = NaN*zeros(size(lon1));

un_up = NaN*zeros(size(Ui));

for i = 1:length(Ui)
    time = cum_doy(i);
    flat = lat(i);
    flon = lon1(i);
    
    if (isnan(flat)==1 | isnan(flon) ==1)
        glat(i) = NaN;
        glon(i) = NaN;
    else 
        glat(i) = max( find( abs(unlat-flat) == min(abs(unlat-flat)) ) );
        glon(i) = max( find( abs(unlon-flon) == min(abs(unlon-flon)) ) );
        
        if (isnan(time)~=1)
            un_up(i) = un(glon(i),glat(i),time);
        end
    end
    
end


vn_up = NaN*zeros(size(Vi));

for i = 1:length(Vi)
    time = cum_doy(i);
    flat = lat(i);
    flon = lon1(i);
    
    if (isnan(flat)==1 | isnan(flon) ==1)
        glat(i) = NaN;
        glon(i) = NaN;
    else 
        glat(i) = max( find( abs(vnlat-flat) == min(abs(vnlat-flat)) ) );
        glon(i) = max( find( abs(vnlon-flon) == min(abs(vnlon-flon)) ) );
        
        if (isnan(time)~=1)
            vn_up(i) = vn(glon(i),glat(i),time);
        end
    end
    
end

%% finding out grid location - cesm %%

idx_extra = find(cum_doy > 15330); 
cum_doy(idx_extra) = [];
Vi(idx_extra) = [];
Ui(idx_extra) = [];
year(idx_extra) = [];
lat(idx_extra) = [];
lon(idx_extra) = [];

lon1 = lon;
idxlon = find(lon < 0);
lon1(idxlon) = lon(idxlon)+360;

glat = NaN*zeros(size(lat));
glon = NaN*zeros(size(lon1));

uc_up = NaN*zeros(size(Ui));

for i = 1:length(Ui)
    time = cum_doy(i);
    flat = lat(i);
    flon = lon1(i);
    
    if (isnan(flat)==1 | isnan(flon) ==1)
        glat(i) = NaN;
        glon(i) = NaN;
    else 
        glat(i) = max( find( abs(uclat-flat) == min(abs(uclat-flat)) ) );
        glon(i) = max( find( abs(uclon-flon) == min(abs(uclon-flon)) ) );
        
        if (isnan(time)~=1)
            uc_up(i) = uc(glon(i),glat(i),time);
        end
    end
    
end


vc_up = NaN*zeros(size(Vi));

for i = 1:length(Vi)
    time = cum_doy(i);
    flat = lat(i);
    flon = lon1(i);
    
    if (isnan(flat)==1 | isnan(flon) ==1)
        glat(i) = NaN;
        glon(i) = NaN;
    else 
        glat(i) = max( find( abs(vclat-flat) == min(abs(vclat-flat)) ) );
        glon(i) = max( find( abs(vclon-flon) == min(abs(vclon-flon)) ) );
        
        if (isnan(time)~=1)
            vc_up(i) = vc(glon(i),glat(i),time);
        end
    end
    
end

%% annual mean - 부표 %%

yearx = year(1:end-1);
yeary = 1979:1:2020;

tt=[];
yy=[];
for j = 1979:1:2020
    idxy = find(yearx == j);
    vlct= sqrt(Ui.^2 + Vi.^2);
    vlct_y=vlct(idxy);
    meanspeed= nanmean(vlct_y); 
    tt=[tt; meanspeed];
end

%% annual mean - u10, v10 : wind data %%

yearx = year(1:end-1);
yeary = 1979:1:2020;

uv10=[];
for j = 1979:1:2020
    idxy = find(yearx == j);
    wd_10= sqrt(u10_up.^2 + v10_up.^2);
    wd_10_y = wd_10(idxy);
    wd_10_m= nanmean(wd_10_y); 
    uv10=[uv10; wd_10_m];
end

%% annual mean - nsidc : ice drift speed data %%

yearx = year(1:end-1);
yeary = 1979:1:2020;

nn=[];
for j = 1979:1:2020
    idxy = find(yearx == j);
    oc_n = sqrt((un_up.*0.01).^2+(vn_up.*0.01).^2);
    oc_n_y = oc_n(idxy);
    oc_n_m= nanmean(oc_n_y); 
    nn=[nn; oc_n_m];
end

%% annual mean - cesm : ice drift speed data %%

yearx = year(1:end-1);
yeary = 1979:1:2020;

cc=[];
for j = 1979:1:2020
    idxy = find(yearx == j);
    oc_c = sqrt(uc_up.^2+vc_up.^2);
    oc_c_y = oc_c(idxy);
    oc_c_m= nanmean(oc_c_y); 
    cc=[cc; oc_c_m];
end

%% plotting all %% 

figure(1)
plot(yeary, tt, '-b')
hold on; grid on;
plot(yeary, cc, 'm')
plot(yeary, nn, '-r')
hold off;

xlabel('year')
ylabel('velocity')
xticks([1979:6:2021])
xlim([1978 2021])
legend('buoys','CESM', 'NSIDC')
title('Annual-mean speed (m/s)')

%% fractional ice drift speed %%

wd_10= sqrt(u10_up.^2 + v10_up.^2); % wind : u10, v10
vlct= sqrt(Ui.^2 + Vi.^2); % ice drift speed : buoys
oc_n = sqrt((un_up.*0.01).^2+(vn_up.*0.01).^2); % ice drift speed : nsidc
oc_c = sqrt(uc_up.^2+vc_up.^2); % ice drift speed : cesm

nan_oc_n = find(oc_n == 0);
oc_n(nan_oc_n) = NaN;

% buoys

yearx = year(1:end-1);
yeary = 1979:1:2020;

byf=[];
for j = 1979:1:2020
    idxy = find(yearx == j);
    byfrac = (vlct(idxy) ./ wd_10(idxy)) .* 100;
    byfrac_y_m= nanmean(byfrac); 
    byf=[byf; byfrac_y_m];
end

% nsidc

yearx = year(1:end-1);
yeary = 1979:1:2020;

nf=[];
for j = 1979:1:2020
    idxy = find(yearx == j);
    nfrac = (oc_n(idxy) ./ wd_10(idxy)) .* 100;
    nfrac_y_m= nanmean(nfrac); 
    nf=[nf; nfrac_y_m];
end

% cesm

yearx = year(1:end-1);
yeary = 1979:1:2020;

cf=[];
for j = 1979:1:2020
    idxy = find(yearx == j);
    cfrac = (oc_c(idxy) ./ wd_10(idxy)) .* 100;
    cfrac_y_m= nanmean(cfrac); 
    cf=[cf; cfrac_y_m];
end

%% plotting fractional ice drift speed %%

figure(2)
plot(yeary, byf, '-b')
hold on; grid on;
plot(yeary, cf, 'm')
plot(yeary, nf, '-r')
hold off;

xlabel('year')
ylabel('percentage')
xticks([1979:6:2021])
xlim([1978 2021])
legend('buoys','CESM', 'NSIDC')
title('fractional ice drift speed (%)')


%% if1 wind speed >= 5 m/s %%

ab = find(wd_10 >= 5);
wd_ab = wd_10(ab);

% buoys

yearxx = year(ab);
yeary = 1979:1:2020;

byf_u5=[];
for j = 1979:1:2020
    idxy = find(yearxx == j);
    byfrac5 = (vlct(idxy) ./ wd_ab(idxy)) .* 100;
    byfrac_y_m5= nanmean(byfrac5); 
    byf_u5=[byf_u5; byfrac_y_m5];
end

% nsidc

yearxx = year(ab);
yeary = 1979:1:2020;

nf_u5=[];
for j = 1979:1:2020
    idxy = find(yearxx == j);
    nfrac5 = (oc_n(idxy) ./ wd_ab(idxy)) .* 100;
    nfrac_y_m5= nanmean(nfrac5); 
    nf_u5=[nf_u5; nfrac_y_m5];
end

% cesm

yearxx = year(ab);
yeary = 1979:1:2020;

cf_u5=[];
for j = 1979:1:2020
    idxy = find(yearxx == j);
    cfrac5 = (oc_c(idxy) ./ wd_ab(idxy)) .* 100;
    cfrac_y_m5= nanmean(cfrac5); 
    cf_u5=[cf_u5; cfrac_y_m5];
end

% plotting

figure(3)
plot(yeary, byf_u5, '-b')
hold on; grid on;
plot(yeary, cf_u5, 'm')
plot(yeary, nf_u5, '-r')
hold off;

xlabel('year')
ylabel('percentage')
xticks([1979:6:2021])
xlim([1978 2021])
legend('buoys','CESM', 'NSIDC')
title('fractional ice drift speed, wdspeed>=5m/s (%)')

%% if2 wind speed <= 5 m/s %%

bl = find(wd_10 <= 5);
wd_bl = wd_10(bl);

% buoys

yearxy = year(bl);
yeary = 1979:1:2020;

byf_b5=[];
for j = 1979:1:2020
    idxy = find(yearxy == j);
    byfrac5 = (vlct(idxy) ./ wd_bl(idxy)) .* 100;
    byfrac_y_m5= nanmean(byfrac5); 
    byf_b5=[byf_b5; byfrac_y_m5];
end

% nsidc

yearxy = year(bl);
yeary = 1979:1:2020;

nf_b5=[];
for j = 1979:1:2020
    idxy = find(yearxy == j);
    nfrac5 = (oc_n(idxy) ./ wd_bl(idxy)) .* 100;
    nfrac_y_m5= nanmean(nfrac5); 
    nf_b5=[nf_b5; nfrac_y_m5];
end

% cesm

yearxy = year(bl);
yeary = 1979:1:2020;

cf_b5=[];
for j = 1979:1:2020
    idxy = find(yearxy == j);
    cfrac5 = (oc_c(idxy) ./ wd_bl(idxy)) .* 100;
    cfrac_y_m5= nanmean(cfrac5); 
    cf_b5=[cf_b5; cfrac_y_m5];
end

% plotting

figure(4)
plot(yeary, byf_b5, '-b')
hold on; grid on;
plot(yeary, cf_b5, 'm')
plot(yeary, nf_b5, '-r')
hold off;

xlabel('year')
ylabel('percentage')
xticks([1979:6:2021])
xlim([1978 2021])
legend('buoys','CESM', 'NSIDC')
title('fractional ice drift speed, wdspeed<=5m/s (%)')