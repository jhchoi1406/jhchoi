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

% ㅂㅏ람 부표, sic, 다 써서 평균 내기
% 부표 지나가는 지점에 66만 시공간 좌표 
% 그림을 잘 그리는 방법 공부해보기
% xticks 크게 하기
% 1. 2012 년 이후로 빨라졌는데 해빙 밀도가 높은 곳에서 빨라졌는지 낮은 곳에서 빨라졌는지를 conc 90퍼 이상인것만 ㄴ뽑아서 똑같은그림을 그리면 가능
% 숙제 2. xticks yticks 크게 하기
% 해빙 밀도가 높은 곳에서 돌아다녔는지 histogram 을 통해 알아보기
% hist(sic_up) 0.1과 0.8사이는 데이터가 없느데, 해빙 밀도가 0.1 과 0.8 사이의 지점이 에어리어가 적어서
% 그렇다. 해빙이 거의 없거나 있는 곳은 90퍼 넘어가는 지역이 많은데, 왜냐면 인공위성으로 관측할 수 있는 해상도가 높지 않아서
% 그렇다. 인공위성으로 볼 수 있는곳은 100km, 150km 지점에서 가능해서 정교하지 못함. 100 150 사이는 꽉 차있는데,
% 10키로 더 자세하게 들어가면 해빙 밀도가 높은것도 있고 낮은것도 있고 섞여있는데 근데 인공위성이 그렇게까지 자세히는 관측을 못해서
% 거의 없ㅇ거나 1에 가깝거나 이 정도밖에 못함. 그래서 histgram이 양극화가 되어 있는 것이다. 해빙 밀도를 자세히 알아야지
% 해빙 속력이 빨라지는 이유를 알 수 잇는데 인공위성 해상도가 낮아서 100-150 큰 정도의 데이터만으로는 sia ice
% speed가 늘어나는 이유를 알기 어렵다.


