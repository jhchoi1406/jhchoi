%% Reading buoy data and calculating ice velocity %% 
clear 
% New data % 
load buoys_1979_2010.mat

ID0 = data(:,1); 
year0 = data(:,2); 
hour0 = data(:,3); 
doy0 = data(:,4); 
lat0 = data(:,5); 
lon0 = data(:,6); 

idx = find(hour0 == 0); 

ID = ID0(idx); 
year = year0(idx); 
hour = hour0(idx); 
doy = doy0(idx); 
lat = lat0(idx); 
lon = lon0(idx); 

idx_lon = find(lon < -180); 
lon(idx_lon) = NaN; 

idx_lat = find(lat<60); 
lat(idx_lat) = NaN; 

idx_jump1 = find( diff(doy) ~= 1);
idx_jump2 = find( diff(ID)~=0 );
idx_jump = find( diff(doy) ~= 1 | diff(ID)~=0 | diff(year)~=0 );

%% importing wind data %%

ncfile='u10_1979-2010_ERA5.nc';

u10=ncread(ncfile,'u10');
rlat = ncread(ncfile, 'lat');
rlon = ncread(ncfile, 'lon');

ncfile='v10_1979-2010_ERA5.nc';

v10=ncread(ncfile,'v10');
vlat = ncread(ncfile, 'lat');
vlon = ncread(ncfile, 'lon');

%% importing SIC data %%

nc_sic='sic_1979-2020_NSIDC.nc';

sic=ncread(nc_sic,'cdr_seaice_conc');
clat = ncread(nc_sic, 'lat');
clon = ncread(nc_sic, 'lon');

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

idx_doy = find(diff(doy) ~= 1);
Ui(idx_doy) = NaN;
Vi(idx_doy) = NaN;

% idxi = find(abs(Ui) > 100);

idxL = find(abs(Ui) > 1);
Ui(idxL) = NaN;

idxL = find(abs(Vi) > 1);
Vi(idxL) = NaN;

% for i=1:length(Ui_up)
%     if (lon_f(i)>0 & lon_b(i)<0)
%         Ui_up(i) = (r(i)/360).*(360-(lon_f(i) - lon_b(i)))/dt;
%     end
%     if (lon_f(i)<0 & lon_b(i)>0)
%         Ui_up(i) = (r(i)/360).*(360-(-lon_f(i) + lon_b(i)))/dt;
%     end
% end
%
% idxi1 = find(abs(Ui_up)>20); 

%% finding out time location %%

day = 365*ones(2010-1979+1,1);
day(2:4:end) = 366;
sumday0 = cumsum(day);
sumday = [0; sumday0];

time_loc = zeros(size(Ui));

for i = 1:length(time_loc)
    time_loc(i) = sumday(year(i)-1979+1)+doy(i);
end

idx_extra = find(time_loc > 11688); 
time_loc(idx_extra) = NaN;


%% finding out grid location %%

lon1 = lon;
idxlon = find(lon < 0);
lon1(idxlon) = lon(idxlon)+360;

glat = NaN*zeros(size(lat));
glon = NaN*zeros(size(lon1));

u10_up = NaN*zeros(size(Ui));

for i = 1:length(Ui)
    time = time_loc(i);
    flat = lat(i);
    flon = lon1(i);
    
    if (isnan(flat)==1 | isnan(flon) ==1)
        glat(i) = NaN;
        glon(i) = NaN;
    else 
        glat(i) = max( find( abs(rlat-flat) == min(abs(rlat-flat)) ) );
        glon(i) = max( find( abs(rlon-flon) == min(abs(rlon-flon)) ) );
        
        if (isnan(time)~=1)
            u10_up(i) = u10(glon(i),glat(i),time);
        end
    end
    
end


v10_up = NaN*zeros(size(Vi));

for i = 1:length(Vi)
    time = time_loc(i);
    flat = lat(i);
    flon = lon1(i);
    
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

for i = 1:length(Vi)
    time = time_loc(i);
    flat = lat(i);
    flon = lon1(i);
    
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

gg=[];
yearx = year(1:end-1)
yeary = 1979:1:2010;
for i = 1979:1:2010;
    idxy = find(yearx == i);
    sicmean = nanmean(sic_up(idxy));
    gg=[gg; sicmean];
end



    


%% wind speed & ocean speed %%

ocspeed = sqrt(Ui.^2 + Vi.^2);
wdspeed = sqrt(u10_up.^2+v10_up.^2);

idxw = find(wdspeed >= 1.0);

wdspeed_idxw = wdspeed(idxw);
ocspeed_idxw = ocspeed(idxw);

frac_speed = 100*(ocspeed_idxw./wdspeed_idxw);

wt = find(90>=doy);
sm = find(273>=doy & doy>=213);

wdspeed_wt = wdspeed(wt);
ocspeed_wt = ocspeed(wt);
wdspeed_sm = wdspeed(sm);
ocspeed_sm = ocspeed(sm);

ww=[];

yearx = year(1:end-1);
yeary = 1979:1:2010;

for i = 1979:1:2010
     idxy = find(yearx(idxw) == i);
     frac_speed0 = 100*((ocspeed_idxw) ./ (wdspeed_idxw));
     frac_speed = frac_speed0(idxy); 
     mean_frac_speed = nanmean(frac_speed);
     ww=[ww; mean_frac_speed];
end

%% plotting %% 

corr(u10_up, Ui, 'rows', 'complete');

figure(1); 
scatter(u10_up, Ui); 

figure(2);
ocspeed = sqrt(Ui.^2 + Vi.^2);
wdspeed = sqrt(u10_up.^2+v10_up.^2);

scatter(wdspeed, ocspeed);

figure(3);

yearx = year(1:end-1);
yeary = 1979:1:2010;

tt=[];
for i = 1979:1:2010 
    idxy = find(yearx == i);
    ocspeedtt= sqrt(Ui(idxy).^2 + Vi(idxy).^2);
    ocmeanspeed= nanmean(ocspeedtt); 
    tt=[tt; ocmeanspeed];
end

plot(yeary, tt);


    

