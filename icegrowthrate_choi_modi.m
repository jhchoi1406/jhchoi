%% Importing data %%
clear;

ncfile1 = 'climsnow_JRA.cice.h.2000-2020.nc';
ncdisp(ncfile1)
ncfile2 = 'hist_JRA.cice.h.2000-2020.nc';
ncdisp(ncfile2)

lat = ncread(ncfile1, 'lat');
lon = ncread(ncfile1, 'lon');
time = ncread(ncfile1, 'time');
Tsfc = ncread(ncfile1, 'Tsfc'); % snow/ice surface temperature
aice = ncread(ncfile1, 'aice'); % ice area (aggregate)
albice = ncread(ncfile1, 'albice'); % bare ice albedo
albsno = ncread(ncfile1, 'albsno'); % snow albedo
hi = ncread(ncfile1, 'hi'); % grid cell mean ice thickness
hs = ncread(ncfile1, 'hs'); % grid cell mean snow thickness
snow = ncread(ncfile1, 'snow'); % snowfall rate (cpl)

%% calculating sea ice growth %%

cp = 1005 ;
cd =  0.0013;
ki = 2.04;
ks = 0.31;
lf = 3.340*10^5;
ls = 2.265*10^6;
tf = 271.3;
rhoi = 930;
rhoa = 1.275;
sigma = 5.67*10^-8;
u = 2.56;
ta = 249.85;
qa = 0.57;
flw = 182.1;
p = 10^5;
hi = [1.0 1.5 2.0];
hs = [0:0.01:0.4];

for i = 1:length(hs)
    for j = 1:length(hi)
       syms ts
       eqa = sigma*ts^4 - flw + rhoa*cp*cd*u*(ts-ta) - (ki*ks*(tf-ts))/(ki*hs(i)+ks*hi(j)) == 0;
       equt(:,i,j) = vpasolve(eqa,ts);
       equtt = equt(2,:,:);
       uplwr = sigma*equtt.^4;
       shf = rhoa*cp*cd*u*(equtt-ta);
       equt_ts(:,i,j) = ki*ks*(tf-equtt(:,i,j))/(ki*hs(i)+ks*hi(j));
       igrate(:,i,j) = equt_ts(:,i,j)/(rhoi*lf);
       igrate_t(:,i,j) = igrate(:,i,j)*60*60*24*30*100;
    end
end

%% sea ice growth plotting %%

figure(1)
plot(hs*100,igrate_t(:,:,1),hs*100,igrate_t(:,:,2),hs*100,igrate_t(:,:,3),'linewidth',2)
newcolors = {'r','k','b'};
colororder(newcolors)
hold on; grid on;

%% winter(11-3) snow thickness % ice growth rate %%    

ncfile = 'GIAF_JRA.cice.h.1979-2020.nc';
ncdisp(ncfile)

lat0 = ncread(ncfile, 'lat'); 
up70 = find(lat0>=70);
lat = lat0(up70);
hi0 = ncread(ncfile, 'hi'); 
hi = hi0(:,up70,:);
hs0 = ncread(ncfile, 'hs');
hs = hs0(:,up70,:);
time = ncread(ncfile, 'time');
lon = ncread(ncfile, 'lon'); 
aice00 = ncread(ncfile, 'aice'); 
aice0 = aice00(:,up70,:);
snow0 = ncread(ncfile, 'snow');
snow = snow0(:,up70,:);
flwdn0 = ncread(ncfile, 'flwdn');
flwdn = flwdn0(:,up70,:);

% 1) 11~3 averaged snow depth
clear avsd
clear aice
years = 1980:1:2020;
for i = 1:length(years)
    avsd(:,:,i) = mean(hs(:,:,12*(i-1)+11:12*i+3), 3, 'omitnan');
    aice(:,:,i) = mean(aice0(:,:,12*(i-1)+11:12*i+3), 3, 'omitnan');
    wtflwdn(:,:,i) = mean(flwdn(:,:,12*(i-1)+11:12*i+3), 3, 'omitnan');
    wtsnow(:,:,i) = mean(snow(:,:,12*(i-1)+11:12*i+3), 3, 'omitnan');    
end
    
% 2) 3 minus 11 sea ice thickness = monthly growth rate
clear misit;
clear hi_yy;
for i = 1:length(years)
    misit(:,:,i) = (hi(:,:,(15+12*(i-1))) - hi(:,:,(11+12*(i-1))))/4;
    hi_yy(:,:,i) = nanmean(hi(:,:,(12*(i-1)+11:12*(i-1)+15)),3);
end


% weighted average %
ap_cosd = ones(size(avsd));
for i = 1:length(lat)
    ap_cosd(:,i,:) = cosd(lat(i));
end

idx_a0 = nan*zeros(size(avsd));
idx_b = nan*zeros(length(lon), length(lat));

yy = length(time)/12;
for i = 1:yy-1
    idx15 = find(aice(:,:,i) >= 0.15);
    idx_b(idx15) = 1;
    idx_a(:,:,i) = idx_b;
end

hs_cosd = avsd.*ap_cosd.*idx_a ;
hi_cosd = misit.*ap_cosd.*idx_a ;
wtflwdn_cosd = wtflwdn.*ap_cosd.*idx_a ;
wtsnow_cosd = wtsnow.*ap_cosd.*idx_a ;

wgtavg_hs = squeeze(sum(sum(hs_cosd,1,'omitnan'),2,'omitnan'))./squeeze(sum(sum(ap_cosd.*idx_a,1,'omitnan'),2,'omitnan'));
wgtavg_hi = squeeze(sum(sum(hi_cosd,1,'omitnan'),2,'omitnan'))./squeeze(sum(sum(ap_cosd.*idx_a,1,'omitnan'),2,'omitnan'));
wgtavg_wtflwdn = squeeze(sum(sum(wtflwdn_cosd,1,'omitnan'),2,'omitnan'))./squeeze(sum(sum(ap_cosd.*idx_a,1,'omitnan'),2,'omitnan'));
wgtavg_wtsnow = squeeze(sum(sum(wtsnow_cosd,1,'omitnan'),2,'omitnan'))./squeeze(sum(sum(ap_cosd.*idx_a,1,'omitnan'),2,'omitnan'));

scatter(wgtavg_hs*100,wgtavg_hi*100,'m')
slr = fitlm(wgtavg_hs*100,wgtavg_hi*100);
plot(slr)
xlabel('Snow depth (cm)','fontsize',13)
ylabel('Ice growth rate (cm / month)','fontsize',13)
title('sensitivity of ice growth to the thickness of snow and ice','fontsize',13)
legend('1.0 m','2.0 m','3.0 m')
hold off;      

wths_mn = nanmean(nanmean(nanmean(avsd)))

snow_flwdn = fitlm(wgtavg_wtflwdn,wgtavg_wtsnow/100);

figure(2)
plot(snow_flwdn)
xlabel('down longwave flux (W/m2)','fontsize',13)
ylabel('snowfall rate (m/day)','fontsize',13)
title('correlation of snowfall and downward longwave radiation','fontsize',13)
grid on;

%% 10 sea ice area %%

yr = 1979:1:2020;
for i = 1:length(yr)
    month10(:,:,i) = aice0(:,:,10+12*(i-1));
end

r = 6400; %km
for i = 1:length(lat)
    grid_area(i) = (2*pi*r*cosd(lat(i))/length(lon)) * (pi*r/length(lat));
end

for i = 1:length(yr)
    for j = 1:length(lat)
        lat_m0(:,j,i) = sum(month10(:,j,i)*grid_area(j),'omitnan');
        lat_m = squeeze(squeeze(sum(lat_m0,2,'omitnan')));
    end
end

figure(3)
plot(yr,lat_m)
xlabel('years','fontsize',13)
ylabel('sea ice area','fontsize',13)
title('sea ice area in October','fontsize',13)
grid on;

        

        












