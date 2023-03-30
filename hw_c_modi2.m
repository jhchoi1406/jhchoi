clear;

ncfile = 'GIAF_JRA.cice.h.1979-2020.nc';
ncdisp(ncfile)

lat0 = ncread(ncfile, 'lat'); 
lon0 = ncread(ncfile, 'lon'); 
hi0 = ncread(ncfile, 'hi'); 
hs0 = ncread(ncfile, 'hs');
time = ncread(ncfile, 'time');
aice0 = ncread(ncfile, 'aice'); 
snow0 = ncread(ncfile, 'snow');
flwdn0 = ncread(ncfile, 'flwdn');
snoice0 = ncread(ncfile, 'snoice');
melts0 = ncread(ncfile, 'melts');
ice_present0 = ncread(ncfile, 'ice_present');
evap0 = ncread(ncfile, 'evap');

LAT1 = 70;
LON1 = 0; LON2 = 360;

idx_lon = find(lon0>=LON1 & lon0<=LON2);
idx_lat = find(lat0>=LAT1);

lon = lon0(idx_lon);
lat = lat0(idx_lat);
hs = hs0(idx_lon,idx_lat,:);
aice = aice0(idx_lon,idx_lat,:);
snow = snow0(idx_lon,idx_lat,:);
snoice = snoice0(idx_lon,idx_lat,:);
melts = melts0(idx_lon,idx_lat,:);
ice_prst = ice_present0(idx_lon,idx_lat,:);
evap = evap0(idx_lon,idx_lat,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

snowA = aice.*(snow - snoice - melts + evap);

ap_cosd = ones(size(snow));
for i = 1:length(lat)
    ap_cosd(:,i,:) = cosd(lat(i));  
end

idxNaN = find(isnan(snowA) == 1);
ap_cosd(idxNaN) = nan;

snowA_cosd = snowA.*ap_cosd;
snowA_cosd(snowA_cosd<0)=0;

avg_snowA = squeeze(nansum(nansum(snowA_cosd,1),2))./squeeze(nansum(nansum(ap_cosd,1),2));

yy = length(time)/12;
for i = 1:yy-1
    sf_ent(i) = nansum(avg_snowA(9+12*(i-1):15+12*(i-1)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Q1

yy_m = 1980:1:2000;

m9 = nanmean(aice(:,:,9:12:9+12*20),3) ;
m10 = nanmean(aice(:,:,10:12:10+12*20),3) ; 
m11 = nanmean(aice(:,:,11:12:11+12*20),3) ; 

aice9 = aice;
aice10 = aice;
aice11 = aice;
aice_et = aice;

for i = 1:yy-1;
    aice9(:,:,9+12*(i-1))= m9(:,:);
    aice10(:,:,10+12*(i-1))= m10(:,:);
    aice11(:,:,11+12*(i-1))= m11(:,:);
    
    aice_et(:,:,9+12*(i-1))= m9(:,:);
    aice_et(:,:,10+12*(i-1))= m10(:,:);
    aice_et(:,:,11+12*(i-1))= m11(:,:);
end

snowA9 = aice9.*(snow - snoice - melts + evap);
snowA10 = aice10.*(snow - snoice - melts + evap);
snowA11 = aice11.*(snow - snoice - melts + evap);
snowA_et = aice_et.*(snow - snoice - melts + evap);

snowA_cosd9 = snowA9.*ap_cosd;
snowA_cosd10 = snowA10.*ap_cosd;
snowA_cosd11 = snowA11.*ap_cosd;
snowA_cosd_et = snowA_et.*ap_cosd;

snowA_cosd9(snowA_cosd9<0)=0;
snowA_cosd10(snowA_cosd10<0)=0;
snowA_cosd11(snowA_cosd11<0)=0;
snowA_cosd_et(snowA_cosd_et<0)=0;

avg_m9 = squeeze(nansum(nansum(snowA_cosd9,1),2))./squeeze(nansum(nansum(ap_cosd,1),2));
avg_m10 = squeeze(nansum(nansum(snowA_cosd10,1),2))./squeeze(nansum(nansum(ap_cosd,1),2));
avg_m11 = squeeze(nansum(nansum(snowA_cosd11,1),2))./squeeze(nansum(nansum(ap_cosd,1),2));
avg_et = squeeze(nansum(nansum(snowA_cosd_et,1),2))./squeeze(nansum(nansum(ap_cosd,1),2));

for i = 1:yy-1
    sf_avg_m9(i) = nansum(avg_m9(9+12*(i-1):15+12*(i-1)));
    sf_avg_m10(i) = nansum(avg_m10(9+12*(i-1):15+12*(i-1)));
    sf_avg_m11(i) = nansum(avg_m11(9+12*(i-1):15+12*(i-1)));
    sf_avg_et(i) = nansum(avg_et(9+12*(i-1):15+12*(i-1)));
end

years = 1980:1:2020;

f1 = figure;
f1.Position = [100 100 1400 600];
subplot(1,4,1)
plot(years,90*sf_ent)
hold on; grid on;
plot(years,90*sf_avg_m9)
hold off;
title('9월 비교','fontsize',18)
xlabel('years','fontsize',15)
ylabel('9-3 accumulated snow fall','fontsize',15)
legend('no edit','edited sep data')

subplot(1,4,2)
plot(years,90*sf_ent)
hold on; grid on;
plot(years,90*sf_avg_m10)
hold off;
title('10월 비교','fontsize',18)
xlabel('years','fontsize',15)
ylabel('9-3 accumulated snow fall','fontsize',15)
legend('no edit','edited oct data')

subplot(1,4,3)
plot(years,90*sf_ent)
hold on; grid on;
plot(years,90*sf_avg_m11)
hold off;
title('11월 비교','fontsize',18)
xlabel('years','fontsize',15)
ylabel('9-3 accumulated snow fall','fontsize',15)
legend('no edit','edited oct data')

subplot(1,4,4)
plot(years,90*sf_ent)
hold on; grid on;
plot(years,90*sf_avg_et)
hold off;
title('가을철 전체 비교','fontsize',18)
xlabel('years','fontsize',15)
ylabel('9-3 accumulated snow fall','fontsize',15)
legend('no edit','edited oct data')

diff_m9 = mean(abs(90*(sf_ent - sf_avg_m9)))
diff_m10 = mean(abs(90*(sf_ent - sf_avg_m10)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Q2
% if1) because of decreasing snow quantity

clear snow_yy
clear snow_12

yy_m = 1980:1:2000;
for i = 1:length(yy_m)
    snow_yy = nanmean(snow(:,:,9+12*(i-1):15+12*(i-1)),3) ;
end

snow_12 = snow;

for i = 12*33+9:12*33+15;
   snow_12(:,:,i)= snow_yy(:,:);
end

snowA_sw = aice.*(snow_12 - snoice - melts + evap);

ap_cosd = ones(size(snow));
for i = 1:length(lat)
    ap_cosd(:,i,:) = cosd(lat(i));  
end

idxNaN = find(isnan(snowA_sw) == 1);
ap_cosd(idxNaN) = nan;

snowA_cosd_sw = snowA_sw.*ap_cosd;
snowA_cosd_sw(snowA_cosd_sw<0)=0;

avg_snowA_sw = squeeze(nansum(nansum(snowA_cosd_sw,1),2))./squeeze(nansum(nansum(ap_cosd,1),2));

yy = length(time)/12;
for i = 1:yy-1
    sf_ent_sw(i) = nansum(avg_snowA_sw(9+12*(i-1):15+12*(i-1)));
end

% if2) because of low aice in autumn

yy_m = 1980:1:2000;
for i = 1:length(yy_m)
    aice_yy = nanmean(aice(:,:,9+12*(i-1):15+12*(i-1)),3) ;
end

aice_12 = aice;

for i = 12*33+9:12*33+15;
   aice_12(:,:,i)= aice_yy(:,:);
end

snowA_ai = aice_12.*(snow - snoice - melts + evap);

ap_cosd = ones(size(snow));
for i = 1:length(lat)
    ap_cosd(:,i,:) = cosd(lat(i));  
end

idxNaN = find(isnan(snowA_ai) == 1);
ap_cosd(idxNaN) = nan;

snowA_cosd_ai = snowA_ai.*ap_cosd;
snowA_cosd_ai(snowA_cosd_ai<0)=0;

avg_snowA_ai = squeeze(nansum(nansum(snowA_cosd_ai,1),2))./squeeze(nansum(nansum(ap_cosd,1),2));

yy = length(time)/12;
for i = 1:yy-1
    sf_ent_ai(i) = nansum(avg_snowA_ai(9+12*(i-1):15+12*(i-1)));
end

figure(2)
plot(years, 90*sf_ent)
hold on; grid on;
plot(years, 90*sf_ent_sw)
hold off;
title('replace 2013 snow to averaged snow','fontsize',18)
xlabel('years','fontsize',15)
ylabel('9-3 accumulated snow fall','fontsize',15)
legend('no edit','edited snow data')

figure(3)
plot(years, 90*sf_ent)
hold on; grid on;
plot(years, 90*sf_ent_ai)
hold off;
title('replace 2013 aice to averaged aice','fontsize',18)
xlabel('years','fontsize',15)
ylabel('9-3 accumulated snow fall','fontsize',15)
legend('no edit','edited aice data')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






























% hs_cosd = hs.*ap_cosd;
% hs_133 = hs_cosd(:,:,12*33+3);
% for i = 1:
% avg = hs_cosd(:,:,12*20+3:12*41+3















