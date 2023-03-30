clear; clc; close all;

ncfile = 'sf_1979-2020_ERA5.nc';
ncdisp(ncfile);
sf = ncread(ncfile, 'sf');
lat = ncread(ncfile, 'lat');
lon = ncread(ncfile, 'lon');

fid = 'N_09_extent_v3.0.csv';
Ncsv = readtable(fid);
extent = Ncsv.extent(1:42);

%% sic data
ncfile3 = 'sic_1979-2020_NSIDC.nc';
%ncdisp(ncfile3);
sic = ncread(ncfile3,'cdr_seaice_conc');
%sic = sic0(:,:,365*(2001-1979)+1:end);

%해빙밀도 연평균
clear sic_yy;
for i = 1:length(sic)/365
    sic_yy(:,:,i) = nanmean(sic(:,:,365*(i-1)+1:365*i),3);
end

%% 6~8월 sea ice 에 내리는 누적 snowfall 과 9월 sea ice extent 의 상관관계 %%
yy = 1979:2020;
for i = 1:length(yy)
    sf_S(:,:,i) = (24*1000*3.3)*nansum(sf(:,:,365*(i-1)+152:365*(i-1)+243),3);
    sf_58(:,:,i) = (24*1000*3.3)*nansum(sf(:,:,365*(i-1)+121:365*(i-1)+243),3);
    sf_78(:,:,i) = (24*1000*3.3)*nansum(sf(:,:,365*(i-1)+182:365*(i-1)+243),3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% area weighting for the seasonal mean data
weight3d = ones(size(sf_S));
for j = 1:length(lat)
    weight3d(:,j,:) = cosd(lat(j));
end
sf_S_weight = sf_S.*weight3d;

%해빙밀도 15%이상 지역 뽑기 
clear idx3d;
idx3d = nan*zeros(size(sf_S));
idx2d = nan*zeros(length(lon), length(lat));

for i = 1:length(yy)
    idx15 = find(sic_yy(:,:,i)>=0.15);
    idx2d(idx15) = 1;
    idx3d(:,:,i) = idx2d;
end

%weighting한 연평균 알베도 중 해빙밀도가 15%이상 되는 지역 
%sf_S_yy_15 = sf_S_yy.*idx3d;
sf_S_weight_15 = sf_S_weight.*idx3d;

%weight average
weight3d_15 = weight3d.*idx3d;
sf_S_weight_15_avg = squeeze(nansum(nansum(sf_S_weight_15,2),1))./squeeze(nansum(nansum(weight3d_15,2),1));

%% Plotting %%
figure(1)
yyaxis left
plot(yy,sf_S_weight_15_avg,'LineWidth',1.0)
ylabel('Jun-Aug accumulated Snowfall (mm)','fontsize',13);
xlim([1978 2021])
grid on;

yyaxis right
plot(yy, extent,'LineWidth',1.0)
ylabel('September sea ice extent (10^{6} km^{2})','fontsize',13)
xlabel('year','fontsize',13)

%% Correlation coefficients %%
corr(squeeze(nanmean(nanmean(sf_S,2),1)), extent)
corr(detrend(squeeze(nanmean(nanmean(sf_S,2),1))), detrend(extent)) % 6~8월
corr(squeeze(nanmean(nanmean(sf_58,2),1)), extent)
corr(detrend(squeeze(nanmean(nanmean(sf_58,2),1))), detrend(extent)) % 5~8월
corr(squeeze(nanmean(nanmean(sf_78,2),1)), extent)
corr(detrend(squeeze(nanmean(nanmean(sf_78,2),1))), detrend(extent)) %7~8월
