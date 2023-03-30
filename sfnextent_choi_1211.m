%% Importing sf & sea ice extent data %%
clear;

ncfile = 'sf_1979-2020_ERA5.nc';
lat = ncread(ncfile, 'lat');
lon = ncread(ncfile, 'lon');
time = ncread(ncfile, 'time');
sf = ncread(ncfile,'sf');

csvfile = readmatrix('N_09_extent_v3.0.csv');
year0 = csvfile(:,1);
year = year0(1:end-1);
extent0= csvfile(:,5);
extent = extent0(1:end-1);


%% calculating snowfall accumulation %%

yy = length(time)/365;
for i = 1:yy
    acsf_s(:,:,i) = sum(sf(:,:,365*(i-1)+152:365*(i-1)+243),3);
    acsf(i) = nanmean(nanmean(acsf_s(:,:,i)));
end

%% plotting snowfall accumulation and sea ice extent by using yyaxis %%

yyaxis left
plot(year,acsf)
yyaxis right
plot(year,extent)

%% correlation coefficient between snowfall accumulation(6,7,8) and sea ice extent(9) %% 

corrcoef(acsf,extent)

%% month = 5,6,7,8 %%

yy = length(time)/365;
for i = 1:yy
    acsf_s_58(:,:,i) = sum(sf(:,:,365*(i-1)+121:365*(i-1)+243),3);
    acsf_58(i) = nanmean(nanmean(acsf_s_58(:,:,i)));
end

corrcoef(acsf_58,extent)

%% month = 7,8 %%

yy = length(time)/365;
for i = 1:yy
    acsf_s_78(:,:,i) = sum(sf(:,:,365*(i-1)+182:365*(i-1)+243),3);
    acsf_78(i) = nanmean(nanmean(acsf_s_78(:,:,i)));
end

corrcoef(acsf_78,extent)











