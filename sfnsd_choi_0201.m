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

%% entire %%

snowA = aice.*snow;

ap_cosd = ones(size(snow));
for i = 1:length(lat)
    ap_cosd(:,i,:) = cosd(lat(i));  
end

idxNaN = find(isnan(snowA) == 1);
ap_cosd(idxNaN) = nan;

snowA_cosd = snowA.*ap_cosd;
hs_cosd = hs.*ap_cosd;

avg_snowA = squeeze(nansum(nansum(snowA_cosd,1),2))./squeeze(nansum(nansum(ap_cosd,1),2));
avg_hs = squeeze(nansum(nansum(hs_cosd,1),2))./squeeze(nansum(nansum(ap_cosd,1),2));

clear sf93;
clear sd3;
yy = length(time)/12;
for i = 1:yy-1
    sf93(i) = sum(avg_snowA(9+12*(i-1):15+12*(i-1)));
    sd3(i) = avg_hs(15+12*(i-1));    
end
corr(sf93', sd3')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear sd3_2D_yy
clear sd3_2D_clim
yy = length(time)/12;
for i = 1:yy-1
    sd3_2D_yy(:,:,i) = nanmean(hs(:,:,15+12*(i-1)),3);    
end

sd3_2D_clim = nanmean(sd3_2D_yy,3);

clear sf93_2D_yy
clear sf93_2D_clim
yy = length(time)/12;
for i = 1:yy-1
    sf93_2D_yy(:,:,i) = mean(nansum(snowA(:,:,9+12*(i-1):15+12*(i-1)),3),3);
end

sf93_2D_clim = nanmean(sf93_2D_yy,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear sf_modi
clear sf_modi_clim
snow_a = aice.*(snow - snoice - melts + evap);
snow_a(snow_a<0)=0;

yy = length(time)/12;
for i = 1:yy-1
    sf_modi(:,:,i) = mean(nansum(snow_a(:,:,9+12*(i-1):15+12*(i-1)),3),3);
end

sf_modi_clim = nanmean(sf_modi,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear p_modi
clear p_clim
ice_prst(find(ice_prst<0.75))=0;
snow_p = ice_prst.*aice.*(snow - snoice - melts + evap);

for i = 1:yy-1
    p_modi(:,:,i) = mean(nansum(snow_p(:,:,9+12*(i-1):15+12*(i-1)),3),3);
end

p_clim = nanmean(p_modi,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plotting %%

years = 1980:1:2020;

figure(1)
yyaxis left
plot(years,sf93*30*3)
title('comparison of snow depth and accumulated snow fall in entire arctic area','fontsize',12)
xlabel('years','fontsize',13)
ylabel('accumulated snow fall during September-March','fontsize',13)
xlim([1975 2025])
grid on;

yyaxis right
plot(years,sd3*100)
ylabel('snow depth in March','fontsize',13)

cmin = -1.5;
cmax = 32.5;

mm = [-40:2:40];

tmp = jet(23);
%tmp = flipud(tmp0);
jet1 = tmp;
% jet1(11:12,:) = 1;
jet1(1:1,:) = 1;
MYMAPs = jet1;

m_proj('stereographic', 'lat', 90, 'radius', 25);

figure(2)
subplot(1,2,1)
m_contourf(lon, lat, 100*sd3_2D_clim', mm); hold on;
caxis([cmin cmax]);
m_grid('linewi',2,'tickdir','out');

axis square
axis off
shading flat

colormap(MYMAPs);
colorbar('ytick', 2*mm, 'fontsize', 14)
m_coast('patch', [0.7 0.7 0.7], 'edgecolor', 'k'); hold off;

subplot(1,2,2)
m_contourf(lon, lat, 90*sf93_2D_clim', mm); hold on;
caxis([cmin cmax]);
m_grid('linewi',2,'tickdir','out');

axis square
axis off
shading flat

colormap(MYMAPs);
colorbar('ytick', 2*mm, 'fontsize', 14)
m_coast('patch', [0.7 0.7 0.7], 'edgecolor', 'k'); hold off;

figure(3)
m_contourf(lon, lat, 90*sf_modi_clim', mm); hold on;
caxis([cmin cmax]);
m_grid('linewi',2,'tickdir','out');

axis square
axis off
shading flat

colormap(MYMAPs);
colorbar('ytick', 2*mm, 'fontsize', 14)
m_coast('patch', [0.7 0.7 0.7], 'edgecolor', 'k'); hold off;

figure(4)
m_contourf(lon, lat, 90*p_clim', mm); hold on;
caxis([cmin cmax]);
m_grid('linewi',2,'tickdir','out');

axis square
axis off
shading flat

colormap(MYMAPs);
colorbar('ytick', 2*mm, 'fontsize', 14)
m_coast('patch', [0.7 0.7 0.7], 'edgecolor', 'k'); hold off;
