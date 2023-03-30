clear;

%% importing data %%

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
evap00 = ncread(ncfile, 'evap');
uvel0 = ncread(ncfile, 'uvel');
vvel0 = ncread(ncfile, 'vvel');
Qref0 = ncread(ncfile, 'Qref');
Tsfc00 = ncread(ncfile, 'Tsfc');
uatm0 = ncread(ncfile, 'uatm');
vatm0 = ncread(ncfile, 'vatm');
flat0 = ncread(ncfile, 'flat');

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
evap0 = evap00(idx_lon,idx_lat,:);
uvel = uvel0(idx_lon,idx_lat,:);
vvel = vvel0(idx_lon,idx_lat,:);
Qref = Qref0(idx_lon,idx_lat,:);
Tsfc0 = Tsfc00(idx_lon,idx_lat,:);
uatm = uatm0(idx_lon,idx_lat,:);
vatm = vatm0(idx_lon,idx_lat,:);
flat = flat0(idx_lon,idx_lat,:);

clear ab_hum
clear s_hum

Tsfc = Tsfc0+273.15 ; %단위 K

ab_hum = (6.1078.*exp(17.2693882.*(Tsfc-273.16)./(Tsfc-35.86))).*100; %단위 Pa
s_hum = (0.622.*ab_hum)./(10^5 - 0.378.*ab_hum) ;  

Hl = 25.04*10^5 * 2.14*10^-3 * sqrt(uatm.^2+vatm.^2) ...
   .* (s_hum- Qref/1000) ; 


idx_a = nan*zeros(size(Hl));
idx_b = nan*zeros(length(lon), length(lat));

for i = 1:length(aice)
    idx9 = find(aice(:,:,i) >= 0.9);
    idx_b(idx9) = 1;
    idx_a(:,:,i) = idx_b;
end

%% simulated latent heat flux - flat %%

flat9 = idx_a.*flat;
for i = 1:length(time)/12
    flat456(:,:,i) = (flat9(:,:,4+12*(i-1))+flat9(:,:,5+12*(i-1))+flat9(:,:,6+12*(i-1)))/3;
end

flat_clim = nanmean(flat456,3);

%% calculated latent heat flux - Hl9 %%

Hl9 = idx_a.*Hl;
for i = 1:length(time)/12
    Hl456(:,:,i) = (Hl9(:,:,4+12*(i-1))+Hl9(:,:,5+12*(i-1))+Hl9(:,:,6+12*(i-1)))/3;
end

Hl_clim = nanmean(Hl456,3);

%% plotting %%

cmin = -20 ;
cmax = -20;

mm = [-20:0.5:20];

tmp = jet(23);
jet1 = tmp;
jet1(12,:) = 1;
MYMAPs = jet1;

m_proj('stereographic', 'lat', 90, 'radius', 25);

figure(1)
m_contourf(lon, lat, -flat_clim', mm); 
caxis([-10.5 10.5]);
m_grid('linewi',2,'tickdir','out');

%axis square
%axis off
%shading flat

colormap(MYMAPs);
colorbar('ytick', 2*mm, 'fontsize', 14)
m_coast('patch', [0.7 0.7 0.7], 'edgecolor', 'k'); hold off;

figure(2)
m_contourf(lon, lat, Hl_clim', mm); hold on;
caxis([-10.5 10.5]);
m_grid('linewi',2,'tickdir','out');

axis square
axis off
shading flat

colormap(MYMAPs);
colorbar('ytick', 2*mm, 'fontsize', 14)
m_coast('patch', [0.7 0.7 0.7], 'edgecolor', 'k'); hold off;
































