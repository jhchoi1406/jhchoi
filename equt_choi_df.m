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

figure;plot(hs*100,igrate_t(:,:,1),hs*100,igrate_t(:,:,2),hs*100,igrate_t(:,:,3),'linewidth',2)
newcolors = {'r','k','b'};
colororder(newcolors)
grid on;
xlabel('Snow depth (cm)','fontsize',13)
ylabel('Ice growth rate (cm / month)','fontsize',13)
title('sensitivity of ice growth to the thickness of snow and ice','fontsize',13)
legend('1.0 m','2.0 m','3.0 m')

%%      
       
        
        
        
        
        
        
        
        
        