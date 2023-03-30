# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 17:42:36 2022

@author: PC
"""

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy
import matplotlib.ticker as mticker
import matplotlib.colors as colors
import matplotlib as mpl
import os
import wget
from netCDF4 import Dataset
import sys
sys.path.append('C:/Users/PC/Desktop/jhc/03120320')
np.set_printoptions(threshold=sys.maxsize)

'''
mm=['01','02','03','04','05','06','07','08','09','10','11','12']
for i in range(1992,1995):
    for j in mm:
        date = str(i)+str(j)
        url_s = "https://n5eil01u.ecs.nsidc.org/PM/NSIDC-0081.002/2021.12.14/NSIDC0081_SEAICE_PS_N25km_20211214_v2.0.nc"
        wget.download(url_s)

inpath = 'C:/Users/PC/Desktop/jhc/rename'
file_list = os.listdir(inpath)
for file in file_list:
    if file.endswith("v04r00.nc"):
        file.split("_")
        newname = str('_'.join(file.split('_')[:5])+'_'+'_'.join(file.split('_')[6:]))
        f_old = str(inpath)+'/'+str(file)
        f_new = str(inpath)+'/'+str(newname)
        os.rename(f_old, f_new)
'''

# NORTH


month = np.arange(1,13)
date01 = [] ; 
box0 = [] ; box1 = [] ; box2 = [] ; box3 = []

for i in range(1981,2022):
    for j in month:
        if j<10:
            ym1 = str(i)+'0'+str(j)
        else:
            ym1 = str(i)+str(j)
        date01.append(ym1)
date=np.array(date01)

for i in range(len(date)):   
    ncfile_path = 'C:/Users/PC/Desktop/jhc/g02202/north/seaice_conc_monthly_nh_{0}_v04r00.nc'.format(date[i])
    ncfile = Dataset(ncfile_path)
    cdr0 = ncfile.variables['cdr_seaice_conc_monthly']
    cdr0 = np.squeeze(cdr0, axis=0)
    box0.append(cdr0)
    
sic_mon_n = np.array(box0)
sic_mon_n[sic_mon_n>250] = np.NAN

sic_n= sic_mon_n.reshape(41,12,448,304)
sic_mn_n = np.mean(sic_n, axis=0)

sic_jan1 = sic_n[0:7,0,:,:]
sic_jan2 = sic_n[8:,0,:,:]
sic_jan = np.mean(np.concatenate((sic_jan1,sic_jan2),axis=0),axis=0)
sic_feb = np.mean(sic_n[:,1,:,:],axis=0)
sic_mar = np.mean(sic_n[:,2,:,:],axis=0)
sic_apr = np.mean(sic_n[:,3,:,:],axis=0)
sic_may = np.mean(sic_n[:,4,:,:],axis=0)
sic_jun = np.mean(sic_n[:,5,:,:],axis=0)
sic_jul1 = sic_n[0:3,6,:,:]
sic_jul2 = sic_n[4:,6,:,:]
sic_jul = np.mean(np.concatenate((sic_jul1,sic_jul2),axis=0),axis=0)
sic_aug = np.mean(sic_n[:,7,:,:],axis=0)
sic_sep = np.mean(sic_n[:,8,:,:],axis=0)
sic_oct = np.mean(sic_n[:,9,:,:],axis=0)
sic_nov = np.mean(sic_n[:,10,:,:],axis=0)
sic_dec1 = sic_n[0:6,11,:,:]
sic_dec2 = sic_n[7:,11,:,:]
sic_dec = np.mean(np.concatenate((sic_dec1,sic_dec2),axis=0),axis=0)


cmap1 = mpl.colors.ListedColormap(['#FFFFFF','#1C1C38','#333266','#3E4993','#3E68AF','#4988BC','#5FA6C7',
                                  '#80C3D2','#B3DFE3','#E7FBFB'])

'''
fig = plt.figure(figsize=(5, 8))
ax = fig.add_subplot(111, projection=ccrs.NorthPolarStereo(central_longitude=-47))     

ax.set_extent([-55, -42, 59, 70], crs=ccrs.PlateCarree())

resol = '50m'  # use data at this scale
land = cartopy.feature.NaturalEarthFeature('physical', 'land', \
    scale=resol, edgecolor='k', facecolor=cfeature.COLORS['land'])
ocean = cartopy.feature.NaturalEarthFeature('physical', 'ocean', \
    scale=resol, edgecolor='none', facecolor=cfeature.COLORS['water'])

# plot sequence is important
ax.add_feature(ocean, linewidth=0.2 )
ax.add_feature(land, facecolor='green')

ax.set_title('central_longitude on -47 $^\circ$W')

# for cartopy version 0.18 only
ax.gridlines(draw_labels=True, linewidth=0.5, linestyle='--', color='black')
'''


plt.figure(1)
#1
ax = plt.subplot(projection=ccrs.NorthPolarStereo(central_longitude=0))

ax.set_extent([-180,180,90,0], ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')

gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = False, linewidth = 0.6, color = 'gray', alpha = 0.5, linestyle = '--')

plt.imshow(sic_jan, cmap = 'jet', transform = ccrs.PlateCarree())

ax.set_title('NSIDC north SICmean Jan', fontsize = 7, fontweight = 'bold', pad=5)



#2
ax = plt.subplot(2,3,2,projection=ccrs.SouthPolarStereo(central_longitude=0))

ax.set_extent([-180,180,-90,-53], ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')

gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = False, linewidth = 0.6, color = 'gray', alpha = 0.5, linestyle = '--')

plotfield = ax.pcolormesh(lon,lat,sic_feb[:-1,:-1], cmap = cmap1, transform = ccrs.PlateCarree())

ax.set_title('NSIDC south SICmean Feb', fontsize = 7, fontweight = 'bold', pad=5)


#3
ax = plt.subplot(2,3,3,projection=ccrs.SouthPolarStereo(central_longitude=0))

ax.set_extent([-180,180,-90,-53], ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')

gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = False, linewidth = 0.6, color = 'gray', alpha = 0.5, linestyle = '--')

plotfield = ax.pcolormesh(lon,lat,sic_mar[:-1,:-1], cmap = cmap1, transform = ccrs.PlateCarree())

ax.set_title('NSIDC south SICmean Mar', fontsize = 7, fontweight = 'bold', pad=5)



#4
ax = plt.subplot(2,3,4,projection=ccrs.SouthPolarStereo(central_longitude=0))

ax.set_extent([-180,180,-90,-53], ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')

gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = False, linewidth = 0.6, color = 'gray', alpha = 0.5, linestyle = '--')

plotfield = ax.pcolormesh(lon,lat,sic_apr[:-1,:-1], cmap = cmap1, transform = ccrs.PlateCarree())

ax.set_title('NSIDC south SICmean Apr', fontsize = 7, fontweight = 'bold', pad=5)



#5
ax = plt.subplot(2,3,5,projection=ccrs.SouthPolarStereo(central_longitude=0))

ax.set_extent([-180,180,-90,-53], ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')

gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = False, linewidth = 0.6, color = 'gray', alpha = 0.5, linestyle = '--')

plotfield = ax.pcolormesh(lon,lat,sic_mar[:-1,:-1], cmap = cmap1, transform = ccrs.PlateCarree())

ax.set_title('NSIDC south SICmean Mar', fontsize = 7, fontweight = 'bold', pad=5)



#6   
ax = plt.subplot(2,3,6,projection=ccrs.SouthPolarStereo(central_longitude=0))

ax.set_extent([-180,180,-90,-53], ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')

gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = False, linewidth = 0.6, color = 'gray', alpha = 0.5, linestyle = '--')

plotfield = ax.pcolormesh(lon,lat,sic_jun[:-1,:-1], cmap = cmap1, transform = ccrs.PlateCarree())

ax.set_title('NSIDC south SICmean Jun', fontsize = 7, fontweight = 'bold', pad=5)


plt.savefig('NSIDC_SIC_sh_jan2jun.png', dpi=500, bbox_inches='tight')
plt.show()

plt.figure(2)
#7
ax = plt.subplot(2,3,1,projection=ccrs.SouthPolarStereo(central_longitude=0))

ax.set_extent([-180,180,-90,-53], ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')

gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = False, linewidth = 0.6, color = 'gray', alpha = 0.5, linestyle = '--')

plotfield = ax.pcolormesh(lon,lat,sic_jul[:-1,:-1], cmap = cmap1, transform = ccrs.PlateCarree())

ax.set_title('NSIDC south SICmean Jul', fontsize = 7, fontweight = 'bold', pad=5)



#2
ax = plt.subplot(2,3,2,projection=ccrs.SouthPolarStereo(central_longitude=0))

ax.set_extent([-180,180,-90,-53], ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')

gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = False, linewidth = 0.6, color = 'gray', alpha = 0.5, linestyle = '--')

plotfield = ax.pcolormesh(lon,lat,sic_aug[:-1,:-1], cmap = cmap1, transform = ccrs.PlateCarree())

ax.set_title('NSIDC south SICmean Aug', fontsize = 7, fontweight = 'bold', pad=5)


#3
ax = plt.subplot(2,3,3,projection=ccrs.SouthPolarStereo(central_longitude=0))

ax.set_extent([-180,180,-90,-53], ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')

gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = False, linewidth = 0.6, color = 'gray', alpha = 0.5, linestyle = '--')

plotfield = ax.pcolormesh(lon,lat,sic_sep[:-1,:-1], cmap = cmap1, transform = ccrs.PlateCarree())

ax.set_title('NSIDC south SICmean Sep', fontsize = 7, fontweight = 'bold', pad=5)



#4
ax = plt.subplot(2,3,4,projection=ccrs.SouthPolarStereo(central_longitude=0))

ax.set_extent([-180,180,-90,-53], ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')

gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = False, linewidth = 0.6, color = 'gray', alpha = 0.5, linestyle = '--')

plotfield = ax.pcolormesh(lon,lat,sic_oct[:-1,:-1], cmap = cmap1, transform = ccrs.PlateCarree())

ax.set_title('NSIDC south SICmean Oct', fontsize = 7, fontweight = 'bold', pad=5)



#5
ax = plt.subplot(2,3,5,projection=ccrs.SouthPolarStereo(central_longitude=0))

ax.set_extent([-180,180,-90,-53], ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')

gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = False, linewidth = 0.6, color = 'gray', alpha = 0.5, linestyle = '--')

plotfield = ax.pcolormesh(lon,lat,sic_nov[:-1,:-1], cmap = cmap1, transform = ccrs.PlateCarree())

ax.set_title('NSIDC south SICmean Nov', fontsize = 7, fontweight = 'bold', pad=5)



#6   
ax = plt.subplot(2,3,6,projection=ccrs.SouthPolarStereo(central_longitude=0))

ax.set_extent([-180,180,-90,-53], ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')

gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = False, linewidth = 0.6, color = 'gray', alpha = 0.5, linestyle = '--')

plotfield = ax.pcolormesh(lon,lat,sic_dec[:-1,:-1], cmap = cmap1, transform = ccrs.PlateCarree())

ax.set_title('NSIDC south SICmean Dec', fontsize = 7, fontweight = 'bold', pad=5)


plt.savefig('NSIDC_SIC_sh_jul2dec.png', dpi=500, bbox_inches='tight')
plt.show()
