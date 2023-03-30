# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 09:24:19 2022

@author: PC
"""
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
import matplotlib.colors as colors
import matplotlib as mpl
import os
import wget
from netCDF4 import Dataset
from matplotlib import animation, rc
import sys
sys.path.append('C:/Users/PC/Desktop/jhc/03120320')
np.set_printoptions(threshold=sys.maxsize)

'''        
mar = ['0312', '0313', '0319', '0320']
for i in range(2008,2011):
    for j in mar:
        date = str(i)+str(j)
        url_s = "ftp://sidads.colorado.edu/DATASETS/NOAA/G02202_V4/south/daily/{0}/seaice_conc_daily_sh_{1}_f17_v04r00.nc".format(i, date)
        wget.download(url_s)

inpath = 'C:/Users/PC/Desktop/jhc/03120320'
file_list = os.listdir(inpath)
for file in file_list:
    if file.endswith("v04r00.nc"):
        file.split("_")
        newname = str('_'.join(file.split('_')[:5])+'_'+'_'.join(file.split('_')[6:]))
        f_old = str(inpath)+'/'+str(file)
        f_new = str(inpath)+'/'+str(newname)
        os.rename(f_old, f_new)
'''
# importing data

nc_grid0 = 'C:/Users/PC/Desktop/jhc/g02202/south/seaice_conc_monthly_icdr_sh_f18_202103_v01r00.nc'
nc_grid = Dataset(nc_grid0)
lat = nc_grid.variables['latitude']
lon = nc_grid.variables['longitude']

day_1213 = ['0312','0313'] ; day_1920 = ['0319','0320']
date_221213 = ['20220312', '20220313'] ; date_221920 = ['20220319','20220320']
date01 = [] ; date02 = []
box0 = [] ; box1 = [] ; box2 = [] ; box3 = []

## 1981-2010

for i in range(1981,2011):
    for j in day_1213:
        ym1 =int(str(i)+str(j))
        date01.append(ym1)       
date_1213=np.array(date01)

for i in range(1981,2011):
    for j in day_1920:
        ym2 =int(str(i)+str(j))
        date02.append(ym2)       
date_1920=np.array(date02)

for i in range(len(date_1213)):   
    ncfile_path = 'C:/Users/PC/Desktop/jhc/03120320/seaice_conc_daily_sh_{0}_v04r00.nc'.format(date_1213[i])
    ncfile = Dataset(ncfile_path)
    cdr0 = ncfile.variables['cdr_seaice_conc']
    cdr0 = np.squeeze(cdr0, axis=0)
    box0.append(cdr0)
sic_1213 = np.array(box0)
sic_1213[sic_1213>250] = np.NAN

for i in range(len(date_1920)):   
    ncfile_path = 'C:/Users/PC/Desktop/jhc/03120320/seaice_conc_daily_sh_{0}_v04r00.nc'.format(date_1920[i])
    ncfile = Dataset(ncfile_path)
    cdr0 = ncfile.variables['cdr_seaice_conc']
    cdr0 = np.squeeze(cdr0, axis=0)
    box1.append(cdr0)
sic_1920 = np.array(box1)
sic_1920[sic_1920>250] = np.NAN  

## 2022   

for i in range(2):   
    ncfile_path = 'C:/Users/PC/Desktop/jhc/03120320/seaice_conc_daily_icdr_sh_{0}_f18_v02r00.nc'.format(date_221213[i])
    ncfile = Dataset(ncfile_path)
    cdr0 = ncfile.variables['cdr_seaice_conc']
    cdr0 = np.squeeze(cdr0, axis=0)
    box2.append(cdr0)
sic_221213 = np.array(box2)
sic_221213[sic_221213>250] = np.NAN

for i in range(2):   
    ncfile_path = 'C:/Users/PC/Desktop/jhc/03120320/seaice_conc_daily_icdr_sh_{0}_f18_v02r00.nc'.format(date_221920[i])
    ncfile = Dataset(ncfile_path)
    cdr0 = ncfile.variables['cdr_seaice_conc']
    cdr0 = np.squeeze(cdr0, axis=0)
    box3.append(cdr0)
sic_221920 = np.array(box3)
sic_221920[sic_221920>250] = np.NAN

# 1) 2022년 1213 평균 sic
mn_221213 = 100*(np.mean(sic_221213, axis=0))

# 2) 2022년 1920 평균 sic

mn_221920 = 100*(np.mean(sic_221920, axis=0))

# 3) 1981-2010년 1213 평균 sic

mn_1213 = 100*(np.mean(sic_1213, axis=0))

# 4) 1981-2010년 1920 평균 sic

mn_1920 = 100*(np.mean(sic_1920, axis=0))

# 5) 1213 sic ano

ano_1213 = mn_221213 - mn_1213

# 6) 1920 sic ano

ano_1920 = mn_221920 - mn_1920

# 7) 1920 평균 sic - 1213 평균 sic

mn_mi = mn_221920 - mn_221213

# 8) 1920 sic ano - 1213 sic ano

ano_mi = ano_1920 - ano_1213

# plotting 1

fig = plt.figure(figsize=(6,5))
ax = plt.subplot(projection=ccrs.SouthPolarStereo(central_longitude=-230))

cmap1 = mpl.colors.ListedColormap(['#FFFFFF','#1C1C38','#333266','#3E4993','#3E68AF','#4988BC','#5FA6C7',
                                  '#80C3D2','#B3DFE3','#E7FBFB'])

bounds = [0,15,25,35,45,55,65,75,85,90,100]
norm = mpl.colors.BoundaryNorm(bounds, cmap1.N)

ax.set_extent([-280, -190, -90, -57], ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')

gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = False, linewidth = 0.6, color = 'gray', alpha = 0.5, linestyle = '--')
gl.xlocator = mticker.FixedLocator([60,90,120,150,180])

gl2 = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.0, linestyle='--')
gl2.rotate_labels = False
gl2.xlocator=mticker.FixedLocator([45,60,75,90,105,120,135,150])
gl2.ylocator=mticker.FixedLocator([-80, -70, -60])
gl2.bottom_labels = gl2.right_labels=False

gl2.xlabel_style = {'size':10, 'color': 'black', 'rotation': 0, 'rotation_mode': 'anchor'}
gl2.ylabel_style = {'size':8, 'color': 'k', 'rotation': 0, 'rotation_mode': 'anchor'}

plotfield = ax.pcolormesh(lon,lat,mn_221213[:-1,:-1], cmap = cmap1, transform = ccrs.PlateCarree())

cbar = fig.colorbar(plotfield, shrink=0.63, ticks=bounds, boundaries=bounds , orientation = 'vertical', spacing='uniform', extendfrac='auto', norm=norm)
cbar.ax.tick_params(direction = 'in', labelsize = 9)
cbar.set_ticks([15,25,35,45,55,65,75,85,90])

ax.set_title('NSIDC SIC 2022-Mar-(12-13)', fontsize = 14, fontweight = 'bold', pad=15)

plt.subplots_adjust(top=1, bottom=0.1, left=0.1, right=0.95, wspace=0.2, hspace=0.)

plt.savefig('NSIDC_SIC_2022_Mar_(12-13).png', dpi=500, bbox_inches='tight')
plt.show()

# plotting 2

fig = plt.figure(figsize=(6,5))
ax = plt.subplot(projection=ccrs.SouthPolarStereo(central_longitude=-230))

ax.set_extent([-280, -190, -90, -57], ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')


gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = False, linewidth = 0.6, color = 'gray', alpha = 0.5, linestyle = '--')
gl.xlocator = mticker.FixedLocator([60,90,120,150,180])

gl2 = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.0, linestyle='--')
gl2.rotate_labels = False
gl2.xlocator=mticker.FixedLocator([45,60,75,90,105,120,135,150])
gl2.ylocator=mticker.FixedLocator([-80, -70, -60])
gl2.bottom_labels = gl2.right_labels=False

gl2.xlabel_style = {'size':10, 'color': 'black', 'rotation': 0, 'rotation_mode': 'anchor'}
gl2.ylabel_style = {'size':8, 'color': 'k', 'rotation': 0, 'rotation_mode': 'anchor'}

plotfield = ax.pcolormesh(lon,lat,mn_221920[:-1,:-1], cmap = cmap1, transform = ccrs.PlateCarree())

cbar = fig.colorbar(plotfield, shrink=0.63, ticks=bounds, boundaries=bounds, orientation = 'vertical', spacing='uniform', extendfrac='auto', norm=norm)
cbar.ax.tick_params(direction = 'in', labelsize = 9)
cbar.set_ticks([15,25,35,45,55,65,75,85,90])

ax.set_title('NSIDC SIC 2022-Mar-(19-20)', fontsize = 14, fontweight = 'bold', pad=15)

plt.subplots_adjust(top=1, bottom=0.1, left=0.1, right=0.95, wspace=0.2, hspace=0.)

plt.savefig('NSIDC_SIC_2022_Mar_(19-20).png', dpi=500, bbox_inches='tight')
plt.show()

# plotting 3

fig = plt.figure(figsize=(6,5))
ax = plt.subplot(projection=ccrs.SouthPolarStereo(central_longitude=-230))

ax.set_extent([-280, -190, -90, -57], ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')


gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = False, linewidth = 0.6, color = 'gray', alpha = 0.5, linestyle = '--')
gl.xlocator = mticker.FixedLocator([60,90,120,150,180])

gl2 = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.0, linestyle='--')
gl2.rotate_labels = False
gl2.xlocator=mticker.FixedLocator([45,60,75,90,105,120,135,150])
gl2.ylocator=mticker.FixedLocator([-80, -70, -60])
gl2.bottom_labels = gl2.right_labels=False

gl2.xlabel_style = {'size':10, 'color': 'black', 'rotation': 0, 'rotation_mode': 'anchor'}
gl2.ylabel_style = {'size':8, 'color': 'k', 'rotation': 0, 'rotation_mode': 'anchor'}

plotfield = ax.pcolormesh(lon,lat,mn_1213[:-1,:-1], cmap = cmap1, transform = ccrs.PlateCarree())

cbar = fig.colorbar(plotfield, shrink=0.63, ticks=bounds, boundaries=bounds, orientation = 'vertical', spacing='uniform', extendfrac='auto', norm=norm)
cbar.ax.tick_params(direction = 'in', labelsize = 9)
cbar.set_ticks([15,25,35,45,55,65,75,85,90])

ax.set_title('NSIDC SIC 8110-Mar-(12-13)', fontsize = 14, fontweight = 'bold', pad=15)

plt.subplots_adjust(top=1, bottom=0.1, left=0.1, right=0.95, wspace=0.2, hspace=0.)

plt.savefig('NSIDC_SIC_8110-Mar-(12-13).png', dpi=500, bbox_inches='tight')
plt.show()

# plotting 4

fig = plt.figure(figsize=(6,5))
ax = plt.subplot(projection=ccrs.SouthPolarStereo(central_longitude=-230))

ax.set_extent([-280, -190, -90, -57], ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')


gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = False, linewidth = 0.6, color = 'gray', alpha = 0.5, linestyle = '--')
gl.xlocator = mticker.FixedLocator([60,90,120,150,180])

gl2 = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.0, linestyle='--')
gl2.rotate_labels = False
gl2.xlocator=mticker.FixedLocator([45,60,75,90,105,120,135,150])
gl2.ylocator=mticker.FixedLocator([-80, -70, -60])
gl2.bottom_labels = gl2.right_labels=False

gl2.xlabel_style = {'size':10, 'color': 'black', 'rotation': 0, 'rotation_mode': 'anchor'}
gl2.ylabel_style = {'size':8, 'color': 'k', 'rotation': 0, 'rotation_mode': 'anchor'}

plotfield = ax.pcolormesh(lon,lat,mn_1920[:-1,:-1], cmap = cmap1, transform = ccrs.PlateCarree())

cbar = fig.colorbar(plotfield, shrink=0.63, ticks=bounds, boundaries=bounds, orientation = 'vertical', spacing='uniform', extendfrac='auto', norm=norm)
cbar.ax.tick_params(direction = 'in', labelsize = 9)
cbar.set_ticks([15,25,35,45,55,65,75,85,90])

ax.set_title('NSIDC SIC 8110-Mar-(19-20)', fontsize = 14, fontweight = 'bold', pad=15)

plt.subplots_adjust(top=1, bottom=0.1, left=0.1, right=0.95, wspace=0.2, hspace=0.)

plt.savefig('NSIDC_SIC_8110-Mar-(19-20).png', dpi=500, bbox_inches='tight')
plt.show()

# plotting 5

fig = plt.figure(figsize=(6,5))
ax = plt.subplot(projection=ccrs.SouthPolarStereo(central_longitude=-230))

cmap2 = mpl.colors.ListedColormap(['#071E46','#072F6B','#08579C','#2171B5','#2171B5','#4292C7','#5AA0CD',
                                  '#78BFD6','#AADCE6','#DBF5FF','#FFFFFF','#FFFFFF','#FFE0E0','#FCBBAA','#FC9272','#FB6A4A','#F03C2B','#CC181E','#CC181E','#A60F14','#780A0F','#5F0000'])

bounds = [-55,-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45,50,55]
norm = mpl.colors.BoundaryNorm(bounds, cmap2.N)

ax.set_extent([-280, -190, -90, -57], ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')

gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = False, linewidth = 0.6, color = 'gray', alpha = 0.5, linestyle = '--')
gl.xlocator = mticker.FixedLocator([60,90,120,150,180])

gl2 = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.0, linestyle='--')
gl2.rotate_labels = False
gl2.xlocator=mticker.FixedLocator([45,60,75,90,105,120,135,150])
gl2.ylocator=mticker.FixedLocator([-80, -70, -60])
gl2.bottom_labels = gl2.right_labels=False

gl2.xlabel_style = {'size':10, 'color': 'black', 'rotation': 0, 'rotation_mode': 'anchor'}
gl2.ylabel_style = {'size':8, 'color': 'k', 'rotation': 0, 'rotation_mode': 'anchor'}

plotfield0 = ax.pcolormesh(lon,lat,ano_1213[:-1,:-1], cmap = cmap2, norm=norm, transform = ccrs.PlateCarree())

cbar = fig.colorbar(plotfield0, shrink=0.63, ticks=bounds, boundaries=bounds, orientation = 'vertical',spacing='uniform', extendfrac='auto', norm=norm)
cbar.ax.tick_params(direction = 'in', labelsize = 9)
cbar.set_ticks([-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45,50])

ax.set_title('NSIDC SIC anomaly (2022 Mar 12-13)', fontsize = 14, fontweight = 'bold', pad=15)

plt.subplots_adjust(top=1, bottom=0.1, left=0.1, right=0.95, wspace=0.2, hspace=0.)

plt.savefig('NSIDC_SIC_anomaly_(2022_Mar_12-13).png', dpi=500, bbox_inches='tight')
plt.show()

# plotting 6

fig = plt.figure(figsize=(6,5))
ax = plt.subplot(projection=ccrs.SouthPolarStereo(central_longitude=-230))

ax.set_extent([-280, -190, -90, -57], ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')


gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = False, linewidth = 0.6, color = 'gray', alpha = 0.5, linestyle = '--')
gl.xlocator = mticker.FixedLocator([60,90,120,150,180])

gl2 = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.0, linestyle='--')
gl2.rotate_labels = False
gl2.xlocator=mticker.FixedLocator([45,60,75,90,105,120,135,150])
gl2.ylocator=mticker.FixedLocator([-80, -70, -60])
gl2.bottom_labels = gl2.right_labels=False

gl2.xlabel_style = {'size':10, 'color': 'black', 'rotation': 0, 'rotation_mode': 'anchor'}
gl2.ylabel_style = {'size':8, 'color': 'k', 'rotation': 0, 'rotation_mode': 'anchor'}

plotfield1 = ax.pcolormesh(lon,lat,ano_1920[:-1,:-1],cmap = cmap2, norm=norm, transform = ccrs.PlateCarree())

cbar = fig.colorbar(plotfield1, shrink=0.63, ticks=bounds, boundaries=bounds, orientation = 'vertical',spacing='uniform', extendfrac='auto', norm=norm)
cbar.ax.tick_params(direction = 'in', labelsize = 9)
cbar.set_ticks([-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45,50])

ax.set_title('NSIDC SIC anomaly (2022 Mar 19-20)', fontsize = 14, fontweight = 'bold', pad=15)

plt.subplots_adjust(top=1, bottom=0.1, left=0.1, right=0.95, wspace=0.2, hspace=0.)

plt.savefig('NSIDC_SIC_anomaly_(2022_Mar_19-20).png', dpi=500, bbox_inches='tight')
plt.show()

# plotting 7

fig = plt.figure(figsize=(6,5))
ax = plt.subplot(projection=ccrs.SouthPolarStereo(central_longitude=-230))

ax.set_extent([-280, -190, -90, -57], ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')


gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = False, linewidth = 0.6, color = 'gray', alpha = 0.5, linestyle = '--')
gl.xlocator = mticker.FixedLocator([60,90,120,150,180])

gl2 = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.0, linestyle='--')
gl2.rotate_labels = False
gl2.xlocator=mticker.FixedLocator([45,60,75,90,105,120,135,150])
gl2.ylocator=mticker.FixedLocator([-80, -70, -60])
gl2.bottom_labels = gl2.right_labels=False

gl2.xlabel_style = {'size':10, 'color': 'black', 'rotation': 0, 'rotation_mode': 'anchor'}
gl2.ylabel_style = {'size':8, 'color': 'k', 'rotation': 0, 'rotation_mode': 'anchor'}

plotfield2 = ax.pcolormesh(lon,lat,mn_mi[:-1,:-1], cmap = cmap2, norm=norm, transform = ccrs.PlateCarree())

cbar = fig.colorbar(plotfield2, shrink=0.63, ticks=bounds, boundaries=bounds, orientation = 'vertical',spacing='uniform', extendfrac='auto', norm=norm)
cbar.ax.tick_params(direction = 'in', labelsize = 9)
cbar.set_ticks([-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45,50])

ax.set_title('a) SIC difference', fontsize = 14, fontweight = 'bold', pad=15)

plt.subplots_adjust(top=1, bottom=0.1, left=0.1, right=0.95, wspace=0.2, hspace=0.)

plt.savefig('NSIDC_SIC_weekly_difference_(2022_Mar_19-20)_minus_(2022_Mar_12-13).png', dpi=500, bbox_inches='tight')
plt.show()

# plotting 8

fig = plt.figure(figsize=(6,5))
ax = plt.subplot(projection=ccrs.SouthPolarStereo(central_longitude=-230))

ax.set_extent([-280, -190, -90, -57], ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')


gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = False, linewidth = 0.6, color = 'gray', alpha = 0.5, linestyle = '--')
gl.xlocator = mticker.FixedLocator([60,90,120,150,180])

gl2 = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.0, linestyle='--')
gl2.rotate_labels = False
gl2.xlocator=mticker.FixedLocator([45,60,75,90,105,120,135,150])
gl2.ylocator=mticker.FixedLocator([-80, -70, -60])
gl2.bottom_labels = gl2.right_labels=False

gl2.xlabel_style = {'size':10, 'color': 'black', 'rotation': 0, 'rotation_mode': 'anchor'}
gl2.ylabel_style = {'size':8, 'color': 'k', 'rotation': 0, 'rotation_mode': 'anchor'}

plotfield3 = ax.pcolormesh(lon,lat,ano_mi[:-1,:-1], cmap = cmap2, norm=norm, transform = ccrs.PlateCarree())

cbar = fig.colorbar(plotfield3, shrink=0.63, ticks=bounds, boundaries=bounds, orientation = 'vertical',spacing='uniform', extendfrac='auto', norm=norm)
cbar.ax.tick_params(direction = 'in', labelsize = 9)
cbar.set_ticks([-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45,50])

ax.set_title('b) SIC anomaly difference', fontsize = 14, fontweight = 'bold', pad=15)

plt.subplots_adjust(top=1, bottom=0.1, left=0.1, right=0.95, wspace=0.2, hspace=0.)

plt.savefig('NSIDC_SIC_weekly_anomaly_difference_(2022_Mar_19-20)_minus_(2022_Mar_12-13).png', dpi=500, bbox_inches='tight')
plt.show()



















