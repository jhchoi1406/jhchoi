# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 10:00:27 2022

@author: PC
"""

# pip --trusted-host pypi.org --trusted-host files.pythonhosted.org install netCDF4

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
import sys
sys.path.append('C:/Users/PC/Desktop/jhc/03120320')
np.set_printoptions(threshold=sys.maxsize)

# data for finding grid location

nc_grid0 = 'C:/Users/PC/Desktop/jhc/g02202/seaice_conc_monthly_icdr_sh_f18_202103_v01r00.nc'
nc_grid = Dataset(nc_grid0)
lat = nc_grid.variables['latitude']
lon = nc_grid.variables['longitude']

mm = ['01','02','03','04','05','06']

# downloading NSIDC data
'''
for i in range(2008,2022):
    for j in mm:
        date = str(i)+str(j)
        url_s = "ftp://sidads.colorado.edu/DATASETS/NOAA/G02202_V4/south/monthly/seaice_conc_monthly_sh_{0}_f17_v04r00.nc".format(date)
        wget.download(url_s)
'''
# renaming data
'''
inpath = 'C:/Users/PC/Desktop/jhc/g02202'
file_list = os.listdir(inpath)
for file in file_list:
    if file.endswith("v04r00.nc"):
        file.split("_")
        newname = str('_'.join(file.split('_')[:5])+'_'+'_'.join(file.split('_')[6:]))
        f_old = str(inpath)+'/'+str(file)
        f_new = str(inpath)+'/'+str(newname)
        os.rename(f_old, f_new)
'''

# G02202

mm = ['01','02','03','04','05','06']

years = list(range(1981,2011))
date0 = []
box0 = []

for i in range(1981,2011):
    for j in mm:
        ym =int(str(i)+str(j))
        date0.append(ym)       
date=np.array(date0)

for i in range(len(date)):   
    ncfile_path = 'C:/Users/PC/Desktop/jhc/g02202/seaice_conc_monthly_sh_{0}_v04r00.nc'.format(date[i])
    ncfile = Dataset(ncfile_path)
    cdr0 = ncfile.variables['cdr_seaice_conc_monthly']
    cdr = np.squeeze(cdr0, axis=0)
    box0.append(cdr)
    
mn = np.array(box0)
mn[mn>250] = np.NAN

mn_8110= mn.reshape(30,6,332,316)
mn_mn = np.mean(mn_8110, axis=0)

# G10016

ano_8110 = np.zeros((6,332,316))
box0_n = []
cdr_t = []

for i in range(6):   
    ncfile_path_n = 'C:/Users/PC/Desktop/jhc/g02202/seaice_conc_monthly_icdr_sh_20220{0}_f18_v02r00.nc'.format(str(i+1))
    ncfile_n = Dataset(ncfile_path_n)
    cdr0_n = ncfile_n.variables['cdr_seaice_conc_monthly']
    cdr_n = np.squeeze(cdr0_n, axis=0)
    cdr_t.append(cdr_n)
    
cdr_tt = np.array(cdr_t)
cdr_tt[cdr_tt>250] = np.NAN

ano_8110 = 100*(cdr_tt - mn_mn)


# plotting

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

plotfield = ax.pcolormesh(lon,lat,ano_8110[2,:-1,:-1], cmap = cmap2, norm=colors.CenteredNorm(), transform = ccrs.PlateCarree())

cbar = fig.colorbar(plotfield, shrink=0.63, ticks=bounds, boundaries=bounds, orientation = 'vertical',spacing='uniform', extendfrac='auto', norm=norm)
cbar.ax.tick_params(direction = 'in', labelsize = 9)
cbar.set_ticks([-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45,50])

ax.set_title('NSIDC SICano 2022-Mar Monthly', fontsize = 14, fontweight = 'bold', pad=15)

plt.subplots_adjust(top=1, bottom=0.1, left=0.1, right=0.95, wspace=0.2, hspace=0.)

plt.savefig('ano_mar.png', dpi=500, bbox_inches='tight')
plt.show()




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

plotfield = ax.pcolormesh(lon,lat,mn_mn[2,:-1,:-1], cmap = cmap1, transform = ccrs.PlateCarree())



ax.set_title('NSIDC SICmean (1981-2010)-Mar Monthly', fontsize = 14, fontweight = 'bold', pad=15)

plt.subplots_adjust(top=1, bottom=0.1, left=0.1, right=0.95, wspace=0.2, hspace=0.)

plt.savefig('mean_mar.png', dpi=500, bbox_inches='tight')
plt.show()




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

plotfield = ax.pcolormesh(lon,lat,cdr_tt[2,:-1,:-1], cmap = cmap1, transform = ccrs.PlateCarree())

ax.set_title('NSIDC SICmean 2022-Mar Monthly', fontsize = 14, fontweight = 'bold', pad=15)

plt.subplots_adjust(top=1, bottom=0.1, left=0.1, right=0.95, wspace=0.2, hspace=0.)

plt.savefig('mean_mar2022.png', dpi=500, bbox_inches='tight')
plt.show()




'''
fig = plt.figure(figsize=(6,5))
ax = plt.subplot(projection=ccrs.SouthPolarStereo(central_longitude=-230))

ax.set_extent([-280, -190, -90, -57], ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='#F7F6F6')


gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = False, linewidth = 0.6, color = 'gray', alpha = 0.5, linestyle = '--')
gl.xlocator = mticker.FixedLocator([60,90,120,150,180])

gl2 = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.0, linestyle='--')
gl2.rotate_labels = False
gl2.xlocator=mticker.FixedLocator([45,60,75,90,105,120,135,150])
gl2.ylocator=mticker.FixedLocator([-80, -70, -60])
gl2.bottom_labels = gl2.right_labels=False

gl2.xlabel_style = {'size':10, 'color': 'black', 'rotation': 0, 'rotation_mode': 'anchor'}
gl2.ylabel_style = {'size':8, 'color': 'k', 'rotation': 0, 'rotation_mode': 'anchor'}

plotfield = ax.pcolormesh(lon,lat,ano_8110[2,:-1,:-1], cmap = 'RdBu_r', norm=colors.CenteredNorm(), transform = ccrs.PlateCarree())

cbar = fig.colorbar(plotfield, shrink=0.63, orientation = 'vertical', extend='both')
cbar.ax.tick_params(direction = 'in', labelsize = 9)

ax.set_title('NSIDC SICano 2022-Mar Monthly', fontsize = 14, fontweight = 'bold', pad=15)

plt.subplots_adjust(top=1, bottom=0.1, left=0.1, right=0.95, wspace=0.2, hspace=0.)

plt.savefig('ano_mar1.png', dpi=100, bbox_inches='tight')
plt.show()
'''
