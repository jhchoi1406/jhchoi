# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 16:30:10 2022

@author: PC
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn import linear_model
import matplotlib.colors 
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
import matplotlib as mpl
import os
import wget
from netCDF4 import Dataset
import sys
sys.path.append('C:/Users/PC/Desktop/jhc/03120320')
np.set_printoptions(threshold=sys.maxsize)

time0 = np.arange(0,34) ; time = time0.reshape(-1,1)
time1 = np.arange(0,32) ; time2 = time1.reshape(-1,1)
month = np.arange(1,13)
date01 = [] ; 
box0 = [] ; box1 = [] ; box2 = [] ; box3 = []

for i in range(1979,2013):
    for j in month:
        if j<10:
            ym1 = str(i)+'0'+str(j)
        else:
            ym1 = str(i)+str(j)
        date01.append(ym1)
date=np.array(date01)

extra0 = Dataset('C:/Users/PC/Desktop/jhc/g02202/north/seaice_conc_monthly_nh_197812_v04r00.nc')
extra = extra0.variables['cdr_seaice_conc_monthly']
extra_cdr = np.array(np.squeeze(extra, axis=0))
extra_cdr[extra_cdr>250] = np.NAN

for i in range(len(date)):   
    ncfile_path = 'C:/Users/PC/Desktop/jhc/g02202/north/seaice_conc_monthly_nh_{0}_v04r00.nc'.format(date[i])
    ncfile = Dataset(ncfile_path)
    cdr0 = ncfile.variables['cdr_seaice_conc_monthly']
    cdr0 = np.squeeze(cdr0, axis=0)
    box0.append(cdr0)
    
sic_mon = np.array(box0)
sic_mon[sic_mon>250] = np.NaN

sic = sic_mon.reshape(34,12,448,304)

# monthly mean

djf_data = np.zeros([34,3,448,304])
djf_data[0,0,:,:] = extra_cdr
djf_data[0,1:3,:,:] = sic[0,0:2,:,:]
sic[9,0,:,:] = sic[9,1,:,:]
sic[8,11,:,:] = sic[9,1,:,:]
for i in range(1,34):
    djf_data[i,0,:,:] = sic[i-1,11,:,:]
    djf_data[i,1:3,:,:] = sic[i,0:2,:,:]
    djf  = 100*np.mean(djf_data, axis=1)

mam_data = np.zeros([34,3,448,304])
for i in range(34):
    mam_data[i,0:3,:,:] = sic[i,2:5,:,:]
    mam = 100*np.mean(mam_data, axis=1)

jja_data = np.zeros([34,3,448,304])
jja_er1 = sic[5,5,:,:]
jja_er2 = sic[5,7,:,:]
jja_er = np.zeros([2,448,304])
jja_er[0,:,:] = jja_er1
jja_er[1,:,:] = jja_er2
jja_er_mn = 100*np.mean(jja_er, axis=0)
for i in range(34):
    jja_data[i,0:3,:,:] = sic[i,5:8,:,:]
    jja = 100*np.mean(jja_data, axis=1)
jja[5,:,:] = jja_er_mn

son_data = np.zeros([34,3,448,304])
for i in range(34):
    son_data[i,0:3,:,:] = sic[i,8:11,:,:]
    son = 100*np.mean(son_data, axis=1)    

# calculating coef by using linear regression

djf_lr = np.zeros([448,304])
for i in range(448):
    for j in range(304):
        if sum(np.isnan(djf[:,i,j])) == 0:
            ols = linear_model.LinearRegression()
            model = ols.fit(time, djf[:,i,j])
            lr = model.coef_
        else:
            lr = np.NaN
        djf_lr[i,j] = lr

mam_lr = np.zeros([448,304])
for i in range(448):
    for j in range(304):
        if sum(np.isnan(mam[:,i,j])) == 0:
            ols = linear_model.LinearRegression()
            model = ols.fit(time, mam[:,i,j])
            lr = model.coef_
        else:
            lr = np.NaN
        mam_lr[i,j] = lr

jja_lr = np.zeros([448,304])
for i in range(448):
    for j in range(304):
        if sum(np.isnan(jja[:,i,j])) == 0:
            ols = linear_model.LinearRegression()
            model = ols.fit(time, jja[:,i,j])
            lr = model.coef_
        else:
            lr = np.NaN
        jja_lr[i,j] = lr

son_lr = np.zeros([448,304])
for i in range(448):
    for j in range(304):
        if sum(np.isnan(son[:,i,j])) == 0:
            ols = linear_model.LinearRegression()
            model = ols.fit(time, son[:,i,j])
            lr = model.coef_
        else:
            lr = np.NaN
        son_lr[i,j] = lr

# plotting

cmap2 = mpl.colors.ListedColormap(['#071E46','#072F6B','#08579C','#2171B5','#4292C7','#5AA0CD',
                                  '#78BFD6','#AADCE6','#DBF5FF','#FFFFFF','#FFFFFF','#FFE0E0','#FCBBAA','#FC9272','#FB6A4A','#F03C2B','#CC181E','#A60F14','#780A0F','#5F0000'])

vmin = -2.5
vmax = 2.5
norm = mpl.colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

fig = plt.figure(figsize=(7,10), dpi=500)

# JJA

ax = plt.subplot(2,2,1,projection=ccrs.NorthPolarStereo())

ax.add_feature(cfeature.OCEAN, facecolor='lightgray')

plotfield = plt.imshow(jja_lr, cmap = cmap2, norm=norm)

ax.set_title('Summer (JJA)', fontsize = 17, pad=5)

cbar = fig.colorbar(plotfield, shrink=0.67, orientation = 'vertical',spacing='uniform', extendfrac='auto', norm=norm)
cbar.ax.tick_params(direction = 'in', labelsize = 12)
cbar.set_label('Trend in % per yr', rotation=270, size=13, labelpad=15)
cbar.set_ticks([-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5])

# SON

ax = plt.subplot(2,2,2,projection=ccrs.NorthPolarStereo())

ax.add_feature(cfeature.OCEAN, facecolor='lightgray')

plotfield = plt.imshow(son_lr, cmap = cmap2, norm=norm)

ax.set_title('Autumn (SON)', fontsize = 17, pad=5)

cbar = fig.colorbar(plotfield, shrink=0.67, orientation = 'vertical',spacing='uniform', extendfrac='auto', norm=norm)
cbar.ax.tick_params(direction = 'in', labelsize = 12)
cbar.set_label('Trend in % per yr', rotation=270, size=13, labelpad=15)
cbar.set_ticks([-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5])

# DJF

ax = plt.subplot(2,2,3,projection=ccrs.NorthPolarStereo())

ax.add_feature(cfeature.OCEAN, facecolor='lightgray')

plotfield = plt.imshow(djf_lr, cmap = cmap2, norm=norm)

ax.set_title('Winter (DJF)', fontsize = 17, pad=5)

cbar = fig.colorbar(plotfield, shrink=0.67, orientation = 'vertical',spacing='uniform', extendfrac='auto', norm=norm)
cbar.ax.tick_params(direction = 'in', labelsize = 12)
cbar.set_label('Trend in % per yr', rotation=270, size=13, labelpad=15)
cbar.set_ticks([-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5])

# MAM

ax = plt.subplot(2,2,4,projection=ccrs.NorthPolarStereo())

ax.add_feature(cfeature.OCEAN, facecolor='lightgray')

plotfield = plt.imshow(mam_lr, cmap = cmap2, norm=norm)

ax.set_title('Spring (MAM)', fontsize = 17, pad=5)

cbar = fig.colorbar(plotfield, shrink=0.67, orientation = 'vertical',spacing='uniform', extendfrac='auto', norm=norm)
cbar.ax.tick_params(direction = 'in', labelsize = 12)
cbar.set_label('Trend in % per yr', rotation=270, size=13, labelpad=15)
cbar.set_ticks([-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5])

fig.suptitle('Sea ice concentration trends during 1979-2012',fontsize = 20, fontweight = 'bold', y=1)

fig.tight_layout()
plt.savefig('NSIDC_SIC_nh_djf2son.png', dpi=500, bbox_inches='tight')
plt.show()