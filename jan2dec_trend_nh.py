# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 16:42:29 2022

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

time0 = np.arange(0,43) ; time = [time0] 
month = np.arange(1,13)
date01 = [] ; 
box0 = [] ; box1 = [] ; box2 = [] ; box3 = []

for i in range(1979,2022):
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
sic_mon_n[sic_mon_n>250] = np.NaN

sic_n0= sic_mon_n.reshape(43,12,448,304)

jan1 = sic_n0[0:9,0,:,:]
jan2 = sic_n0[10:,0,:,:]
jan = 100*np.concatenate((jan1,jan2),axis=0)
feb = 100*sic_n0[:,1,:,:] 
mar = 100*sic_n0[:,2,:,:]
apr = 100*sic_n0[:,3,:,:]
may = 100*sic_n0[:,4,:,:]
jun = 100*sic_n0[:,5,:,:]
jul1 = sic_n0[0:5,6,:,:]
jul2 = sic_n0[6:,6,:,:]
jul = 100*np.concatenate((jul1,jul2),axis=0)
aug = 100*sic_n0[:,7,:,:]
sep = 100*sic_n0[:,8,:,:]
oct = 100*sic_n0[:,9,:,:]
nov = 100*sic_n0[:,10,:,:]
dec1 = sic_n0[0:8,11,:,:]
dec2 = sic_n0[9:,11,:,:]
dec = 100*np.concatenate((dec1,dec2),axis=0)

#1

lr_jan = np.zeros([448,304])
time0_ab = np.arange(0,42)
time_ab = time0_ab.reshape(-1,1)
for i in range(448):
    for j in range(304):
        if sum(np.isnan(jan[:,i,j])) == 0:
            ols = linear_model.LinearRegression()
            model = ols.fit(time_ab, jan[:,i,j])
            lr = model.coef_
        else:
            lr = np.NaN
        lr_jan[i,j] = lr

#2

lr_feb = np.zeros([448,304])
time0 = np.arange(0,43)
time = time0.reshape(-1,1)
for i in range(448):
    for j in range(304):
        if sum(np.isnan(feb[:,i,j])) == 0:
            ols = linear_model.LinearRegression()
            model = ols.fit(time, feb[:,i,j])
            lr = model.coef_
        else:
            lr = np.NaN
        lr_feb[i,j] = lr

#3

lr_mar = np.zeros([448,304])
time0 = np.arange(0,43)
time = time0.reshape(-1,1)
for i in range(448):
    for j in range(304):
        if sum(np.isnan(mar[:,i,j])) == 0:
            ols = linear_model.LinearRegression()
            model = ols.fit(time, mar[:,i,j])
            lr = model.coef_
        else:
            lr = np.NaN
        lr_mar[i,j] = lr

#4

lr_apr = np.zeros([448,304])
time0 = np.arange(0,43)
time = time0.reshape(-1,1)
for i in range(448):
    for j in range(304):
        if sum(np.isnan(apr[:,i,j])) == 0:
            ols = linear_model.LinearRegression()
            model = ols.fit(time, apr[:,i,j])
            lr = model.coef_
        else:
            lr = np.NaN
        lr_apr[i,j] = lr

#5

lr_may = np.zeros([448,304])
time0 = np.arange(0,43)
time = time0.reshape(-1,1)
for i in range(448):
    for j in range(304):
        if sum(np.isnan(may[:,i,j])) == 0:
            ols = linear_model.LinearRegression()
            model = ols.fit(time, may[:,i,j])
            lr = model.coef_
        else:
            lr = np.NaN
        lr_may[i,j] = lr

#6

lr_jun = np.zeros([448,304])
time0 = np.arange(0,43)
time = time0.reshape(-1,1)
for i in range(448):
    for j in range(304):
        if sum(np.isnan(jun[:,i,j])) == 0:
            ols = linear_model.LinearRegression()
            model = ols.fit(time, jun[:,i,j])
            lr = model.coef_
        else:
            lr = np.NaN
        lr_jun[i,j] = lr

#7


lr_jul = np.zeros([448,304])
time0_ab = np.arange(0,42)
time_ab = time0_ab.reshape(-1,1)
for i in range(448):
    for j in range(304):
        if sum(np.isnan(jul[:,i,j])) == 0:
            ols = linear_model.LinearRegression()
            model = ols.fit(time_ab, jul[:,i,j])
            lr = model.coef_
        else:
            lr = np.NaN
        lr_jul[i,j] = lr
        
#8

lr_aug = np.zeros([448,304])
time0 = np.arange(0,43)
time = time0.reshape(-1,1)
for i in range(448):
    for j in range(304):
        if sum(np.isnan(aug[:,i,j])) == 0:
            ols = linear_model.LinearRegression()
            model = ols.fit(time, aug[:,i,j])
            lr = model.coef_
        else:
            lr = np.NaN
        lr_aug[i,j] = lr
        
#9

lr_sep = np.zeros([448,304])
time0 = np.arange(0,43)
time = time0.reshape(-1,1)
for i in range(448):
    for j in range(304):
        if sum(np.isnan(sep[:,i,j])) == 0:
            ols = linear_model.LinearRegression()
            model = ols.fit(time, sep[:,i,j])
            lr = model.coef_
        else:
            lr = np.NaN
        lr_sep[i,j] = lr

#10

lr_oct = np.zeros([448,304])
time0 = np.arange(0,43)
time = time0.reshape(-1,1)
for i in range(448):
    for j in range(304):
        if sum(np.isnan(oct[:,i,j])) == 0:
            ols = linear_model.LinearRegression()
            model = ols.fit(time, oct[:,i,j])
            lr = model.coef_
        else:
            lr = np.NaN
        lr_oct[i,j] = lr
        
#11

lr_nov = np.zeros([448,304])
time0 = np.arange(0,43)
time = time0.reshape(-1,1)
for i in range(448):
    for j in range(304):
        if sum(np.isnan(nov[:,i,j])) == 0:
            ols = linear_model.LinearRegression()
            model = ols.fit(time, nov[:,i,j])
            lr = model.coef_
        else:
            lr = np.NaN
        lr_nov[i,j] = lr

#12

lr_dec = np.zeros([448,304])
time0_ab = np.arange(0,42)
time_ab = time0_ab.reshape(-1,1)
for i in range(448):
    for j in range(304):
        if sum(np.isnan(dec[:,i,j])) == 0:
            ols = linear_model.LinearRegression()
            model = ols.fit(time_ab, dec[:,i,j])
            lr = model.coef_
        else:
            lr = np.NaN
        lr_dec[i,j] = lr

cmap2 = mpl.colors.ListedColormap(['#071E46','#072F6B','#08579C','#2171B5','#4292C7','#5AA0CD',
                                  '#78BFD6','#AADCE6','#DBF5FF','#FFFFFF','#FFFFFF','#FFE0E0','#FCBBAA','#FC9272','#FB6A4A','#F03C2B','#CC181E','#A60F14','#780A0F','#5F0000'])

vmin = -2.5
vmax = 2.5
norm = mpl.colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

plt.figure(1)


fig = plt.figure(figsize=(12,10), dpi=100)


#1
ax = plt.subplot(2,3,1,projection=ccrs.NorthPolarStereo())

ax.add_feature(cfeature.OCEAN, facecolor='lightgray')


plotfield = plt.imshow(lr_jan, cmap = cmap2, norm=norm)

ax.set_title('January', fontsize=15)

cbar = fig.colorbar(plotfield, shrink=0.75, orientation = 'vertical',spacing='uniform', extendfrac='auto', norm=norm)
cbar.ax.tick_params(direction = 'in', labelsize = 12)
cbar.set_label('Trend in % per yr', rotation=270, size=13, labelpad=15)
cbar.set_ticks([-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5])

#2
ax = plt.subplot(2,3,2,projection=ccrs.NorthPolarStereo())

ax.add_feature(cfeature.OCEAN, facecolor='lightgray')


plotfield = plt.imshow(lr_feb, cmap = cmap2, norm=norm)

ax.set_title('February', fontsize=15)

cbar = fig.colorbar(plotfield, shrink=0.75, orientation = 'vertical',spacing='uniform', extendfrac='auto', norm=norm)
cbar.ax.tick_params(direction = 'in', labelsize = 12)
cbar.set_label('Trend in % per yr', rotation=270, size=13, labelpad=15)
cbar.set_ticks([-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5])

#3
ax = plt.subplot(2,3,3,projection=ccrs.NorthPolarStereo())

ax.add_feature(cfeature.OCEAN, facecolor='lightgray')


plotfield = plt.imshow(lr_mar, cmap = cmap2, norm=norm)

ax.set_title('March', fontsize=15)

cbar = fig.colorbar(plotfield, shrink=0.75, orientation = 'vertical',spacing='uniform', extendfrac='auto', norm=norm)
cbar.ax.tick_params(direction = 'in', labelsize = 12)
cbar.set_label('Trend in % per yr', rotation=270, size=13, labelpad=15)
cbar.set_ticks([-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5])

#4
ax = plt.subplot(2,3,4,projection=ccrs.NorthPolarStereo())

ax.add_feature(cfeature.OCEAN, facecolor='lightgray')


plotfield = plt.imshow(lr_apr, cmap = cmap2, norm=norm)

ax.set_title('April', fontsize=15)

cbar = fig.colorbar(plotfield, shrink=0.75, orientation = 'vertical',spacing='uniform', extendfrac='auto', norm=norm)
cbar.ax.tick_params(direction = 'in', labelsize = 12)
cbar.set_label('Trend in % per yr', rotation=270, size=13, labelpad=15)
cbar.set_ticks([-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5])

#5
ax = plt.subplot(2,3,5,projection=ccrs.NorthPolarStereo())

ax.add_feature(cfeature.OCEAN, facecolor='lightgray')


plotfield = plt.imshow(lr_may, cmap = cmap2, norm=norm)

ax.set_title('May', fontsize=15)

cbar = fig.colorbar(plotfield, shrink=0.75, orientation = 'vertical',spacing='uniform', extendfrac='auto', norm=norm)
cbar.ax.tick_params(direction = 'in', labelsize = 12)
cbar.set_label('Trend in % per yr', rotation=270, size=13, labelpad=15)
cbar.set_ticks([-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5])

#6
ax = plt.subplot(2,3,6,projection=ccrs.NorthPolarStereo())

ax.add_feature(cfeature.OCEAN, facecolor='lightgray')


plotfield = plt.imshow(lr_jun, cmap = cmap2, norm=norm)

ax.set_title('June', fontsize=15)

cbar = fig.colorbar(plotfield, shrink=0.75, orientation = 'vertical',spacing='uniform', extendfrac='auto', norm=norm)
cbar.ax.tick_params(direction = 'in', labelsize = 12)
cbar.set_label('Trend in % per yr', rotation=270, size=13, labelpad=15)
cbar.set_ticks([-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5])

fig.suptitle('Sea ice concentration trends during 1979-2021',fontsize = 17, fontweight = 'bold', y=1)

fig.tight_layout()
plt.savefig('NSIDC_SIC_nh_jan2jun_tr.png', dpi=500, bbox_inches='tight')
plt.show()








plt.figure(2)


fig = plt.figure(figsize=(12,10), dpi=100)


#7
ax = plt.subplot(2,3,1,projection=ccrs.NorthPolarStereo())

ax.add_feature(cfeature.OCEAN, facecolor='lightgray')


plotfield = plt.imshow(lr_jul, cmap = cmap2, norm=norm)

ax.set_title('July', fontsize=15)

cbar = fig.colorbar(plotfield, shrink=0.75, orientation = 'vertical',spacing='uniform', extendfrac='auto', norm=norm)
cbar.ax.tick_params(direction = 'in', labelsize = 12)
cbar.set_label('Trend in % per yr', rotation=270, size=13, labelpad=15)
cbar.set_ticks([-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5])

#8
ax = plt.subplot(2,3,2,projection=ccrs.NorthPolarStereo())

ax.add_feature(cfeature.OCEAN, facecolor='lightgray')


plotfield = plt.imshow(lr_aug, cmap = cmap2, norm=norm)

ax.set_title('August', fontsize=15)

cbar = fig.colorbar(plotfield, shrink=0.75, orientation = 'vertical',spacing='uniform', extendfrac='auto', norm=norm)
cbar.ax.tick_params(direction = 'in', labelsize = 12)
cbar.set_label('Trend in % per yr', rotation=270, size=13, labelpad=15)
cbar.set_ticks([-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5])

#9
ax = plt.subplot(2,3,3,projection=ccrs.NorthPolarStereo())

ax.add_feature(cfeature.OCEAN, facecolor='lightgray')


plotfield = plt.imshow(lr_sep, cmap = cmap2, norm=norm)

ax.set_title('September', fontsize=15)

cbar = fig.colorbar(plotfield, shrink=0.75, orientation = 'vertical',spacing='uniform', extendfrac='auto', norm=norm)
cbar.ax.tick_params(direction = 'in', labelsize = 12)
cbar.set_label('Trend in % per yr', rotation=270, size=13, labelpad=15)
cbar.set_ticks([-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5])

#10
ax = plt.subplot(2,3,4,projection=ccrs.NorthPolarStereo())

ax.add_feature(cfeature.OCEAN, facecolor='lightgray')


plotfield = plt.imshow(lr_oct, cmap = cmap2, norm=norm)

ax.set_title('October', fontsize=15)

cbar = fig.colorbar(plotfield, shrink=0.75, orientation = 'vertical',spacing='uniform', extendfrac='auto', norm=norm)
cbar.ax.tick_params(direction = 'in', labelsize = 12)
cbar.set_label('Trend in % per yr', rotation=270, size=13, labelpad=15)
cbar.set_ticks([-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5])

#11
ax = plt.subplot(2,3,5,projection=ccrs.NorthPolarStereo())

ax.add_feature(cfeature.OCEAN, facecolor='lightgray')


plotfield = plt.imshow(lr_nov, cmap = cmap2, norm=norm)

ax.set_title('November', fontsize=15)

cbar = fig.colorbar(plotfield, shrink=0.75, orientation = 'vertical',spacing='uniform', extendfrac='auto', norm=norm)
cbar.ax.tick_params(direction = 'in', labelsize = 12)
cbar.set_label('Trend in % per yr', rotation=270, size=13, labelpad=15)
cbar.set_ticks([-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5])

#12
ax = plt.subplot(2,3,6,projection=ccrs.NorthPolarStereo())

ax.add_feature(cfeature.OCEAN, facecolor='lightgray')


plotfield = plt.imshow(lr_dec, cmap = cmap2, norm=norm)

ax.set_title('December', fontsize=15)

cbar = fig.colorbar(plotfield, shrink=0.75, orientation = 'vertical',spacing='uniform', extendfrac='auto', norm=norm)
cbar.ax.tick_params(direction = 'in', labelsize = 12)
cbar.set_label('Trend in % per yr', rotation=270, size=13, labelpad=15)
cbar.set_ticks([-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5])

fig.suptitle('Sea ice concentration trends during 1979-2021',fontsize = 17, fontweight = 'bold', y=1)

fig.tight_layout()
plt.savefig('NSIDC_SIC_nh_jul2dec_tr.png', dpi=500, bbox_inches='tight')
plt.show()

