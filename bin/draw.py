'''
Author       : linfd 3039562364@qq.com
Date         : 2024-10-21 23:19:39
LastEditTime : 2024-10-22 17:32:17
FilePath     : \computeRelevantGravityParameters\bin\draw.py
Description  : 数据可视化
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.colors as mcolors

def draw():
    data = np.loadtxt('out/T.txt')
    longitudes = data[:, 1]
    latitudes = data[:, 0]
    values = data[:, 2]

    central_longitude = (longitudes.min() + longitudes.max()) / 2
    central_latitude = (latitudes.min() + latitudes.max()) / 2
    
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Orthographic(central_longitude, central_latitude))
    
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linestyle=':')

    # 设置地图范围，只显示数据区域
    ax.set_extent([longitudes.min(), longitudes.max(), latitudes.min(), latitudes.max()], crs=ccrs.PlateCarree())
    
    lon_grid = np.linspace(longitudes.min(), longitudes.max(), 100)
    lat_grid = np.linspace(latitudes.min(), latitudes.max(), 100)
    lon_grid, lat_grid = np.meshgrid(lon_grid, lat_grid)
    
    grid_z = griddata((longitudes, latitudes), values, (lon_grid, lat_grid), method='cubic')
    
    cmap = mcolors.LinearSegmentedColormap.from_list('custom_cmap', ['green', 'yellow', 'red'])
    levels = np.linspace(grid_z.min(), grid_z.max(), 200)
    
    contour = ax.contourf(lon_grid, lat_grid, grid_z, transform=ccrs.PlateCarree(), cmap=cmap, levels=levels, extend='both')
    plt.colorbar(contour, ax=ax, orientation='vertical', shrink=0.5, label='Disturbance Potential T (m$^2$/s$^2$)')
    
    plt.title('Disturbance Potential T Distribution', fontsize=16)

    # 显示网格线
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, linestyle='--', color='gray', alpha=0.5)
    gl.xformatter = LongitudeFormatter(zero_direction_label=True)
    gl.yformatter = LatitudeFormatter()
    gl.top_labels = False
    gl.right_labels = False
    
    plt.show()

if __name__ == '__main__':
    draw()
