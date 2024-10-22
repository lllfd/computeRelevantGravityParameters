'''
Author       : linfd 3039562364@qq.com
Date         : 2024-10-21 23:19:39
LastEditTime : 2024-10-22 10:58:06
FilePath     : \computeRelevantGravityParameters\bin\draw.py
Description  : 
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

def draw():
    data = np.loadtxt('out/DeltaG.txt')
    longitudes = data[:, 1]
    latitudes = data[:, 0]
    values = data[:, 2]
    lon_grid = np.linspace(longitudes.min(), longitudes.max(), 100)
    lat_grid = np.linspace(latitudes.min(), latitudes.max(), 100)
    lon_grid, lat_grid = np.meshgrid(lon_grid, lat_grid)
    grid_z = griddata((longitudes, latitudes), values, (lon_grid, lat_grid), method='cubic')
    plt.figure(figsize=(10, 8))
    im = plt.pcolormesh(lon_grid, lat_grid, grid_z, cmap='hot', shading='auto')
    plt.colorbar(im)  
    plt.grid(True)
    plt.title('Colored 2D Image with Values')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.show()

if __name__ == '__main__':
    draw()
