#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: nastaran
Email: nastaran.mrad@gmail.com
This code read the phantom data, create masks and plot the data and masks
"""

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt


#%%
def create_circular_mask(h, w, center1=None, center2 = None, center3 = None, radius=None):
    """
     *** input: the width and height of the image, center and the radius for drawing circular masks
     *** output: masks for blood, lipid, and calcium 
    """
    if center1 is None: # use the middle of the image
        center = (int(w/2), int(h/2))
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center1 = np.sqrt((X - center1[0])**2 + (Y-center1[1])**2)
    dist_from_center2 = np.sqrt((X - center2[0])**2 + (Y-center2[1])**2)
    dist_from_center3 = np.sqrt((X - center3[0])**2 + (Y-center3[1])**2)

    mask1 = dist_from_center1 <= radius
    mask2 = dist_from_center2 <= radius
    mask3 = dist_from_center3 <= radius
    
    return mask1, mask2, mask3

#%%
def create_final_mask(m1,m2,m3):
    """
     put three masks in an image and make them with
     *** input: created masks, i.e., the output of create_circular_mask function
     *** output: binary mask
    """
    invert_m1 = 1-m1
    invert_m2 = 1-m2
    invert_m3 = 1-m3
       
    masks = np.multiply((np.multiply(invert_m1, invert_m2)),invert_m3)
    finalMask = 1 - masks
    return finalMask

#%%
def color_masks(finalMask, baseDIR, dataType, dataName):
    """
    *** input: binary mask
    *** output: RGB mask
    
    """
    for i in range(finalMask.shape[0]):
        for j in range(finalMask.shape[1]):
            if 26<=j<=42 and finalMask[i][j] == 1:
                finalMask[i][j] = 255
            if 86<=j<=102 and finalMask[i][j] == 1:
                finalMask[i][j] = 125
            elif 56<=j<=72 and finalMask[i][j] == 1:
                finalMask[i][j] = 62
    sio.savemat(baseDIR + 'mask_'+dataType + dataName + '.mat',{'mask':finalMask})
    plt.figure()
    plt.imshow(finalMask)
    return finalMask
#%%    
def readData(baseDIR, dataType, dataName):
    """
        read data and plot     
    """
    matContent = sio.loadmat(baseDIR + dataType + dataName + '.mat')    
    data = matContent[dataType]
    data_cropped = data[200:460,0:260,120]  
    ## plot the data
    plt.figure()
    plt.style.use('grayscale')
    plt.imshow(X = 20*np.log10((data_cropped)),vmin = 70, vmax = 110)
  
#%%   
if __name__ == '__main__':
    data_type = 'US'
    data_name = '_blood_down_L12_3V_20190318_nir'
    baseDIR = 'data path'
    
    readData(baseDIR, data_type, data_name)
    mask1, mask2, mask3 = create_circular_mask(260, 260, center1 = [34, 148], center2 = [94, 148], center3 = [64, 88], radius = 8)
    final_mask = create_final_mask(mask1, mask2, mask3)
    color_masks(final_mask, baseDIR, data_type, data_name)


"""
plot the circle
fig , ax = plt.subplots()
ax.set(xlim = (0,260), ylim= (260,0))

#plt.axis([0, 150, 0, 260])
c = plt.Circle((34, 148), 10, color='r')
ax.add_artist(c)

"""

    