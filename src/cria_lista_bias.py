#!/usr/bin/env python
# coding: utf-8

# Este código recebe uma lista dos modos de operação
# da câmera escritos em JSON e cria uma lista de imagens
# de bias com os mesmos modos de operação.

#20/02/2020. Denis Varise Bernardes.

##import json
##
##file_dir = r'C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Experimento_monocromador\date_12_3'
##file_name = '\log.txt'
##arq = open(file_dir + '\\' + file_name, 'r')
##lines = arq.read().splitlines()
##arq.close()
##
##arq = open(file_dir + '\\' +  'bias_list.txt', 'w')
##for line in lines:
##    line = json.loads(line)
##    line['t_exp'] = 1e-5
##    line['obs_type'] = 'bias'
##    json.dump(line, arq, sort_keys=True)
##    arq.write('\n')
##arq.close()



from astropy.io import fits
from sys import exit
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

def profiles(image, x_max, y_max): 
    x = np.take(image, x_max, axis=0)    
    y = np.take(image, y_max, axis=1) 
##    plt.plot(x)
##    plt.show(),exit()
    return x, y #these are the horizontal and vertical profiles through the star's centroid


def interpolate_width(axis, max_star_flux):
    half_max = max_star_flux/2
    x = np.linspace(0, len(axis), len(axis))

    # Do the interpolation
    spline = UnivariateSpline(x, axis-half_max, s=0)
    r1, r2 = spline.roots()

    return r2-r1 #this is the FWHM along the specified axis

file_dir = r'C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Experimento_monocromador\date_12_3'
file_name = 'EM_10MHz_texp1.fits'
bias_name = 'bias.fits'
img_data = fits.getdata(file_dir + '\\' + file_name)[0].astype(float)
img_bias = np.asarray(fits.getdata(file_dir + '\\' + bias_name)[0]).astype(float)
img_data -= img_bias
##plt.imshow(img_data, origin='lower left')
##plt.show(),exit()
img_shape = img_data.shape
working_mask = np.ones(img_shape,bool)
ym, xm = np.indices(img_shape, dtype='float32')
x,y = 540, 528
r = np.sqrt((xm - x)**2 + (ym - y)**2)
mask = (r < 100) * working_mask
plt.imshow(mask, origin = 'lower left')
plt.show()
star_flux = img_data[np.where(mask)]
max_star_flux = max(star_flux)
x_max, y_max = np.where(img_data==max_star_flux)
x_max, y_max  = x_max[0], y_max[0]

#img_data = img_data[x_max-50:x_max+50, y_max-50:y_max+50]
print(x_max, y_max)

horizontal, vertical = profiles(img_data, x_max, y_max)
fwhm_x = interpolate_width(horizontal, max_star_flux)
fwhm_y = interpolate_width(vertical, max_star_flux)
print(fwhm_x, fwhm_y)
