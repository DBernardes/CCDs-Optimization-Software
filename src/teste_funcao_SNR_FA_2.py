
lista_modos = [
{'em_mode': 1, 'em_gain': 2, 'hss':  1, 'preamp': 1, 'binn': 1, 'sub_img': 1024, 'min_t_exp': 80, 'max_t_exp': 100}
]


import SNR_Calculation_Bib as snrc
import Optimize_Acquisition_Rate_Bib as oar
import Optimize_Signal_Noise_Ratio_Bib as osnr
import Acquisition_Rate_Calculation_Bib as arc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import Modos_Operacao_Bib as mob
import collections
import json
from sys import exit
from hyperopt import hp, tpe, rand, Trials, fmin
from hyperopt.pyll.stochastic import sample
from hyperopt.pyll import scope
from copy import copy
import gc, os
import random as rd
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
from math import ceil



def function_fa(parameters = []):    
    t_exp = parameters[0]
    em_mode = parameters[1]
    em_gain = parameters[2]
    hss = parameters[3]
    preamp = parameters[4]
    binn = parameters[5]
    sub_img = parameters[6]
    CCD_temp = parameters[7]
    sky_flux = parameters[8]
    star_flux = parameters[9]
    n_pixels_star = parameters[10]
    serial_number = parameters[11]    
    max_fa = parameters[14]
    min_fa = parameters[15]        

    ARC = arc.AcquisitionRateCalc()
    ARC.write_operation_mode(em_mode, hss, binn, sub_img, t_exp)
    ARC.seleciona_t_corte()
    ARC.calc_acquisition_rate()
    acq_rate = ARC.acquisition_rate      
    norm_acq_rate = (acq_rate-min_fa)/(max_fa - min_fa)
    return norm_acq_rate

def function_snr(parameters = []):    
    t_exp = parameters[0]
    em_mode = parameters[1]
    em_gain = parameters[2]
    hss = parameters[3]
    preamp = parameters[4]
    binn = parameters[5]
    sub_img = parameters[6]
    ccd_temp = parameters[7]
    sky_flux = parameters[8]
    star_flux = parameters[9]
    n_pix_star = parameters[10]
    serial_number = parameters[11]
    max_snr = parameters[12]
    min_snr = parameters[13]
    snr_target = parameters[16]    
    SNRC = snrc.SignalToNoiseRatioCalc(t_exp = t_exp,
                              em_mode = em_mode,
                              em_gain = em_gain,
                              hss = hss,
                              preamp = preamp,
                              binn = binn,
                              ccd_temp = ccd_temp,
                              sky_flux = sky_flux,
                              star_flux = star_flux,
                              n_pix_star = n_pix_star,
                              serial_number = serial_number)    
    SNRC.calc_RN()
    SNRC.calc_DC()
    SNRC.calc_SNR()
    snr = SNRC.get_SNR()    
    if snr < snr_target: snr = 0   
    norm_snr = (snr - min_snr)/(max_snr - min_snr)
    return norm_snr


def function(parameters = []):
    return (function_snr(parameters) * function_fa(parameters))


aux = 10**((20 - 2)/2.5)
star_flux = 56122.295000000006*aux
sky_flux = 12.298897076737294
n_pix = 305

mode = lista_modos[0]
ARC = arc.AcquisitionRateCalc()
ARC.write_operation_mode(em_mode = mode['em_mode'], hss = mode['hss'], binn = mode['binn'], sub_img = mode['sub_img'], t_exp = mode['min_t_exp'])
ARC.seleciona_t_corte()
ARC.calc_acquisition_rate()
max_fa = float(ARC.acquisition_rate)
ARC.write_operation_mode(em_mode = mode['em_mode'], hss = mode['hss'], binn = mode['binn'], sub_img = mode['sub_img'], t_exp = mode['max_t_exp'])
ARC.seleciona_t_corte()
ARC.calc_acquisition_rate()
min_fa = float(ARC.acquisition_rate)




SNRC = snrc.SignalToNoiseRatioCalc(em_mode = mode['em_mode'], hss = mode['hss'], binn = mode['binn'], em_gain = mode['em_gain'], preamp = mode['preamp'], t_exp = mode['max_t_exp'], ccd_temp = -60, sky_flux = sky_flux , star_flux = star_flux, n_pix_star = n_pix, serial_number = 9916)
SNRC.set_gain_value()
SNRC.calc_RN()
SNRC.calc_DC()
SNRC.calc_SNR()
max_snr = SNRC.get_SNR()
SNRC = snrc.SignalToNoiseRatioCalc(em_mode = mode['em_mode'], hss = mode['hss'], binn = mode['binn'], em_gain = mode['em_gain'], preamp = mode['preamp'], t_exp = mode['min_t_exp'], ccd_temp = -60, sky_flux = sky_flux , star_flux = star_flux, n_pix_star = n_pix, serial_number = 9916)
SNRC.set_gain_value()
SNRC.calc_RN()
SNRC.calc_DC()
SNRC.calc_SNR()
min_snr = SNRC.get_SNR()





max_output = 0
func_output = [[],[],[],[]]
n_repeat = ceil(1/len(lista_modos))        
if n_repeat < 1: n_repeat = 1        
for i in range(n_repeat):
    for mode in lista_modos:
        for t_exp in np.linspace(mode['min_t_exp'],mode['max_t_exp'],20):
            em_gain = 1            
            if mode['em_mode'] == 1:
                em_gain = mode['em_gain']#rd.uniform(2,mode['em_gain'])            
            new_mode = [t_exp, mode['em_mode'], em_gain, mode['hss'], mode['preamp'], mode['binn'], mode['sub_img'], -60, sky_flux, star_flux, n_pix, 9916, max_snr, min_snr, max_fa , min_fa, 0]        
            output = function(new_mode)            
            func_output[0].append(output)                
            func_output[1].append(t_exp)
            func_output[2].append(function_snr(new_mode))
            func_output[3].append(function_fa(new_mode))



##x1 = np.linspace(0.9,1,10)
##x2 = np.linspace(1,0.8,10)
##plt.plot(x1*x2)
##plt.show(),exit()

plt.figure(figsize = (12,5))
plt.subplot(131)
plt.plot(func_output[1], func_output[2], marker = 'o', color='b')
plt.xlabel('Tempo de exposição (s)')
plt.ylabel('SNR')
plt.xlim(mode['min_t_exp']-0.1, mode['max_t_exp']+0.1)

plt.subplot(132)
plt.plot(func_output[1], func_output[3],  marker = 'o', color='b')
plt.xlabel('Tempo de exposição (s)')
plt.ylabel('FA')
plt.xlim(mode['min_t_exp']-0.1, mode['max_t_exp']+0.1)

plt.subplot(133)
plt.plot(func_output[1], func_output[0], marker =  'o', color='b')
plt.xlabel('Tempo de exposição (s)')
plt.ylabel('SNR x FA')
plt.xlim(mode['min_t_exp']-0.1, mode['max_t_exp']+0.1)

plt.show()








