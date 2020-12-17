
lista_modos = [
{'em_mode': 1, 'em_gain': 2, 'hss':  1, 'preamp': 1, 'binn': 1, 'sub_img': 1024, 'min_t_exp': 40, 'max_t_exp': 100}
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
    mean_fa = parameters[13]        

    ARC = arc.AcquisitionRateCalc()
    ARC.write_operation_mode(em_mode, hss, binn, sub_img, t_exp)
    ARC.seleciona_t_corte()
    ARC.calc_acquisition_rate()
    acq_rate = ARC.acquisition_rate      
    norm_acq_rate = acq_rate/mean_fa    
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
    mean_snr = parameters[12]    
    snr_target = parameters[14]    
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
    norm_snr = snr/mean_snr    
    return norm_snr


def function(parameters = []):
    return function_snr(parameters)*function_fa(parameters)

max_output = 0
func_output = [[],[],[],[]]
n_repeat = ceil(1/len(lista_modos))        
if n_repeat < 1: n_repeat = 1        
for i in range(n_repeat):
    for mode in lista_modos:
        for t_exp in np.linspace(mode['min_t_exp'],mode['max_t_exp'],20):
            em_gain = 1            
            if mode['em_mode'] == 1:
                em_gain = 50#rd.uniform(2,mode['em_gain'])            
            new_mode = [t_exp, mode['em_mode'], em_gain, mode['hss'], mode['preamp'], mode['binn'], mode['sub_img'], -60, 22.85988123378231, 3099.5717251070682, 2788, 9916, 85.0300915874527, 0.025 , 0]        
            output = function(new_mode)            
            func_output[0].append(output)                
            func_output[1].append(t_exp)
            func_output[2].append(function_snr(new_mode))
            func_output[3].append(function_fa(new_mode))

#print(np.mean(func_output[0])), exit()
##top_10_idx = np.argsort(func_output[0])[-11:]
##top_10_modes = [func_output[3][i] for i in top_10_idx]
##top_10_output = [func_output[0][i] for i in top_10_idx]
##
##
##for i in range(len(top_10_output)):
##    index = len(top_10_output) - i - 1
##    mode = top_10_modes[index]
##    em_mode = 'Conv'
##    if mode['em_mode'] == 1: em_mode = 'EM'
##    print(em_mode, func_output[2][index], mode['hss'], mode['preamp'], mode['binn'], mode['sub_img'], func_output[1][index],top_10_output[index])
##

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

exit()

fig = plt.figure()
list_fake2Dlines = []
list_labels = []
ax = fig.add_subplot((111), projection='3d')
ax.scatter(func_output[1], func_output[2], func_output[0], c='blue', marker='o', alpha=0.5)
fake2Dline1 = mpl.lines.Line2D([0],[0], linestyle="none", c='blue', marker = 'o')
list_fake2Dlines.append(fake2Dline1)
list_labels.append(r'30 MHz')
ax.set_xlabel('Tempo de exposição (s)')
ax.set_ylabel('Ganho EM')
ax.set_zlabel('SNR*FA')
ax.legend(list_fake2Dlines, list_labels, numpoints = 1, loc='upper left')
plt.show()

##list_operations_modes = lista_modos
##n_repeat = 1
##losses_SNR = []
##losses_FA = []
##for i in range(n_repeat):
##    for mode in list_operations_modes:            
##        em_mode = mode['em_mode']
##        hss = mode['hss']
##        binn = mode['binn']
##        sub_img = mode['sub_img']        
##        preamp = mode['preamp']
##
##        for t_exp in np.random.uniform(mode['min_t_exp'], mode['max_t_exp'],50):
##            em_gain = 0
##            SNRC = snrc.SignalToNoiseRatioCalc(t_exp, em_mode, em_gain, hss, preamp, binn, -60, 22.85988123378231, 3099.5717251070682, 2788, 9916)
##            SNRC.set_gain_value()
##            SNRC.calc_RN()
##            SNRC.calc_DC()
##            SNRC.calc_SNR()
##            losses_SNR.append(SNRC.get_SNR())           
##            
##            #print(em_mode, hss, binn, sub_img, t_exp)
##            ARC = arc.AcquisitionRateCalc()
##            ARC.write_operation_mode(em_mode, hss, binn, sub_img, t_exp)
##            ARC.seleciona_t_corte()
##            ARC.calc_acquisition_rate()
##            losses_FA.append(float(ARC.return_acquisition_rate()))
##mean_snr =np.mean(losses_SNR)
##mean_fa = np.mean(losses_FA)
###print(losses_SNR, losses_FA)
##print(mean_snr, mean_fa),exit()






