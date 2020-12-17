#!/usr/bin/env python
# coding: utf-8

# Esta biblioteca calcula a FA dos CCDs em função da magnitude de um objeto para uma SNR fixa

#25/10/2019. Denis Varise Bernardes.


import os
from sys import exit
import SNR_Calculation_Bib as snrc
import Photon_Flux_Calc_Bib as pfc
import json
import numpy as np
import matplotlib.pyplot as plt
import Optimize_Acquisition_Rate_Bib as oar
import Optimize_Signal_Noise_Ratio_Bib as osnr
import locale


class profile_FA_MAG:

    def __init__(self, snr, obj_magnitude, FA, si_modes, bin_modes):
        self.snr = snr
        self.obj_mag = obj_magnitude
        self.si_modes = si_modes
        self.bin_modes = bin_modes
        self.FA = FA
        self.star_flux = 0
        self.hats24_magnitude = 12.25                
        self.hats24_flux = 56122.295000000006 #e-/s
        self.sky_flux = 12.298897076737294 #e-/pix/s
        self.n_pixels = 305
        self.bias_level = 500


    def calc_star_sky_flux_from_magnitude(self):                
        aux = 10**((self.hats24_magnitude - self.obj_mag)/2.5)
        self.star_flux = self.hats24_flux*aux


    def get_best_mode(self):
        OAR =  oar.OptimizeAcquisitionRate(self.FA, self.si_modes, self.bin_modes)
        OAR.determine_operation_modes()        
        obj_lista_modos_FA = OAR.read_MOB_obj()        
        #OAR.print_MOB_list(),exit()

        OSNR = osnr.OptSignalNoiseRation(self.snr, 9916,  -60,  self.n_pixels,  self.sky_flux,  self.star_flux,  self.bias_level)         
        OSNR.write_MOB_obj(obj_lista_modos_FA)       
        OSNR.determine_operation_modes_minimun_SNR()                                                           
        obj_list_of_modes = OSNR.read_MOB_obj()              
        #OSNR.print_MOB_list(),exit()

        OAR.write_MOB_obj(obj_list_of_modes)
        self.best_mode = OAR.determine_fastest_operation_mode()
        

    def get_best_FA(self):
        return self.best_mode['max_acq_rate']
        

  
magnitudes = np.linspace(5, 20, 15)
colors = ['b', 'r', 'g', 'k']
i=0
for snr in [1,10,100,1000]:
    new_mags = []
    max_fas = []  
    for mag in magnitudes:
        PFAMAG = profile_FA_MAG(snr=snr, obj_magnitude=mag, FA = 0.01, si_modes=[256, 512], bin_modes=[2])
        PFAMAG.calc_star_sky_flux_from_magnitude()
        PFAMAG.get_best_mode()
        try:
            max_fa = PFAMAG.get_best_FA()
            max_fas.append(max_fa)
            new_mags.append(mag)
        except: break
    plt.semilogy(new_mags, max_fas, '-', c=colors[i], label='SNR = '+str(snr))

    new_mags = []
    max_fas = []  
    for mag in magnitudes:
        PFAMAG = profile_FA_MAG(snr=snr, obj_magnitude=mag, FA = 0.01, si_modes=[1024], bin_modes=[1])
        PFAMAG.calc_star_sky_flux_from_magnitude()
        PFAMAG.get_best_mode()
        try:
            max_fa = PFAMAG.get_best_FA()
            max_fas.append(max_fa)
            new_mags.append(mag)
        except: break
    plt.semilogy(new_mags, max_fas, '--', c=colors[i])      
    i+=1

fontsize = 15
locale.setlocale(locale.LC_NUMERIC, "de_DE")
plt.rcdefaults()
plt.rcParams['axes.formatter.use_locale'] = True
plt.rcParams.update({'font.size': fontsize})
plt.xlabel('Magnitude I', fontsize = fontsize)
plt.ylabel('Acquisition Rate (fps)', fontsize = fontsize)
plt.legend(loc='lower left')
plt.xticks(rotation=0, fontsize=fontsize)
plt.yticks(rotation=0, fontsize=fontsize)
plt.show()

