#!/usr/bin/env python
#coding: utf-8
#Denis Varise Bernardes.
#12/12/2019.


import SNR_Calc as snrc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import random as rd
import Modos_Operacao_Bib as MOB
import collections
import os, json
from sys import exit
from hyperopt import hp, tpe, rand, Trials, fmin
from hyperopt.pyll.stochastic import sample
from hyperopt.pyll import scope
from copy import copy



class OptSignalNoiseRation:

    def __init__(self, snr_target, serial_number, ccd_temp, n_pix_star, sky_flux, star_flux, bias_level):        
        self.MOB = MOB.ModosOperacao()
        self.space = []
        self.best_mode = {}        
        self.list_all_modes  = []
        self.filtered_list = []
        self.best_sub_img = []
        self.new_list = []        
        self.hss = [[],[]]
        self.binn = [[],[]]
        self.sub_img = []        
        self.ccd_temp = ccd_temp
        self.serial_number = serial_number
        self.gain = 0 
        self.dark_noise = 0
        self.set_dc() 
        self.snr_target = snr_target
        
        self.sky_flux = sky_flux #e-/pix/s                
        self.star_flux = star_flux #e-/s        
        self.n_pix_star = n_pix_star
        self.bias_level = bias_level
       
        
    def write_mode_to_MOB_class(self, em_mode, em_gain, hss, preamp, binn, sub_img, t_exp):
        #Write a mode in the list of the modes class
        self.MOB.write_mode(em_mode, em_gain, hss, preamp, binn, sub_img, t_exp)

   
    def print_MOB_list(self):        
        lista = self.MOB.get_list_of_modes()
        for modo in lista:
            print(modo)

    def write_MOB_obj(self, obj):
        #Write a object with the list of modes in the class
        self.MOB = obj

    def read_MOB_obj(self):
        #Reads the object with the list of modes
        return self.MOB


    def set_gain(self, em_mode, hss, preamp):
        #Selectes the CCD gain based on the provided operation mode
        gain = 0
        if em_mode == 1:
            if hss == 30:
                if preamp == 1:
                    gain = 17.2
                if preamp == 2:
                    gain = 5.27
            if hss == 20:
                if preamp == 1:
                    gain = 16.4
                if preamp == 2:
                    gain = 4.39
            if hss == 10:
                if preamp == 1:
                    gain = 16.0
                if preamp == 2:
                    gain = 3.96
            if hss == 1:
                if preamp == 1:
                    gain = 15.9
                if preamp == 2:
                    gain = 3.88
        else:
            if hss == 1:
                if preamp == 1:
                    gain = 3.37
                if preamp == 2:
                    gain = 0.8
            if hss == 0.1:
                if preamp == 1:
                    gain = 3.35
                if preamp == 2:
                    gain = 0.8
        self.gain = gain  


    def set_dc(self):
        #Calculates the dark current based ib the CCD temperature.
        #These used equations model the dark current of the SPARC4 CCDs,
        #and can be found in: D V Bernardes et al 2018 PASP 130 095002
        T = self.ccd_temp
        if self.serial_number == 9914:
            self.dark_noise = 24.66*np.exp(0.0015*T**2+0.29*T) 
        if self.serial_number == 9915:
            self.dark_noise = 35.26*np.exp(0.0019*T**2+0.31*T)
        if self.serial_number == 9916:
            self.dark_noise = 9.67*np.exp(0.0012*T**2+0.25*T)
        if self.serial_number == 9917:
            self.dark_noise = 5.92*np.exp(0.0005*T**2+0.18*T)


    def calc_max_em_gain(self, max_t_exp, min_t_exp):
        #Calculation of the maximum EM gain allowed as a function of the CCD operating mode.
        #This calculation takes into account the maximum amount of 100 photons per pixel for which
        #the EM mode is better than the Conventional one.
        sky_flux = self.sky_flux
        star_flux = self.star_flux
        n_pix = self.n_pix_star
        dn = self.dark_noise
        gain = self.gain
        max_fotons = 100
        bias = self.bias_level
        max_ADU = (2**16)*0.8
        
        aux = (sky_flux + star_flux/n_pix + dn) * max_t_exp / gain     
        max_em_gain = (max_ADU - bias)/aux

        #if the photons/pix ir bigger than 100 photons, it is calculated the EM gain and the
        #exposure time values that accomplish this limit
        if aux > max_fotons:
            max_em_gain = (max_ADU - bias)/(max_fotons/gain)
            max_t_exp = max_fotons/(sky_flux + star_flux/n_pix + dn)
        #If the max exposure time found in the previous step is smaller than the min exposure time,
        #this mode is discarded
        if max_t_exp < min_t_exp:
            max_em_gain = 0
            max_t_exp = 0       
        return max_em_gain, max_t_exp


        

    def determine_operation_modes_minimun_SNR(self):
        #This function determines the operating modes that meet a minimum SNR.
        #For each mode, the maximum allowed EM gain is calculated.
        #This gain is used to calculate the minimum exposure time allowed to reach the SNR.
        #The selected modes are passed to the MOB object mode list

        #iterates each mode of the list of selected mode
        for mode in self.MOB.get_list_of_modes():            
            em_mode = mode['em_mode']
            hss = mode['hss']
            binn = mode['binn']
            max_t_exp = mode['max_t_exp']            
            sub_img = mode['sub_img']
            for preamp in [1,2]:                
                max_em_gain = 0
                if em_mode == 1:
                    #calculates the CCD gain
                    self.set_gain(em_mode, hss, preamp)
                    #calculates the EM gain
                    max_em_gain, max_t_exp = self.calc_max_em_gain(mode['max_t_exp'], mode['min_t_exp'])
                    #if the returned EM gain is 0, this mode is discarded
                    if max_em_gain == 0:
                        print('The number of photons per pixel is above 100. This mode was rejected.')                        
                        continue
                    #limits the maximum EM gain to 300x
                    if max_em_gain > 300: max_em_gain = 300
                #Starts the SNR library
                SNRC = snrc.SignalToNoiseRatioCalc(max_t_exp,
                                                   em_mode,
                                                   max_em_gain,
                                                   hss,
                                                   preamp,
                                                   binn,
                                                   self.ccd_temp,
                                                   self.sky_flux,
                                                   self.star_flux,
                                                   self.n_pix_star,
                                                   self.serial_number)
                SNRC.set_gain_value()
                SNRC.calc_RN()
                SNRC.calc_DC()
                #Calculates the minimum exposure time needed to achieve the provided SNR
                min_t_exp = SNRC.calc_minimun_texp_provided_SNR(self.snr_target)
                #limits the minimum exposure time in 0.00001 s
                if min_t_exp < 1e-5: min_t_exp = 1e-5
                #Add the selected mode in a new list of modes in the MOB class in the dictionary form
                if min_t_exp <= max_t_exp:
                    dic = {'em_mode':em_mode,
                           'em_gain':max_em_gain,
                           'hss':hss,
                           'preamp':preamp,
                           'binn':binn,
                           'sub_img':sub_img,
                           'min_t_exp':min_t_exp,
                           'max_t_exp':max_t_exp}                                        
                    self.filtered_list.append(dic)                
        self.MOB.write_list_of_modes(self.filtered_list)
        
        

    def duplicate_list_of_modes_for_PA12(self):
        #Creates a list of allowed modes. In this step, repeated sub_img modes are discarded.
        #However, it is still necessary to remove the modes with overlapping maximum texp.
        #Each iteration receives one of the preamp values: 1 or 2.
        self.list_all_modes = self.MOB.get_list_of_modes()           
        for preamp in [1,2]: 
            for mode in self.list_all_modes:
                new_mode = {'em_mode':mode['em_mode'], 'hss':mode['hss'], 'preamp':preamp, 'binn':mode['binn'], 'max_t_exp':mode['max_t_exp'], 'min_t_exp':mode['min_t_exp']}
                if new_mode not in self.filtered_list: self.filtered_list.append(new_mode)        
        
        


    def remove_repeat_modes(self):
        #This function eliminates repeated modes
        #which have overlapping values of exposure time.
        new_list = []
        for i in range(len(self.filtered_list)-1):          
            mode_before = self.filtered_list[i]            
            mode_after = self.filtered_list[i+1]           
            if (mode_before['em_mode'] == mode_after['em_mode']):
                if(mode_before['hss'] == mode_after['hss']):
                    if (mode_before['binn'] == mode_after['binn']):
                        if(mode_before['preamp'] == mode_after['preamp']):
                            if (mode_before['max_t_exp'] < mode_after['max_t_exp']):
                                new_list.append(mode_before)
        for mode in new_list:
            self.filtered_list.remove(mode)        


    def calc_min_snr(self):
        #This function calculates the maximum and minimum
        #values of the SNR for the selected operation modes.
        max_snr = 0
        min_snr = 1e5
        best_mode = {}
        #iterates each modes of the list of selected modes
        for mode in self.MOB.get_list_of_modes():            
            em_mode = mode['em_mode']
            hss = mode['hss']
            binn = mode['binn']            
            preamp = mode['preamp']
            min_t_exp = mode['min_t_exp']
            max_em_gain = 0
            #if the EM mode is conventional, EM gain = 1
            min_em_gain = 1                        
            if mode['em_mode'] == 1:
                #if the EM mode is EM, EM gain = 2
                min_em_gain = 2
                #calculates the CCD gain
                self.set_gain(em_mode, hss, preamp)
                #Calculates the maximum EM gain
                max_em_gain, max_t_exp = self.calc_max_em_gain(mode['max_t_exp'], mode['min_t_exp'])
                #if the returned EM gain is 0, this mode is discarded
                if max_em_gain == 0:
                    print('The number of photons per pixel is above 100. This mode was rejected.')
                    continue
                #limites the EM gain in 300x
                if max_em_gain > 300: max_em_gain = 300
            #starts the class of the SNR calculation
            SNRC = snrc.SignalToNoiseRatioCalc(min_t_exp,
                                               em_mode,
                                               max_em_gain,
                                               hss,
                                               preamp,
                                               binn,
                                               self.ccd_temp,
                                               self.sky_flux,
                                               self.star_flux,
                                               self.n_pix_star,
                                               self.serial_number)
            SNRC.set_gain_value()
            SNRC.calc_RN()
            SNRC.calc_DC()
            SNRC.calc_SNR()
            #Calculates the SNR for the respective mode
            snr = SNRC.get_SNR()
            #if the obtained SNR value is smaller than the current smalest SNR
            #this new mode is selected as the new minimum
            if snr < min_snr: min_snr = snr            
        self.best_mode = best_mode        
        return min_snr
        


    def calc_best_mode(self):
        #Calculates the best value of the SNR based on the selected operation modes.
        best_snr = 0
        best_mode = {}
        #This function eliminates repeated modes
        #which have overlapping values of exposure time.
        self.remove_repeat_modes()
        #iterates each modes of the list of selected modes        
        for mode in self.filtered_list:            
            em_mode = mode['em_mode']
            hss = mode['hss']
            binn = mode['binn']            
            preamp = mode['preamp']
            max_em_gain = 0
            max_t_exp = mode['max_t_exp']
            if mode['em_mode'] == 1:
                #Calculates the CCD gain
                self.set_gain(em_mode, hss, preamp)
                #Calculates the maximim EM gain
                max_em_gain, max_t_exp = self.calc_max_em_gain(mode['max_t_exp'], mode['min_t_exp'])
                #if the returned EM gain is 0, this mode is discarded
                if max_em_gain == 0:
                    print('The number of photons per pixel is above 100. This mode was rejected.')
                    continue
                #limits the maximum EM gain in 300x
                if max_em_gain > 300: max_em_gain = 300
            #Starts the class of the SNR calculation            
            SNRC = snrc.SignalToNoiseRatioCalc(max_t_exp,
                                               em_mode,
                                               max_em_gain,
                                               hss, preamp,
                                               binn,
                                               self.ccd_temp,
                                               self.sky_flux,
                                               self.star_flux,
                                               self.n_pix_star,
                                               self.serial_number)
            SNRC.set_gain_value()
            SNRC.calc_RN()
            SNRC.calc_DC()
            SNRC.calc_SNR()
            #Calculates the maximum SNR
            snr = SNRC.get_SNR()            
            #If the calculated SNR is greater than the current best SNR,
            #this mode is selected as the new best mode
            if snr > best_snr:                
                best_snr = snr
                best_mode = mode
                best_mode['em_gain'] = max_em_gain
                best_mode['t_exp'] = max_t_exp
                best_mode['SNR'] = snr               
                del best_mode['max_t_exp']
                del best_mode['min_t_exp']
        self.best_mode = best_mode        
        return best_snr
      

 

    def print_filtered_list_of_modes(self):
        for mode in self.filtered_list:
            print(mode)


            
    def find_sub_img_best_mode(self):
        #There are some cases where the optimum mode has more than one sub_image mode.
        #Thus, this function was create to deal with them.   
        for mode in self.list_all_modes:            
            if mode['em_mode'] == self.best_mode['em_mode']:
                if mode['hss'] == self.best_mode['hss']:
                    if mode['binn'] == self.best_mode['binn']:                        
                        if self.best_mode['t_exp'] <= mode['max_t_exp']:
                            self.best_sub_img.append(mode['sub_img'])       
        
    

    def print_best_mode(self):                 
        if self.best_mode['em_mode'] == 1:
            print('\nEM Mode')
            print('-------')
            print('EM gain: ', self.best_mode['em_gain'])
        else:
            print('\nConventional Mode')
            print('-----------------')            
        print('Exposure time (s): ', self.best_mode['t_exp'])            
        print('HSS: ', self.best_mode['hss'])
        print('Preamp: ', self.best_mode['preamp'])
        print('Binning: ', self.best_mode['binn'])        
        print('Sub image: ', max(self.best_sub_img))
        print('\nBest SNR: ', self.best_mode['SNR'])
        if self.snr_target > self.best_mode['SNR']: print('\nIt was not possible to reach the provided SNR.')
        


    def export_optimal_setup(self, img_directory, file_base_name, star_radius, obj_coords):
        #This functions exports the obtainde best mode to a .txt file
        #To acocmplish this, it is used the json format
        dic={}        
        if self.best_mode['em_mode'] == 1:
            dic['em_mode'] = 'EM'            
        else:
            dic['em_mode'] = 'CONV'
        dic['em_gain'] = self.best_mode['em_gain']
        dic['t_exp'] = self.best_mode['t_exp']
        dic['hss'] = self.best_mode['hss']
        dic['preamp'] = self.best_mode['preamp']
        dic['bin'] = self.best_mode['binn']
        dic['sub_img'] = max(self.best_sub_img)
        dic['output'] = self.best_mode['SNR']
        dic['obs_type'] = 'object'
        dic['img_name'] = file_base_name + '_OPTMODE.fits'
        dic['star_radius'] = star_radius
        try: dic['obj_coords'] = '(%i,%i)'%(obj_coords[0], obj_coords[1])
        except:1
        dic['kinetic_series_length'] = 1
        
        file_name = img_directory + file_base_name + '_OPTSETUP.txt'
        with open(file_name, 'w') as arq:
            json.dump(dic, arq, indent = 4, sort_keys=True)
            arq.close()






        
