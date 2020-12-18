#!/usr/bin/env python
# coding: utf-8

# Esta biblioteca gerencia o processo de otimização dos modos de
# operação dos CCDs. 

#10/01/2020. Denis Varise Bernardes.

import Optimize_SNR as osnr
import Optimize_AR as oar
import Optimize_SNR_and_AR as osnrar
import SNR_Calc as snrc
import AR_Calc as arc
import Modos_Operacao_Bib as mob
import FWHM
from sys import exit
from useful_functions import get_obs_setup
from math import ceil
from hyperopt import rand
import os
import numpy as np
import random as rd
from copy import copy




class Optimize_Operation_Mode:

    def __init__(self, img_dir, algorithm):
        snr, acq_rate, obj_magnitude, sub_img_modes, binn_modes, serial_number, ccd_temp, max_evals, file_base_name, export_arq, export_loss, export_bias, use_pre_img, img_name, obj_coords, bias_name, sky_radius = get_obs_setup(img_dir)
        self.algorithm = algorithm
        self.snr = snr
        self.acq_rate = acq_rate
        self.obj_magnitude = obj_magnitude
        self.max_evals = max_evals        
        self.file_base_name = file_base_name
        self.export_arq = export_arq
        self.export_loss = export_loss        
        self.export_bias = export_bias
        
        self.sub_img_modes = sub_img_modes
        self.binn_modes = binn_modes
        self.serial_number = serial_number
        self.ccd_temp = ccd_temp        
        
        self.sky_flux = 12.298897076737294 #e-/pix/s
        self.hats24_flux = 56122.295000000006 #e-/s
        self.star_flux = 0
        self.n_pix_star = 305
        self.hats24_magnitude = 12.25
        
        self.use_pre_img = use_pre_img
        self.img_dir = img_dir
        self.img_name = img_name       
        self.obj_coords = obj_coords
        self.bias_name = bias_name 
        self.sky_radius = sky_radius
        self.star_radius = 0
        self.bias_level = 500        

        self.MOB = []

       

    def verify_provides_modes(self):
        #Verifies if the provided sub_img and bin modes are allowed        
        for sub_img in self.sub_img_modes:
            if sub_img not in [1024,512,256]:
                print('\nInvalid sub-image mode! [%i]'%sub_img)                
                exit()
        for binn in self.binn_modes:
            if binn not in [1,2]:
                print('\nInvalid binning mode! [%i]'%binn)                
                exit()



    def calc_star_flux(self):
        #Select the method for the star flux calculation bases on the information provided to the code.
        if self.use_pre_img == 'y': self.calc_star_sky_flux_from_preimg()
        if self.use_pre_img == 'n': self.calc_star_sky_flux_from_magnitude()

        

    def calc_star_sky_flux_from_preimg(self):
        #Calculates the star flux based on the pre-image data.
        #For more information, access: https://github.com/DBernardes/FWHM
        FWHM_obj = FWHM.fwhm(img_name = self.img_dir + '\\' +  self.img_name,                     
                             xy_star = self.obj_coords,
                             sky_radius = self.sky_radius,
                             bias_name = self.img_dir + '\\' +  self.bias_name)                             
        FWHM_obj.read_star_img()
        FWHM_obj.get_max_count()
        FWHM_obj.set_centroid()
        fwhm, star_radius, x, y = FWHM_obj.calc_FWHM()
        FWHM_obj.read_bias_img()
        FWHM_obj.calc_dark_current()
        FWHM_obj.read_exp_time()
        FWHM_obj.read_em_gain()
        FWHM_obj.calc_star_sky_flux()
        snr, rn, sky_flux, star_flux, n_pixels, bias_level = FWHM_obj.calc_SNR()
        self.sky_flux = sky_flux
        self.star_flux = star_flux
        self.n_pix_star = n_pixels
        self.obj_coords = [x,y]
        self.star_radius = star_radius
        self.bias_level = bias_level        

        
    def calc_star_sky_flux_from_magnitude(self):
        #Calculates, through the Pogson equation, the star flux based on the magnitude provided to the software. 
        #This function uses the star flux and magnitude values of the HATS24 star as reference values.       
        aux = 10**((self.hats24_magnitude - self.obj_magnitude)/2.5)
        self.star_flux = self.hats24_flux*aux
        


    def optimize(self, opt_param):
        #Choose the optimization method based on the option provided to the software.
        if opt_param == 1: self.Optimize_SNR()
        if opt_param == 2: self.Optimize_FA()
        if opt_param == 3: self.Optimize_SNR_and_AR()

        

    def Optimize_SNR(self):        
        #Starts the object that determines which modes meet the acquisition rate limit.        
        OAR =  oar.OptimizeAcquisitionRate(acquisition_rate = self.acq_rate,
                                           sub_img_modes=self.sub_img_modes,
                                           binn_modes=self.binn_modes)
        #Determines those modes which accomplish the acquisition rate requirement
        OAR.determine_operation_modes()        
        obj_lista_modos = OAR.read_MOB_obj()
        if obj_lista_modos.get_list_of_modes() == []:
            print('\n\tNo mode meets the requirements provided.')
            print('\tStopping execution.')
            exit()
            
        #Creates the class object that performs the SNR optimization
        OSNR = osnr.OptSignalNoiseRation(snr_target = self.snr,
                                         serial_number = self.serial_number,
                                         ccd_temp = self.ccd_temp,
                                         n_pix_star = self.n_pix_star,
                                         sky_flux = self.sky_flux,
                                         star_flux = self.star_flux,
                                         bias_level = self.bias_level)
        #Write in the class the list of the modes selected in the previous step
        OSNR.write_MOB_obj(obj_lista_modos)
        #Duplicates the list of modes for the PREAMP options 1 and 2.        
        OSNR.duplicate_list_of_modes_for_PA12()
        #Removes repeat modes
        OSNR.remove_repeat_modes()       
        #Select the mode with the largest SNR
        best_snr = OSNR.calc_best_mode()
        #Seeks for the largest allowd sub-image option    
        OSNR.find_sub_img_best_mode()        
        #Prints the best mode
        OSNR.print_best_mode()
        #Exports the optimum mode to an external .txt file
        if 's' in self.export_arq:
           OSNR.export_optimal_setup(self.img_dir, self.file_base_name, self.star_radius, self.obj_coords) 

  

    def Optimize_FA(self):
        repeat = True
        fa_target = self.acq_rate
        while repeat == True:
            #Starts the object for the determination of the modes that accomplish the minimum acquisition rate            
            OAR =  oar.OptimizeAcquisitionRate(acquisition_rate = self.acq_rate,
                                               sub_img_modes=self.sub_img_modes,
                                               binn_modes=self.binn_modes)            
            #Determines which modes meet the provided bin and sub-images options            
            OAR.determine_operation_modes()
            #Returns the object with the selected operation modes            
            obj_lista_modos_FA = OAR.read_MOB_obj()                                         
            #-------------------------------------------------------------------------                   
            #Starts the object for the determination of the modes that accomplish the minimum SNR     
            OSNR = osnr.OptSignalNoiseRation(serial_number = self.serial_number,
                                             snr_target = self.snr,
                                             ccd_temp = self.ccd_temp,
                                             n_pix_star = self.n_pix_star,
                                             sky_flux = self.sky_flux,
                                             star_flux = self.star_flux,
                                             bias_level = self.bias_level)
            #Write the list of the modes selected in the previous step            
            OSNR.write_MOB_obj(obj_lista_modos_FA)
            #Determines which modes meet the minimum provided SNR             
            OSNR.determine_operation_modes_minimun_SNR()
            #Returns the list of selected modes            
            obj_list_of_modes = OSNR.read_MOB_obj()
            #-------------------------------------------------------------------------                               
            #If there is no mode that meets the requirements provided,
            #the minimum acquisition rate is multiplied by 0.8, and the previous steps are performed again
            if obj_list_of_modes.get_list_of_modes() == []:                              
                self.acq_rate *= 0.8
                repeat = True                                        
            else: repeat = False          
        #Write the list of selected operation modes
        OAR.write_MOB_obj(obj_list_of_modes)              
        #Determines the operation mode with the best acquisition rate
        best_mode = OAR.determine_fastest_operation_mode()
        #Prints on the screen the best mode
        OAR.print_best_mode()        
        if fa_target > self.acq_rate:
            print('\nUsed FA= ', round(self.acq_rate,2), 'Hz')
            print('It was not possible to reach the desirable FA')
        #Exports the optimum mode to an external .txt file
        if 's' in self.export_arq:
           OAR.export_optimal_setup(self.img_dir, self.file_base_name, self.star_radius, self.obj_coords, self.acq_rate)
       



    def Optimize_SNR_and_AR(self):
        #Starts the object for the determination of the modes that accomplish the minimum acquisition rate                   
        OAR =  oar.OptimizeAcquisitionRate(acquisition_rate = self.acq_rate, sub_img_modes=self.sub_img_modes, binn_modes=self.binn_modes)    
        #Determines which modes meet the provided bin and sub-images options    
        OAR.determine_operation_modes()
        #Returns the object with the selected operation modes   
        obj_lista_modos_FA = OAR.read_MOB_obj()       
        #-------------------------------------------------------------------------   
        repeat = True
        initial_snr = self.snr
        while repeat == True:            
            #Starts the object for the determination of the modes that accomplish the minimum SNR
            OSNR = osnr.OptSignalNoiseRation(serial_number = self.serial_number,
                                             snr_target = self.snr,
                                             ccd_temp = self.ccd_temp,
                                             n_pix_star = self.n_pix_star,
                                             sky_flux = self.sky_flux,
                                             star_flux = self.star_flux,
                                             bias_level = self.bias_level)
            #Write the list of selected operation modes in the previous step
            OSNR.write_MOB_obj(copy(obj_lista_modos_FA))            
            #Selects those modes that accomplish the minimum SNR
            OSNR.determine_operation_modes_minimun_SNR()        
            #Returns the list of selected modes
            obj_list_of_modes = OSNR.read_MOB_obj()                  
            #If there is no mode that meets the requirements provided,
            #the minimum SNR is multiplied by 0.8, and the previous steps are performed again
            if obj_list_of_modes.get_list_of_modes() == []:                              
                self.snr *= 0.8
                repeat = True                                        
            else: repeat = False
        #-------------------------------------------------------------------------   
        #Starts the object the optimiza both the SNR and the acquisition rante        
        OSNRAR = osnrar.Opt_SignalNoiseRatio_AcquisitionRate(snr_target = self.snr,
                                                           acq_rate = self.acq_rate,
                                                           sub_img_modes = self.sub_img_modes,
                                                           binn_modes = self.binn_modes,
                                                           serial_number = self.serial_number,
                                                           ccd_temp = self.ccd_temp,
                                                           n_pix_star = self.n_pix_star,
                                                           sky_flux = self.sky_flux,
                                                           star_flux = self.star_flux,
                                                           bias_level = self.bias_level)
        #Write the selected modes in the previous steps        
        OSNRAR. write_MOB_obj(obj_list_of_modes)        
        #Calculates the mean and standard deviation of the SNR and acquisition rate.
        #However, it can happen that the SNR is not reached if the chosen texp and em_gain are very small.
        #In these cases, the function discards the respective value
        OSNRAR.SNR_FA_ranges(allow_every_mode = 'n')        
        #Creates the space of states in the hyperopt library format
        OSNRAR.create_space(allow_every_mode = 'n')       
        #Runs the optimzation
        OSNRAR.run_bayesian_opt(max_evals=self.max_evals, algorithm = self.algorithm)
        #Prints the best mode on the screen
        OSNRAR.print_best_mode()
        if initial_snr > self.snr:
            print('\nUsed SNR= ', round(self.snr,2))
            print('It was not possible to reach the desirable SNR')

        #Exports the optimzation iterations to a .txt file
        if 's' in self.export_loss:            
            OSNRAR.creat_log_parameters(self.img_dir, self.file_base_name)
        #Exports the optimal mode to a .txt file
        if 's' in self.export_arq:
            OSNRAR.export_optimal_setup(self.img_dir, self.file_base_name, self.star_radius, self.obj_coords, self.snr)
        



    
