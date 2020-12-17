#!/usr/bin/env python
# coding: utf-8

# Esta biblioteca calcula a curva da SNR para os modos
# ótimo e do usuário.

#25/10/2019. Denis Varise Bernardes.


import os
from sys import exit
import SNR_Calculation_Bib as snrc
import Photon_Flux_Calc_Bib as pfc
import json
import numpy as np
import matplotlib.pyplot as plt
#Locale settings
import locale

class SNR_profile:

    def __init__(self, dir_path, plot_ndots = 10):
        self.dir_path = dir_path
        self.opt_SNR_file = ''
        self.opt_FA_file = ''
        self.star_flux = 0
        self.sky_flux = 0
        self.n_pixels = 0
        self.snr = 0
        self.dc = 0
        self.rn = 0
        self.opt_SNR_profile = []
        self.opt_FA_profile = []
        self.user_SNR_profile = []
        self.plot_ndots = plot_ndots
        

    def get_opt_files_names(self):
        files = os.listdir(self.dir_path)        
        for file in files:
            if '_SNR_OPTSETUP' in file:
                self.opt_SNR_file = file
            if '_FA_OPTSETUP' in  file:
                self.opt_FA_file = file
            if '_BOTH_OPTSETUP' in file:
                self.opt_both_file = file

    def read_operation_mode_files(self):
        arq = open(self.dir_path + '\\' + self.opt_SNR_file, 'r')
        self.opt_SNR_mode = json.load(arq)        
        arq.close()

        arq = open(self.dir_path + '\\' + self.opt_FA_file, 'r')
        self.opt_FA_mode = json.load(arq)      
        arq.close()

        arq = open(self.dir_path + '\\' + self.opt_both_file, 'r')
        self.opt_BOTH_mode = json.load(arq)
        arq.close()

    
        

    def calc_Photon_Flux(self):
        #Pega as informações da observação para o cálculo do fluxo da estrela        
        setup_file = 'observation_setup.txt'
        arq = open(self.dir_path + '\\' + setup_file, 'r')
        lines = arq.read().splitlines()
        arq.close()
        pre_img_name = ''
        bias_img_name = ''
        xy_star = ()
        star_radius = 0        
        for line in lines:                        
            if 'Nome da imagem do objeto' in line: pre_img_name = line.split('=')[1].replace(' ','')  
            if '.gz' not in pre_img_name: pre_img_name+='.gz'
            if 'Nome da imagem de bias' in line: bias_img_name = line.split('=')[1].replace(' ','')  
            if '.gz' not in bias_img_name: bias_img_name +='.gz'
            if 'Coordenadas' in line:
                line = line.split('=')[1].split(',')
                xy_star = (int(line[0]), int(line[1]))
            if 'Raio maximo do ceu' in line: star_radius = int(line.split('=')[1])
            if 'Relacao sinal-ruido' in line: user_snr = float(line.split('=')[1])                     
        
        bias_img_name = self.dir_path + '\\' + bias_img_name
        pre_img_name = self.dir_path + '\\' + pre_img_name
        PFC = pfc.PhotonFluxCalc(img_name = pre_img_name, bias_name = bias_img_name, xy_star = xy_star , sky_radius = star_radius, ccd_serial=9917)        
        PFC.read_bias_img()        
        PFC.read_img_get_info()        
        PFC.get_max_count()        
        PFC.set_centroid()        
        PFC.calc_FWHM()        
        PFC.calc_star_sky_flux()        
        PFC.calc_SNR()
        snr, fwhm, dc, rn, bias, sky_flux, star_flux, n_pixels, star_radius, x, y, t_exp, binn, em_mode = PFC.get_results()
        self.user_mode = {
            'star_flux': star_flux,
            'sky_flux': sky_flux,
            'n_pixels': n_pixels,
            'snr': snr,
            'dc': dc,
            'rn': rn,
            't_exp': t_exp,
            'binn': binn,
            'em_mode':em_mode}
        


    def get_opt_mode_SNR_profile(self):
        opt_SNR_mode = self.opt_SNR_mode
        user_mode = self.user_mode
        max_em_gain  = opt_SNR_mode['em_gain']
        max_t_exp    = opt_SNR_mode['t_exp']
        #min_t_exp   = opt_FA_mode['t_exp']
        em_mode      = opt_SNR_mode['em_mode']
        hss          = opt_SNR_mode['hss']
        binn         = opt_SNR_mode['bin']            
        preamp       = opt_SNR_mode['preamp']
        sub_img      = opt_SNR_mode['sub_img']
        star         = user_mode['star_flux']
        sky          = user_mode['sky_flux']
        n_pix        = user_mode['n_pixels']
    
        opt_SNR_profile = []
        #t_exp_range = np.linspace(min_t_exp*0.95, max_t_exp, self.plot_ndots)
        t_exp_range = np.linspace(1e-5, max_t_exp, self.plot_ndots)
        for t_exp in t_exp_range:
            if em_mode ==  'CONV':em_mode=0
            if em_mode == 'EM': em_mode=1            
            SNRC = snrc.SignalToNoiseRatioCalc(t_exp = t_exp, em_mode = em_mode, em_gain = max_em_gain, hss = hss, preamp = preamp, binn = binn, ccd_temp = -60, sky_flux = sky , star_flux = star, n_pix_star = n_pix, serial_number = 9917)
            SNRC.set_gain_value()
            SNRC.calc_RN()
            SNRC.calc_DC()
            SNRC.calc_SNR()
            snr = SNRC.get_SNR()
            #print(snr)
            opt_SNR_profile.append(snr)
        self.opt_SNR_profile = opt_SNR_profile


    def get_opt_mode_FA_profile(self):
        opt_FA_mode = self.opt_FA_mode
        user_mode   = self.user_mode
        max_em_gain = opt_FA_mode['em_gain']
        max_t_exp   = self.opt_SNR_mode['t_exp']
        min_t_exp   = opt_FA_mode['t_exp']
        em_mode     = opt_FA_mode['em_mode']
        hss         = opt_FA_mode['hss']
        binn        = opt_FA_mode['bin']            
        preamp      = opt_FA_mode['preamp']
        sub_img     = opt_FA_mode['sub_img']
        star         = user_mode['star_flux']
        sky          = user_mode['sky_flux']
        n_pix        = user_mode['n_pixels']
        
        opt_FA_profile = []
        #t_exp_range = np.linspace(min_t_exp*0.95, max_t_exp, self.plot_ndots)
        t_exp_range = np.linspace(1e-5, max_t_exp, self.plot_ndots)
        for t_exp in t_exp_range:
            if em_mode ==  'CONV':em_mode=0
            if em_mode == 'EM': em_mode=1
            SNRC = snrc.SignalToNoiseRatioCalc(t_exp = t_exp, em_mode = em_mode, em_gain = max_em_gain, hss = hss, preamp = preamp, binn = binn, ccd_temp = -60, sky_flux = sky , star_flux = star, n_pix_star = n_pix, serial_number = 9917)
            SNRC.set_gain_value()
            SNRC.calc_RN()
            SNRC.calc_DC()
            SNRC.calc_SNR()
            snr = SNRC.get_SNR()
            opt_FA_profile.append(snr)
        self.opt_FA_profile = opt_FA_profile


    def get_opt_mode_BOTH_profile(self):
        opt_BOTH_mode = self.opt_BOTH_mode        
        user_mode = self.user_mode
        max_em_gain      = opt_BOTH_mode['em_gain']
        max_t_exp        = self.opt_SNR_mode['t_exp']
        #min_t_exp       = opt_BOTH_mode['t_exp']
        em_mode          = opt_BOTH_mode['em_mode']
        hss              = opt_BOTH_mode['hss']
        binn             = opt_BOTH_mode['bin']            
        preamp           = opt_BOTH_mode['preamp']
        sub_img          = opt_BOTH_mode['sub_img']
        star         = user_mode['star_flux']
        sky          = user_mode['sky_flux']
        n_pix        = user_mode['n_pixels']  
    
        opt_BOTH_profile = []
        #t_exp_range = np.linspace(min_t_exp*0.95, max_t_exp, self.plot_ndots)
        t_exp_range = np.linspace(1e-5, max_t_exp, self.plot_ndots)
        for t_exp in t_exp_range:
            if em_mode ==  'CONV':em_mode=0
            if em_mode == 'EM': em_mode=1            
            SNRC = snrc.SignalToNoiseRatioCalc(t_exp = t_exp, em_mode = em_mode, em_gain = max_em_gain, hss = hss, preamp = preamp, binn = binn, ccd_temp = -60, sky_flux = sky , star_flux = star, n_pix_star = n_pix, serial_number = 9917)
            SNRC.set_gain_value()
            SNRC.calc_RN()
            SNRC.calc_DC()
            SNRC.calc_SNR()
            snr = SNRC.get_SNR()            
            opt_BOTH_profile.append(snr)
        self.opt_BOTH_profile = opt_BOTH_profile
        

    def get_user_mode_SNR_profile(self):
        user_mode = self.user_mode
        max_t_exp = self.opt_SNR_mode['t_exp']       
        min_t_exp = self.opt_FA_mode['t_exp']
        em_gain   = 1
        n_pix     = user_mode['n_pixels']
        binn      = user_mode['binn']
        dc        = user_mode['dc']
        rn        = user_mode['rn']
        star      = user_mode['star_flux']
        sky       = user_mode['sky_flux']
        nf        = 1    
        if user_mode['em_mode'] == 1: nf=1.4        

        user_SNR_profile = []        
        t_exp_range = np.linspace(1e-5, max_t_exp, self.plot_ndots)
        for t_exp in t_exp_range:
            #print(snr, sky, star, n_pix),exit()            
            aux = np.sqrt(star * t_exp * nf**2 + n_pix * ( (rn/em_gain/binn)**2 + (sky + dc)*t_exp * nf**2))        
            SNR = (star*t_exp) / aux            
            user_SNR_profile.append(SNR)
        self.user_SNR_profile = user_SNR_profile


    def calc_SNR_opt_mode_with_user_t_exp(self):
        opt_SNR_mode = self.opt_SNR_mode
        user_mode = self.user_mode
        em_gain      = opt_SNR_mode['em_gain']       
        em_mode      = opt_SNR_mode['em_mode']
        hss          = opt_SNR_mode['hss']
        binn         = opt_SNR_mode['bin']            
        preamp       = opt_SNR_mode['preamp']
        t_exp        = self.user_mode['t_exp']
        star         = user_mode['star_flux']
        sky          = user_mode['sky_flux']
        n_pix        = user_mode['n_pixels']
        SNRC = snrc.SignalToNoiseRatioCalc(t_exp = t_exp, em_mode = em_mode, em_gain = em_gain, hss = hss, preamp = preamp, binn = binn, ccd_temp = -60, sky_flux = sky , star_flux = star, n_pix_star = n_pix, serial_number = 9917)
        SNRC.set_gain_value()
        SNRC.calc_RN()
        SNRC.calc_DC()
        SNRC.calc_SNR()
        SNR = SNRC.get_SNR()        
        return SNR   


    def get_modes(self):          
        return self.opt_SNR_mode, self.opt_FA_mode, self.opt_BOTH_mode, self.user_mode
    

    def get_profiles(self):
        return self.opt_SNR_profile, self.opt_FA_profile, self.opt_BOTH_profile, self.user_SNR_profile
       

dir_path = r'C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Imagens OPD'
directory = '17nov24'
sub_directory = dir_path + '\\' + directory
##directories = os.listdir(dir_path)
##for directory in directories:
##    if os.path.isdir(dir_path + '\\' + directory):
##        if 'Noites nao usadas' not in directory:            
##            if '17mar' in directory: directory+= '\\arsco'
##            if '17jul28' in directory : continue
##            if 'ZIP' in directory : continue
##            if 'Graphs' in directory: continue
##            sub_directory = dir_path + '\\' + directory            

SNRP = SNR_profile(sub_directory, 10)
SNRP.get_opt_files_names()
SNRP.read_operation_mode_files()
SNRP.calc_Photon_Flux()
SNRP.get_opt_mode_SNR_profile()
SNRP.get_opt_mode_FA_profile()
SNRP.get_opt_mode_BOTH_profile()
SNRP.get_user_mode_SNR_profile()
opt_snr_user_t_exp = SNRP.calc_SNR_opt_mode_with_user_t_exp()
opt_SNR_mode, opt_FA_mode, opt_BOTH_mode, user_mode = SNRP.get_modes()
opt_SNR_profile, opt_FA_profile, opt_BOTH_profile, user_SNR_profile = SNRP.get_profiles()
#opt_snr_user_t_exp*=0.8
fontsize = 14
locale.setlocale(locale.LC_NUMERIC, "de_DE")
plt.rcdefaults()
plt.rcParams['axes.formatter.use_locale'] = True
plt.rcParams.update({'font.size': fontsize})
plt.figure(figsize=(7,5))

t_exp_opt = np.linspace(1e-5, opt_SNR_mode['t_exp'], 10)
t_exp_user = np.linspace(1e-5, user_mode['t_exp'], 10)

#----------------------FA profile------------------------------------
plt.plot(t_exp_opt, opt_FA_profile, '-', c='r')#, label='Modo FA ótimo')
plt.plot(opt_FA_mode['t_exp']*1.035, user_mode['snr'], 'o', c='r')

#----------------------BOTH profile------------------------------------
plt.plot(t_exp_opt, opt_BOTH_profile, '-', c='r')#, label='Modo SNR x FA ótimo')
plt.plot(opt_BOTH_mode['t_exp'], opt_BOTH_mode['SNR'], 'o', c='r')

#----------------------SNR profile--------------------------------------
plt.plot(t_exp_opt, opt_SNR_profile, '-', c='r', label='Modo SNR ótimo')
plt.plot(opt_SNR_mode['t_exp'], opt_SNR_mode['output'], 'o', c='r')

#------------------------User Profile---------------------------------
plt.plot(t_exp_opt, user_SNR_profile, '-', c='b', label='Modo do usuário')
plt.plot(user_mode['t_exp'], user_mode['snr'], 'o', c='b')
plt.plot(user_mode['t_exp'], opt_snr_user_t_exp*0.98, 'o', c='r')


plt.hlines(user_mode['snr'], 0, user_mode['t_exp'], colors='k', linestyles='--', linewidth=0.9)
plt.vlines(opt_SNR_mode['t_exp'], 0, opt_SNR_mode['output'], colors='k', linestyles='--', linewidth=0.9)
plt.vlines(user_mode['t_exp'], 0, opt_snr_user_t_exp, colors='k', linestyles='--', linewidth=0.9)


plt.xlim(opt_FA_mode['t_exp']*0.6, opt_SNR_mode['t_exp']*1.01)
plt.ylim(user_mode['snr']*0.6, opt_SNR_mode['output']*1.1)
plt.xlabel('Tempo de exposição (s)', fontsize = fontsize)
plt.ylabel('SNR', fontsize = fontsize)
plt.legend(loc='upper left', fontsize = fontsize)
if 'arsco' in directory: directory = directory.split('\\')[0]
plt.show()
#plt.savefig(dir_path + '\Graphs' + '\\' + directory + '.png')
plt.close()


##fontsize = 14
##plt.figure(figsize=(7,5))
##plt.annotate('FA ótima', xy=(0.5, 0.5), fontsize = fontsize)
##plt.annotate('SNR ótima', xy=(0.5, 0.4), fontsize = fontsize)
##plt.annotate('SNR x FA ótima', xy=(0.5, 0.3), fontsize = fontsize)
##plt.annotate('Modo do usuário', xy=(0.5, 0.2), fontsize = fontsize)
##plt.savefig(dir_path + '\Graphs' + '\\' + directory + '_labels' + '.png')

            
