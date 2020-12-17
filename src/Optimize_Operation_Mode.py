#!/usr/bin/env python
# coding: utf-8

# Esta biblioteca gerencia o processo de otimização dos modos de
# operação dos CCDs. 

#10/01/2020. Denis Varise Bernardes.

import Optimize_Signal_Noise_Ratio_Bib as osnr
import Optimize_Acquisition_Rate_Bib as oar
import Optimize_Signal_Noise_Ratio_and_Acquisition_Rate_Bib as osnrar
import SNR_Calculation_Bib as snrc
import Acquisition_Rate_Calculation_Bib as arc
import Modos_Operacao_Bib as mob
import Photon_Flux_Calc_Bib as pfc
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
        snr, acq_rate, obj_magnitude, sub_img_modes, binn_modes, serial_number, ccd_temp, max_evals, file_base_name, export_arq, export_loss, export_bias, use_pre_img, pre_img_name, obj_coords, bias_img_name, sky_radius = get_obs_setup(img_dir)
        #print(snr, acq_rate, obj_magnitude, sub_img_modes, binn_modes, serial_number, ccd_temp, n_pix_star, max_evals)
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
        self.img_directory = img_dir+ '\\'
        self.pre_img_name = pre_img_name
        #self.pre_img_name = image_name
        self.obj_coords = obj_coords
        self.bias_img_name = bias_img_name 
        self.sky_radius = sky_radius
        self.star_radius = 0
        self.bias_level = 500

        self.MOB = []

        #os.chdir(img_directory) #muda para o diretório das imagens
        

    def verify_provides_modes(self):
        for sub_img in self.sub_img_modes:
            if sub_img not in [1024,512,256]:
                print('\nModo sub-image inválido! [%i]'%sub_img)
                #print(sub_img)
                exit()
        for binn in self.binn_modes:
            if binn not in [1,2]:
                print('\nModo binning inválido! [%i]'%binn)
                #print(binn)
                exit()

    def calc_star_flux(self):
        if self.use_pre_img == 's': self.calc_star_sky_flux_from_preimg()
        if self.use_pre_img == 'n': self.calc_star_sky_flux_from_magnitude()        

    def calc_star_sky_flux_from_preimg(self):
        #Inicia o objeto para a redução da pré-imagem
        pre_img_name = self.img_directory + self.pre_img_name
        bias_img_name = self.img_directory + self.bias_img_name        
        PFC = pfc.PhotonFluxCalc(img_name = pre_img_name, bias_name = bias_img_name, xy_star = self.obj_coords , sky_radius = self.sky_radius, ccd_serial=self.serial_number)
        #Lê a imagem de bias fornecida
        PFC.read_bias_img()
        #Lê as informações do header da imagem
        PFC.read_img_get_info()
        #Procura pelo pixel com valor máximo
        PFC.get_max_count()
        #Configura o centroide da estrela com base no valor máximo de contagem encontrado
        PFC.set_centroid()
        #Calcula a FWHM
        PFC.calc_FWHM()
        #Calcula o fluxo da estrela e do céu com base no raio encontrado
        PFC.calc_star_sky_flux()
        #Calcula a SNR 
        PFC.calc_SNR()
        snr, fwhm, dc, rn, bias, sky_flux, star_flux, n_pixels, star_radius, x, y, t_exp, binn, em_mode = PFC.get_results()        
        self.sky_flux = sky_flux
        self.star_flux = star_flux
        self.n_pix_star = n_pixels
        self.obj_coords = [x,y]
        self.star_radius = star_radius
        self.bias_level = bias
##        print(snr)
##        print(self.star_flux)
##        print(self.sky_flux)
##        print(n_pixels)        
##        exit()


    def calc_star_sky_flux_from_magnitude(self):
        #Com base na magnitude fornecido pelo usuário, é calculado o fluxo de fótons do objeto
        # em relação do fluxo do hats24.
        aux = 10**((self.hats24_magnitude - self.obj_magnitude)/2.5)
        self.star_flux = self.hats24_flux*aux
        #print(self.star_flux),exit()


    def optimize(self, fixar_param):
        if fixar_param == 1: self.Optimize_SNR_provided_FA()
        if fixar_param == 2: self.Optimize_FA_provided_SNR()
        if fixar_param == 3: self.Optimize_SNR_and_AR()

    def Optimize_SNR_provided_FA(self):
        repeat = True
##        fa_target = self.acq_rate
##        while repeat == True:          
        # Inicializa o objeto para a determinacao da frequencia de aquisicao
        OAR =  oar.OptimizeAcquisitionRate(acquisition_rate = self.acq_rate,
                                           sub_img_modes=self.sub_img_modes,
                                           binn_modes=self.binn_modes)       
        #Determina quais modos se encaixam nos parametros passados
        OAR.determine_operation_modes()
        #Retorna um objeto contendo os modos de operacao selecionados
        obj_lista_modos = OAR.read_MOB_obj()
        #Cria a lista de modos à partir do objeto MOB        
        #OAR.print_MOB_list(), exit()

##            if obj_list_of_modes.get_list_of_modes() == []:                              
##                    self.acq_rate *= 0.8
##                    repeat = True                                        
##                else: repeat = False

        if obj_lista_modos.get_list_of_modes() == []:
            print('\n\tNenhum modo atende aos requisitos fornecidos.')
            print('\tStopping execution.')
            exit()        
        #cria o objeto da classe que executa o metodo de otimizacao
        OSNR = osnr.OptSignalNoiseRation(snr_target = self.snr,
                                         serial_number = self.serial_number,
                                         ccd_temp = self.ccd_temp,
                                         n_pix_star = self.n_pix_star,
                                         sky_flux = self.sky_flux,
                                         star_flux = self.star_flux,
                                         bias_level = self.bias_level)

        #escreve na classe a lista dos modos que atendem ao requisito da Freq.
        OSNR.write_MOB_obj(obj_lista_modos)       
        #imprime na tela os modos de operação selecionados pela biblioteca acquisition_frequency
        #OSNR.print_MOB_list(),exit()
        #duplica a lista dos modos selecionados para preamp 1 e 2
        OSNR.duplicate_list_of_modes_for_PA12()
        OSNR.remove_repeat_modes()
        #OSNR.create_space()
        #OSNR.run_bayesian_opt(max_evals=self.max_evals, algorithm = self.algorithm)
        
        #Realiza um loop calculando a melhor SNR para cada modo
        # O modo com a maior SNR será selecionado
        best_snr = OSNR.calc_best_mode()
        #encontra o maior sub_img para o modo selecionado        
        OSNR.find_sub_img_best_mode()        
        #printa na tela o melhor modo
        OSNR.print_best_mode()
        
        if 's' in self.export_arq:
           OSNR.export_optimal_setup(self.img_directory, self.file_base_name, self.star_radius, self.obj_coords)

        
##        #cria o dominio da funcao de interesse
##        OSNR.create_space(allow_every_mode = 'n')        
##        #imprime a lista filtrada dos modos de operação. Não deve aparecer modo redundante para o caso do SNR
##        OSNR.print_filtered_list_of_modes(),exit()
##        #executa o metodo de otimizacao bayesiano        
##        OSNR.run_bayesian_opt(max_evals=self.max_evals, algorithm = self.algorithm)        
##        if 's' in self.export_loss:            
##            OSNR.creat_log_parameters(self.img_directory, self.file_base_name)        
##        if 's' in self.export_bias:
##           OSNR.create_bias_list(self.img_directory, self.file_base_name)



    def Optimize_FA_provided_SNR(self):
        repeat = True
        fa_target = self.acq_rate
        while repeat == True:  
            # Inicializa o objeto para a determinacao da frequencia de aquisicao ótima
            OAR =  oar.OptimizeAcquisitionRate(acquisition_rate = self.acq_rate,
                                               sub_img_modes=self.sub_img_modes,
                                               binn_modes=self.binn_modes)
            
            #Determina quais modos se encaixam nos parametros passados (bin e sub-img)
            OAR.determine_operation_modes()
            #Retorna um objeto contendo os modos de operacao selecionados
            obj_lista_modos_FA = OAR.read_MOB_obj()
            #Cria a lista de modos à partir do objeto MOB        
            #OAR.print_MOB_list(),exit()

                   
            #----------------------------------------------------------------------------
            #Faço essa inverção porque preciso o texp máximo de cada modo que atinge a FA
            #mínima. Com base nisso, calculo o valor do ganho EM máximo correspondente.                     
            #cria o objeto da classe que encontra os modos de SNR permitidos
            OSNR = osnr.OptSignalNoiseRation(serial_number = self.serial_number,
                                             snr_target = self.snr,
                                             ccd_temp = self.ccd_temp,
                                             n_pix_star = self.n_pix_star,
                                             sky_flux = self.sky_flux,
                                             star_flux = self.star_flux,
                                             bias_level = self.bias_level)
            #escreve na classe a lista dos modos que atendem ao requisito da Freq.
            OSNR.write_MOB_obj(obj_lista_modos_FA)
            #OAR.print_MOB_list(),exit()
            # determina os modos de operação que atendem ao SNR mínimo fornecido
            OSNR.determine_operation_modes_minimun_SNR()                                                   
            # Lê o objeto Modos de Operação contendo a lista dos modos selecionados
            obj_list_of_modes = OSNR.read_MOB_obj()      
            # Imprime a lista de objetos na tela
            #OSNR.print_MOB_list()

            if obj_list_of_modes.get_list_of_modes() == []:                              
                self.acq_rate *= 0.8
                repeat = True                                        
            else: repeat = False
        #----------------------------------------------------------------------------
        

        list_of_modes = [        
        {'em_mode': 1, 'em_gain': 2, 'hss': 10, 'preamp': 1, 'binn': 1, 'sub_img': 512, 'min_t_exp': 0.03, 'max_t_exp': 15},
        {'em_mode': 1, 'em_gain': 2, 'hss': 10, 'preamp': 1, 'binn': 1, 'sub_img': 256, 'min_t_exp': 0.03, 'max_t_exp': 15}        
        ]    
        #obj_list_of_modes.write_list_of_modes(list_of_modes)

        
        # Escreve dentro da classe o objeto contendo os modos de operação selecionados pela OSNR
        OAR.write_MOB_obj(obj_list_of_modes)       
        #Imprime na tela os modos de operação fornecidos
        #OAR.print_MOB_list(),exit()
        #Determina o modo de operação mais rápido contido na lista do objeto MOB
        best_mode = OAR.determine_fastest_operation_mode()
        #Imprime na tela o melhor modo
        OAR.print_best_mode()        
        if fa_target > self.acq_rate:
            print('\nUsed FA= ', round(self.acq_rate,2), 'Hz')
            print('It was not possible to reach the desirable FA')

        if 's' in self.export_arq:
           OAR.export_optimal_setup(self.img_directory, self.file_base_name, self.star_radius, self.obj_coords, self.acq_rate)
       

    def Optimize_SNR_and_AR(self):
        # Inicializa o objeto para a determinacao da frequencia de aquisicao
        OAR =  oar.OptimizeAcquisitionRate(acquisition_rate = self.acq_rate, sub_img_modes=self.sub_img_modes, binn_modes=self.binn_modes)    
        #Determina quais modos se encaixam nos parametros passados
        OAR.determine_operation_modes()
        #Retorna um objeto contendo os modos de operacao selecionados
        obj_lista_modos_FA = OAR.read_MOB_obj()
        #Cria a lista de modos à partir do objeto MOB        
        #OAR.print_MOB_list(),exit()
           
        repeat = True
        initial_snr = self.snr
        while repeat == True:            
            #cria o objeto da classe que encontra os modos de SNR permitidos
            OSNR = osnr.OptSignalNoiseRation(serial_number = self.serial_number,
                                             snr_target = self.snr,
                                             ccd_temp = self.ccd_temp,
                                             n_pix_star = self.n_pix_star,
                                             sky_flux = self.sky_flux,
                                             star_flux = self.star_flux,
                                             bias_level = self.bias_level)
            #escreve na classe a lista dos modos que atendem ao requisito da Freq.
            OSNR.write_MOB_obj(copy(obj_lista_modos_FA))
            #OSNR.print_MOB_list()
            # determina os modos de operação que atendem ao SNR mínimo fornecido
            OSNR.determine_operation_modes_minimun_SNR()        
            # Lê o objeto Modos de Operação contendo a lista dos modos selecionados
            obj_list_of_modes = OSNR.read_MOB_obj()       
            # Imprime a lista de objetos na tela
            #OSNR.print_MOB_list()
            
            if obj_list_of_modes.get_list_of_modes() == []:                              
                self.snr *= 0.8
                repeat = True                                        
            else: repeat = False
        
        #Cria o objeto que otimiza a SNR e a FA ao mesmo tempo
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
        

        
        list_of_modes = [        
        {'em_mode': 1, 'em_gain': 2, 'hss':  1, 'preamp': 1, 'binn': 1, 'sub_img': 1024, 'min_t_exp': 10, 'max_t_exp': 15}           
        ]    
        #obj_list_of_modes.write_list_of_modes(list_of_modes)

        #Escreve na classe o obj com a lista de modos selecionados pela FA e SNR
        OSNRAR. write_MOB_obj(obj_list_of_modes)
        #OSNRAR.print_MOB_list(),exit()
        # Realiza várias iterações para o cálculo da média e desvio padrão da SNR e FA
        # Possui um problema. A SNR pode não ser atingida caso o texp e o em_gain escolhidos
        # sejam muito pequenos. Quando isso ocorre, a função retorna valor zero (descarta)
        OSNRAR.SNR_FA_ranges(allow_every_mode = 'n')        
        #Cria o espaço de estados no formato do MOB.
        OSNRAR.create_space(allow_every_mode = 'n')       
        # Roda o MOB para a quantidade de iterações fornecida
        OSNRAR.run_bayesian_opt(max_evals=self.max_evals, algorithm = self.algorithm)
        #printa na tela o melhor modo
        OSNRAR.print_best_mode()
        #cria um arquivo contendo o log dos valores do ruido obtidos em cada iteracao
        #OSNRAR.creat_log_loss()

        if initial_snr > self.snr:
            print('\nUsed SNR= ', round(self.snr,2))
            print('It was not possible to reach the desirable SNR')

        if 's' in self.export_loss:            
            OSNRAR.creat_log_parameters(self.img_directory, self.file_base_name)        
        if 's' in self.export_arq:
            OSNRAR.export_optimal_setup(self.img_directory, self.file_base_name, self.star_radius, self.obj_coords, self.snr)
        if 's' in self.export_bias:
            OSNRAR.create_bias_list(self.img_directory, self.file_base_name)

        #Esta função irá avaliar todo o espaço de estados para verificar
        # se realmente o modo ótimo obtido é o melhor de todos.
        #OSNRAR.get_mean_random_SNR_FA(self.img_directory, self.file_base_name)



    def create_space_of_states_with_all_modes(self):
        self.MOB = mob.ModosOperacao()
        vector_em = [0, 1]        
        for em_mode in vector_em:
            vector_hss = [1, 10, 20, 30]
            if em_mode == 0: vector_hss = [0.1, 1]
            
            for hss in vector_hss:
                for preamp in [1,2]:
                    for binn in [1,2]:
                        for sub_img in [1024, 512, 256]:
                            self.MOB.write_mode(em_mode, 300, hss, preamp, binn, sub_img, 1/self.acq_rate, 1e-5)
##        lista = self.MOB.get_list_of_modes()
##        for mode in lista:print(mode)
##        exit()


    def create_space_of_states_provided_FA(self):
        # Inicializa o objeto para a determinacao da frequencia de aquisicao
        OAR =  oar.OptimizeAcquisitionRate(acquisition_rate = self.acq_rate,
                                           sub_img_modes=self.sub_img_modes,
                                           binn_modes=self.binn_modes)       
        #Determina quais modos se encaixam nos parametros passados
        OAR.determine_operation_modes()
        #Retorna um objeto contendo os modos de operacao selecionados
        self.MOB = OAR.read_MOB_obj()
        #Cria a lista de modos à partir do objeto MOB        
        #OAR.print_MOB_list(), exit()


    def get_mean_radom_SNR(self):
        losses_SNR = []
        list_operations_modes = self.MOB.get_list_of_modes()
        n_repeat = ceil(500/len(list_operations_modes))        
        if n_repeat < 1: n_repeat = 1        
        for i in range(n_repeat):
            for mode in list_operations_modes:                
                em_mode = mode['em_mode']
                hss = mode['hss']
                binn = mode['binn']
                sub_img = mode['sub_img']
                preamp = mode['preamp']#rd.choice([1,2])
                max_em_gain = mode['em_gain']
                em_gain = 0
                t_exp = rd.uniform(mode['min_t_exp'], mode['max_t_exp'])
                if em_mode == 1:
                    gain = self.set_gain(em_mode, hss, preamp)
                    max_em_gain, max_t_exp = self.calc_max_em_gain(mode['max_t_exp'], mode['min_t_exp'], gain)
                    if max_em_gain > 300: max_em_gain = 300
                    if max_em_gain == 0 : continue
                    t_exp = rd.uniform(mode['min_t_exp'], max_t_exp)
                    em_gain = rd.uniform(2,max_em_gain)
##                    print(em_mode, hss, binn, sub_img, preamp, max_t_exp, max_em_gain)
##                print(em_mode, hss, binn, sub_img, preamp, t_exp, max_em_gain)
                SNRC = snrc.SignalToNoiseRatioCalc(t_exp, em_mode, em_gain, hss, preamp, binn, self.ccd_temp, self.sky_flux, self.star_flux, self.n_pix_star, self.serial_number)
                SNRC.set_gain_value()
                SNRC.calc_RN()
                SNRC.calc_DC()
                SNRC.calc_SNR()
                losses_SNR.append(SNRC.get_SNR())                                    
               
        mean_snr = np.mean(losses_SNR)
        print('\n mean: ', mean_snr, '\n')



    def get_mean_radom_FA(self):
        losses_FA = []
        list_operations_modes = self.MOB.get_list_of_modes()
        n_repeat = ceil(500/len(list_operations_modes))        
        if n_repeat < 1: n_repeat = 1       
        for i in range(n_repeat):
            for mode in list_operations_modes:                
                em_mode = mode['em_mode']
                hss = mode['hss']
                binn = mode['binn']
                sub_img = mode['sub_img']
                preamp = rd.choice([1,2])
                t_exp = rd.uniform(mode['min_t_exp'], mode['max_t_exp'])
                max_t_exp = mode['max_t_exp']
                max_em_gain = 0
                if em_mode == 1:
                    gain = self.set_gain(em_mode, hss, preamp)
                    max_em_gain, max_t_exp = self.calc_max_em_gain(mode['max_t_exp'], mode['min_t_exp'], gain)
                    if max_em_gain > 300: max_em_gain = 300
                    if max_em_gain == 0 : continue
                    em_gain = rd.uniform(2,max_em_gain)                
                SNRC = snrc.SignalToNoiseRatioCalc(max_t_exp, em_mode, max_em_gain, hss, preamp, binn, self.ccd_temp, self.sky_flux, self.star_flux, self.n_pix_star, self.serial_number)
                SNRC.set_gain_value()
                SNRC.calc_RN()
                SNRC.calc_DC()                
                min_t_exp = SNRC.calc_minimun_texp_provided_SNR(self.snr)
                if max_t_exp < min_t_exp:continue
                t_exp = rd.uniform(min_t_exp, max_t_exp)                
                ARC = arc.AcquisitionRateCalc()
                ARC.write_operation_mode(em_mode, hss, binn, sub_img, t_exp)
                ARC.seleciona_t_corte()
                ARC.calc_acquisition_rate()
                losses_FA.append(float(ARC.return_acquisition_rate()))                                   
               
        mean_fa = np.mean(losses_FA)
        print('\n mean: ', mean_fa, '\n')




    def get_mean_random_SNR_FA(self):
        # Inicializa o objeto para a determinacao da frequencia de aquisicao
        OAR =  oar.OptimizeAcquisitionRate(acquisition_rate = self.acq_rate, sub_img_modes=self.sub_img_modes, binn_modes=self.binn_modes)    
        #Determina quais modos se encaixam nos parametros passados
        OAR.determine_operation_modes()
        #Retorna um objeto contendo os modos de operacao selecionados
        obj_lista_modos_FA = OAR.read_MOB_obj()
        #Cria a lista de modos à partir do objeto MOB        
        #OAR.print_MOB_list(),exit()      

        #cria o objeto da classe que encontra os modos de SNR permitidos
        OSNR = osnr.OptSignalNoiseRation(serial_number = self.serial_number,
                                         snr_target = self.snr,
                                         ccd_temp = self.ccd_temp,
                                         n_pix_star = self.n_pix_star,
                                         sky_flux = self.sky_flux,
                                         star_flux = self.star_flux,
                                         bias_level = self.bias_level)
        #escreve na classe a lista dos modos que atendem ao requisito da Freq.
        OSNR.write_MOB_obj(obj_lista_modos_FA)        
        #OSNR.print_MOB_list(),exit()
        # determina os modos de operação que atendem ao SNR mínimo fornecido
        OSNR.determine_operation_modes_minimun_SNR()        
        # Lê o objeto Modos de Operação contendo a lista dos modos selecionados
        obj_list_of_modes = OSNR.read_MOB_obj()       
        # Imprime a lista de objetos na tela
        #OSNR.print_MOB_list(),exit()
        
        if obj_list_of_modes.get_list_of_modes() == []:
            print('\n\tNenhum modo atende aos requisitos fornecidos.')
            print('\tStopping execution.')
            exit()


        #------------------------------------------------------------------------------------------------------------------------

        losses_FA = []
        losses_SNR = []
        list_operations_modes = obj_list_of_modes.get_list_of_modes()
        n_repeat = ceil(50/len(list_operations_modes))        
        if n_repeat < 1: n_repeat = 1       
        for i in range(n_repeat):
            for mode in list_operations_modes:                
                em_mode = mode['em_mode']
                hss = mode['hss']
                binn = mode['binn']
                sub_img = mode['sub_img']
                preamp = rd.choice([1,2])
                t_exp = rd.uniform(mode['min_t_exp'], mode['max_t_exp'])                
                em_gain = 0
                if em_mode == 1:                    
                    em_gain = rd.uniform(2,mode['em_gain'])
                    
                ARC = arc.AcquisitionRateCalc()
                ARC.write_operation_mode(em_mode, hss, binn, sub_img, t_exp)
                ARC.seleciona_t_corte()
                ARC.calc_acquisition_rate()
                losses_FA.append(float(ARC.return_acquisition_rate()))
                
                SNRC = snrc.SignalToNoiseRatioCalc(t_exp, em_mode, em_gain, hss, preamp, binn, self.ccd_temp, self.sky_flux, self.star_flux, self.n_pix_star, self.serial_number)
                SNRC.set_gain_value()
                SNRC.calc_RN()
                SNRC.calc_DC()
                SNRC.calc_SNR()
                losses_SNR.append(SNRC.get_SNR())   
               
        mean_fa = np.mean(losses_FA)
        mean_snr = np.mean(losses_SNR)
        #print(mean_fa, mean_snr),exit()

        losses_SNR_FA = []
        list_operations_modes = obj_list_of_modes.get_list_of_modes()
        n_repeat = ceil(500/len(list_operations_modes))        
        if n_repeat < 1: n_repeat = 1       
        #for i in range(n_repeat):
        for mode in list_operations_modes:
            for t_exp in np.linspace(mode['min_t_exp'], mode['max_t_exp'], 10):
                em_mode = mode['em_mode']
                hss = mode['hss']
                binn = mode['binn']
                sub_img = mode['sub_img']
                preamp = rd.choice([1,2])
                #t_exp = rd.uniform(mode['min_t_exp'], mode['max_t_exp'])                
                em_gain = 0
                if em_mode == 1:                    
                    em_gain = rd.uniform(2,mode['em_gain'])
                    
                ARC = arc.AcquisitionRateCalc()
                ARC.write_operation_mode(em_mode, hss, binn, sub_img, t_exp)
                ARC.seleciona_t_corte()
                ARC.calc_acquisition_rate()
                FA = float(ARC.return_acquisition_rate())
                
                #SNRC = snrc.SignalToNoiseRatioCalc(t_exp, em_mode, em_gain, hss, preamp, binn, self.ccd_temp, self.sky_flux, self.star_flux, self.n_pix_star, self.serial_number)
                SNRC = snrc.SignalToNoiseRatioCalc(t_exp, em_mode, 300, hss, preamp, binn, self.ccd_temp, self.sky_flux, self.star_flux, self.n_pix_star, self.serial_number)
                SNRC.set_gain_value()
                SNRC.calc_RN()
                SNRC.calc_DC()
                SNRC.calc_SNR()
                SNR = SNRC.get_SNR()

                loss_snr_fa = FA/mean_fa * SNR/mean_snr
                print(mode, ';', t_exp, ';', loss_snr_fa)
                #print(em_mode, em_gain, hss, preamp, binn, sub_img, t_exp)

            #losses_SNR_FA.append(loss_snr_fa)

##        #print('\n', losses_SNR_FA)
##        mean_fa_snr = np.mean(losses_SNR_FA)
##        print('\nmean = ', mean_fa_snr)
        

               
                 


    def calc_max_em_gain(self, max_t_exp, min_t_exp, gain):
        #Cálculo do ganho EM máximo permitido para cada modo. O cálculo é baseado na quantidade máxima de 100 fótons
        # por pixel do CCD para a qual o modo EM é melhor que o Convencional.
        # Esta função recebe o t_exp máximo de cada modo. Contudo, dentro de um mesmo modo, o ganho EM poderia ser maior
        # quando a iteração escolher um t_exp menor. Talvez, esta seja uma limitação da biblioteca. Ainda não achei solução        
        max_fotons = 100
        dc = 0.0007
        bias = self.bias_level
        max_ADU = (2**16)*0.8        
        aux = (self.sky_flux + 1000000 + dc)*max_t_exp        
        max_em_gain = (max_ADU - bias)/(aux/gain)        
        if aux > max_fotons:
            max_em_gain = (max_ADU - bias)/(max_fotons/gain)
            max_t_exp = max_fotons/(self.sky_flux + self.star_flux/self.n_pix_star + dc)        
        if max_t_exp < min_t_exp:
            max_em_gain = 0
            t_exp = 0    
        return max_em_gain, max_t_exp


    def set_gain(self, em_mode, hss, preamp):
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
        return gain  
