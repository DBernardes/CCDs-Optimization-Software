##def function(img_dir, opt_mode):
##    #!/usr/bin/env python
##    # coding: utf-8
##
##    # Este código para a determinacao do
##    # modo de operacao de uma camera iXon Ultra 888 que otimize o ruido
##    # de leitura, a taxa de aquisicao, ou ambos.
##
##    #25/10/2019. Denis Varise Bernardes.
##
##    #img_dir =r'C:\Users\observer\Desktop\Imagens_CCD'
##
##    
##    from sys import exit
##    import Optimize_Operation_Mode as OOM
##    #Cria a classe que otimiza o modo de operação
##    OOM = OOM.Optimize_Operation_Mode(img_dir)
##    #verifica se os modos fornecidos de sub_img e binning estao corretos
##    OOM.verify_provides_modes()
##    #Esta opção serve para o cálculo do fluxo da estrela
##    OOM.calc_star_flux()
##    #Fixar(1-SNR, 2-Freq. Acq., 3-Ambos)
##    OOM.optimize(opt_mode)


### make a table of Gaussian sources
from astropy.table import Table
from astropy.io.fits import writeto, getdata, Header
##from photutils.datasets import make_gaussian_prf_sources_image
##from photutils.datasets import make_noise_image
##import matplotlib.pyplot as plt
##from sys import exit
##
##table = Table()
##table['amplitude'] = [100]
##table['x_0'] = [100]
##table['y_0'] = [100]
##table['sigma'] = [10]
##
##shape = (200, 200)
##image = make_gaussian_prf_sources_image(shape, table)
##image += make_noise_image(shape, distribution='gaussian', mean=1167., stddev=3.37)
##hdr = Header()
##                       
##hdr['NAXIS1']  =(1024, 'length of data axis 1')
##hdr['NAXIS2']  =(1024, 'length of data axis 2')                          
##hdr['EXTEND']  = ('T', 'FITS dataset may contain extensions')                                             
##hdr['COMMENT'] = 'and Astrophysics, volume 376, page 359; bibcode: 2001A&A...376..3'
##hdr['ACQMODE'] = ('Single  ', 'Acquisition Mode')
##hdr['READMODE']= ('Image   ', 'Readout Mode')
##hdr['IMGRECT'] = ('1, 1024,1024, 1', 'Image Format')
##hdr['HBIN']    = ('1       ', 'Horizontal Binning')
##hdr['VBIN']    = ('1       ','Vertical Binning')
##hdr['TRIGGER'] = ('Internal','Trigger Mode')
##hdr['EXPOSURE']= ('20,00000', 'Total Exposure Time')
##hdr['TEMP']    = ('-80,1   ', 'Temperature')
##hdr['READTIME']= (1.0E-006 ,'Pixel readout time ')
##hdr['VSHIFT']  = ('4.33E-06', 'Vertical Shift Speed')
##hdr['OUTPTAMP']= ('Conventional', 'Output Amplifier')
##hdr['PREAMP']  = ('1x', 'Pre Amplifier Gain')
##hdr['SERNO']   = ('9917', 'Serial Number')
##hdr['DATE']   = ('2017-07-14T00:00:58', 'File Creation Date (YYYY-MM-HHThh:mm:ss)')
##hdr['FRAME']   = ('2017-07-14T00:00:58.642', 'Start of Frame Exposure')
##hdr['IMAGE']  = ('hats-24_I_transito_001', 'Nome do arquivo')
##writeto('image.fits', image, overwrite=True, header=hdr)

#------------------------------------------------------------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
#Esta funcao calcula o SNR da estrela



##star = np.linspace(1,1000,1000)
##def calc_SNR_before(rn, em_gain, NF=1.41):
##    t_exp = 1
##    n_pix = 305
##    dc = 0.1    
##    sky = 12.298897076737294 #e-/pix/s
##    star = 100
##    
##    aux = np.sqrt(star*t_exp*NF**2*em_gain + n_pix * (rn**2 + (sky + dc)*t_exp*NF**2*em_gain))        
##    SNR = (star*t_exp*em_gain) / aux    
##    return SNR
##
##def calc_SNR_after(rn, em_gain, star, NF=1.41):
##    t_exp = 1
##    n_pix = 305
##    dc = 0.1    
##    sky = 12.298897076737294 #e-/pix/s       
##    
##    aux = np.sqrt(star*t_exp*NF**2 + n_pix * (rn**2/em_gain**2 + (sky + dc)*t_exp*NF**2 ))        
##    SNR = (star*t_exp) / aux    
##    return SNR


##rn = 209
##plt.subplot(121)
##SNR = calc_SNR_before(rn, em_gain, NF=1)
##plt.plot(em_gain,SNR, '-', c = 'b', label='before')
##plt.xlabel('EM gain')
##plt.ylabel('SNR')
##plt.title('before')
##plt.subplot(122)
##SNR = calc_SNR_after(rn, em_gain)
##plt.plot(em_gain,SNR, '-', c = 'r', label='after')
##plt.title('after')
##plt.xlabel('EM gain')
##plt.ylabel('SNR')
##plt.show()

##rn = 12.20
##em_gain = 300
##SNR = calc_SNR_after(rn, em_gain, star)
##plt.plot(star,SNR, '-', c = 'b', label='EMCCD')
##
##rn = 3.46
##em_gain = 1
##SNR = calc_SNR_after(rn, em_gain, star, NF=1)
##plt.plot(star,SNR, '-', c = 'r', label='CCD ideal')
##
##plt.legend()
##plt.xlabel('photons number per pix')
##plt.ylabel('SNR')
##plt.show()


##def calc_SNR(em_gain, rn, star, nf = 1.41):
##    #Esta funcao calcula o SNR da estrela
##    t_exp = 1  
##    
##    n_pix = 305
##    dc = 0.1     
##    sky = 12    
##    #print(star, t_exp, n_pix,rn,sky ,dc),exit()
##    aux = np.sqrt(star * t_exp * nf**2 + n_pix * (rn**2/em_gain**2 + (sky + dc)*t_exp * nf**2))        
##    SNR = (star*t_exp) / aux
##    return SNR
##
##star = np.linspace(1,1000,100)
##em_gain = 1000
##rn = 12.20
##SNR = calc_SNR(em_gain, rn, star)
##plt.loglog(star, SNR, c='b')
##em_gain = 1
##rn = 3.47
##SNR = calc_SNR(em_gain, rn, star, nf=1)
##plt.loglog(star, SNR, c='r')
##plt.show()


##from collections import OrderedDict
##from astropy.modeling.models import Moffat2D
##from photutils.datasets import make_model_sources_image
##import numpy as np
##
##model = Moffat2D()
##n_sources = 1
##shape = (200, 200)
##table = Table()
##table['amplitude'] = [15000]
##table['x_0'] = [100]
##table['y_0'] = [100]
##table['sigma'] = [1]
##data = make_model_sources_image(shape, model, table)
##writeto('image.fits', data, overwrite=True)

#------------------------------------------------------------------------------------------------------------------------------------
##
##from sys import exit
##from photutils.datasets import make_gaussian_prf_sources_image, make_gaussian_sources_image
##from photutils.datasets import make_noise_image
##from astropy.table import Table
##from astropy.io.fits import writeto
##import numpy as np
##import matplotlib.pyplot as plt
##
##def calc_SNR(t_exp, gain, rn, em_gain, radius, nf, binn):
##    
##    
##    shape = (200,200)
##    y, x = np.indices(shape, dtype='float32')
##    star_amp = (2000)* t_exp * em_gain * binn**2 / gain
##    
##    table = Table()
##    table['amplitude'] = [star_amp]
##    table['x_mean'] = [100]
##    table['y_mean'] = [100]
##    table['x_stddev'] = [3/binn]
##    table['y_stddev'] = [3/binn]
##    table['theta'] = np.radians(np.array([0]))
##    star = make_gaussian_sources_image(shape, table)
##    
##    working_mask = np.ones(shape,bool)
##    r = np.sqrt((x - 100)**2 + (y - 100)**2)       
##    mask = (r < radius/binn) * working_mask   
##    n_pix = len(np.where(mask)[0])
##    star = sum(star[np.where(mask)])*gain/em_gain
##    sky = 12.298897076737294
##    dc = 0.0007
##    
##    print(gain, star, sky, n_pix, dc, t_exp, rn, em_gain, nf),exit()
##    
##    image_noise = np.sqrt(star * nf**2 + n_pix * ( (rn/em_gain)**2 + (sky + dc) * nf**2 * t_exp * binn**2))    
##    SNR = star/image_noise
##    return SNR
##
##
##t_exp = [20, 20, 20, 20, 1, 1, 1, 1, 1, 1, 1, 1]
##gain = [3.37, 3.37, 0.8, 0.8, 15.9, 15.9, 3.38, 3.38, 3.38, 16, 16.4, 17.2]
##rn = [6.67, 6.94, 4.76, 4.79, 25.18, 34.72, 15.14, 12.71, 16, 85.60, 164.41, 292.40]
##em_gain = [1, 1, 1, 1, 20, 20, 20, 20, 25, 20, 20, 20]
##radius = [21.92, 21.94, 21.94, 21.95, 22.03, 22.02, 21.95, 21.92, 21.93, 21.96, 21.93, 21.92]
##nf = [1, 1, 1, 1, 1.41, 1.41, 1.41, 1.41, 1.41, 1.41, 1.41, 1.41]
##binn = [1,2,1,2,1,2,1,2,2,1,1,1]
###print(len(t_exp), len(gain), len(rn),len(em_gain), len(radius),len( nf))
##for i in range(1):
##    i=11
##    SNR = calc_SNR(t_exp[i], gain[i], rn[i], em_gain[i], radius[i], nf[i], binn[i])
##    print(SNR)

#------------------------------------------------------------------------------------------------------------------------------------

##import Photon_Flux_Calc_Bib as pfc
##import os
##from sys import exit
##import math
##
##dir_path = r'C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Testes do codigo\Teste setup ao redor do otimo\SNR'
##bias_name = 'bias_final_002.fits'
##lista = os.listdir(dir_path)
##lista.remove(bias_name)
##
##new_list = []
##for name in lista:
##    if '.fits' in name :new_list.append(name)
##
##for img in new_list:
##    if 'BIAS' not in img:
##        #Inicia o objeto para a redução da pré-imagem
##        pre_img_name = dir_path + '\\' + img
##        bias_img_name = dir_path + '\\' + img.split('.fits')[0] + '_BIAS.fits'        
##        PFC = pfc.PhotonFluxCalc(img_name = pre_img_name, bias_name = bias_img_name, xy_star = (100,100) , sky_radius = 15, ccd_serial=9917)
##        #Lê a imagem de bias fornecida
##        PFC.read_bias_img()
##        #Lê as informações do header da imagem
##        PFC.read_img_get_info()
##        #Procura pelo pixel com valor máximo
##        PFC.get_max_count()
##        #Configura o centroide da estrela com base no valor máximo de contagem encontrado
##        PFC.set_centroid()
##        #Calcula a FWHM
##        PFC.calc_FWHM()
##        #Calcula o fluxo da estrela e do céu com base no raio encontrado
##        PFC.calc_star_sky_flux()
##        #Calcula a SNR 
##        PFC.calc_SNR()
##        snr, fwhm, dc, rn, bias, sky_flux, star_flux, n_pixels, star_radius, x, y, t_exp = PFC.get_results()    
##        print(img, snr, star_radius)
##



#------------------------------------------------------------------------------------------------------------------------------------

##import Photon_Flux_Calc_Bib as pfc
##import os
##from sys import exit
##import math
##import astropy.io.fits as fits
##import numpy as np
##
##dir_path = r'C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Imagens OPD\17mar24\arsco'
##bias_name = 'bias_010.fits.gz'
##xy_star = (224,227)
##tag = 'b_ar'
##sky_radius = 15
##t_exp = 2.0
##
##lista = os.listdir(dir_path)
##lista.remove(bias_name)
##
##new_list = []
##for name in lista:
##    if '.fits' in name :new_list.append(name)
##
##dic = {}
##dic['img_name'] = []
##dic['snr'] = []
##dic['t_exp'] = []
##snr_total = 0
##i=0
##for img in new_list:
##    if tag not in img or '.fits' not in img: continue    
##    try:
##        #Inicia o objeto para a redução da pré-imagem
##        pre_img_name = dir_path + '\\' + img
##        header = fits.getheader(pre_img_name)
##        exposure = header['exposure'].split(',')
##        
##        if float(exposure[0]+'.'+exposure[1]) != t_exp : continue
##        
##        #bias_img_name = dir_path + '\\' + img.split('.fits')[0]+'_BIAS.fits'
##        bias_img_name = dir_path + '\\' + bias_name
##        PFC = pfc.PhotonFluxCalc(img_name = pre_img_name, bias_name = bias_img_name, xy_star = xy_star , sky_radius = sky_radius, ccd_serial=9917)
##        #Lê a imagem de bias fornecida
##        PFC.read_bias_img()
##        #Lê as informações do header da imagem
##        PFC.read_img_get_info()
##        #Procura pelo pixel com valor máximo
##        PFC.get_max_count()
##        #Configura o centroide da estrela com base no valor máximo de contagem encontrado
##        PFC.set_centroid()
##        #Calcula a FWHM
##        PFC.calc_FWHM()
##        #Calcula o fluxo da estrela e do céu com base no raio encontrado
##        PFC.calc_star_sky_flux()
##        #Calcula a SNR 
##        PFC.calc_SNR()
##        snr, fwhm, dc, rn, bias, sky_flux, star_flux, n_pixels, star_radius, x, y, t_exp = PFC.get_results()        
##        xy_star = (x,y)
##        if math.isnan(snr) == True: continue
##        snr_total += snr
##        i+=1
##
##       
##
##    except:continue
##        
##print(i,snr_total)
        
##        dic['img_name'].append(img)
##        dic['snr'].append(snr)
##        dic['t_exp'].append(t_exp)
    

##top_10_idx = np.argsort(dic['snr'])[-20:]
##top_10_names = [dic['img_name'][i] for i in top_10_idx]
##top_10_snr = [dic['snr'][i] for i in top_10_idx]
##top_10_t_exp = [dic['t_exp'][i] for i in top_10_idx]
##len_top10 = len(top_10_snr)
##for i in range(len_top10):            
##    print(top_10_names[len_top10-i-1], top_10_snr[len_top10-i-1], top_10_t_exp[len_top10-i-1])
##
#------------------------------------------------------------------------------------------------------------------------------------

##import json
##import Acquisition_Rate_Calculation_Bib as arc
##import SNR_Calculation_Bib as snrc
##
##dir_path = r'C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Testes do codigo\Teste setup ao redor do otimo\FA'
##arq = open(dir_path + '\\' + 'modos ao redor.txt', 'r')
##lines = arq.read().splitlines()
##
##dic = {}
##dic['mode'] = []
##dic['fa'] = []
##dic['t_exp'] = []
##for line in lines:
##    mode = json.loads(line)
##    min_t_exp = mode['min_t_exp']
##    em_mode = mode['em_mode']
##    hss = mode['hss']
##    binn = mode['binn']
##    sub_img = mode['sub_img']
##    preamp = mode['preamp']
##    max_em_gain = mode['em_gain']
##    ccd_temp = -60
##    sky_flux = 12.188074909054697 #e-/pix/s
##    star_flux = 113201.32209025856 #e-/s
##    n_pix_star = 1481
##    serial_number = 9917
##
##    
##    ARC = arc.AcquisitionRateCalc()        
##    ARC.write_operation_mode(mode['em_mode'], mode['hss'], mode['binn'],  mode['sub_img'], mode['min_t_exp']) 
##    ARC.seleciona_t_corte()
##    ARC.calc_acquisition_rate()    
##    max_acq_rate = float(ARC.return_acquisition_rate())
##  
##    dic['mode'].append(mode)
##    dic['fa'].append(max_acq_rate)  
##
####    SNRC = snrc.SignalToNoiseRatioCalc(min_t_exp, em_mode, max_em_gain, hss, preamp, binn, ccd_temp, sky_flux, star_flux, n_pix_star, serial_number)
####    SNRC.set_gain_value()
####    SNRC.calc_RN()
####    SNRC.calc_DC()
####    SNRC.calc_SNR()
####    snr = SNRC.get_SNR()
##
##    #print(em_mode, max_em_gain, hss, preamp, binn, sub_img, min_t_exp, max_acq_rate)
##
##top_10_idx = np.argsort(dic['fa'])[-11:]
##top_10_names = [dic['mode'][i] for i in top_10_idx]
##top_10_snr = [dic['fa'][i] for i in top_10_idx]
##len_top10 = len(top_10_snr)
##for i in range(len_top10):
##    s='1'
##    mode = top_10_names[len_top10-i-1]
##    if mode['em_mode'] == 1: s = '2'
##    if mode['hss'] == 0.1: s+='1'
##    if mode['hss'] == 1: s+='2'
##    if mode['hss'] == 10: s+='4'
##    if mode['hss'] == 20: s+='5'
##    if mode['hss'] == 30: s+='6'
##    s+= str(mode['preamp'])+str(mode['binn'])
##    if mode['sub_img'] == 256: s+='1'
##    if mode['sub_img'] == 512: s+='2'
##    if mode['sub_img'] == 1024: s+='3'
##    
##    print(s, mode['min_t_exp'], mode['em_gain'], top_10_snr[len_top10-i-1])


#-------------------------------------------------------------------------------------------------------------------------------
##import os
##import json
##import numpy as np
##from sys import exit
##dir_path = r'C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Testes do codigo\Teste setup ao redor do otimo\SNR'
##file_name = 'INITIAL_SETUPS.txt'
##
##arq = open(dir_path + '\\' + file_name, 'r')
##arq = arq.read().splitlines()
##
##for line in arq:
##    mode = json.loads(line)
##    em_mode = 'Conv'
##    em_gain = 1
##    if mode['em_mode'] == 1:
##        em_mode='EM'
##        em_gain = mode['em_gain']
##    print(em_mode, em_gain, mode['hss'], mode['preamp'], mode['bin'], mode['t_exp'])
  

#-------------------------------------------------------------------------------------------------------------------------------
##from sys import exit
##from math import ceil
##import json
##import Acquisition_Rate_Calculation_Bib as arc
##import random as rd
##import SNR_Calculation_Bib as snrc
##import numpy as np
##import os
##
##dir_path = r'C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Testes do codigo\Teste performance'
##file_name = 'CONV_1MHz_PA1_B1_TEXP20_G1_S3000_RAND_output_all_modes.txt'
##
##
##
##arq = open(dir_path + '\\' + file_name, 'r')
##lines = arq.read().splitlines()
##arq.close()
##dic = {}
##dic['output']=[]
##dic['mode']=[]
##dic['t_exp']=[]
##for line in lines:
##    line = line.split(';')
##    mode = json.loads(line[0])      
##    dic['mode'].append(mode)
##    dic['t_exp'].append(line[1]) 
##    dic['output'].append(float(line[2]))
##
##print(np.mean(dic['output'])),exit()
##
##top_10_idx = np.argsort(dic['snr'])[-11:]
##top_10_modes = [dic['mode'][i] for i in top_10_idx]
##top_10_output = [dic['snr'][i] for i in top_10_idx]
##top_10_t_exp = [dic['t_exp'][i] for i in top_10_idx]
##
##
##for i in range(len(top_10_output)):
##    index = len(top_10_output) - i - 1
##    mode = top_10_modes[index]
##    s='1'    
##    if mode['em_mode'] == 1: s = '2'
##    if mode['hss'] == 0.1: s+='1'
##    if mode['hss'] == 1: s+='2'
##    if mode['hss'] == 10: s+='4'
##    if mode['hss'] == 20: s+='5'
##    if mode['hss'] == 30: s+='6'
##    s+= str(mode['preamp'])+str(mode['binn'])
##    if mode['sub_img'] == 256: s+='1'
##    if mode['sub_img'] == 512: s+='2'
##    if mode['sub_img'] == 1024: s+='3'
##     
##    print(s, top_10_t_exp[index], mode['em_gain'], top_10_output[index])


#-------------------------------------------------------------------------------------------------------------------------------

##import astropy.io.fits as fits
##import os
##path = r'C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Imagens OPD\18mar06'
##lista_imagens = os.listdir(path)
##
##t_exp = 0
##desirable_t_exp = 180.0
##i=0
##for image in lista_imagens:
##    if 'OPD' in image and '.fits' in image:
##        header = fits.getheader(path + '\\' + image)                
##        #print(image, header['AIRMASS'], header['exposure'])
##        string_t_exp = header['exposure'].split(',')
##        string_t_exp = float(string_t_exp[0] + '.' + string_t_exp[1])        
##        if desirable_t_exp == string_t_exp:
##            i+=1
##        
##print(i)
        
        


#-------------------------------------------------------------------------------------------------------------------------------


##import os
##import json
##
##path = r'C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Testes do codigo\Teste performance\SNR_FA'
##
##lista = os.listdir(path)
##
##for file in lista:
##    if 'OPTSETUP' in file:
##        arq = open(path + '\\' + file, 'r')
##        mode = arq.read()
##        arq.close()
##        mode = json.loads(mode)        
##        print(file, mode['output'])



#-------------------------------------------------------------------------------------------------------------------------------




##path = r'C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Imagens OPD'
##import os
##import Acquisition_Rate_Calculation_Bib as arc
##import json
##from sys import exit
##
##files = os.listdir(path)
##lista = []
##for file in files:
##    if os.path.isdir(path+'\\'+file):
##        if 'Noites nao usadas' not in file:
##            #print(file)
##            if '17mar' in file: file+= '\\arsco'
##            sub_dir = path+'\\'+file
##            sub_files = os.listdir(sub_dir)
##            for sub_file in sub_files:
##                if 'SNR_FA_OPTSETUP' in sub_file or 'FA_SNR_OPTSETUP' in sub_file:
##                    name = path+'\\'+file + '\\' + sub_file
##                    arq = open(name,'r')
##                    arq = arq.read()                    
##                    mode = json.loads(arq)
##                    em_mode = mode['em_mode']
##                    em_gain = mode['em_gain']
##                    hss = mode['hss']
##                    preamp = mode['preamp']
##                    binn = mode['bin']
##                    sub_img = mode['sub_img']
##                    t_exp = mode['t_exp']
##
##                    ARC = arc.AcquisitionRateCalc()
##                    ARC.write_operation_mode(em_mode = em_mode, hss = hss, binn = binn, sub_img = sub_img, t_exp = t_exp)
##                    ARC.seleciona_t_corte()
##                    ARC.calc_acquisition_rate()
##                    FA = float(ARC.acquisition_rate)
##                    #FA = str(FA).split('.')
##                    SNR = round(mode['output']*(mode['mean_fa']*mode['mean_snr'])/FA,2)
##                    SNR = str(SNR).split('.')
##                    print(file, SNR[0]+','+SNR[1], FA)
##                    
##            
##exit()
    

##for i in range(len(modes)):
##    s = str(modes[i])
##    t_exp = t_exps[i]
##    em_gain = em_gains[i]
##    em_mode = 0
##    if s[0] == '2': em_mode = 1
##    hss = 0.1
##    if s[1] == '2': hss = 1
##    if s[1] == '4': hss = 10
##    if s[1] == '5': hss = 20
##    if s[1] == '6': hss = 30
##    preamp = int(s[2])
##    binn = int(s[3])    
##    sub_img = 256
##    if s[4] == '2': sub_img = 512
##    if s[4] == '3': sub_img = 1024
##
##    ARC = arc.AcquisitionRateCalc()
##    ARC.write_operation_mode(em_mode = em_mode, hss = hss, binn = binn, sub_img = sub_img, t_exp = t_exp)
##    ARC.seleciona_t_corte()
##    ARC.calc_acquisition_rate()
##    FA = round(float(ARC.acquisition_rate),3)
##    FA = str(FA).split('.')
##    print((FA[0]+','+FA[1],2))


#-------------------------------------------------------------------------------------------------------------------------------

##import SNR_Calculation_Bib as snrc
##import Photon_Flux_Calc_Bib as pfc
##
###modos = [22123,22123, 22123, 12113, 22113,  22123, 12113, 12213, 12213, 12213, 22123]  
##snrs = [83.16, 109.32, 187.08, 1594.81, 214.84, 1270.18, 7195.54, 54.9, 34.56, 37.65, 135.95, 106.15]
##
##path = r'C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Imagens OPD'
##
##import os
##import Acquisition_Rate_Calculation_Bib as arc
##import json
##from sys import exit
##i=0
##files = os.listdir(path)
##for file in files:
##    if os.path.isdir(path+'\\'+file):
##        if 'Noites nao usadas' not in file:            
##            if '17mar' in file: file+= '\\arsco'
##            if '17jul28' in file : continue
##            sub_dir = path+'\\'+file
##            sub_files = os.listdir(sub_dir)
##            for sub_file in sub_files:
##                if 'observation_setup.txt' in sub_file:                    
##                    arq = open(sub_dir + '\\' + 'observation_setup.txt', 'r')
##                    lines = arq.read().splitlines()
##                    arq.close()
##                    file_name = ''
##                    bias_name = ''
##                    obj_coords = ()
##                    star_radius = 0
##                    file_dir = sub_dir
##                    for line in lines:                        
##                        if 'Nome da imagem do objeto' in line: file_name = line.split('=')[1].replace(' ','')  
##                        if '.gz' not in file_name: file_name+='.gz'
##                        if 'Nome da imagem de bias' in line: bias_name = line.split('=')[1].replace(' ','')  
##                        if '.gz' not in bias_name: bias_name +='.gz'
##                        if 'Coordenadas' in line:
##                            line = line.split('=')[1].split(',')
##                            obj_coords = (int(line[0]), int(line[1]))
##                        if 'Raio maximo do ceu' in line: star_radius = int(line.split('=')[1])
##                    PFC = pfc.PhotonFluxCalc(img_name = file_dir + '\\' +  file_name, bias_name = file_dir + '\\' +  bias_name, xy_star = obj_coords, sky_radius=10, ccd_serial=9917)
##                    PFC.read_bias_img()
##                    PFC.read_img_get_info()
##                    PFC.get_max_count()
##                    PFC.set_centroid()
##                    PFC.calc_FWHM() 
##                    PFC.calc_star_sky_flux()
##                    minimum_t_exp = PFC.calc_minimun_texp_provided_SNR(snrs[i])
##                    PFC.calc_SNR()
##                    snr, fwhm, dc, rn, bias, sky_flux, star_flux, n_pixels, star_radius, x, y, t_exp= PFC.get_results()             
##                    print(file, snrs[i], snr)
##                    i+=1


#-------------------------------------------------------------------------------------------------------------------------------
##import os
##from sys import exit
##import SNR_Calculation_Bib as snrc
##import Photon_Flux_Calc_Bib as pfc
##import json
##import numpy
##import matplotlib.pyplot as plt
##
##
##dir_path = r'C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Imagens OPD\17mar22\arsco'
##files = os.listdir(dir_path)
##opt_SNR_file = ''
##opt_FA_file = ''
##for file in files:
##    if '_SNR_OPTSETUP' in file:
##        opt_SNR_file = file
##    if '_FA_OPTSETUP' in  file:
##        opt_FA_file = file
##
###Pega as informações da observação para o cálculo do fluxo da estrela        
##setup_file = 'observation_setup.txt'
##arq = open(dir_path + '\\' + setup_file, 'r')
##lines = arq.read().splitlines()
##arq.close()
##pre_img_name = ''
##bias_img_name = ''
##xy_star = ()
##star_radius = 0
##user_snr = 0
##for line in lines:                        
##    if 'Nome da imagem do objeto' in line: pre_img_name = line.split('=')[1].replace(' ','')  
##    if '.gz' not in pre_img_name: pre_img_name+='.gz'
##    if 'Nome da imagem de bias' in line: bias_img_name = line.split('=')[1].replace(' ','')  
##    if '.gz' not in bias_img_name: bias_img_name +='.gz'
##    if 'Coordenadas' in line:
##        line = line.split('=')[1].split(',')
##        xy_star = (int(line[0]), int(line[1]))
##    if 'Raio maximo do ceu' in line: star_radius = int(line.split('=')[1])
##    if 'Relacao sinal-ruido' in line: user_snr = float(line.split('=')[1])
##
##
##
###pre_img_name = 'J0957+0805_0174.fits.gz'
##             
###bias_img_name = dir_path + '\\' + img.split('.fits')[0]+'_BIAS.fits'
##bias_img_name = dir_path + '\\' + bias_img_name
##pre_img_name = dir_path + '\\' + pre_img_name
##PFC = pfc.PhotonFluxCalc(img_name = pre_img_name, bias_name = bias_img_name, xy_star = xy_star , sky_radius = star_radius, ccd_serial=9917)
###Lê a imagem de bias fornecida
##PFC.read_bias_img()
###Lê as informações do header da imagem
##PFC.read_img_get_info()
###Procura pelo pixel com valor máximo
##PFC.get_max_count()
###Configura o centroide da estrela com base no valor máximo de contagem encontrado
##PFC.set_centroid()
###Calcula a FWHM
##PFC.calc_FWHM()
###Calcula o fluxo da estrela e do céu com base no raio encontrado
##PFC.calc_star_sky_flux()
###Calcula a SNR 
##PFC.calc_SNR()
##snr, fwhm, dc, rn, bias, sky_flux, star_flux, n_pixels, star_radius, x, y, t_exp, binn, em_mode = PFC.get_results()
##user_mode = {'user_snr':user_snr , 'dc':dc , 'rn':rn , 'bias':bias , 'sky_flux':sky_flux , 'star_flux':star_flux , 'n_pixels':n_pixels , 't_exp':t_exp, 'binn':binn, 'em_mode':em_mode}
##
##
##def calc_SNR(em_mode, max_em_gain, max_t_exp, hss, binn, preamp, star_flux, sky_flux,n_pixels):
##    if em_mode ==  'CONV':em_mode=0
##    if em_mode == 'EM': em_mode=1
##    SNRC = snrc.SignalToNoiseRatioCalc(t_exp = t_exp, em_mode = em_mode, em_gain = max_em_gain, hss = hss, preamp = preamp, binn = binn, ccd_temp = -70, sky_flux = sky_flux , star_flux = star_flux, n_pix_star = n_pixels, serial_number = 9917)
##    SNRC.set_gain_value()
##    SNRC.calc_RN()
##    SNRC.calc_DC()
##    SNRC.calc_SNR()
##    snr = SNRC.get_SNR()
##    return snr
##
##
##
##   
##arq = open(dir_path + '\\' + opt_SNR_file, 'r')
##opt_mode = json.load(arq)
##arq.close()
##max_em_gain = opt_mode['em_gain']
##max_t_exp = opt_mode['t_exp']
##em_mode = opt_mode['em_mode']
##hss = opt_mode['hss']
##binn = opt_mode['bin']            
##preamp = opt_mode['preamp']
##sub_img = opt_mode['sub_img']
##snr_opt_mode = opt_mode['output']
##snr_list = []
##t_exp_range = np.linspace(1e-5, max_t_exp, 10)
##for t_exp in t_exp_range:
##        snr = calc_SNR(em_mode, max_em_gain, max_t_exp, hss, binn, preamp, star_flux, sky_flux, n_pixels)
##        snr_list.append(snr)     
##plt.plot(t_exp_range, snr_list, '-', c='r', label='Modo ótimo')
##plt.plot(max_t_exp, snr_opt_mode, 'o', c='r')
##
##
##arq = open(dir_path + '\\' + opt_FA_file, 'r')
##opt_mode = json.load(arq)
##arq.close()
##max_t_exp = opt_mode['t_exp']
##snr_opt_mode = user_mode['user_snr']
##plt.plot(max_t_exp, snr_opt_mode, 'o', c='r')
##
##                    
##
##
##max_t_exp = user_mode['t_exp']
##em_gain = 1
##n_pix = user_mode['n_pixels']
##binn = user_mode['binn']
##dc = user_mode['dc']
##rn = user_mode['rn']
##star = user_mode['star_flux']
##sky = user_mode['sky_flux']    
##nf = 1    
##if user_mode['em_mode']==1:nf=1.4
##snr_user = user_mode['user_snr']
##
##
##snr_list = []
##t_exp_range = np.linspace(1e-5, max_t_exp, 10)
##for t_exp in t_exp_range:   
##    aux = np.sqrt(star * t_exp * nf**2 + n_pix * ( (rn/em_gain/binn)**2 + (sky + dc)*t_exp * nf**2))        
##    SNR = (star*t_exp) / aux
##    snr_list.append(SNR)
##plt.plot(t_exp_range, snr_list, '-', c='b', label='Modo do usuário')
##plt.plot(max_t_exp, snr_user, 'o', c='b')
##
##
##plt.xlabel('Tempo de exposição (s)')
##plt.ylabel('SNR')
##plt.legend()
##plt.show()

#-----------------------------------------------------------------------------------------------------------------------------       

##import os
##from sys import exit
##import SNR_Calculation_Bib as snrc
##import Photon_Flux_Calc_Bib as pfc
##import json
##import numpy
##import matplotlib.pyplot as plt
##import Acquisition_Rate_Calculation_Bib as arc
##
##dir_path = r'C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Imagens OPD\17abr14'
##files = os.listdir(dir_path)
##opt_mode_file = ''
##for file in files:
##    if '_FA_OPTSETUP' in file:
##        opt_mode_file = file
##setup_file = 'observation_setup.txt'
##
##arq = open(dir_path + '\\' + setup_file, 'r')
##lines = arq.read().splitlines()
##arq.close()
##pre_img_name = ''
##bias_img_name = ''
##xy_star = ()
##star_radius = 0
##for line in lines:                        
##    if 'Nome da imagem do objeto' in line: pre_img_name = line.split('=')[1].replace(' ','')  
##    if '.gz' not in pre_img_name: pre_img_name+='.gz'
##    if 'Nome da imagem de bias' in line: bias_img_name = line.split('=')[1].replace(' ','')  
##    if '.gz' not in bias_img_name: bias_img_name +='.gz'
##    if 'Coordenadas' in line:
##        line = line.split('=')[1].split(',')
##        xy_star = (int(line[0]), int(line[1]))
##    if 'Raio maximo do ceu' in line: star_radius = int(line.split('=')[1])
##
##             
###bias_img_name = dir_path + '\\' + img.split('.fits')[0]+'_BIAS.fits'
##bias_img_name = dir_path + '\\' + bias_img_name
##pre_img_name = dir_path + '\\' + pre_img_name
##PFC = pfc.PhotonFluxCalc(img_name = pre_img_name, bias_name = bias_img_name, xy_star = xy_star , sky_radius = star_radius, ccd_serial=9917)
###Lê a imagem de bias fornecida
##PFC.read_bias_img()
###Lê as informações do header da imagem
##PFC.read_img_get_info()
###Procura pelo pixel com valor máximo
##PFC.get_max_count()
###Configura o centroide da estrela com base no valor máximo de contagem encontrado
##PFC.set_centroid()
###Calcula a FWHM
##PFC.calc_FWHM()
###Calcula o fluxo da estrela e do céu com base no raio encontrado
##PFC.calc_star_sky_flux()
###Calcula a SNR 
##PFC.calc_SNR()
##snr, fwhm, dc, rn, bias, sky_flux, star_flux, n_pixels, star_radius, x, y, t_exp = PFC.get_results()
##
##
##
##def calc_FA(em_mode, hss, binn, sub_img, t_exp):
##    ARC = arc.AcquisitionRateCalc()
##    ARC.write_operation_mode(em_mode = em_mode, hss = hss, binn = binn, sub_img = sub_img, t_exp = t_exp)
##    ARC.seleciona_t_corte()
##    ARC.calc_acquisition_rate()
##    FA = float(ARC.acquisition_rate)
##    return FA
##
##MOP_file = 'modos_opt_FA.txt'
##arq = open(dir_path + '\\' + MOP_file, 'r')
##lines = arq.read().splitlines()
##
##for line in lines:
##    line = line.split(';')
##    mode = json.loads(line[0])    
##    max_t_exp = mode['max_t_exp']
##    min_t_exp = mode['min_t_exp']
##    em_mode = mode['em_mode']
##    hss = mode['hss']
##    binn = mode['bin']            
##    preamp = mode['preamp']
##    sub_img = mode['sub_img']
##    snr_list = []
##    fa_lista = []
##    t_exp_range = np.linspace(min_t_exp, max_t_exp, 10)
##    for t_exp in t_exp_range:       
##        fa = calc_FA(em_mode, hss, binn, sub_img, t_exp)
##        fa_lista.append(fa)
##    plt.plot(t_exp_range, fa_lista, 'o', c='k', alpha=0.6)
##
##   
##arq = open(dir_path + '\\' + opt_mode_file, 'r')
##opt_mode = json.load(arq)
##arq.close()
##max_em_gain = opt_mode['em_gain']
##max_t_exp = 40
##min_t_exp = opt_mode['t_exp']
##em_mode = opt_mode['em_mode']
##hss = opt_mode['hss']
##binn = opt_mode['bin']            
##preamp = opt_mode['preamp']
##sub_img = opt_mode['sub_img']
##snr_list = []
##fa_lista = []
##t_exp_range = np.linspace(min_t_exp, max_t_exp, 10)
##for t_exp in t_exp_range:        
##        fa = calc_FA(em_mode, hss, binn, sub_img, t_exp)        
##        fa_lista.append(fa)
##plt.plot(t_exp_range, fa_lista, 'o', c='r', label='opt_mode')
##
##
##arq = open(dir_path + '\\' + opt_mode_file, 'r')
##opt_mode = json.load(arq)
##arq.close()                      
##user_mode = {'em_mode':1, 'em_gain':2, 'hss':1, 'preamp':1, 'bin':2 ,'sub_image':1024, 'max_t_exp':40, 'min_t_exp':1e-5}
##max_em_gain = user_mode['em_gain']
##max_t_exp = user_mode['max_t_exp']
##min_t_exp = user_mode['min_t_exp']
##em_mode = user_mode['em_mode']
##hss = user_mode['hss']
##binn = user_mode['bin']            
##preamp = user_mode['preamp']
##sub_img = mode['sub_img']
##snr_list = []
##fa_lista = []
##t_exp_range = np.linspace(min_t_exp, max_t_exp, 10)
##for t_exp in t_exp_range:        
##        fa = calc_FA(em_mode, hss, binn, sub_img, t_exp)
##        fa_lista.append(fa)
###plt.plot(t_exp_range, fa_lista, '-', c='b', label='user_mode')
##plt.xlabel('Tempo de exposição')
##plt.ylabel('SNR')
##plt.legend()
##plt.show()
##

#--------------------------------------------------------------------------------------------------------------------------------

import SNR_Calculation_Bib as snrc
import Photon_Flux_Calc_Bib as pfc
import os
import Acquisition_Rate_Calculation_Bib as arc
import json
from sys import exit
import astropy.io.fits as fits
import numpy as np

def calc_quadratic_equation(a, b, c):
        delta = (b**2) - (4*a*c)
        x1 = (-b-np.sqrt(delta))/(2*a)
        x2 = (-b+np.sqrt(delta))/(2*a)
        return [x1,x2]

path = r'C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Imagens OPD'
i=0
directories = os.listdir(path)
for directory in directories:
    if os.path.isdir(path + '\\' + directory):
        if 'Noites nao usadas' not in directory:            
            if '17mar' in directory: directory+= '\\arsco'
            if '17jul28' in directory : continue
            if 'ZIP' in directory : continue
            #if '17ago08' not in directory:continue
            sub_directory = path + '\\' + directory
            sub_files_names = os.listdir(sub_directory)
            opt_mode = {}
            pre_img_name = ''
            bias_img_name = ''
            xy_star = ()
            star_radius = 0            
            for file in sub_files_names:
                if 'observation_setup.txt' in file:                    
                    arq = open(sub_directory + '\\' + 'observation_setup.txt', 'r')
                    lines = arq.read().splitlines()                             
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
                if 'BOTH_OPTSETUP.txt' in file:                    
                    arq = open(sub_directory + '\\' + file, 'r')
                    opt_mode = json.load(arq)
                    arq.close()
            
            #bias_img_name = dir_path + '\\' + img.split('.fits')[0]+'_BIAS.fits'
            bias_img_name = sub_directory + '\\' + bias_img_name
            pre_img_name = sub_directory + '\\' + pre_img_name
            PFC = pfc.PhotonFluxCalc(img_name = pre_img_name, bias_name = bias_img_name, xy_star = xy_star , sky_radius = star_radius, ccd_serial=9917)
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
            
            image, hdr = fits.getdata(pre_img_name, header=True)            
            em_gain = 1
            n_pix = n_pixels
            dc = dc
            #rn = np.std(image)*float(hdr['gain'])
            star = star_flux
            sky = sky_flux
            nf = 1
            if 'Electron Multiplying' in hdr['OUTPTAMP']: nf = 1.41           
            
            binn = int(hdr['vbin'])
            if opt_mode['SNR'] > snr:
                snr = opt_mode['SNR']            
            
            
            a = star**2
            b = snr**2 * nf**2 * (star+n_pix*(sky+dc))
            c = snr**2*n_pix*(rn/em_gain/binn)**2        
            minimun_t_exps = calc_quadratic_equation(a, -b, -c)       
            for t_exp in minimun_t_exps:
                if t_exp <= 0: minimun_t_exps.remove(t_exp)                 
            t_exp = min(minimun_t_exps)
            print(sub_directory.split('\\')[-1], t_exp, snr)
#-----------------------------------------------------------------------------------------------------------------------------------

##import matplotlib.pyplot as plt
##import numpy as np
##
##x = np.linspace(1,10,3)
##color = ['b','r','k']
##i=0
##for a in [1,2,3]:
##    y = a*x
##    plt.semilogy(x,y, '-', c=color[i])
##    i+=1
##
##plt.show()
#-----------------------------------------------------------------------------------------------------------------------------------
##data = '2017-03-25T04:29:08.946'
##seconds1 = data.split('T')[1]
##seconds1 = seconds1.split(':')
##seconds1 = int(seconds1[0])*3600 + int(seconds1[1])*60 + float(seconds1[2])
##data = '2017-03-25T07:49:38.030'
##seconds2 = data.split('T')[1]
##seconds2 = seconds2.split(':')
##seconds2 = int(seconds2[0])*3600 + int(seconds2[1])*60 + float(seconds2[2])
##
##print(seconds2 + 2 - seconds1)
