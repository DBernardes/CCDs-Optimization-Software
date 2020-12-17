#!/usr/bin/env python
# coding: utf-8

# Este código para a determinacao do
# modo de operacao de uma camera iXon Ultra 888 que otimize o ruido
# de leitura, a taxa de aquisicao, ou ambos.

#25/10/2019. Denis Varise Bernardes.


from sys import exit
import Optimize_Operation_Mode as oom
from hyperopt import tpe, rand
import os


img_dir =r'C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Testes do codigo\Teste_funcao_aquisicao'
bias_file = 'bias_final_002.fits'
lista = os.listdir(img_dir)
lista.remove(bias_file)

new_list = []
for name in lista:
    if '.fits' in name :new_list.append(name)

#for name in new_list:print(name)
#exit()
    

for image_name in new_list:
    if 'BIAS' not in image_name:
        image_name = image_name.split('.fits')[0]
        print('\n'+image_name)
        s = '''

Arquivo contendo a configuracao da observacao a ser otimizada
-------------------------------------------------------------

Relacao sinal-ruido = 100
Taxa de aquisicao (Hz) = 2
Magnitude do objeto = 8 
Modos de sub-imagem = 1024,512,256
Binning dos pixeis = 1,2
Numero serial do CCD = 9917
Temperatura do CCD (oC) = -70 
Numero de iteracao do MOB = 170
Nome do experimento = %s_RAND
Exportar configuracao para txt (s/n)? = s
Exportar lista de imagens de loss (s/n)? = s
Exportar lista de imagens de bias (s/n)? = n



Configuracoes da pre-imagem
----------------------------

Utilizar pre-imagem (s/n) = s
Nome da imagem do objeto = %s
Coordenadas do objeto (x,y) = 101,101
Nome da imagem de bias = %s_BIAS
Raio maximo do ceu = 20       
    '''%(image_name, image_name, image_name)
        arq =  open(img_dir + '\\' + 'observation_setup.txt', 'w')    
        arq.write(s)
        arq.close()
        #Cria a classe que otimiza o modo de operação
        OOM = oom.Optimize_Operation_Mode(img_dir, algorithm = rand.suggest)
        #verifica se os modos fornecidos de sub_img e binning estao corretos
        OOM.verify_provides_modes()
        #Esta opção serve para o cálculo do fluxo da estrela
        OOM.calc_star_flux()
        #OPT(1-SNR, 2-FA, 3-Ambos)
        OOM.optimize(1)

        #OOM.create_space_of_states_provided_FA()
        #OOM.create_space_of_states_with_all_modes()
        #OOM.get_mean_radom_SNR()
        #OOM.get_mean_radom_FA()
        #OOM.get_mean_random_SNR_FA()

    


img_dir =r'C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Testes do codigo\Teste_funcao_aquisicao'

coords = [(508,369),(484,475),(501,231),(789,813),(691,934),(569,984),(515,1001),(845,83),(742,485)]

##i=1
##for coord in coords:
##    s = '''
##Arquivo contendo a configuracao da observacao a ser otimizada
##-------------------------------------------------------------
##
##Relacao sinal-ruido = 50
##Taxa de aquisicao (Hz) = 2
##Magnitude do objeto = 12.6 
##Modos de sub-imagem = 1024,512,256
##Binings dos pixieis = 1,2
##Numero serial do CCD = 9917
##Temperatura do CCD (oC) = -70 
##Numero de iteracao do MOB = 100
##Nome do experimento = RAND_%i
##Exportar configuracao para txt (s/n)? = s
##Exportar lista de imagens de loss (s/n)? = s
##Exportar lista de imagens de bias (s/n)? = n
##
##
##
##Configuracoes da pre-imagem
##----------------------------
##
##Utilizar pre-imagem (s/n) = s 
##Nome da imagem do objeto = hats-24_I_transito_001
##Coordenadas do objeto (x,y) = %i,%i
##Nome da imagem de bias = bias_final_002.fits
##Raio maximo do ceu = 20
##'''%(i, coord[0], coord[1])    
##
##    arq =  open(img_dir + '\\' + 'observation_setup.txt', 'w')    
##    arq.write(s)
##    arq.close()
##
##    #Cria a classe que otimiza o modo de operação
##    OOM = oom.Optimize_Operation_Mode(img_dir, algorithm = rand.suggest)
##    #verifica se os modos fornecidos de sub_img e binning estao corretos
##    OOM.verify_provides_modes()
##    #Esta opção serve para o cálculo do fluxo da estrela
##    OOM.calc_star_flux()
##    #Fixar(1-FA, 2-SNR, 3-Ambos)
##    OOM.optimize(3)
##    i+=1



##for coord in coords:
##    s = '''
##Arquivo contendo a configuracao da observacao a ser otimizada
##-------------------------------------------------------------
##
##Relacao sinal-ruido = 100
##Taxa de aquisicao (Hz) = 2
##Magnitude do objeto = 12.6 
##Modos de sub-imagem = 1024,512,256
##Binings dos pixieis = 1,2
##Numero serial do CCD = 9917
##Temperatura do CCD (oC) = -70 
##Numero de iteracao do MOB = 100
##Nome do experimento = hats-24_I_transito_001
##Exportar configuracao para txt (s/n)? = s
##Exportar lista de imagens de loss (s/n)? = n
##Exportar lista de imagens de bias (s/n)? = n
##
##
##
##Configuracoes da pre-imagem
##----------------------------
##
##Utilizar pre-imagem (s/n) = s 
##Nome da imagem do objeto = hats-24_I_transito_001
##Coordenadas do objeto (x,y) = %i,%i
##Nome da imagem de bias = bias_final_002.fits
##Raio maximo do ceu = 20
##'''%(coord[0], coord[1])
##    arq =  open(img_dir + '\\' + 'observation_setup.txt', 'w')    
##    arq.write(s)
##    arq.close()
##    #Cria a classe que otimiza o modo de operação
##    OOM = oom.Optimize_Operation_Mode(img_dir, algorithm = tpe.suggest)
##    #verifica se os modos fornecidos de sub_img e binning estao corretos
##    OOM.verify_provides_modes()
##    #Esta opção serve para o cálculo do fluxo da estrela
##    OOM.calc_star_flux()
##    #OOM.create_space_of_states_provided_FA()
##    OOM.create_space_of_states_with_all_modes()
##    #OOM.get_mean_radom_SNR()
##    OOM.get_mean_random_SNR_FA()
##


##lista = [
##C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Imagens OPD\17abr14,
##C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Imagens OPD\17abr15,
##C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Imagens OPD\17abr16,
##C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Imagens OPD\17ago08,
##C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Imagens OPD\17jul02,
##C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Imagens OPD\17jul13,
##C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Imagens OPD\17jul28,
##C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Imagens OPD\17jun03,
##C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Imagens OPD\17mar22,
##C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Imagens OPD\17mar23,
##C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Imagens OPD\17mar24,
##C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Imagens OPD\17nov24,
##C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Imagens OPD\18mar06]
##
##
##
##
