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


img_dir =r'C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Testes do codigo'

#Cria a classe que otimiza o modo de operação
OOM = oom.Optimize_Operation_Mode(img_dir, algorithm = tpe.suggest)
#verifica se os modos fornecidos de sub_img e binning estao corretos
OOM.verify_provides_modes()
#Esta opção serve para o cálculo do fluxo da estrela
OOM.calc_star_flux()
#Optimize(1-SNR, 2-FA, 3-Ambos)
OOM.optimize(3)

#OOM.create_space_of_states_provided_FA()
#OOM.create_space_of_states_with_all_modes()
#OOM.get_mean_radom_SNR()
#OOM.get_mean_radom_FA()
#OOM.get_mean_random_SNR_FA()


