#!/usr/bin/env python
# coding: utf-8
#25/10/2019. Denis Varise Bernardes.


from sys import exit
import Optimize_Operation_Mode as oom
from hyperopt import tpe, rand
import os


img_dir =r'C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Testes do codigo'
#Create the optimization class
OOM = oom.Optimize_Operation_Mode(img_dir, algorithm = tpe.suggest)
#Check the provided sub-image and binning modes
OOM.verify_provides_modes()
#Calculate the star flux
OOM.calc_star_flux()
#Optimize(1-SNR, 2-FA, 3-Both)
OOM.optimize(3)




