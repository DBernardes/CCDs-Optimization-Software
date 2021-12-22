#!/usr/bin/env python
# coding: utf-8
# 25/10/2019. Denis Varise Bernardes.


import os
from pathlib import Path
from sys import exit

from hyperopt import rand, tpe

import Optimize_Operation_Mode as oom

img_dir = os.path.dirname(os.getcwd()) + "\example\\"
# Create the optimization class
OOM = oom.Optimize_Operation_Mode(img_dir, algorithm=tpe.suggest)
# Check the provided sub-image and binning modes
OOM.verify_provides_modes()
# Calculate the star flux
OOM.calc_star_flux()
# Optimize(1-SNR, 2-FA, 3-Both)
OOM.optimize(3)
