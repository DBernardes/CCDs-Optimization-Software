#!/usr/bin/env python
# coding: utf-8
# 25/10/2019. Denis Varise Bernardes.


import os

from hyperopt import rand, tpe

import Optimize_Operation_Mode as oom

img_dir = os.path.join(os.getcwd(), "example")
OOM = oom.Optimize_Operation_Mode(img_dir, algorithm=tpe.suggest)
OOM.verify_provides_modes()
OOM.calc_star_flux()
# Optimize(1-SNR, 2-FA, 3-Both)
OOM.optimize(1)
