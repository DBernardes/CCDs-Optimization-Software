#!/usr/bin/env python
# coding: utf-8
# 10/01/2020. Denis Varise Bernardes.

import json
import os
from sys import exit

from hyperopt import rand, tpe


def create_input_file(img_dir):
    """Create input file"""

    file_parameters = {
        "snr": 100,
        "acq_rate": 1,
        "mag": 12,
        "t_exp": [1e-5, 84600],
        "em_mode": ["EM", "Conv"],
        "readout_rate": [0.1, 1, 10, 20, 30],
        "preamp": [1, 2],
        "sub_img": [1024, 512, 256],
        "bin": [1, 2],
        "temperature": -70,
        "serial_number": 9917,
        "use_pre_img": True,
        "pre_img_name": "",
        "path": "",
        "star_coords": [0, 0],
        "bias_img_name": "",
        "sky_radius": 20,
        "export_setup_file": True,
        "export_iterations_file": False,
        "export_bias_file": False,
        "exp_name": "",
        "iteration_number": 10,
        "optimization_algorithm": tpe.suggest,
    }
    file_path = os.path.join(img_dir, "observation_setup.txt")
    with open(file_path, "w") as arq:
        json.dump(file_parameters, arq)
        arq.close()


def get_obs_setup(img_dir):
    # This function reads the info of the observation

    # Check if there is the 'observation_setup' in the provided directory
    # Otherwise, this function will create one
    create_arq_obs_setup(img_dir)

    snr = 0
    acq_rate = 0
    obj_magnitude = 0
    sub_img_modes = 0
    binn_modes = 0
    serial_number = 0
    ccd_temp = 0
    n_pix_star = 0
    max_evals = 0
    fixar_param = 0
    export_arq = ""
    export_loss = ""
    export_bias = ""

    use_pre_img = ""
    pre_img_name = ""
    obj_coords = ()
    bias_img_name = ""
    sky_radius = 0

    # opens the file and reads the info of each line
    arq = open(img_dir + "\\" + "observation_setup.txt", "r")
    lines = arq.read().splitlines()
    for line in lines:
        line = line.split("=")
        if "Signal-to-Noise ratio" in line[0]:
            snr = float(line[1])
        if "Acquisition rate" in line[0]:
            acq_rate = float(line[1])
        if "Object magnitude" in line[0]:
            obj_magnitude = float(line[1])
        if "Sub-image modes" in line[0]:
            sub_img_modes = [int(sub_img) for sub_img in line[1].split(",")]
        if "Binning modes" in line[0]:
            binn_modes = [int(binn) for binn in line[1].split(",")]
        if "CCD serial number" in line[0]:
            serial_number = int(line[1])
        if "CCD temperature" in line[0]:
            ccd_temp = float(line[1])
        if "Iterations number of the optimization" in line[0]:
            max_evals = int(line[1])
        if "Experiment name" in line[0]:
            file_base_name = line[1].replace(" ", "")
        if "Export the opt mode to a txt file" in line[0]:
            export_arq = line[1].replace(" ", "")
        if "Export the list of the iterations" in line[0]:
            export_loss = line[1].replace(" ", "")
        if "Export the list of the bias images" in line[0]:
            export_bias = line[1].replace(" ", "")

        if "Use a pre-image" in line[0]:
            use_pre_img = line[1].strip().replace(" ", "")
        if use_pre_img == "y":
            if "Image name" in line[0]:
                pre_img_name = line[1].strip().replace(" ", "")
                if ".fits" not in pre_img_name:
                    pre_img_name += ".fits"
            if "Object coordinates" in line[0]:
                obj_coords = (int(line[1].split(",")[0]), int(line[1].split(",")[1]))
            if "Bias image name" in line[0]:
                bias_img_name = line[1].strip().replace(" ", "")
                if ".fits" not in bias_img_name:
                    bias_img_name += ".fits"
            if "Maximum sky radius" in line[0]:
                sky_radius = int(line[1])

    return (
        snr,
        acq_rate,
        obj_magnitude,
        sub_img_modes,
        binn_modes,
        serial_number,
        ccd_temp,
        max_evals,
        file_base_name,
        export_arq,
        export_loss,
        export_bias,
        use_pre_img,
        pre_img_name,
        obj_coords,
        bias_img_name,
        sky_radius,
    )
