#!/usr/bin/env python
# coding: utf-8
# 10/01/2020. Denis Varise Bernardes.

from sys import exit


def create_arq_obs_setup(img_dir):
    # This function tries to open the file with the info of the observation.
    # If this file does not exists, this function creates it.
    try:
        arq = open(img_dir + "\\" + "observation_setup.txt", "r")
    except:
        print("Fill in the file with the observation configuration.")
        s = """File with the observation setup to be used in optimization.
-------------------------------------------------------------

Signal-to-Noise ratio = 100
Acquisition rate (Hz) = 1
Object magnitude = 12.25 
Sub-image modes = 1024,512,256
Binning modes = 1,2
CCD serial number = 9917
CCD temperature (oC) = -70 
Iterations number of the optimization = 170
Experiment name = EXP1
Export the opt mode to a txt file (y/n) = y
Export the list of the iterations (y/n) = y
Export the list of the bias images (y/n) = n



Pre-image configuration
----------------------------
Use a pre-image (y/n) = n
Image name = 
Object coordinates (x,y) = 
Bias image name = 
Maximum sky radius = 20
        """
        arq = open(img_dir + "\\" + "observation_setup.txt", "w")
        arq.write(s)
        arq.close()
        exit()


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
