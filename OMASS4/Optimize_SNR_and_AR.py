#!/usr/bin/env python
# coding: utf-8
# Denis Varise Bernardes.
# 15/01/2020.


import collections
import gc
import json
import os
import random as rd
from copy import copy
from math import ceil
from sys import exit

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from hyperopt import Trials, fmin, hp, rand, tpe
from hyperopt.pyll import scope
from hyperopt.pyll.stochastic import sample

import AR_Calc as arc
import Modos_Operacao_Bib as mob
import Optimize_AR as oar
import Optimize_SNR as osnr
import SNR_Calc as snrc


def function_fa(parameters=[]):
    t_exp = parameters[0]
    em_mode = parameters[1]
    em_gain = parameters[2]
    hss = parameters[3]
    preamp = parameters[4]
    binn = parameters[5]
    sub_img = parameters[6]
    CCD_temp = parameters[7]
    sky_flux = parameters[8]
    star_flux = parameters[9]
    n_pixels_star = parameters[10]
    serial_number = parameters[11]
    max_fa = parameters[14]
    min_fa = parameters[15]

    ARC = arc.AcquisitionRateCalc()
    ARC.write_operation_mode(em_mode, hss, binn, sub_img, t_exp)
    ARC.seleciona_t_corte()
    ARC.calc_acquisition_rate()
    acq_rate = ARC.acquisition_rate
    norm_acq_rate = (acq_rate - min_fa) / (max_fa - min_fa)
    return norm_acq_rate


def function_snr(parameters=[]):
    t_exp = parameters[0]
    em_mode = parameters[1]
    em_gain = parameters[2]
    hss = parameters[3]
    preamp = parameters[4]
    binn = parameters[5]
    sub_img = parameters[6]
    ccd_temp = parameters[7]
    sky_flux = parameters[8]
    star_flux = parameters[9]
    n_pix_star = parameters[10]
    serial_number = parameters[11]
    max_snr = parameters[12]
    min_snr = parameters[13]
    snr_target = parameters[16]
    SNRC = snrc.SignalToNoiseRatioCalc(
        t_exp=t_exp,
        em_mode=em_mode,
        em_gain=em_gain,
        hss=hss,
        preamp=preamp,
        binn=binn,
        ccd_temp=ccd_temp,
        sky_flux=sky_flux,
        star_flux=star_flux,
        n_pix_star=n_pix_star,
        serial_number=serial_number,
    )
    SNRC.calc_RN()
    SNRC.calc_DC()
    SNRC.calc_SNR()
    snr = SNRC.get_SNR()
    if snr < snr_target:
        snr = 0
        min_snr = 0
    norm_snr = (snr - min_snr) / (max_snr - min_snr)
    return norm_snr


def function(parameters=[]):
    # This is the function that will be optimized. It is composed by the
    # normalized values of the SNR and the acquisition rate. These values, in turn,
    # are calculated by using their respective libraries. For this reason, it is needed to provide
    # to this funciton a list with the operation mode of the CCD, the star and sky flux, and the
    # maximum and minimum values of the SNR and acquisition rate. If the calculates values by the functions
    # are smaller than their target values, the functions return zero.
    return -function_snr(parameters) * function_fa(parameters)


# -------------------------------------------------------------------


class Opt_SignalNoiseRatio_AcquisitionRate:
    def __init__(
        self,
        snr_target,
        acq_rate,
        sub_img_modes,
        binn_modes,
        serial_number,
        ccd_temp,
        n_pix_star,
        sky_flux,
        star_flux,
        bias_level,
    ):
        self.hss = [[], []]
        self.binn = [[], []]
        self.sub_img = []
        self.ccd_temp = ccd_temp
        self.serial_number = serial_number
        self.gain = 0
        self.dark_noise = 0
        self.set_dc()
        self.binn_modes = binn_modes
        self.sub_img_modes = sub_img_modes

        self.sky_flux = sky_flux  # e-/pix/s
        self.star_flux = star_flux  # e-/s
        self.n_pix_star = n_pix_star
        self.bias_level = bias_level

        self.best_mode = {}
        self.space = []
        self.list_all_modes = []
        self.best_sub_img = []
        self.acq_rate_target = acq_rate
        self.snr_target = snr_target
        self.losses_SNR = []
        self.losses_FA = []
        self.MOB = mob.ModosOperacao()
        self.max_snr = 0
        self.min_snr = 0
        self.max_fa = 0
        self.min_fa = 0
        self.space_all_modes = []
        self.new_list = []

    def set_gain(self, em_mode, hss, preamp):
        # Select the value of the CCD gain
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
        self.gain = gain

    def set_dc(self):
        # Calculates the dark current based ib the CCD temperature.
        # These used equations model the dark current of the SPARC4 CCDs,
        # and can be found in: D V Bernardes et al 2018 PASP 130 095002
        T = self.ccd_temp
        if self.serial_number == 9914:
            self.dark_noise = 24.66 * np.exp(0.0015 * T ** 2 + 0.29 * T)
        if self.serial_number == 9915:
            self.dark_noise = 35.26 * np.exp(0.0019 * T ** 2 + 0.31 * T)
        if self.serial_number == 9916:
            self.dark_noise = 9.67 * np.exp(0.0012 * T ** 2 + 0.25 * T)
        if self.serial_number == 9917:
            self.dark_noise = 5.92 * np.exp(0.0005 * T ** 2 + 0.18 * T)

    def write_MOB_obj(self, obj):
        # Write a object with the list of modes in the class
        self.MOB = obj

    def print_MOB_list(self):
        lista = self.MOB.get_list_of_modes()
        for modo in lista:
            print(modo)

    def SNR_FA_ranges(self):
        # This function calculates hte minimum and maximum values of the
        # SNR and the AR for the normalization of the cost function.

        # Starts the class of the SNR optimization
        OSNR = osnr.OptSignalNoiseRation(
            snr_target=self.snr_target,
            serial_number=self.serial_number,
            ccd_temp=self.ccd_temp,
            n_pix_star=self.n_pix_star,
            sky_flux=self.sky_flux,
            star_flux=self.star_flux,
            bias_level=self.bias_level,
        )
        OSNR.write_MOB_obj(copy(self.MOB))
        OSNR.duplicate_list_of_modes_for_PA12()
        OSNR.remove_repeat_modes()
        # Calculates the maximum SNR
        max_snr = OSNR.calc_best_mode()
        # Calculates the minimum SNR
        min_snr = OSNR.calc_min_snr()
        # Starts the class of the acquisition rate optimization
        OAR = oar.OptimizeAcquisitionRate(
            acquisition_rate=self.acq_rate_target,
            sub_img_modes=self.sub_img_modes,
            binn_modes=self.binn_modes,
        )
        OAR.write_MOB_obj(copy(self.MOB))
        # Calculates the minimum acquisition rate
        min_fa = OAR.determine_min_acquisition_rate()
        # Calculates the maximum acquisition rate
        best_mode, max_fa = OAR.determine_fastest_operation_mode()
        self.max_snr = max_snr
        self.min_snr = min_snr
        self.max_fa = max_fa
        self.min_fa = min_fa

    def calc_max_em_gain(self, max_t_exp, min_t_exp):
        # Calculation of the maximum EM gain allowed as a function of the CCD operating mode.
        # This calculation takes into account the maximum amount of 100 photons per pixel for which
        # the EM mode is better than the Conventional one.
        sky_flux = self.sky_flux
        star_flux = self.star_flux
        n_pix = self.n_pix_star
        dn = self.dark_noise
        gain = self.gain
        max_fotons = 100
        bias = self.bias_level
        max_ADU = (2 ** 16) * 0.8

        aux = (sky_flux + star_flux / n_pix + dn) * max_t_exp / gain
        max_em_gain = (max_ADU - bias) / aux

        # if the photons/pix ir bigger than 100 photons, it is calculated the EM gain and the
        # exposure time values that accomplish this limit
        if aux > max_fotons:
            max_em_gain = (max_ADU - bias) / (max_fotons / gain)
            max_t_exp = max_fotons / (sky_flux + star_flux / n_pix + dn)
        # If the max exposure time found in the previous step is smaller than the min exposure time,
        # this mode is discarded
        if max_t_exp < min_t_exp:
            max_em_gain = 0
            max_t_exp = 0
        return max_em_gain, max_t_exp

    def create_space(self):
        # This function is responsible to get the list of selected modes and
        # convert it in the hyperopt library space of states
        i = 0
        self.new_list = []
        space_all_modes = []
        # iterates each mode of the list of selected modes
        for mode in self.MOB.get_list_of_modes():
            max_em_gain = 0
            max_t_exp = mode["max_t_exp"]
            min_t_exp = mode["min_t_exp"]
            # the exposure time is uniform between its maximum and minimu vales
            t_exp = hp.uniform("t_exp_" + str(i), min_t_exp, max_t_exp)
            # the EM mode is a discreet choice between 0 and 1.
            em_mode = hp.choice("em_mode_" + str(i), [mode["em_mode"]])
            # if the EM mode is conventional, the em gain is fixed in 1
            em_gain = hp.choice("em_gain_" + str(i), [max_em_gain])
            # if the EM mode is EM...
            if mode["em_mode"] == 1:
                # Sets the CCD gain
                self.set_gain(mode["em_mode"], mode["hss"], mode["preamp"])
                # Calculates the maximum EM gain
                max_em_gain, max_t_exp = self.calc_max_em_gain(max_t_exp, min_t_exp)
                # Recreates the space of states for the exposure time
                t_exp = hp.uniform("t_exp_" + str(i), min_t_exp, max_t_exp)
                # if the EM gain is zero, the mode is discarded
                if max_em_gain == 0:
                    print(
                        "The number of photons per pixel is above 100. This mode was rejected."
                    )
                    continue
                # limits the maximum EM gain in 300x
                if max_em_gain > 300:
                    max_em_gain = 300
                # Sets the EM gain as uniform within its maximum and minimum values
                em_gain = hp.uniform("em_gain_" + str(i), 2, max_em_gain)
            # the HSS is a discreet choice between its allowed values
            hss = hp.choice("hss_" + str(i), [mode["hss"]])
            # the PREAMP is a discreet choice between 1 and 2
            preamp = hp.choice("preamp_" + str(i), [mode["preamp"]])
            # the BINNING is a discreet choice between 1 and 2
            binn = hp.choice("binn_" + str(i), [mode["binn"]])
            # the Sub-Image is a discreet choice between its allowed values
            sub_img = hp.choice("sub_img_" + str(i), [mode["sub_img"]])
            # the CCD temperature is fixed for the povided value
            ccd_temp = hp.choice("ccd_temp_" + str(i), [self.ccd_temp])
            # the sky flux is fixed for the povided value
            sky_flux = hp.choice("sky_flux_" + str(i), [self.sky_flux])
            # the star flux is fixed for the povided value
            star_flux = hp.choice("obj_flux_" + str(i), [self.star_flux])
            # the number of pixels is fixed for the povided value
            n_pix_star = hp.choice("n_pix_star_" + str(i), [self.n_pix_star])
            # the serial number is fixed for the povided value
            serial_number = hp.choice("serial_number" + str(i), [self.serial_number])
            # creates a list with the CCD operation mode in the hyperopt space of states
            new_mode = [
                t_exp,
                em_mode,
                em_gain,
                hss,
                preamp,
                binn,
                sub_img,
                ccd_temp,
                sky_flux,
                star_flux,
                n_pix_star,
                serial_number,
            ]
            # Appends this mode to a list
            self.space_all_modes.append(new_mode)
            # This is a mirror list of the hyperopt list for later data management
            self.new_list.append(mode)
            i += 1
        self.space = hp.choice("operation_mode", self.space_all_modes)

    def run_bayesian_opt(self, max_evals, algorithm=tpe.suggest):
        gc.collect()
        # Create the algorithm
        tpe_algo = algorithm
        # Create a trials object
        self.tpe_trials = Trials()
        # Normalization parameters
        max_snr = self.max_snr
        min_snr = self.min_snr
        max_fa = self.max_fa
        min_fa = self.min_fa
        # Run evals with the tpe algorithm
        best_mode = fmin(
            fn=function,
            space=self.space + [max_snr, min_snr, max_fa, min_fa, self.snr_target],
            algo=tpe_algo,
            trials=self.tpe_trials,
            max_evals=max_evals,
            rstate=np.random.RandomState(50),
        )
        # Gets the best mode
        index_list_modes = best_mode["operation_mode"]
        chosen_mode = self.new_list[index_list_modes]
        # Recovers the mode parameters through the mirror list
        t_exp = best_mode["t_exp_" + str(index_list_modes)]
        em_mode = chosen_mode["em_mode"]
        em_gain = best_mode["em_gain_" + str(index_list_modes)]
        hss = chosen_mode["hss"]
        preamp = chosen_mode["preamp"]
        binn = chosen_mode["binn"]
        sub_img = chosen_mode["sub_img"]
        # Write the best mode as a dictionary format
        self.best_mode = {
            "t_exp": t_exp,
            "em_mode": em_mode,
            "em_gain": em_gain,
            "hss": hss,
            "preamp": preamp,
            "binn": binn,
            "sub_img": sub_img,
            "SNR*FA": -self.tpe_trials.best_trial["result"]["loss"],
        }

    def print_best_mode(self):
        if self.best_mode["em_mode"] == 1:
            print("\nEM Mode")
            print("-------")
            print("EM gain: ", self.best_mode["em_gain"])
        else:
            print("\nConventional Mode")
            print("-----------------")
        print("Exposure time (s): ", self.best_mode["t_exp"])
        print("HSS: ", self.best_mode["hss"])
        print("Preamp: ", self.best_mode["preamp"])
        print("Binning: ", self.best_mode["binn"])
        print("Sub image: ", self.best_mode["sub_img"])
        print("\nBest SNR*FA: ", self.best_mode["SNR*FA"])

    def creat_log_parameters(self, path, file_base_name):
        # In this function, it is created a file with the log of the used parameters of each hyperopt iteration.

        # Get the hyperopt trials
        opt_log = self.tpe_trials.idxs_vals[1]
        op_modes_list = opt_log["operation_mode"]
        dic_em_gain = {}
        dic_t_exp = {}
        # sorts the iterations by the index and uses it to get the EM gain and exposure time parameters
        # in the same order of the trials of the operation modes
        for item, count in collections.Counter(op_modes_list).items():
            # Get the current EM gain value by the hyperopt log
            dic_em_gain[str(item)] = opt_log["em_gain_" + str(item)]
            # Get the current exposure time value by the hyperopt log
            dic_t_exp[str(item)] = opt_log["t_exp_" + str(item)]

        dic = {}
        # opens a .txt file to record the hyperopt trials
        arq = open(path + file_base_name + "_LOG.txt", "w")
        # iterates the trials
        for i in range(len(op_modes_list)):
            # get the index of the trial
            item = op_modes_list[i]
            # get the current operation mode by the mirror list
            mode = self.new_list[item]
            # sets the dictionary EM mode as conventional
            dic["em_mode"] = "CONV"
            # if the EM mode is EM, sets the dictionary EM mode as EM
            if mode["em_mode"] == 1:
                dic["em_mode"] = "EM"
            # Pops the EM gain values obtained in the previous step
            dic["em_gain"] = int(dic_em_gain[str(item)].pop(0))
            # gets the HHS from the mirror list
            dic["hss"] = int(mode["hss"])
            if mode["hss"] < 1:
                dic["hss"] = float(mode["hss"])
            # gets the PREAMP from the mirror list
            dic["preamp"] = int(mode["preamp"])
            # gets the BINNING from the mirror list
            dic["bin"] = int(mode["binn"])
            # gets the SUB-IMAGE from the mirror list
            dic["sub_img"] = self.best_mode["sub_img"]
            # sets the kinetic series length as 1 (for posterior data treatment)
            dic["kinetic_series_length"] = 1
            # sets the observation type as object
            dic["obs_type"] = "object"
            # sets the star radius as zero
            dic["star_radius"] = 0
            # sets the coordinates of the object as ''
            dic["obj_coords"] = ""
            num = str(i)
            while len(num) < 3:
                num = "0" + num
            # sets the name of the used pre-image
            dic["img_name"] = file_base_name + "_O_" + num + ".fits"
            # Pops the exposure time values obtained in the previous step
            dic["t_exp"] = dic_t_exp[str(item)].pop(0)
            # sets the values of the cost function
            dic["output"] = -self.tpe_trials.results[i]["loss"]
            # write the hyperopt trial to the .txt file in the json format
            json.dump(dic, arq, sort_keys=True)
            arq.write("\n")
        arq.close()

    def create_bias_list(self, path, file_base_name):
        # This function creates a list of bias images with the same operation mode
        # used in the hyperopt trials. This list is created to be used
        # in an automatic image acquisition software.
        # This function works by reading the written modes by the creat_log_parameters() function,
        # and replacing the exposure time value by 0.00001 s.
        arq = open(path + file_base_name + "_LOG.txt", "r")
        lines = arq.read().splitlines()
        arq.close()
        arq = open(path + file_base_name + "_BIAS.txt", "w")
        for line in lines:
            line = json.loads(line)
            line["t_exp"] = 1e-5
            line["obs_type"] = "bias"
            line["img_name"] = line["img_name"].replace("_O_", "_B_")
            json.dump(line, arq)
            arq.write("\n")
        arq.close()

    def export_optimal_setup(
        self, img_directory, file_base_name, star_radius, obj_coords, snr_target
    ):
        # This function exports the optimal mode obtained by the hyperopt library to a .txt file.
        em_mode = self.best_mode["em_mode"]
        hss = self.best_mode["hss"]
        binn = self.best_mode["binn"]
        sub_img = self.best_mode["sub_img"]
        t_exp = self.best_mode["t_exp"]
        em_gain = self.best_mode["em_gain"]
        preamp = self.best_mode["preamp"]
        # Start the SNR calculation library
        SNRC = snrc.SignalToNoiseRatioCalc(
            t_exp,
            em_mode,
            em_gain,
            hss,
            preamp,
            binn,
            self.ccd_temp,
            self.sky_flux,
            self.star_flux,
            self.n_pix_star,
            self.serial_number,
        )
        SNRC.set_gain_value()
        SNRC.calc_RN()
        SNRC.calc_DC()
        SNRC.calc_SNR()
        # Calculates the SNR value for the respective optimum mode.
        # This step is needed because it is not possible recover
        # the maximum SNR used for the normalization
        snr = SNRC.get_SNR()
        # Starts the acquisition rate calculation class
        ARC = arc.AcquisitionRateCalc()
        ARC.write_operation_mode(em_mode, hss, binn, sub_img, t_exp)
        ARC.seleciona_t_corte()
        # Calculates the acquisition rate value.
        # This step is needed because it is not possible recover
        # the maximum acquisition rate used for the normalization
        ARC.calc_acquisition_rate()
        fa = float(ARC.return_acquisition_rate())

        dic = {}
        # Creates a dictionary with the optimum mode
        if self.best_mode["em_mode"] == 1:
            dic["em_mode"] = "EM"
        else:
            dic["em_mode"] = "CONV"
        dic["em_gain"] = int(self.best_mode["em_gain"])
        dic["t_exp"] = self.best_mode["t_exp"]
        dic["hss"] = self.best_mode["hss"]
        dic["preamp"] = self.best_mode["preamp"]
        dic["bin"] = self.best_mode["binn"]
        dic["sub_img"] = self.best_mode["sub_img"]
        dic["output"] = self.best_mode["SNR*FA"]
        dic["obs_type"] = "object"
        dic["img_name"] = file_base_name + "_OPTMODE.fits"
        dic["star_radius"] = star_radius
        try:
            dic["obj_coords"] = "(%i,%i)" % (obj_coords[0], obj_coords[1])
        except:
            1
        dic["kinetic_series_length"] = 1
        dic["max_snr"] = self.max_snr
        dic["min_snr"] = self.min_snr
        dic["max_fa"] = self.max_fa
        dic["min_fa"] = self.min_fa
        dic["snr_target"] = snr_target
        dic["SNR"] = snr
        dic["FA"] = fa
        file_name = img_directory + file_base_name + "_OPTSETUP.txt"
        # Write the optimum mode to a .txt file in the json format
        with open(file_name, "w") as arq:
            json.dump(dic, arq, indent=4, sort_keys=True)
            arq.close()
