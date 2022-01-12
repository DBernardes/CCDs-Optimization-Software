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

import numpy as np
from Acq_Rate_Calculation import Acquisition_Rate_Calculation
from hyperopt import Trials, fmin, hp, rand, tpe
from hyperopt.pyll import scope
from Operation_Modes import Operation_Modes
from Opt_Acquisition_Rate import Optimize_Acquisition_Rate
from Opt_SNR import Optimize_SNR
from SNR_Calculation import SNR_Calculation


def function_fa(parameters=[]):
    t_exp = parameters[0]
    em_mode = parameters[1]
    em_gain = parameters[2]
    readout_rate = parameters[3]
    preamp = parameters[4]
    binn = parameters[5]
    sub_img = parameters[6]
    temperature = parameters[7]
    sky_flux = parameters[8]
    star_flux = parameters[9]
    n_pixels_star = parameters[10]
    serial_number = parameters[11]
    max_fa = parameters[14]
    min_fa = parameters[15]

    mode = {
        "t_exp": t_exp,
        "em_mode": em_mode,
        "em_gain": em_gain,
        "readout_rate": readout_rate,
        "preamp": preamp,
        "bin": binn,
        "sub_img": sub_img,
    }
    ar_calc = Acquisition_Rate_Calculation()
    ar_calc.write_operation_mode(mode)
    ar_calc.seleciona_t_corte()
    ar_calc.calc_acquisition_rate()
    acq_rate = ar_calc.acquisition_rate
    norm_acq_rate = (acq_rate - min_fa) / (max_fa - min_fa)
    return norm_acq_rate


def function_snr(parameters=[]):
    t_exp = parameters[0]
    em_mode = parameters[1]
    em_gain = parameters[2]
    readout_rate = parameters[3]
    preamp = parameters[4]
    binn = parameters[5]
    sub_img = parameters[6]
    temperature = parameters[7]
    sky_flux = parameters[8]
    star_flux = parameters[9]
    n_pix_star = parameters[10]
    serial_number = parameters[11]
    max_snr = parameters[12]
    min_snr = parameters[13]
    snr_target = parameters[16]
    mode = {
        "t_exp": t_exp,
        "em_mode": em_mode,
        "em_gain": em_gain,
        "readout_rate": readout_rate,
        "preamp": preamp,
        "bin": binn,
        "sub_img": sub_img,
    }
    snr_calc = SNR_Calculation(
        mode,
        temperature=temperature,
        sky_flux=sky_flux,
        star_flux=star_flux,
        n_pix_star=n_pix_star,
        serial_number=serial_number,
    )
    snr_calc.calc_RN()
    snr_calc.calc_DC()
    snr_calc.calc_SNR()
    snr = snr_calc.get_SNR()
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


class Opt_SNR_AR:
    def __init__(
        self,
        snr_target,
        acq_rate,
        temperature,
        serial_number,
        n_pix_star,
        sky_flux,
        star_flux,
        bias_level,
    ):
        self.hss = [[], []]
        self.binn = [[], []]
        self.sub_img = []
        self.temperature = temperature
        self.serial_number = serial_number
        self.gain = 0
        self.dark_noise = 0
        self.set_dc()

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
        self.MOB = Operation_Modes()
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
        T = self.temperature
        if self.serial_number == 9914:
            self.dark_noise = 24.66 * np.exp(0.0015 * T ** 2 + 0.29 * T)
        if self.serial_number == 9915:
            self.dark_noise = 35.26 * np.exp(0.0019 * T ** 2 + 0.31 * T)
        if self.serial_number == 9916:
            self.dark_noise = 9.67 * np.exp(0.0012 * T ** 2 + 0.25 * T)
        if self.serial_number == 9917:
            self.dark_noise = 5.92 * np.exp(0.0005 * T ** 2 + 0.18 * T)

    def write_operation_modes(self, _list):
        # Write a object with the list of modes in the class
        self.operation_modes = _list

    def SNR_FA_ranges(self):
        opt_snr = Optimize_SNR(
            snr_target=self.snr_target,
            serial_number=self.serial_number,
            temperature=self.temperature,
            n_pix_star=self.n_pix_star,
            sky_flux=self.sky_flux,
            star_flux=self.star_flux,
            bias_level=self.bias_level,
        )

        opt_snr.write_operation_modes(self.operation_modes)
        max_snr = opt_snr.calc_best_mode()
        min_snr = opt_snr.calc_min_snr()

        opt_ar = Optimize_Acquisition_Rate(
            acquisition_rate=self.acq_rate_target,
        )
        opt_ar.write_operation_modes(self.operation_modes)
        min_fa = opt_ar.determine_min_acquisition_rate()
        best_mode, max_fa = opt_ar.determine_fastest_operation_mode()
        self.max_snr = max_snr
        self.min_snr = min_snr
        self.max_fa = max_fa
        self.min_fa = min_fa

    def calc_max_em_gain(self, max_t_exp, min_t_exp, max_em_gain):
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
        new_max_em_gain = (max_ADU - bias) / aux

        # if the photons/pix ir bigger than 100 photons, it is calculated the EM gain and the
        # exposure time values that accomplish this limit
        if aux > max_fotons:
            new_max_em_gain = (max_ADU - bias) / (max_fotons / gain)
            max_t_exp = max_fotons / (sky_flux + star_flux / n_pix + dn)
        if new_max_em_gain < max_em_gain:
            max_em_gain = new_max_em_gain
        # If the max exposure time found in the previous step is smaller than the min exposure time,
        # this mode is discarded
        if max_t_exp < min_t_exp:
            max_em_gain = 0
            max_t_exp = 0
        return max_em_gain, max_t_exp

    def create_space(self):
        i = 0
        self.new_list = []
        space_all_modes = []
        for mode in self.operation_modes:
            new_mode = copy(mode)
            max_em_gain = max(mode["em_gain"])
            min_em_gain = min(mode["em_gain"])
            max_t_exp = max(mode["t_exp"])
            min_t_exp = min(mode["t_exp"])
            t_exp = hp.uniform("t_exp_" + str(i), min_t_exp, max_t_exp)
            em_mode = hp.choice("em_mode_" + str(i), [mode["em_mode"]])
            em_gain = hp.choice("em_gain_" + str(i), [max_em_gain])
            if mode["em_mode"] == "EM":
                em_gain = hp.uniform("em_gain_" + str(i), min_em_gain, max_em_gain)
            readout_rate = hp.choice("hss_" + str(i), [mode["readout_rate"]])
            preamp = hp.choice("preamp_" + str(i), [mode["preamp"]])
            binn = hp.choice("binn_" + str(i), [mode["bin"]])
            sub_img = hp.choice("sub_img_" + str(i), [mode["sub_img"]])
            temperature = hp.choice("temperature_" + str(i), [self.temperature])
            sky_flux = hp.choice("sky_flux_" + str(i), [self.sky_flux])
            star_flux = hp.choice("obj_flux_" + str(i), [self.star_flux])
            n_pix_star = hp.choice("n_pix_star_" + str(i), [self.n_pix_star])
            serial_number = hp.choice("serial_number" + str(i), [self.serial_number])
            new_mode = [
                t_exp,
                em_mode,
                em_gain,
                readout_rate,
                preamp,
                binn,
                sub_img,
                temperature,
                sky_flux,
                star_flux,
                n_pix_star,
                serial_number,
            ]
            self.space_all_modes.append(new_mode)
            self.new_list.append(mode)
            i += 1
        self.space = hp.choice("operation_mode", self.space_all_modes)

    def run_bayesian_opt(self, max_evals, algorithm=tpe.suggest):
        gc.collect()
        tpe_algo = algorithm
        self.tpe_trials = Trials()
        max_snr = self.max_snr
        min_snr = self.min_snr
        max_fa = self.max_fa
        min_fa = self.min_fa
        best_mode = fmin(
            fn=function,
            space=self.space + [max_snr, min_snr, max_fa, min_fa, self.snr_target],
            algo=tpe_algo,
            trials=self.tpe_trials,
            max_evals=max_evals,
        )
        index_list_modes = best_mode["operation_mode"]
        chosen_mode = self.new_list[index_list_modes]
        t_exp = best_mode["t_exp_" + str(index_list_modes)]
        em_mode = chosen_mode["em_mode"]
        em_gain = best_mode["em_gain_" + str(index_list_modes)]
        readout_rate = chosen_mode["readout_rate"]
        preamp = chosen_mode["preamp"]
        binn = chosen_mode["bin"]
        sub_img = chosen_mode["sub_img"]
        self.best_mode = {
            "t_exp": t_exp,
            "em_mode": em_mode,
            "em_gain": em_gain,
            "readout_rate": readout_rate,
            "preamp": preamp,
            "bin": binn,
            "sub_img": sub_img,
            "SNR*FA": -self.tpe_trials.best_trial["result"]["loss"],
        }

    def print_best_mode(self):
        if self.best_mode["em_mode"] == "EM":
            print("\nEM Mode")
            print("-------")
            print("EM gain: ", self.best_mode["em_gain"])
        else:
            print("\nConventional Mode")
            print("-----------------")
        print("Exposure time (s): ", self.best_mode["t_exp"])
        print("Readout rate: ", self.best_mode["readout_rate"])
        print("Preamp: ", self.best_mode["preamp"])
        print("Binning: ", self.best_mode["bin"])
        print("Sub image: ", self.best_mode["sub_img"])
        print("\nBest SNR*FA: ", self.best_mode["SNR*FA"])

    def creat_log_parameters(self, path, file_base_name):
        opt_log = self.tpe_trials.idxs_vals[1]
        op_modes_list = opt_log["operation_mode"]
        dic_em_gain = {}
        dic_t_exp = {}
        for item, count in collections.Counter(op_modes_list).items():
            dic_em_gain[str(item)] = opt_log["em_gain_" + str(item)]
            dic_t_exp[str(item)] = opt_log["t_exp_" + str(item)]

        dic = {}
        arq = open(path + file_base_name + "_LOG.json", "w")
        for i in range(len(op_modes_list)):
            item = op_modes_list[i]
            mode = self.new_list[item]
            dic["em_gain"] = int(dic_em_gain[str(item)].pop(0))
            dic["readout_rate"] = int(mode["readout_rate"])
            if mode["readout_rate"] < 1:
                dic["readout_rate"] = float(mode["readout_rate"])
            dic["preamp"] = int(mode["preamp"])
            dic["bin"] = int(mode["bin"])
            dic["sub_img"] = self.best_mode["sub_img"]
            dic["kinetic_series_length"] = 1
            dic["obs_type"] = "object"
            dic["star_radius"] = 0
            dic["obj_coords"] = ""
            num = str(i)
            while len(num) < 3:
                num = "0" + num
            dic["img_name"] = file_base_name + "_O_" + num + ".fits"
            dic["t_exp"] = dic_t_exp[str(item)].pop(0)
            dic["output"] = -self.tpe_trials.results[i]["loss"]
            json.dump(dic, arq, sort_keys=True)
            arq.write("\n")
        arq.close()

    def create_bias_list(self, path, file_base_name):
        """This function creates a list of bias images with the same operation mode
        used in the hyperopt trials. This list is created to be used
        in an automatic image acquisition software.
        This function works by reading the written modes by the creat_log_parameters() function,
        and replacing the exposure time value by 0.00001 s"""

        arq = open(path + file_base_name + "_LOG.json", "r")
        lines = arq.read().splitlines()
        arq.close()
        arq = open(path + file_base_name + "_BIAS.json", "w")
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
        SNRC = SNR_Calculation(
            self.best_mode,
            self.temperature,
            self.sky_flux,
            self.star_flux,
            self.n_pix_star,
            self.serial_number,
        )
        SNRC.calc_SNR()
        snr = SNRC.get_SNR()
        ARC = Acquisition_Rate_Calculation()
        ARC.write_operation_mode(self.best_mode)
        ARC.seleciona_t_corte()
        ARC.calc_acquisition_rate()
        acq_rate = float(ARC.return_acquisition_rate())

        self.best_mode["star_radius"] = star_radius
        self.best_mode["obj_coords"] = [int(coord) for coord in obj_coords]
        self.best_mode["max_snr"] = self.max_snr
        self.best_mode["min_snr"] = self.min_snr
        self.best_mode["max_acq_rate"] = self.max_fa
        self.best_mode["min_acq_rate"] = self.min_fa
        self.best_mode["snr_target"] = snr_target
        self.best_mode["SNR"] = snr
        self.best_mode["ACQ_RATE"] = acq_rate

        suffix = "_OPTSETUP.json"
        if file_base_name == "":
            suffix = "OPTSETUP.json"
        file_name = os.path.join(img_directory, file_base_name + suffix)
        with open(file_name, "w") as arq:
            json.dump(self.best_mode, arq, indent=4, sort_keys=True)
            arq.close()
