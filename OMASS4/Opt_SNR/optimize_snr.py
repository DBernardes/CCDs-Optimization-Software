#!/usr/bin/env python
# coding: utf-8
# Denis Varise Bernardes.
# 12/12/2019.


import json
import os
from copy import copy
from sys import exit

import numpy as np
from hyperopt.pyll import scope
from hyperopt.pyll.stochastic import sample
from Operation_Modes import Operation_Modes
from SNR_Calculation import SNR_Calculation


class Optimize_SNR:
    """Optimize the Signal to Noise Ratio Class."""

    _MAX_FOTONS = 100
    _MAX_ADU = (2 ** 16) * 0.8

    def __init__(
        self,
        snr_target,
        serial_number,
        temperature,
        n_pix_star,
        sky_flux,
        star_flux,
        bias_level,
    ):

        self.temperature = temperature
        self.serial_number = serial_number
        self.snr_target = snr_target
        self.sky_flux = sky_flux  # e-/pix/s
        self.star_flux = star_flux  # e-/s
        self.n_pix_star = n_pix_star
        self.bias_level = bias_level
        self._set_dc()

    def write_operation_modes(self, _list):
        """Write the operation modes list into the class."""
        self.operation_modes = _list
        self._calc_max_em_gain()

    def read_operation_modes(self):
        """Read the list of operation modes."""
        return self.operation_modes

    def _set_gain(self, mode):
        em_mode = mode["em_mode"]
        readout_rate = mode["readout_rate"]
        preamp = mode["preamp"]
        gain = 0
        if em_mode == "EM":
            if readout_rate == 30:
                if preamp == 1:
                    gain = 17.2
                if preamp == 2:
                    gain = 5.27
            if readout_rate == 20:
                if preamp == 1:
                    gain = 16.4
                if preamp == 2:
                    gain = 4.39
            if readout_rate == 10:
                if preamp == 1:
                    gain = 16.0
                if preamp == 2:
                    gain = 3.96
            if readout_rate == 1:
                if preamp == 1:
                    gain = 15.9
                if preamp == 2:
                    gain = 3.88
        else:
            if readout_rate == 1:
                if preamp == 1:
                    gain = 3.37
                if preamp == 2:
                    gain = 0.8
            if readout_rate == 0.1:
                if preamp == 1:
                    gain = 3.35
                if preamp == 2:
                    gain = 0.8
        self.gain = gain

    def _set_dc(self):
        T = self.temperature
        if self.serial_number == 9914:
            self.dark_noise = 24.66 * np.exp(0.0015 * T ** 2 + 0.29 * T)
        if self.serial_number == 9915:
            self.dark_noise = 35.26 * np.exp(0.0019 * T ** 2 + 0.31 * T)
        if self.serial_number == 9916:
            self.dark_noise = 9.67 * np.exp(0.0012 * T ** 2 + 0.25 * T)
        if self.serial_number == 9917:
            self.dark_noise = 5.92 * np.exp(0.0005 * T ** 2 + 0.18 * T)

    def _calc_max_em_gain(self):
        new_list = []
        for mode in self.operation_modes:
            if mode["em_mode"] == "EM":
                max_t_exp = max(mode["t_exp"])
                min_t_exp = min(mode["t_exp"])
                max_em_gain = max(mode["em_gain"])
                min_em_gain = min(mode["em_gain"])
                self._set_gain(mode)
                total_fotons = (
                    self.sky_flux + self.star_flux / self.n_pix_star + self.dark_noise
                ) * max_t_exp
                if total_fotons > self._MAX_FOTONS:
                    total_fotons = self._MAX_FOTONS
                    max_t_exp = self._MAX_FOTONS / (
                        self.sky_flux
                        + self.star_flux / self.n_pix_star
                        + self.dark_noise
                    )
                em_gain = (self._MAX_ADU - self.bias_level) / (total_fotons / self.gain)
                if em_gain < max_em_gain:
                    max_em_gain = em_gain
                if max_t_exp < min_t_exp:
                    continue
                mode["t_exp"] = [min_t_exp, max_t_exp]
                mode["em_gain"] = [min_em_gain, max_em_gain]
            new_list.append(mode)
        self.operation_modes = new_list

    def select_operation_modes_minimun_snr(self):
        """Select the operation modes that accomplish the SNR value.

        This function determines the operating modes that meet a minimum SNR.
        For each mode, the maximum allowed EM gain is calculated.
        This gain is used to calculate the minimum exposure time allowed to reach the SNR.
        The selected modes are passed to the MOB object mode list
        """
        new_list = []
        for mode in self.operation_modes:
            min_t_exp = min(mode["t_exp"])
            max_t_exp = max(mode["t_exp"])
            mode["t_exp"] = min_t_exp
            min_em_gain = min(mode["em_gain"])
            max_em_gain = max(mode["em_gain"])
            mode["em_gain"] = min_em_gain

            snr_calc = SNR_Calculation(
                mode,
                self.temperature,
                self.sky_flux,
                self.star_flux,
                self.n_pix_star,
                self.serial_number,
            )
            t_exp = snr_calc.calc_minimun_texp_provided_snr(self.snr_target)
            if t_exp > min_t_exp:
                min_t_exp = t_exp
            if min_t_exp <= max_t_exp:
                mode["t_exp"] = [min_t_exp, max_t_exp]
                mode["em_gain"] = [min_em_gain, max_em_gain]
                if mode["em_mode"] == "Conv":
                    mode["em_gain"] = [1]
                new_list.append(mode)
        self.operation_modes = new_list

    def calc_min_snr(self):
        min_snr = 1e5
        best_mode = {}
        for mode in self.operation_modes:
            new_mode = copy(mode)
            new_mode["t_exp"] = max(new_mode["t_exp"])
            new_mode["em_gain"] = max(new_mode["em_gain"])
            snr_calc = SNR_Calculation(
                new_mode,
                self.temperature,
                self.sky_flux,
                self.star_flux,
                self.n_pix_star,
                self.serial_number,
            )
            snr_calc.calc_SNR()
            snr = snr_calc.get_SNR()
            if snr < min_snr:
                min_snr = snr
        self.best_mode = best_mode
        return min_snr

    def calc_best_mode(self):
        best_snr = 0
        best_modes = []
        for mode in self.operation_modes:
            new_mode = copy(mode)
            new_mode["t_exp"] = max(new_mode["t_exp"])
            new_mode["em_gain"] = max(new_mode["em_gain"])
            snr_calc = SNR_Calculation(
                new_mode,
                self.temperature,
                self.sky_flux,
                self.star_flux,
                self.n_pix_star,
                self.serial_number,
            )
            snr_calc.calc_SNR()
            snr = snr_calc.get_SNR()
            if snr > best_snr:
                best_snr = snr
                best_modes = [new_mode]
            elif snr == best_snr:
                best_modes.append(new_mode)
        self.best_modes = best_modes
        self.find_largest_sub_img()
        self.best_snr = best_snr

        return best_snr

    def find_largest_sub_img(self):
        best_mode = {}
        largest_sub_img = 0
        for mode in self.best_modes:
            if mode["sub_img"] > largest_sub_img:
                best_mode = mode
        self.best_mode = mode

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
        print("\nBest SNR: ", self.best_snr)
        if self.snr_target > self.best_snr:
            print("\nThe provided SNR could not be reach.")

    def export_optimal_setup(
        self, img_directory, file_base_name, star_radius, obj_coords
    ):
        # This functions exports the obtainde best mode to a .txt file
        # To acocmplish this, it is used the json format
        self.best_mode["best_snr"] = self.best_snr
        self.best_mode["star_radius"] = star_radius

        self.best_mode["obj_coords"] = "(%i,%i)" % (obj_coords[0], obj_coords[1])

        if file_base_name != "":
            file_base_name += "_"

        file_name = os.path.join(img_directory, file_base_name + "OPTSETUP.json")
        with open(file_name, "w") as arq:
            json.dump(self.best_mode, arq, indent=4, sort_keys=True)
            arq.close()
