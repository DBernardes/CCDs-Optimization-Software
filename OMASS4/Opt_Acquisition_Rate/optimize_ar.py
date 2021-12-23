#!/usr/bin/env python
# coding: utf-8
# 24/10/2019. Denis Varise Bernardes.

import json
import os
from sys import exit

from Acq_Rate_Calculation import Acquisition_Rate_Calculation
from Operation_Modes import Operation_Modes


class Optimize_Acquisition_Rate:
    def __init__(self, acquisition_rate, ccd_operation_mode):
        self.MOB = Operation_Modes()
        self.ARC = Acquisition_Rate_Calculation()
        self.acquisition_rate = acquisition_rate
        self.ccd_operation_mode = ccd_operation_mode

    def write_mode_to_MOB_class(
        self, em_mode, em_gain, hss, preamp, binn, sub_img, t_exp
    ):
        # Write the provided mode in the class
        self.MOB.write_mode(em_mode, em_gain, hss, preamp, binn, sub_img, t_exp)

    def print_MOB_list(self):
        # Prints the list of the available modes in the class
        lista = self.MOB.get_list_of_modes()
        for modo in lista:
            print(modo)

    def write_MOB_obj(self, obj):
        # Write a object with the list of operation modes
        self.MOB = obj

    def read_MOB_obj(self):
        # Reads the object with the list of operation modes
        return self.MOB

    def determine_operation_modes(self):
        # This functions determines those operation modes that
        # accomplish the minimum acquisition rate requirement.
        for em_mode in self.ccd_operation_mode["em_mode"]:
            for hss in self.ccd_operation_mode["readout_rate"]:
                if em_mode == "Conv" and hss not in [0.1, 1]:
                    continue
                if em_mode == "EM" and hss == 0.1:
                    continue
                for binn in self.ccd_operation_mode["bin"]:
                    for sub_img in self.ccd_operation_mode["sub_img"]:
                        # Starts the class for the acquisition rate calculation
                        self.ARC.write_operation_mode(
                            em_mode,
                            hss,
                            binn,
                            sub_img,
                            min(self.ccd_operation_mode["t_exp"]),
                        )
                        # Selects the cut-off exposure time
                        self.ARC.seleciona_t_corte()
                        # Calculates tha acquisition rate for the provided mode
                        self.ARC.calc_acquisition_rate()
                        # If the acquisition rate is greater than the provided limit,
                        # this mode is selected
                        if self.ARC.acquisition_rate >= self.acquisition_rate:
                            # Calculates the maximum exposure time that still acomplish the acquisition rate requirement
                            max_t_exp = (
                                self.ARC.calc_texp_provided_acquisition_frequency(
                                    self.acquisition_rate
                                )
                            )
                            if max_t_exp > max(self.ccd_operation_mode["t_exp"]):
                                max_t_exp = max(self.ccd_operation_mode["t_exp"])
                            self.write_mode_to_MOB_class(
                                em_mode, 0, hss, 0, binn, sub_img, max_t_exp
                            )

    def determine_min_acquisition_rate(self):
        # Given a list of modes,
        # this functions calculates the highest and the smallest acquisition rates
        min_acq_rate = 1e3
        # Iterates each mode of the list of selected modes
        for mode in self.MOB.get_list_of_modes():
            # Starts the class for the acquisition rate calculation
            self.ARC.write_operation_mode(
                mode["em_mode"],
                mode["hss"],
                mode["binn"],
                mode["sub_img"],
                mode["max_t_exp"],
            )
            # Selects the cut-off exposure time
            self.ARC.seleciona_t_corte()
            # Calculates tha acquisition rate for the provided mode
            self.ARC.calc_acquisition_rate()
            # If the acquisition rate is smaller than the current smallest acquisition rate,
            # this value is selected
            if self.ARC.acquisition_rate < min_acq_rate:
                min_acq_rate = float(self.ARC.return_acquisition_rate())
        return min_acq_rate

    def determine_fastest_operation_mode(self):
        # Given a list of modes, this functions determine the mode with the highest acquisition rate
        max_acq_rate = 0
        # Iterates each mode of the list of selected modes
        for mode in self.MOB.get_list_of_modes():
            # Starts the class for the acquisition rate calculation
            self.ARC.write_operation_mode(
                mode["em_mode"],
                mode["hss"],
                mode["binn"],
                mode["sub_img"],
                mode["min_t_exp"],
            )
            # Selects the cut-off exposure time
            self.ARC.seleciona_t_corte()
            # Calculates tha acquisition rate for the provided mode
            self.ARC.calc_acquisition_rate()
            # If the acquisition rate is greater than the current largest acquisition rate,
            # this mode is selected
            if self.ARC.acquisition_rate > max_acq_rate:
                max_acq_rate = float(self.ARC.return_acquisition_rate())
                self.best_mode = mode
                self.best_mode["max_acq_rate"] = max_acq_rate
            # These ifs are needed for those cases where there are several
            # sub-images options for the same optimum mode
            if self.ARC.acquisition_rate == max_acq_rate:
                if mode["hss"] == self.best_mode["hss"]:
                    if mode["em_mode"] == self.best_mode["em_mode"]:
                        if mode["binn"] == self.best_mode["binn"]:
                            if mode["sub_img"] > self.best_mode["sub_img"]:
                                self.best_mode["sub_img"] = mode["sub_img"]
        return self.best_mode, max_acq_rate

    def print_best_mode(self):
        if self.best_mode["em_mode"] == 1:
            print("\nEM Mode")
            print("-------")
            print("EM gain: ", self.best_mode["em_gain"])
        else:
            print("\nConventional Mode")
            print("-----------------")
        print("Exposure time (s): ", self.best_mode["min_t_exp"])
        print("HSS: ", self.best_mode["hss"])
        print("Preamp: ", self.best_mode["preamp"])
        print("Binning: ", self.best_mode["binn"])
        print("Sub image: ", self.best_mode["sub_img"])
        print("\nBest Acquisition Rate: ", self.best_mode["max_acq_rate"])
        if self.best_mode["max_acq_rate"] < self.acquisition_rate:
            print("\nIt was not possible to reach the acquisition rate")

    def export_optimal_setup(
        self, img_directory, file_base_name, star_radius, obj_coords, FA
    ):
        # Exports the best obtained operation mode to a .txt file
        # To accomplish this, it is used the json format
        dic = {}
        if self.best_mode["em_mode"] == 1:
            dic["em_mode"] = "EM"
        else:
            dic["em_mode"] = "CONV"
        dic["em_gain"] = self.best_mode["em_gain"]
        dic["t_exp"] = self.best_mode["min_t_exp"]
        dic["hss"] = self.best_mode["hss"]
        dic["preamp"] = self.best_mode["preamp"]
        dic["bin"] = self.best_mode["binn"]
        dic["sub_img"] = self.best_mode["sub_img"]
        dic["output"] = float(self.best_mode["max_acq_rate"])
        dic["obs_type"] = "object"
        dic["img_name"] = file_base_name + "_OPTMODE.fits"
        dic["star_radius"] = star_radius
        try:
            dic["obj_coords"] = "(%i,%i)" % (obj_coords[0], obj_coords[1])
        except:
            1
        dic["FA"] = FA

        file_name = os.path.join(img_directory, file_base_name + "_OPTSETUP.txt")
        with open(file_name, "w") as arq:
            json.dump(dic, arq, indent=4, sort_keys=True)
            arq.close()
