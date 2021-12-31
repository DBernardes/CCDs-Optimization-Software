#!/usr/bin/env python
# coding: utf-8
# 24/10/2019. Denis Varise Bernardes.

import json
import os
from copy import copy
from sys import exit

from Acq_Rate_Calculation import Acquisition_Rate_Calculation
from Operation_Modes import Operation_Modes


class Optimize_Acquisition_Rate:
    """Optimiza Acquisition Rate class."""

    def __init__(self, acquisition_rate):
        self.arc = Acquisition_Rate_Calculation()
        self.acquisition_rate = acquisition_rate
        self.operation_modes = []

    def write_operation_modes(self, _list):
        """Write the list of operation modes."""

        self.operation_modes = _list

    def read_operation_modes(self):
        """Read the list of operation modes."""

        return self.operation_modes

    def select_operation_modes(self):
        """Select the operation modes.

        This functions determines those operation modes that
        accomplish the minimum acquisition rate requirement.
        """
        new_list = []
        for mode in self.operation_modes:
            min_t_exp = min(mode["t_exp"])
            max_t_exp = max(mode["t_exp"])
            mode["t_exp"] = min_t_exp
            self.arc.write_operation_mode(mode)
            self.arc.seleciona_t_corte()
            acq_rate = self.arc.calc_acquisition_rate()
            if acq_rate >= self.acquisition_rate:
                new_max_t_exp = self.arc.calculate_maximum_t_exp(self.acquisition_rate)
                if new_max_t_exp < max_t_exp:
                    max_t_exp = new_max_t_exp
                mode["t_exp"] = [min_t_exp, max_t_exp]
                new_list.append(mode)
        self.operation_modes = new_list

    def determine_min_acquisition_rate(self):
        min_acq_rate = 1e3
        for mode in self.operation_modes:
            new_mode = copy(mode)
            new_mode["t_exp"] = max(new_mode["t_exp"])
            self.arc.write_operation_mode(new_mode)
            self.arc.seleciona_t_corte()
            self.arc.calc_acquisition_rate()
            if self.arc.acquisition_rate < min_acq_rate:
                min_acq_rate = float(self.arc.return_acquisition_rate())
        return min_acq_rate

    def determine_fastest_operation_mode(self):
        max_acq_rate = 0
        best_modes = []
        for mode in self.operation_modes:
            new_mode = copy(mode)
            new_mode["t_exp"] = min(new_mode["t_exp"])
            self.arc.write_operation_mode(new_mode)
            self.arc.seleciona_t_corte()
            self.arc.calc_acquisition_rate()
            if self.arc.acquisition_rate > max_acq_rate:
                max_acq_rate = float(self.arc.return_acquisition_rate())
                best_modes = [new_mode]
            if self.arc.acquisition_rate == max_acq_rate:
                best_modes.append(new_mode)
            self.best_mode = self.find_largest_sub_img(best_modes)
            self.best_mode["max_acq_rate"] = max_acq_rate

        return self.best_mode, max_acq_rate

    def find_largest_sub_img(self, best_modes):
        largest_sub_img = 0
        best_mode = {}
        for mode in best_modes:
            if mode["sub_img"] > largest_sub_img:
                best_mode = mode
                largest_sub_img = mode["sub_img"]
        return best_mode

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
        print("\nBest Acquisition Rate: ", self.best_mode["max_acq_rate"])
        if self.best_mode["max_acq_rate"] < self.acquisition_rate:
            print("\nIt was not possible to reach the acquisition rate")

    def export_optimal_setup(
        self, img_directory, file_base_name, star_radius, obj_coords, FA
    ):
        self.best_mode["star_radius"] = star_radius
        self.best_mode["obj_coords"] = "(%i,%i)" % (obj_coords[0], obj_coords[1])
        self.best_mode["FA"] = FA
        if file_base_name != "":
            file_base_name += "_"
        file_name = os.path.join(img_directory, file_base_name + "OPTSETUP.txt")
        with open(file_name, "w") as arq:
            json.dump(self.best_mode, arq, indent=4, sort_keys=True)
            arq.close()
