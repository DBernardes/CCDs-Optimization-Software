#!/usr/bin/env python
# coding: utf-8
# Denis Varise Bernardes.
# 08/10/2019.

import math
import os
from sys import exit

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d


class Acquisition_Rate_Calculation:
    """Acquisition Rate Calculation class."""

    def __init__(self):
        self.acquisition_rate = 0
        self.t_corte = 0
        self.tab_valores_readout = os.path.join(
            "OMASS4",
            "Acq_Rate_Calculation",
            "Spreadsheets",
            "Readout_Values.csv",
        )
        self.max_t_exp = 0

    def write_operation_mode(self, operation_mode):
        """Write the provided operation mode to the class."""
        self.t_exp = operation_mode["t_exp"]
        self.em_mode = operation_mode["em_mode"]
        self.readout_rate = operation_mode["readout_rate"]
        self.binn = operation_mode["bin"]
        self.sub_img = operation_mode["sub_img"]

    def create_spreadsheet_name(self):
        """Define the spreadshett name.

        This function mounts the string of the spreadsheet name the to calculate the acquisition
        rate of the CCD."""

        tab_name = "X" + str(self.sub_img) + "B" + str(self.binn) + "HSS"
        if self.readout_rate == 0.1:
            tab_name += "01"
        else:
            tab_name += str(self.readout_rate)
        return os.path.join(
            "OMASS4", "Acq_Rate_Calculation", "Spreadsheets", tab_name + ".csv"
        )

    def calc_acquisition_rate_texp_greater_tcorte(self):
        # If the exposure time ir greater than the critial redout time,
        # this function should be used
        self.acquisition_rate = 1 / self.t_exp

    def calc_acquisition_rate_texp_lower_tcorte(self):
        # If the exposure time ir smaller than the critial redout time,
        # this function should be used

        # Reads the spreadsheet
        tab_name = self.create_spreadsheet_name()
        columns = pd.read_csv(tab_name, dtype=np.float64)

        # Get the column with the exposure times
        t_exp_column = columns["TEXP (s)"]
        t_exp_column = self.filter_list(t_exp_column)
        # Get the column with the acquisition rate
        acq_rate_column = columns["FREQ (fps)"]
        acq_rate_column = self.filter_list(acq_rate_column)
        # Interpolates the data
        f = interp1d(t_exp_column, acq_rate_column)
        # Calculates the acquisition rate
        self.acquisition_rate = f(self.t_exp)

    def filter_list(self, lista):
        # this function receives a pandas spreadsheet column, and
        # removes annoying values
        for i in range(len(lista)):
            if math.isnan(lista[i]):
                lista = lista[0:i]
                break
            else:
                pass
        return lista

    def seleciona_t_corte(self):
        # This function selects the times that CCD spends to read the image.
        # These times were obtained in the work presented in
        # https://repositorio.unifei.edu.br/jspui/handle/123456789/2201
        indice_tab_t_corte = 0
        if self.readout_rate == 30:
            if self.sub_img == 1024:
                if self.binn == 2:
                    indice_tab_t_corte = 1
                if self.binn == 1:
                    indice_tab_t_corte = 2
            if self.sub_img == 512:
                if self.binn == 2:
                    indice_tab_t_corte = 3
                if self.binn == 1:
                    indice_tab_t_corte = 4
            if self.sub_img == 256:
                if self.binn == 2:
                    indice_tab_t_corte = 5
                if self.binn == 1:
                    indice_tab_t_corte = 6
        if self.readout_rate == 20:
            if self.sub_img == 1024:
                if self.binn == 2:
                    indice_tab_t_corte = 7
                if self.binn == 1:
                    indice_tab_t_corte = 8
            if self.sub_img == 512:
                if self.binn == 2:
                    indice_tab_t_corte = 9
                if self.binn == 1:
                    indice_tab_t_corte = 10
            if self.sub_img == 256:
                if self.binn == 2:
                    indice_tab_t_corte = 11
                if self.binn == 1:
                    indice_tab_t_corte = 12
        if self.readout_rate == 10:
            if self.sub_img == 1024:
                if self.binn == 2:
                    indice_tab_t_corte = 13
                if self.binn == 1:
                    indice_tab_t_corte = 14
            if self.sub_img == 512:
                if self.binn == 2:
                    indice_tab_t_corte = 15
                if self.binn == 1:
                    indice_tab_t_corte = 16
            if self.sub_img == 256:
                if self.binn == 2:
                    indice_tab_t_corte = 17
                if self.binn == 1:
                    indice_tab_t_corte = 18
        if self.readout_rate == 1:
            if self.sub_img == 1024:
                if self.binn == 2:
                    indice_tab_t_corte = 19
                if self.binn == 1:
                    indice_tab_t_corte = 20
            if self.sub_img == 512:
                if self.binn == 2:
                    indice_tab_t_corte = 21
                if self.binn == 1:
                    indice_tab_t_corte = 22
            if self.sub_img == 256:
                if self.binn == 2:
                    indice_tab_t_corte = 23
                if self.binn == 1:
                    indice_tab_t_corte = 24
        if self.readout_rate == 0.1:
            if self.sub_img == 1024:
                if self.binn == 2:
                    indice_tab_t_corte = 25
                if self.binn == 1:
                    indice_tab_t_corte = 26
            if self.sub_img == 512:
                if self.binn == 2:
                    indice_tab_t_corte = 27
                if self.binn == 1:
                    indice_tab_t_corte = 28
            if self.sub_img == 256:
                if self.binn == 2:
                    indice_tab_t_corte = 29
                if self.binn == 1:
                    indice_tab_t_corte = 30
        # Reads the spreadsheet
        columns = pd.read_csv(self.tab_valores_readout)
        # Get the readout values
        readout_times = [float(value) for value in columns["Tempo critico"][1:]]
        # Select the readout value based on the index obtained in previous step.
        self.t_corte = readout_times[indice_tab_t_corte - 1]

    def calc_acquisition_rate(self):
        # For t_exp < t_c: interpolates the spreadsheet data
        if self.t_exp < self.t_corte:
            self.calc_acquisition_rate_texp_lower_tcorte()
        # For t_exp > t_c: calculates 1/t_exp
        if self.t_exp >= self.t_corte:
            self.calc_acquisition_rate_texp_greater_tcorte()

        return self.acquisition_rate

    def calculate_maximum_t_exp(self, target_acquisition_rate):
        acq_rate = target_acquisition_rate
        freq_corte = 1 / self.t_corte
        if acq_rate <= freq_corte:
            self.max_t_exp = 1 / acq_rate
        if acq_rate > freq_corte:
            # Name of the spreadsheet with the texp vs acquisition rate curve
            tab_name = self.create_spreadsheet_name()
            # Reads the spreadsheet
            columns = pd.read_csv(tab_name, dtype=np.float64)
            # Get the exposure time column
            t_exp_column = columns["TEXP (s)"]
            t_exp_column = self.filter_list(t_exp_column)
            # Get the acquisition rate column
            acq_rate_column = columns["FREQ (fps)"]
            acq_rate_column = self.filter_list(acq_rate_column)
            # Calculates the exposure time given the acquisition rate
            f = interp1d(acq_rate_column, t_exp_column)
            # Calculates the exposure time
            self.max_t_exp = float(f(acq_rate))
        return self.max_t_exp

    def return_acquisition_rate(self):
        return self.acquisition_rate
