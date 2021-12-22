#!/usr/bin/env python
# coding: utf-8
# Denis Varise Bernardes.
# 08/10/2019.

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d


class ReadNoiseCalc:
    def __init__(self):
        self.noise = 0

    def write_operation_mode(self, em_mode, em_gain, hss, preamp, binn):
        # Write the operation mode to the class
        self.em_mode = em_mode
        self.em_gain = em_gain
        self.hss = hss
        self.preamp = preamp
        self.binn = binn

    def get_operation_mode(self):
        print("em_mode = ", self.em_mode)
        print("em_gain = ", self.em_gain)
        print("hss = ", self.hss)
        print("preamp = ", self.preamp)
        print("binn = ", self.binn)

    def calc_read_noise(self):
        # For the conventional mode, it is returned the table's read noise  value
        if self.em_mode == 0:
            self.calc_read_noise_conventional()
        # Para o modo EM, calcula em funcao do ganho_em
        if self.em_mode == 1:
            self.calc_read_noise_em()

    def calc_read_noise_conventional(self):
        # Based on the operation mode, this function selects the
        # respective read noise value
        indice_tab = 0
        if self.hss == 1:
            if self.preamp == 1:
                if self.binn == 1:
                    indice_tab = 17
                if self.binn == 2:
                    indice_tab = 18
            if self.preamp == 2:
                if self.binn == 1:
                    indice_tab = 19
                if self.binn == 2:
                    indice_tab = 20
        if self.hss == 0.1:
            if self.preamp == 1:
                if self.binn == 1:
                    indice_tab = 21
                if self.binn == 2:
                    indice_tab = 22
            if self.preamp == 2:
                if self.binn == 1:
                    indice_tab = 23
                if self.binn == 2:
                    indice_tab = 24
        # Reads the spreadsheet
        df = pd.read_excel("Spreadsheets/Read_Noise_Values.xlsm")
        columns = pd.DataFrame(df)
        column_noise = columns["Noise"]
        self.noise = column_noise[indice_tab]

    def calc_read_noise_em(self):
        # This funciton calculates the read noise for a given value of the EM gain
        hss = self.hss
        # Reads the spreadsheet
        tab_name = (
            "Spreadsheets/RN_PA"
            + str(int(self.preamp))
            + "B"
            + str(int(self.binn))
            + "HSS"
            + str(int(hss))
            + ".xlsm"
        )
        df = pd.read_excel(tab_name)
        columns = pd.DataFrame(df)
        column_noise = columns["Noise (e-)"][0:11]
        column_em_gain = columns["EM Gain"][0:11]
        # It is done an interpolation for the read noise and EM gain values found in the spreadsheet.
        f = interp1d(column_em_gain, column_noise)
        # Calculates the read noise
        self.noise = f(self.em_gain)
