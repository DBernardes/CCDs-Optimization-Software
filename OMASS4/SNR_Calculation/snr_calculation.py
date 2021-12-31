#!/usr/bin/env python
# coding: utf-8
# Denis Varise Bernardes.
# 26/11/2019.

from math import exp, sqrt
from sys import exit

import numpy as np
from Read_Noise_Calculation import Read_Noise_Calculation


class SNR_Calculation:
    def __init__(
        self,
        mode,
        temperature,
        sky_flux,
        star_flux,
        n_pix_star,
        serial_number,
    ):
        self.t_exp = mode["t_exp"]
        self.em_mode = mode["em_mode"]
        self.readout_rate = mode["readout_rate"]
        self.preamp = mode["preamp"]
        self.binn = mode["bin"]
        self.temperature = temperature
        self.sky_flux = sky_flux
        self.star_flux = star_flux
        self.n_pix_star = n_pix_star
        self.serial_number = serial_number

        self.dark_noise = 0
        self.SNR = 0
        self.gain = 0
        self.read_noise = 0
        self.noise_factor = 1
        self.em_gain = 1
        if self.em_mode == "EM":
            self.noise_factor = 1.41
            self.em_gain = mode["em_gain"]

        self.set_gain_value()
        self.calc_RN()
        self.calc_DC()

    def calc_SNR(self):

        # This function calculates the SNR of the star
        t_exp = self.t_exp
        em_gain = self.em_gain
        n_pix = self.n_pix_star
        binn = self.binn
        dc = self.dark_noise
        rn = self.read_noise
        star = self.star_flux
        sky = self.sky_flux
        nf = self.noise_factor

        aux = np.sqrt(
            star * t_exp * nf ** 2
            + n_pix * ((rn / em_gain / binn) ** 2 + (sky + dc) * t_exp * nf ** 2)
        )
        self.SNR = (star * t_exp) / aux

    def calc_minimun_texp_provided_snr(self, snr):
        # This function calculates the minimum exposure time needed
        # to achieve a provided value of SNR
        em_gain = self.em_gain
        n_pix = self.n_pix_star
        dc = self.dark_noise
        rn = self.read_noise
        star = self.star_flux
        sky = self.sky_flux
        nf = self.noise_factor
        binn = self.binn

        # equation: a * texp^2 - b* texp - c
        # This equation was obtained by isolating the exposure time in the SNR equation
        a = star ** 2
        b = snr ** 2 * nf ** 2 * (star + n_pix * (sky + dc))
        c = snr ** 2 * n_pix * (rn / em_gain / binn) ** 2
        # Calculates the roots of the equation
        minimun_t_exps = self.calc_quadratic_equation(a, -b, -c)
        # Then, it is selected the minimum non negative exposure time
        for t_exp in minimun_t_exps:
            if t_exp <= 0:
                minimun_t_exps.remove(t_exp)
        return min(minimun_t_exps)

    def calc_quadratic_equation(self, a, b, c):
        # Calculates the roots of a quadratic equaiton
        delta = (b ** 2) - (4 * a * c)
        x1 = (-b - sqrt(delta)) / (2 * a)
        x2 = (-b + sqrt(delta)) / (2 * a)
        return [x1, x2]

    def calc_DC(self):
        # Calculates the dark current based ib the CCD temperature.
        # These used equations model the dark current of the SPARC4 CCDs,
        # and can be found in: D V Bernardes et al 2018 PASP 130 095002
        T = self.temperature
        if self.serial_number == 9914:
            self.dark_noise = 24.66 * exp(0.0015 * T ** 2 + 0.29 * T)
        if self.serial_number == 9915:
            self.dark_noise = 35.26 * exp(0.0019 * T ** 2 + 0.31 * T)
        if self.serial_number == 9916:
            self.dark_noise = 9.67 * exp(0.0012 * T ** 2 + 0.25 * T)
        if self.serial_number == 9917:
            self.dark_noise = 5.92 * exp(0.0005 * T ** 2 + 0.18 * T)

    def calc_RN(self):
        # Calculates the read noise using the read noise library
        RN = Read_Noise_Calculation()
        RN.write_operation_mode(
            self.em_mode, self.em_gain, self.readout_rate, self.preamp, self.binn
        )
        RN.calc_read_noise()
        self.read_noise = float(RN.noise)

    def set_gain_value(self):
        # Sets the CCD gain
        em_mode = self.em_mode
        readoout_rate = self.readout_rate
        preamp = self.preamp
        gain = 0
        if em_mode == "EM":
            if readoout_rate == 30:
                if preamp == 1:
                    gain = 17.2
                if preamp == 2:
                    gain = 5.27
            if readoout_rate == 20:
                if preamp == 1:
                    gain = 16.4
                if preamp == 2:
                    gain = 4.39
            if readoout_rate == 10:
                if preamp == 1:
                    gain = 16.0
                if preamp == 2:
                    gain = 3.96
            if readoout_rate == 1:
                if preamp == 1:
                    gain = 15.9
                if preamp == 2:
                    gain = 3.88
        else:
            if readoout_rate == 1:
                if preamp == 1:
                    gain = 3.37
                if preamp == 2:
                    gain = 0.8
            if readoout_rate == 0.1:
                if preamp == 1:
                    gain = 3.35
                if preamp == 2:
                    gain = 0.8
        self.gain = gain

    def get_SNR(self):
        return self.SNR

    def print_noise_info(self):
        print("\n------ Noise Info ------ ")
        print("Star Flux: ", self.star_flux, "photons")
        print("Sky Flux: ", self.sky_flux, "photons/pix")
        print("Read Noise: ", self.read_noise, "e-/pix")
        print("Dark Noise: ", self.dark_noise, "e-/pix")
