#!/usr/bin/env python
# coding: utf-8
# 10/01/2020. Denis Varise Bernardes.

import json
import os
from copy import copy
from sys import exit
from webbrowser import Opera

from FWHM import FWHM
from hyperopt import rand, tpe
from hyperopt.pyll.base import repeat
from Opt_Acquisition_Rate import Optimize_Acquisition_Rate
from Opt_SNR import Optimize_SNR
from Opt_SNR_AR import Opt_SNR_AR


class Optimize_Camera:
    """Optimize Camera Class.

    Parameters
    ----------

    input_file_path: str, optional
        The path of the json file with the observation setup.

    Yields
    ------
        A json file with the optimal operation mode of the CCD.


    """

    _INPUT_FILE_NAME = "observation_setup.json"
    _SKY_FLUX = 12.298897076737294  # e-/pix/s
    _HATS24_FLUX = 56122.295000000006  # e-/s
    _N_PIX_STAR = 305
    _HATS24_MAG = 12.25

    def __init__(self, input_file_path=""):

        self._BIAS_LEVEL = 500
        self._STAR_RADIUS = 0
        self._STAR_COORDS = [0, 0]
        self.operation_modes = []
        self.input_file_path = input_file_path
        self._read_input_file()
        self._verify_ccd_parameters()
        self._create_list_operation_modes()

        if self.file_parameters["optimization_algorithm"] == "TPE":
            self.optimization_algorithm = tpe.suggest
        elif self.file_parameters["optimization_algorithm"] == "RAND":
            self.optimization_algorithm = rand.suggest
        else:
            raise ValueError(
                f"Unknow optimization algorithm: {self.file_parameters['optimization_algorithm']}"
            )

        if self.file_parameters["use_pre_img"]:
            self._calc_star_sky_flux_from_pre_img()
        else:
            self._calc_star_sky_flux_from_magnitude()

    def _read_input_file(self):
        try:
            file = os.path.join(self.input_file_path, self._INPUT_FILE_NAME)
            with open(file, "r") as file:
                self.file_parameters = json.load(file)
                file.close()
        except Exception:
            self._create_input_file()
            print("The configuration file was not found. A standard file was created.")
            exit()

    def _create_input_file(self):
        file_parameters = {
            "snr": 50,
            "acq_rate": 1,
            "mag": 12,
            "t_exp": [1e-5, 84600],
            "em_mode": ["EM", "Conv"],
            "em_gain": [2, 300],
            "readout_rate": [0.1, 1, 10, 20, 30],
            "preamp": [1, 2],
            "sub_img": [1024, 512, 256],
            "bin": [1, 2],
            "temperature": -70,
            "serial_number": 9917,
            "use_pre_img": False,
            "pre_img_name": "",
            "star_coords": [0, 0],
            "bias_img_name": "",
            "sky_radius": 20,
            "export_setup_file": True,
            "export_iterations_file": False,
            "export_bias_file": False,
            "exp_name": "",
            "iterations_number": 10,
            "optimization_algorithm": "TPE",
            "optimization_parameter": "SNR",
        }

        file_path = os.path.join(self.input_file_path, self._INPUT_FILE_NAME)
        with open(file_path, "w") as arq:
            json.dump(file_parameters, arq, indent=4)
            arq.close()

    def _verify_ccd_parameters(self):

        em_mode = self.file_parameters["em_mode"]
        em_gain = self.file_parameters["em_gain"]
        readout_rate = self.file_parameters["readout_rate"]
        preamp = self.file_parameters["preamp"]
        t_exp = self.file_parameters["t_exp"]
        temperature = self.file_parameters["temperature"]
        sub_img = self.file_parameters["sub_img"]
        binn = self.file_parameters["bin"]

        for value in em_mode:
            if value not in ["EM", "Conv"]:
                raise ValueError(f"Invalid value for the EM mode: {em_mode}")

        for value in em_gain:
            if value > 300 or value < 2:
                raise ValueError(f"The em gain out of range [2, 300]: {em_gain}")

        for value in preamp:
            if value not in [1, 2]:
                raise ValueError(f"Invalid value for the pre-amplification: {preamp}")

        for value in readout_rate:
            if value not in [0.1, 1, 10, 20, 30]:
                raise ValueError(f"Invalid value for the readout rate: {readout_rate}")
            if value == 0.1 and "Conv" not in em_mode:
                raise ValueError(
                    "The EM mode does not have the readout rate of 0.1 MHz."
                )
            if value in [10, 20, 30] and "EM" not in em_mode:
                raise ValueError(
                    f"The conventional mode does not have the readout rate of {value} MHz."
                )

        for value in binn:
            if value not in [1, 2]:
                raise ValueError(f"Invalid value for the binning: {bin}")

        for value in sub_img:
            if value not in [1024, 512, 256]:
                raise ValueError(f"Invalid value for the binning: {sub_img}")

        for value in t_exp:
            if value < 1e-5 or value > 84600:
                raise ValueError(
                    f"The exposure time out of range [1e-5, 84600]: {t_exp}"
                )

        if temperature < -80 or temperature > 15:
            raise ValueError(f"The temperature out of range [-80, 15]: {temperature}")

    def _create_list_operation_modes(self):
        em_modes = self.file_parameters["em_mode"]
        readout_rates = self.file_parameters["readout_rate"]
        preamps = self.file_parameters["preamp"]
        t_exps = self.file_parameters["t_exp"]
        sub_imgs = self.file_parameters["sub_img"]
        binns = self.file_parameters["bin"]

        for em_mode in em_modes:
            for readout_rate in readout_rates:
                if em_mode == "Conv" and readout_rate not in [0.1, 1]:
                    continue
                if em_mode == "EM" and readout_rate == 0.1:
                    continue
                em_gains = self.file_parameters["em_gain"]
                if em_mode == "Conv":
                    em_gains = [1]
                for binn in binns:
                    for sub_img in sub_imgs:
                        for preamp in preamps:
                            operation_mode = {
                                "em_mode": em_mode,
                                "em_gain": em_gains,
                                "readout_rate": readout_rate,
                                "preamp": preamp,
                                "bin": binn,
                                "sub_img": sub_img,
                                "t_exp": t_exps,
                            }

                            self.operation_modes.append(operation_mode)

    def _calc_star_sky_flux_from_pre_img(self):
        """Calculates the star flux based on the pre-image data."""
        file_path = os.path.join(
            self.input_file_path, self.file_parameters["pre_img_name"]
        )
        star_coords = self.file_parameters["star_coords"]
        sky_radius = self.file_parameters["sky_radius"]
        bias_img_path = os.path.join(
            self.input_file_path, self.file_parameters["bias_img_name"]
        )

        fwhm_obj = FWHM(
            img_name=file_path,
            xy_star=star_coords,
            sky_radius=sky_radius,
            bias_name=bias_img_path,
        )
        fwhm_obj.read_star_img()
        fwhm_obj.get_max_count()
        fwhm_obj.set_centroid()
        fwhm, star_radius, x, y = fwhm_obj.calc_FWHM()
        fwhm_obj.read_bias_img()
        fwhm_obj.calc_dark_current()
        fwhm_obj.read_exp_time()
        fwhm_obj.read_em_gain()
        fwhm_obj.calc_star_sky_flux()
        snr, rn, sky_flux, star_flux, n_pixels, bias_level = fwhm_obj.calc_SNR()
        self._SKY_FLUX = sky_flux
        self.star_flux = star_flux
        self._N_PIX_STAR = n_pixels
        self._STAR_COORDS = [x, y]
        self._STAR_RADIUS = star_radius
        self._BIAS_LEVEL = bias_level

    def _calc_star_sky_flux_from_magnitude(self):
        # Calculates, through the Pogson equation, the star flux based on the magnitude provided to the software.
        # This function uses the star flux and magnitude values of the HATS24 star as reference values.
        aux = 10 ** ((self._HATS24_MAG - self.file_parameters["mag"]) / 2.5)
        self.star_flux = self._HATS24_FLUX * aux

    def optimize(self):
        """Optimize the camera."""

        opt_param = self.file_parameters["optimization_parameter"]

        if opt_param == "SNR":
            self._optimize_snr()
        elif opt_param == "ACQ_RATE":
            self._optimize_ar()
        elif opt_param == "BOTH":
            self._optimize_snr_ar()
        else:
            raise ValueError(f"Unknow optimization parameter: {opt_param}.")

    def _optimize_snr(self):
        repeat = True

        while repeat == True:
            operation_modes = [copy(mode) for mode in self.operation_modes]
            opt_ar = Optimize_Acquisition_Rate(
                acquisition_rate=self.file_parameters["acq_rate"],
            )

            opt_ar.write_operation_modes(operation_modes)
            opt_ar.select_operation_modes()
            operation_modes = opt_ar.read_operation_modes()

            if operation_modes == []:
                self.file_parameters["acq_rate"] *= 0.8
                print("\nNo operation mode meets the provided requirements.")
                print(
                    f"The optimization will be repeated with 80 % of the acquisition rate: {self.file_parameters['acq_rate']:.2f} fps."
                )
                repeat = True
            else:
                repeat = False

        opt_snr = Optimize_SNR(
            snr_target=self.file_parameters["snr"],
            serial_number=self.file_parameters["serial_number"],
            temperature=self.file_parameters["temperature"],
            n_pix_star=self._N_PIX_STAR,
            sky_flux=self._SKY_FLUX,
            star_flux=self.star_flux,
            bias_level=self._BIAS_LEVEL,
        )
        opt_snr.write_operation_modes(operation_modes)
        best_snr = opt_snr.calc_best_mode()
        opt_snr.print_best_mode()
        if self.file_parameters["export_setup_file"]:
            opt_snr.export_optimal_setup(
                self.input_file_path,
                self.file_parameters["exp_name"],
                self._STAR_RADIUS,
                self._STAR_COORDS,
            )

    def _optimize_ar(self):
        repeat = True
        while repeat == True:

            opt_ar = Optimize_Acquisition_Rate(
                acquisition_rate=self.file_parameters["acq_rate"],
            )
            operation_modes = [copy(mode) for mode in self.operation_modes]
            opt_ar.write_operation_modes(operation_modes)
            opt_ar.select_operation_modes()
            operation_modes = opt_ar.read_operation_modes()

            if operation_modes == []:
                self.file_parameters["acq_rate"] *= 0.8
                print("\nNo operation mode meets the provided requirements.")
                print(
                    f"The optimization will be repeated with 80 % of the acquisition rate: {self.file_parameters['acq_rate']:.2f} fps."
                )
                repeat = True
            else:
                repeat = False
                self.operation_modes = operation_modes

        repeat = True
        while repeat == True:
            operation_modes = [copy(mode) for mode in self.operation_modes]
            opt_snr = Optimize_SNR(
                serial_number=self.file_parameters["serial_number"],
                snr_target=self.file_parameters["snr"],
                temperature=self.file_parameters["temperature"],
                n_pix_star=self._N_PIX_STAR,
                sky_flux=self._SKY_FLUX,
                star_flux=self.star_flux,
                bias_level=self._BIAS_LEVEL,
            )

            opt_snr.write_operation_modes(operation_modes)
            opt_snr.select_operation_modes_minimun_snr()
            operation_modes = opt_snr.read_operation_modes()

            if operation_modes == []:
                self.file_parameters["snr"] *= 0.8
                print("\nNo operation mode meets the provided requirements.")
                print(
                    f"The optimization will be repeated with 80 % of the SNR: {self.file_parameters['snr']:.2f}."
                )
                repeat = True
            else:
                repeat = False
                self.operation_modes = operation_modes

        opt_ar.write_operation_modes(self.operation_modes)
        best_mode = opt_ar.determine_fastest_operation_mode()
        opt_ar.print_best_mode()
        if self.file_parameters["export_setup_file"]:
            opt_ar.export_optimal_setup(
                self.input_file_path,
                self.file_parameters["exp_name"],
                self._STAR_RADIUS,
                self._STAR_COORDS,
                self.file_parameters["acq_rate"],
            )

    def _optimize_snr_ar(self):
        repeat = True
        while repeat == True:
            opt_ar = Optimize_Acquisition_Rate(
                acquisition_rate=self.file_parameters["acq_rate"],
            )
            operation_modes = [copy(mode) for mode in self.operation_modes]
            opt_ar.write_operation_modes(operation_modes)
            opt_ar.select_operation_modes()
            operation_modes = opt_ar.read_operation_modes()

            if operation_modes == []:
                self.file_parameters["acq_rate"] *= 0.9
                print("\nNo operation mode meets the provided requirements.")
                print(
                    f"The optimization will be repeated with 90 % of the acquisition rate: {self.file_parameters['acq_rate']:.2f} fps.\n "
                )
                repeat = True
            else:
                repeat = False
                self.operation_modes = operation_modes

        repeat = True
        while repeat == True:
            opt_snr = Optimize_SNR(
                serial_number=self.file_parameters["serial_number"],
                snr_target=self.file_parameters["snr"],
                temperature=self.file_parameters["temperature"],
                n_pix_star=self._N_PIX_STAR,
                sky_flux=self._SKY_FLUX,
                star_flux=self.star_flux,
                bias_level=self._BIAS_LEVEL,
            )
            operation_modes = [copy(mode) for mode in self.operation_modes]
            opt_snr.write_operation_modes(operation_modes)
            opt_snr.select_operation_modes_minimun_snr()
            operation_modes = opt_snr.read_operation_modes()

            if operation_modes == []:
                self.file_parameters["snr"] *= 0.9
                print("\nNo operation mode meets the provided requirements.")
                print(
                    f"The optimization will be repeated with 90 % of the SNR {self.file_parameters['snr']:.2f}.\n "
                )
                repeat = True
            else:
                repeat = False
                self.operation_modes = operation_modes
        # -------------------------------------------------------------------------

        opt_snr_ar = Opt_SNR_AR(
            snr_target=self.file_parameters["snr"],
            acq_rate=self.file_parameters["acq_rate"],
            serial_number=self.file_parameters["serial_number"],
            temperature=self.file_parameters["temperature"],
            n_pix_star=self._N_PIX_STAR,
            sky_flux=self._SKY_FLUX,
            star_flux=self.star_flux,
            bias_level=self._BIAS_LEVEL,
        )

        opt_snr_ar.write_operation_modes(self.operation_modes)
        opt_snr_ar.SNR_FA_ranges()
        opt_snr_ar.create_space()
        opt_snr_ar.run_bayesian_opt(
            max_evals=self.file_parameters["iterations_number"],
            algorithm=self.optimization_algorithm,
        )
        opt_snr_ar.print_best_mode()
        if self.file_parameters["export_iterations_file"]:
            opt_snr_ar.creat_log_parameters(
                self.input_file_path, self.file_parameters["exp_name"]
            )
        if self.file_parameters["export_setup_file"]:
            opt_snr_ar.export_optimal_setup(
                self.input_file_path,
                self.file_parameters["exp_name"],
                self._STAR_RADIUS,
                self._STAR_COORDS,
                self.file_parameters["snr"],
            )
