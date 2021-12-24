from SNR_Calculation import SNR_Calculation

mode = {
    "em_mode": "Conv",
    "em_gain": [1],
    "readout_rate": 1,
    "preamp": 1,
    "bin": 1,
    "sub_img": 1024,
    "t_exp": 1,
}

_SKY_FLUX = 12.298897076737294  # e-/pix/s
_HATS24_FLUX = 56122.295000000006  # e-/s
_N_PIX_STAR = 305
_HATS24_MAG = 12.25
mag = 12
aux = 10 ** ((_HATS24_MAG - mag) / 2.5)
star_flux = _HATS24_FLUX * aux


snr_calc = SNR_Calculation(mode, -70, _SKY_FLUX, star_flux, _N_PIX_STAR, 9917)
snr_calc.calc_DC()
snr_calc.calc_RN()
snr_calc.calc_SNR()
snr = snr_calc.get_SNR()
print(snr)
