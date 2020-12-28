# Software of the Optimization Method of the Electron Multiplying Charge Coupled Device of the Acquisition System of the SPARC4

## Introduction

The Astrophysics Division of the *Instituto Nacional de Pesquisas Espaciais* (INPE) in collaboration with the *Laboratório Nacional de Astrofísica* (LNA) is developing a new astronomical instrument, the Simultaneous Polarimeter and Rapid Camera in Four Bands ([SPARC4](https://www.spiedigitallibrary.org/conference-proceedings-of-spie/8446/844626/Concept-of-SPARC4--a-simultaneous-polarimeter-and-rapid-camera/10.1117/12.924976.full?casa_token=7b-hbhyqIMoAAAAA%3a99lzc7LW-gGeFuEs1N_7ZGdcFS1EiapC3jbzEYyrWT3PDiUP4RXPDEiR9IdfuRvDY7pPetsPx88&SSO=1)). SPARC4 will be installed on the 1.6 m Perkin-Elmer telescope at Observatório Pico dos Dias (OPD), Brazil, and it will allow image acquisition in the four Sloan Digital Sky Survey (SDSS) photometric bands: g, r, i and z. For the acquisition in each band (channel), there is a dedicated iXon Ultra EMCCD, produced by Andor Technology. These devices have an optical window and coating optimized for the spectral range in which they were designed to operate. These cameras also have frame transfer and electron-multiplying capabilities, allowing acquisition rates (AR) of up to 26 fps full-frame (1024 x 1024 pixels) even on faint astronomical objects, which requires high sensitivity for short exposure times.   

The quality of photometric measurements in astronomical observations can be quantified by the signal-to-noise ratio (SNR) and the acquisition rate (AR) provided by these camewras. Either the SNR or AR or both can change depending on the configuration of the operational mode of the CCD. Therefore, an optimal selection of the operational mode for each CCD is important to obtain the best performance of the instrument. 

To solve this problem, we present the Optimization Method of the Electron Multiplying Charge Coupled Device (EMCCDs) of the Acquisition System of the SPARC4 (OMASS4). The OMASS4 uses as figures of merit the signal-to-noise ratio (SNR) and the acquisition rate (AR) as a function of the operation mode of the CCDs. Three different modes of optimization are included in the OMASS4:  (1) optimization of SNR only; (2) optimization of AR only; and (3) optimization of both SNR and AR simultaneously. The first two modes calculate an analytical minimization of the cost function whereas the third mode uses the bayesian optimization method (BOM) to determine the optimum mode of operation. We apply the OMASS4 to find the optimum mode for observations obtained at the Pico dos Dias Observatory, Brazil, and compare the delivered modes of operation and its performance with the ones adopted by the observer. If the OMASS4 had been used as a tool to optimize the CCDs in all of these nights, it would be possible to improve their efficiency in 97.17 %, 65.08 %, and 77.66 % for the optimization modes 1, 2, and 3, respectively. For more information, access: [Otimização dos modos de operação do sistema de aquisição do instrumento SPARC4](https://repositorio.unifei.edu.br/jspui/handle/123456789/2201). This repo presents the software developed to implement the OMASS4.

## Software Description

The software to implement the OMASS4 was developed using Python Language 3.7.4 to determine the optimum operation mode of the SPARC4 CCDs. It is structured into three parts: the initialization, the star flux calculation, and the CCD optimization. For the initialization step, it requires to provide to the software all the information related to the astronomical object, i.e.: an image of the object, its (x,y) coordinates, the maximum star radius, a bias image with the respective used CCD operation mode, the SNR, the AR, allowed SI and Bin modes, CCD temperature, and the iterations number of the BOM. Then, the star flux is calculated for the optimal star radius given by the full width at half maximum (FWHM) parameter. The OMASS4 uses a set of packages to calculate the SNR and the AR values according to the CCD operation mode. The code developed to execute the BOM is based on the library provided by [Koehrsen](https://github.com/WillKoehrsen/hyperparameter-optimization). The used algorithm to model the objective function of the BOM is the tree structured Parzen estimator (TPE). The SNR package operation is based on the methodology presented by [Bernardes, D. V. (2020)](https://repositorio.unifei.edu.br/jspui/handle/123456789/2201). The star flux, sky flux, and the number of star pixels are obtained thorugh a pre-image of the object. The DC noise is calculated according to the model presented by [Bernardes et. al (2018)](https://iopscience.iop.org/article/10.1088/1538-3873/aacb1e/meta?casa_token=QzaY5kK_Yp8AAAAA:Qz_wlI6tq2WMi4sRF-tLvw-S2RwkmkF1_N8i7mReLYSUgim4dqp3yceqyLmlbrgUHt5TTzGYcrnYW_9ttxnfrw) for the four SPARC4 cameras. The read noise is obtained through the characterization presented in [Bernardes, D. V. (2020)](https://repositorio.unifei.edu.br/jspui/handle/123456789/2201). The G value is obtained through the camera datasheet.

The performance of the EM mode is better than the conventional mode until a max value of 100 photons per pixel. So, the t<sub>exp</sub> of each EM mode is limited to accomplish this requirement. Also, the maximum value allowed for the amplification EM gain G<sub>em</sub> is 300x. Values larger than 300x would deteriorate the device. Furthermore, the G<sub>em</sub> must be such that the CCD will not saturate. For this reason, the maximum EM gain allowed was arbitrarily configured to provide a signal up to 80 \% of the pixel well depth. For an image with 16 bits per pixel, this value is 2<sup>16</sup> x 0.8 = 52429 analogical to digital unit (ADU). Given that a pixel value is composed by the star, sky and, dark current signals, and the bias level, the maximum value for the G<sub>em</sub> is

<p align="center">
<img src="https://latex.codecogs.com/svg.latex?G_{\rm&space;em}&space;=&space;\frac{(52429&space;-&space;B)&space;\times&space;G}{(S/n_{\rm&space;p}&space;&plus;&space;S_{\rm&space;sky}&space;&plus;&space;S_{\rm&space;dc})}." title="G_{\rm em} = \frac{(52429 - B) \times G}{(S/n_{\rm p} + S_{\rm sky} + S_{\rm dc})}." />
</p>

where S, S<sub>dc</sub>, and S<sub>sky</sub> represent the photons number of the star, the mean thermoelectrons, and the mean of the photons number of the sky for the acquired image, respectively; n<sub>p</sub> is the number of pixels considered to calculate the S value, G is the gain of the CCD in e-/ADU, and B is the bias level in ADU. 

The AR package operation is based on the characterization presented of the CCDs presented by [Bernardes, D. V. (2020)](https://repositorio.unifei.edu.br/jspui/handle/123456789/2201). For each mode, the AR value will be calculated through interpolation if t<sub>exp</sub> < t<sub>c</sub>, and it is the inverse of the t<sub>exp</sub>, if t<sub>exp</sub> equals or greater than t<sub>c</sub>, where the t<sub>c</sub> is the time spent by the camera to read one image, as a function of the CCD operation mode.

Therefore, the OMASS4 was implemented using the aforementioned packages, being applied to three different optimization modes: optimize SNR (mode 1), optimize AR (mode 2), and optimize both SNR and AR (mode 3). 

* Mode 1: in this mode, the SNR is optimized, keeping the AR fixed. First, it is selected those modes that accomplish the AR requirement. Then, it is calculated the SNR value for each selected mode, using the maximum values for the t<sub>exp</sub> and G<sub>em</sub>. The optimum mode is given by that one with the highest SNR.
    
* Mode 2: in this mode, the AR is optimized, keeping the SNR fixed. Initially, for each mode, it is calculated the minimum t<sub>exp</sub> value that accomplish the SNR requirement, for the maximum G<sub>em</sub> allowed. For this calculation, it is considered the values of the star flux s = S/t<sub>exp</sub> in photons/s, the sky flux s<sub>sky</sub> = S<sub>sky</sub>/t<sub>exp</sub>, in photons/pixel/s, and the dark current s<sub>dc</sub> = S<sub>dc</sub>/t<sub>exp</sub>, in e-/pixel/s . So, the equation for the SNR of the star can be written as follows


<p align="center">
<img src="https://latex.codecogs.com/svg.latex?S&space;=&space;\frac{s&space;\times&space;t_{\rm&space;exp}}{\{&space;s&space;\;&space;t_{\rm&space;exp}&space;\;&space;N_{\rm&space;F}^2&space;&plus;&space;\\&space;n_{\rm&space;p}&space;[\&space;(s_{\rm&space;sky}&space;&plus;&space;s_{\rm&space;dc})&space;\;&space;t_{\rm&space;exp}&space;\;&space;N_{\rm&space;F}^2&space;&plus;&space;\\&space;(\sigma_{\rm&space;ADU}&space;\;&space;G/G_{\rm&space;em})^2&space;]\&space;\}^{1/2}}." title="S = \frac{s \times t_{\rm exp}}{\{ s \; t_{\rm exp} \; N_{\rm F}^2 + \\ n_{\rm p} [\ (s_{\rm sky} + s_{\rm dc}) \; t_{\rm exp} \; N_{\rm F}^2 + \\ (\sigma_{\rm ADU} \; G/G_{\rm em})^2 ]\ \}^{1/2}}." />
</p>
    
    
where &sigma;<sub>ADU</sub> represent the counts' distribution of the acquired image. N<sub>F</sub> is the noise factor and represents and extra noise added to the image because of the use of the EM amplifier. For an Andor EMCCD, N<sub>F</sub> = 1.41. Rearranging the terms of the equation above and isolating t<sub>exp</sub>,
    
<p align="center">
<img src="https://latex.codecogs.com/svg.latex?s^2&space;\;&space;t_{\rm&space;exp}^2&space;-&space;S^2&space;\;&space;N_{\rm&space;F}^2&space;\;&space;[\&space;s&space;&plus;&space;n_{\rm&space;p}&space;(s_{\rm&space;sky}&space;&plus;&space;s_{\rm&space;dc})&space;]\&space;\;&space;t_{\rm&space;exp}&space;-&space;S^2&space;\;&space;n_{\rm&space;p}&space;\;&space;\sigma_{\rm&space;r}^2&space;=&space;0" title="s^2 \; t_{\rm exp}^2 - S^2 \; N_{\rm F}^2 \; [\ s + n_{\rm p} (s_{\rm sky} + s_{\rm dc}) ]\ \; t_{\rm exp} - S^2 \; n_{\rm p} \; \sigma_{\rm r}^2 = 0" />
</p>
    
The minimum t<sub>exp</sub> of the equation above is given by its smallest non-negative root. Therefore, the optimum mode is given through the calculation of the AR of the selected modes for the minimum t<sub>exp</sub>.
    
* Mode 3: in this mode, both SNR and AR are optimized. Initially, it is selected those modes which accomplish the SNR and AR at the same time. The resulting list of modes is used to create the space of states of the BOM. Then, it is calculated the maximum values S<sup>M</sup> and A<sup>M</sup> and the minimum values S<sup>m</sup> and A<sup>m</sup> of the SNR and AR, respectively. They are used in normalization of both parameters into the range between 0 and 1. So, the function to be optimized is given by the multiplication of the normalized signal to noise ratio S<sub>NR</sub> and acquisition rate A values for each operation mode:

<p align="center">
<img src="https://latex.codecogs.com/svg.latex?f&space;=&space;\frac{{S}_{\rm&space;NR}&space;-&space;{S}^{\rm&space;m}}{{S}^{\rm&space;M}&space;-&space;{S}^{\rm&space;m}}&space;\times&space;\frac{{A}&space;-&space;{A}^{\rm&space;m}}{{A}^{\rm&space;M}&space;-&space;{A}^{\rm&space;m}}." title="f = \frac{{S} - {S}^{\rm m}}{{S}^{\rm M} - {S}^{\rm m}} \times \frac{{A} - {A}^{\rm m}}{{A}^{\rm M} - {A}^{\rm m}}." />
</p>

Therefore, the optimum mode for the CCD will be given by the set of parameters obtained through the BOM that maximizes the function given by the equation above. Figure below presents the SNR x AR values obtained as a function of the t<sub>exp</sub>, G<sub>em</sub> and readout rate of the CCD over the BOM iterations. Through this figure, it is possible to see a maximum point for the readout rate of 1 MHz.

<p align="center">   
    <img src="https://github.com/DBernardes/OMASS4/blob/main/iteracoes_MOB_ingles.png" />
</p>

## Running the OMASS4

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites
There are some packages that need to be installed before running the software.

* [astropy](https://www.astropy.org/)
* [hyperopt](https://github.com/WillKoehrsen/hyperparameter-optimization)
* [numpy](https://numpy.org/)
* [pandas](https://pandas.pydata.org/)
* [matplotlib](https://matplotlib.org/)
* [collections](https://docs.python.org/3/library/collections.html)
* [json](https://www.w3schools.com/python/python_json.asp)
* [xlrd](https://xlrd.readthedocs.io/en/latest/)
* [Photutils](https://photutils.readthedocs.io/en/stable/)
* [Scipy](https://www.scipy.org/)

To install these packages it is suggested to use the pip command as follows
```
pip install <package_name>
```

### Installing
Clone this repo using ``` git clone https://github.com/DBernardes/OMASS4.git ```

## Running the tests

To run a simple test, there is an image created artificially in the example directory. If you run the \_\_main\_\_.py file, the OMASS4 will be executed over this image. You can choose between the options to optimize the SNR, the acquisition rate, or both parameters providing the option 1, 2, or 3 for the optimize function, respectively. Also, you can choose to use or not the pre-image available changing the (y/n) parameter in the observation_setup.txt file. When the execution is done, the optimum mode will be printed on the screen, and a .txt file with the resulting information will be created in the image directory.
## Authors and Contact

* **Denis Bernardes**: 

email: denis.bernardes099@gmail.com 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details


