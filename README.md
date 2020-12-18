# Optimization Method of the Electron Multiplying Charge Coupled Device of the Acquisition System of the SPARC4

This repo presents the Optimization Method of the Electron Multiplying Charge Coupled Device (EMCCDs) of the Acquisition System of the [SPARC4](https://www.spiedigitallibrary.org/conference-proceedings-of-spie/8446/844626/Concept-of-SPARC4--a-simultaneous-polarimeter-and-rapid-camera/10.1117/12.924976.full?casa_token=7b-hbhyqIMoAAAAA%3a99lzc7LW-gGeFuEs1N_7ZGdcFS1EiapC3jbzEYyrWT3PDiUP4RXPDEiR9IdfuRvDY7pPetsPx88&SSO=1) (OMASS4). The OMASS4 uses as figures of merit the signal-to-noiseratio (SNR) and the acquisition rate (AR) as a function of the operation mode of the CCDs. Three different modes of optimization are included in the OSASS4:  (1) optimization of SNR only; (2) optimization of AR only; and (3) optimization of both SNR and AR simultaneously. The first two modes calculate an analytical minimization of the costfunction whereas the third mode uses the bayesian optimization method (BOM) to determine the optimum mode of operation. We apply the OMASS4 to find the optimum mode for observations obtained at the Pico dos Dias Observatory, Brazil, and compare the delivered modes of operation and its performance with the ones adopted by the observer. If the OSASS4 had been used as a tool to optimize the CCDs in all of these nights, it would be possible to improve their efficiency in 97.17 %, 65.08 %, and 77.66 % for the optimization modes 1, 2, and 3, respectively.

For more information, access: [Otimização dos modos de operação do sistema de aquisição do instrumento SPARC4](https://repositorio.unifei.edu.br/jspui/handle/123456789/2201)

## Software Description

This section presents the software developed to implemente the OMASS4. This software was developed using Python Language 3.7.4 to determine the optimum operation mode of the SPARC4 CCDs. It was structured into three parts: the initialization, the star flux calculation, and the CCD optimization. For the initialization step, it requires to provide to the software all the information related to the astronomical object, i.e.: an image of the object, its x,y coordinates, the maximum star radius, a bias image with the respective used CCD operation mode, the SNR, the AR, allowed SI and Bin modes, CCD temperature, and the iterations number of the BOM. Then, the star flux is calculated for the optimal star radius given by the FWHM parameter. The OMASS4 uses a set of packages to calculate the SNR and the AR values according to the CCD operation mode. The code developed to execute the BOM is based on the library provided by [Koehrsen](https://github.com/WillKoehrsen/hyperparameter-optimization). The used algorithm to model the objective function of the BOM is the TPE. The SNR package operation is based on the methodology presented by [Bernardes, D. V. (2020)](https://repositorio.unifei.edu.br/jspui/handle/123456789/2201). The star flux, sky flux, and the number of star pixels are obtained thorugh a pre-image of the object. The DC noise is calculated according to the model presented by [Bernardes et. al (2018)](https://iopscience.iop.org/article/10.1088/1538-3873/aacb1e/meta?casa_token=QzaY5kK_Yp8AAAAA:Qz_wlI6tq2WMi4sRF-tLvw-S2RwkmkF1_N8i7mReLYSUgim4dqp3yceqyLmlbrgUHt5TTzGYcrnYW_9ttxnfrw) for the four SPARC4 cameras. The read noise is obtained through the characterization presented in [Bernardes, D. V. (2020)](https://repositorio.unifei.edu.br/jspui/handle/123456789/2201). The G value is obtained through the camera datasheet.

The performance of the EM mode is better than the conventional mode until a max value of 100 photons per pixel. So, the t<sub>exp</sub> of each EM mode is limited to accomplish this requirement. Also, the maximum value allowed for the G<sub>em</sub> is 300x. Values larger than 300x would deteriorate the device. Furthermore, the G<sub>em</sub> must be such that the CCD will not saturate. For this reason, the maximum EM gain allowed was arbitrarily configured to provide a signal up to 80 \% of the pixel well depth. For an image with 16 bits per pixel, this value is 2<sup>16</sup> \times 0.8 = 52429 analogical to digital unit (ADU). Given that a pixel value is composed by the star, sky and, dark current signals, and the bias level, the maximum value for the G<sub>em</sub> is

![imagem](https://latex.codecogs.com/svg.latex?G_%7B%5Crm%20em%7D%20%3D%20%5Cfrac%7B%2852429%20-%20B%29%20%5Ctimes%20G%7D%7B%28S/n_%7B%5Crm%20p%7D%20&plus;%20S_%7B%5Crm%20sky%7D%20&plus;%20S_%7B%5Crm%20dc%7D%29%7D.)

The AR package operation is based on the characterization presented of the CCDs presented by [Bernardes, D. V. (2020)](https://repositorio.unifei.edu.br/jspui/handle/123456789/2201). For each mode, the AR value will be calculated through interpolation if t<sub>exp</sub> < t<sub>c</sub>, and it is the inverse of the t<sub>exp</sub>, if t<sub>exp</sub> equals or greater than t<sub>c</sub>, where the t<sub>c</sub> is the time spent by the camera to read one image, as a function of the CCD operation mode.

Therefore, the OMASS4 was implemented using the aforementioned packages, being applied to three different optimization modes: optimize SNR (mode 1), optimize AR (mode 2), and optimize both SNR and AR (mode 3). 

* Mode 1: in this mode, the SNR is optimized, keeping the AR fixed. First, it is selected those modes that accomplish the AR requirement. Then, it is calculated the SNR value for each selected mode, using the maximum values for the t<sub>exp</sub> and G<sub>em</sub>. The optimum mode is given by that one with the highest SNR.
    
* Mode 2: in this mode, the AR is optimized, keeping the SNR fixed. Initially, for each mode, it is calculated the minimum t<sub>exp</sub> value that accomplish the SNR requirement, for the maximum G<sub>em</sub> allowed. For this calculation, it is considered the values of the star flux s in photons/s; the sky flux s<sub>sky</sub>, in photons/pixel/s; and the dark current s<sub>dc</sub>, in e-/pixel/s . So, the equation for the SNR of the star can be written as follows
    
![imagem](https://latex.codecogs.com/svg.latex?%5Cmathcal%7BS%7D%20%3D%20%5Cfrac%7Bs%20%5Ctimes%20t_%7B%5Crm%20exp%7D%7D%7B%5C%7B%20s%20%5C%3B%20t_%7B%5Crm%20exp%7D%20%5C%3B%20N_%7B%5Crm%20F%7D%5E2%20&plus;%20%5C%5C%20n_%7B%5Crm%20p%7D%20%5B%5C%20%28s_%7B%5Crm%20sky%7D%20&plus;%20s_%7B%5Crm%20dc%7D%29%20%5C%3B%20t_%7B%5Crm%20exp%7D%20%5C%3B%20N_%7B%5Crm%20F%7D%5E2%20&plus;%20%5C%5C%20%28%5Csigma_%7B%5Crm%20ADU%7D%20%5C%3B%20G/G_%7B%5Crm%20em%7D%29%5E2%20%5D%5C%20%5C%7D%5E%7B1/2%7D%7D.)
 
    
 Rearranging the terms of the equation above and isolating t<sub>exp</sub>,
    
![imagem](https://latex.codecogs.com/svg.latex?s%5E2%20%5C%3B%20t_%7B%5Crm%20exp%7D%5E2%20-%20%5Cmathcal%7BS%7D%5E2%20%5C%3B%20N_%7B%5Crm%20F%7D%5E2%20%5C%3B%20%5B%5C%20s%20&plus;%20n_%7B%5Crm%20p%7D%20%28s_%7B%5Crm%20sky%7D%20&plus;%20s_%7B%5Crm%20dc%7D%29%20%5D%5C%20%5C%3B%20t_%7B%5Crm%20exp%7D%20-%20%5Cmathcal%7BS%7D%5E2%20%5C%3B%20n_%7B%5Crm%20p%7D%20%5C%3B%20%5Csigma_%7B%5Crm%20r%7D%5E2%20%3D%200)
    
The minimum t<sub>exp</sub> of the equation above is given by its smallest non-negative root. Therefore, the optimum mode is given through the calculation of the AR of the selected modes for the minimum t<sub>exp</sub>.
    
* Mode 3: in this mode, both SNR and AR are optimized. Initially, it is selected those modes which accomplish the SNR and AR at the same time. The resulting list of modes is used to create the space of states of the BOM. Then, it is calculated the maximum values S<sup>M</sup> and A<sup>M</sup> and the minimum values S<sup>m</sup> and A<sup>m</sup> of the SNR and AR, respectively. They are used in normalization of both parameters into the range between 0 and 1. So, the function to be optimized is given by the multiplication of the normalized SNR and acquisition rate A values for each operation mode:

![imagem](https://latex.codecogs.com/svg.latex?f%20%3D%20%5Cfrac%7B%5Cmathcal%7BS%7D%20-%20%5Cmathcal%7BS%7D%5E%7B%5Crm%20m%7D%7D%7B%5Cmathcal%7BS%7D%5E%7B%5Crm%20M%7D%20-%20%5Cmathcal%7BS%7D%5E%7B%5Crm%20m%7D%7D%20%5Ctimes%20%5Cfrac%7B%5Cmathcal%7BA%7D%20-%20%5Cmathcal%7BA%7D%5E%7B%5Crm%20m%7D%7D%7B%5Cmathcal%7BA%7D%5E%7B%5Crm%20M%7D%20-%20%5Cmathcal%7BA%7D%5E%7B%5Crm%20m%7D%7D.)
    
Therefore, the optimum mode for the CCD will be given by the set of parameters obtained through the BOM that maximizes the function given by the equation above.

## Running the OMASS4

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites
There are some packages that need to be installed before running the software.

* [Astropy](https://www.astropy.org/)
* [Photutils](https://photutils.readthedocs.io/en/stable/)
* [Collections](https://docs.python.org/3/library/collections.html)
* [JSON](https://www.w3schools.com/python/python_json.asp)
* [Pandas](https://pandas.pydata.org/)
* [xlrd](https://xlrd.readthedocs.io/en/latest/)
* [Scipy](https://www.scipy.org/)

To install these packages it is suggested to use the pip command as follows
```
pip install <package_name>
```

### Installing
Clone this repo using ``` git clone https://github.com/DBernardes/OMASS4.git ```

## Running the tests

To run a simple test, you only need to execute the run.py file and the image would be created in your current directory. The run.py file will provide to the AIG the basic information for its execution, that is the star flux, in photons/s; the sky flux, in photons/pix/s, the standard deviation of the Gaussian, in pixels, and the operation mode of the CCD. In particular, the CCD operation mode should be a python dictionary with the control parameters used to configure the acquisition of the SPARC4 cameras. They are the Electron Multiplying Mode (em_mode), the Electron Multiplying Gain (em_gain), the Pre-amplification (preamp), the Horizontal Shift Speed (hss), the Pixels Binning (bin), and the Exposure Time (texp). Below, it is presented the accepted values for each parameter previously described.

- em_mode: 0 or 1
- em_gain: from 2 to 300
- preamp: 1 or 2
- hss: 0.1, 1, 10, 20, and 30
- bin: 1 or 2
- texp: greater or equal than 1e-5

Beyond the paramaters presented before, there are a set of optional paramaters. They are the CCD temperature (ccd_temp), the CCD serial number (serial_number), the image bias level (bias_level), and the directory where the image should be saved (image_dir)

- ccd_temp: from 0 ºC to -70 ºC
- serial_number: 9914, 9915, 9916, or 9917
- bias_level: integer and greater or equal than 1
- image_dir: string

## Authors and Contact

* **Denis Bernardes**: 

email: denis.bernardes099@gmail.com 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

