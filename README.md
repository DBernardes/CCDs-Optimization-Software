# CCDs Optimization Software

This repo presents the Optimization Method for the Electron Multiplying Charge Coupled Device (EMCCDs) of the Acquisition System of the SPARC4 (OMASS4). The OMASS4 was implemented in the python 3.7.4 language, and it found the optimum operation mode of the CCDs for some scientific requirements, based on the signal-to-noiseratio (SNR) and the acquisition rate (AR) as a function of the operation mode of the CCDs as figures of merit. The OMASS4 is structured into three parts: the initialization, the star flux calculation, and the CCD optimization. For the initialization step, it requires to provide to the software all the information related to the astronomical object, i.e.: an image of the object, its x,y coordinates, the maximum star radius, a bias image with the respective used CCD operation mode, the SNR, the AR, allowed SI and Bin modes, CCD temperature, and the iterations number of the BOM. Then, the star flux is calculated for the optimal star radius given by the full width at half maximum (FWHM) parameter, , as described in ref{sec:SNR}. The OMASS4 uses a set of packages to calculate the SNR and the AR values according to the CCD operation mode. The code developed to execute the BOM is based on the library provided by (hyperopt_koehrsen). The used algorithm to model the objective function of the BOM is the tree strucutured Parzen estimator (TPE). The SNR package operation is based on the methodology presented in ref{sec:SNR}. The star flux, sky flux, and the number of star pixels are obtained in the previous step. The DC noise is calculated according to the model presented by citenum{DenisVB} for the four SPARC4 cameras. The read noise is obtained through the characterization presented in \ref{subsec:read_noise_char}. For conventional modes, it is used the values in Table \ref{tab:rn_conv}. For EM modes, it is made an interpolation of the data presented in Fig. \ref{fig:EM_x_RN}. The G value is obtained through the camera datasheet \cite{iXonDataSheet}.

The performance of the EM mode is better than the conventional mode until a max value of 100 photons per pixel \cite{Andor_max_emgain}. So, the t<sub>exp</sub> of each EM mode is limited to accomplish this requirement. Also, the maximum value allowed for the G<sub>em</sub> is 300x. Values larger than 300x would deteriorate the device. Furthermore, the G<sub>em</sub> must be such that the CCD will not saturate. For this reason, the maximum EM gain allowed was arbitrarily configured to provide a signal up to 80 % of the pixel well depth. For an image with 16 bits per pixel, this value is 2<sup>16</sup> x 0.8 = 52429 analogical to digital unit (ADU). Given that a pixel value is composed by the star, sky and, dark current signals, and the bias level, the maximum value for the G<sub>em</sub> is

\begin{equation}
    G_{\rm em} = \frac{(52429 - B) \times G}{(S/n_{\rm p} + S_{\rm sky} + S_{\rm dc})}.
\label{eq:max_em_gain}    
\end{equation}

The AR package operation is based on the characterization presented in Sec.\ref{subsec:frequency_acq_char}. Be the value of the t<sub>c</sub> given by the Table \ref{tab:tempo_critico_all_modes}. For each mode, the AR value will be calculated through interpolation if t<sub>exp</sub> < t<sub>c</sub>, and it is the inverse of the t<sub>exp</sub>, if t_<sub>exp</sub> maior ou igual a t<sub>c</sub>.

Therefore, the OMASS4 was implemented using the aforementioned packages, being applied to three different optimization modes: optimize SNR (mode 1), optimize AR (mode 2), and optimize both SNR and AR (mode 3). 

[Mode 1:] in this mode, the SNR is optimized, keeping the AR fixed. First, it is selected those modes that accomplish the AR requirement. Then, it is calculated the SNR value for each selected mode, using the maximum values for the t<sub>exp</sub> and G<sub>em</sub>. The optimum mode is given by that one with the highest SNR.
    
item[Mode 2:] in this mode, the AR is optimized, keeping the SNR fixed. Initially, for each mode, it is calculated the minimum t<sub>exp</sub> value that accomplish the SNR requirement, for the maximum G<sub>em</sub> allowed. For this calculation, it is considered the values s = S/t<sub>exp</sub>, in photons/s of the star; the s<sub>sky</sub> = S<sub>sky</sub>/t<sub>exp</sub>, in photons/pixel/s of the sky; and the $s<sub>dc</sub> = S<sub>dc</sub>/t<sub>exp</sub>, in e-/pixel/s of the dark current. So, the Eq. \ref{eq:SNR_2} can be rewritten as follows
    
    \begin{equation}
             \mathcal{S} = \frac{s \times t_{\rm exp}}{\{ s \; t_{\rm exp} \; N_{\rm F}^2  + \\
            n_{\rm p} [\ (s_{\rm sky} +  s_{\rm dc}) \; t_{\rm exp} \; N_{\rm F}^2 +  \\
            (\sigma_{\rm ADU} \; G/G_{\rm em})^2 ]\ \}^{1/2}}.
        \label{eq:SNR_3_a}    
    \end{equation}    
      
Rearranging the terms of the Eq. \ref{eq:SNR_3_a} and isolating $t_{\rm exp}$,
    
    \begin{equation}
             s^2 \; t_{\rm exp}^2 -\\
             \mathcal{S}^2 \; N_{\rm F}^2 \; [\ s + n_{\rm p} (s_{\rm sky} + s_{\rm dc}) ]\ \; t_{\rm exp}  -\\
             \mathcal{S}^2 \; n_{\rm p} \; \sigma_{\rm r}^2 = 0 
    \label{eq:texp_min}        
    \end{equation}
    
The minimum t<sub>exp</sub> of Eq. \ref{eq:texp_min} is given by its smallest non-negative root. Therefore, the optimum mode is given through the calculation of the AR of the selected modes for the minimum t<sub>exp</sub>.
    
[Mode 3:] in this mode, both SNR and AR are optimized. Initially, it is selected those modes which accomplish the SNR and AR at the same time. The resulting list of modes is used to create the space of states of the BOM. Then, it is calculated the maximum values \mathcal{S}<sup>M</sup> and \mathcal{A}<sup>M</sup> and the minimum values \mathcal{S}<sup>m</sup> and \mathcal{A}<sup>m</sup> of the SNR and AR, respectively. They are used in normalization of both parameters into the range between 0 and 1. So, the function to be optimized is given by the multiplication of the normalized SNR and acquisition rate \mathcal{A} values for each operation mode:
    
    \begin{equation}
         f = \frac{\mathcal{S} - \mathcal{S}^{\rm m}}{\mathcal{S}^{\rm M} - \mathcal{S}^{\rm m}} \times \frac{\mathcal{A} - \mathcal{A}^{\rm m}}{\mathcal{A}^{\rm M} - \mathcal{A}^{\rm m}}.
    \label{eq:funcao_custo_SNR_FA}        
    \end{equation}
    
Therefore, the optimum mode for the CCD will be given by the set of parameters obtained through the BOM that maximizes the function given by Eq. \ref{eq:funcao_custo_SNR_FA}.

