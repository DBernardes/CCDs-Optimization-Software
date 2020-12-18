# CCDs Optimization Software

This repo presents the software developed for the optimmizaiton of the Charge Coupled Devices (CCDs) of the Simultanoues Polarimeter and Rapid Camera in Four Bands (SPARC4). This software was developed by using the python 3.7.4 language, and it is structured into three parts: the initialization, the star flux calculation, and the CCD optimization. For the initialization step, it requires to provide to the software all the information related to the astronomical object, i.e.: an image of the object, its $x,y$ coordinates, the maximum star radius, a bias image with the respective used CCD operation mode, the SNR, the AR, allowed SI and Bin modes, CCD temperature, and the iterations number of the BOM. Then, the star flux is calculated for the optimal star radius given by the FWHM parameter, as described in Sec.~\ref{sec:SNR}. The OMASS4 uses a set of packages to calculate the SNR and the AR values according to the CCD operation mode. The code developed to execute the BOM is based on the library provided by Ref. \citenum{hyperopt_koehrsen}. The used algorithm to model the objective function of the BOM is the TPE. The SNR package operation is based on the methodology presented in Sec.~\ref{sec:SNR}. The star flux, sky flux, and the number of star pixels are obtained in the previous step. The DC noise is calculated according to the model presented by Ref.~\citenum{DenisVB} for the four SPARC4 cameras. The read noise is obtained through the characterization presented in Sec.~\ref{subsec:read_noise_char}. For conventional modes, it is used the values in Table \ref{tab:rn_conv}. For EM modes, it is made an interpolation of the data presented in Fig. \ref{fig:EM_x_RN}. The $G$ value is obtained through the camera datasheet \cite{iXonDataSheet}.

The performance of the EM mode is better than the conventional mode until a max value of 100 photons per pixel \cite{Andor_max_emgain}. So, the $t_{\rm exp}$ of each EM mode is limited to accomplish this requirement. Also, the maximum value allowed for the $G_{\rm em}$ is 300x. Values larger than 300x would deteriorate the device~\cite{AndorSolis}. Furthermore, the $G_{\rm em}$ must be such that the CCD will not saturate. For this reason, the maximum EM gain allowed was arbitrarily configured to provide a signal up to 80 \% of the pixel well depth. For an image with 16 bits per pixel, this value is $2^{16} \times 0.8 \approx 52429$ ADU. Given that a pixel value is composed by the star, sky and, dark current signals, and the bias level, the maximum value for the $G_{\rm em}$ is

\begin{equation}
    G_{\rm em} = \frac{(52429 - B) \times G}{(S/n_{\rm p} + S_{\rm sky} + S_{\rm dc})}.
\label{eq:max_em_gain}    
\end{equation}

The AR package operation is based on the characterization presented in Sec. \ref{subsec:frequency_acq_char}. Be the value of the $t_{\rm c}$ given by the Table \ref{tab:tempo_critico_all_modes}. For each mode, the AR value will be calculated through interpolation if $t_{\rm exp} < t_{\rm c}$, and it is the inverse of the $t_{\rm exp}$, if $t_{\rm exp} \geq t_{\rm c}$.

Therefore, the OMASS4 was implemented using the aforementioned packages, being applied to three different optimization modes: optimize SNR (mode~1), optimize AR (mode 2), and optimize both SNR and AR (mode 3). 

\begin{description}
    \item[Mode 1:] in this mode, the SNR is optimized, keeping the AR fixed. First, it is selected those modes that accomplish the AR requirement. Then, it is calculated the SNR value for each selected mode, using the maximum values for the $t_{\rm exp}$ and $G_{\rm em}$. The optimum mode is given by that one with the highest SNR.
    
    \item[Mode 2:] in this mode, the AR is optimized, keeping the SNR fixed. Initially, for each mode, it is calculated the minimum $t_{\rm exp}$ value that accomplish the SNR requirement, for the maximum $G_{\rm em}$ allowed. For this calculation, it is considered the values $s = S/t_{\rm exp}$, in photons/s of the star; the $s_{\rm sky} = S_{\rm sky}/t_{\rm exp}$, in photons/pixel/s of the sky; and the $s_{\rm dc} = S_{\rm dc}/t_{\rm exp}$, in e-/pixel/s of the dark current. So, the Eq. \ref{eq:SNR_2} can be rewritten as follows
    
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
    
     The minimum $t_{\rm exp}$ of Eq. \ref{eq:texp_min} is given by its smallest non-negative root. Therefore, the optimum mode is given through the calculation of the AR of the selected modes for the minimum $t_{\rm exp}$.
    
    \item[Mode 3:] in this mode, both SNR and AR are optimized. Initially, it is selected those modes which accomplish the SNR and AR at the same time. The resulting list of modes is used to create the space of states of the BOM. Then, it is calculated the maximum values $\mathcal{S}^{\rm M}$ and $\mathcal{A}^{\rm M}$ and the minimum values $\mathcal{S}^{\rm m}$ and $\mathcal{A}^{\rm m}$ of the SNR and AR, respectively. They are used in normalization of both parameters into the range between 0 and 1. So, the function to be optimized is given by the multiplication of the normalized SNR and acquisition rate $\mathcal{A}$ values for each operation mode:
    
    \begin{equation}
         f = \frac{\mathcal{S} - \mathcal{S}^{\rm m}}{\mathcal{S}^{\rm M} - \mathcal{S}^{\rm m}} \times \frac{\mathcal{A} - \mathcal{A}^{\rm m}}{\mathcal{A}^{\rm M} - \mathcal{A}^{\rm m}}.
    \label{eq:funcao_custo_SNR_FA}        
    \end{equation}
    
    Therefore, the optimum mode for the CCD will be given by the set of parameters obtained through the BOM that maximizes the function given by Eq. \ref{eq:funcao_custo_SNR_FA}.
\end{description}

