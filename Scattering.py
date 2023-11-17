import os
import numpy as np
from numpy import pi, sqrt, arctan, sin, cos, roots, size, real, imag, sqrt, exp
import const
from scipy.interpolate import interp1d
from scipy.signal import periodogram, welch
from gwinc.ifo.noises import ifo_power



frequency = np.linspace(1e-2,100,4096)

def loginterp(x, y, f):
    return np.power(10, np.interp(np.log10(f),np.log10(x),np.log10(y)))


def spectrum_rayleigh_horizontal(ifo):
    frequency = np.linspace(1e-2,100,4096)
    dataSosEnattos = np.loadtxt("acoustic_spectra/rwave_spectrum_SosEnattos.txt")
    dataTerziet = np.loadtxt("acoustic_spectra/rwave_spectrum_Terziet.txt")
    
    if ifo.Seismic.Site=='ET':
        rayleighwave=(10**(
            0.5 * (
                np.log10(gwinc.noise.seismic.seismic_ground_NLNM(frequency))
                + np.log10(gwinc.noise.seismic.seismic_ground_NHNM(frequency))
                )
            )
            )**2
    elif ifo.Seismic.Site=='SosEnattos':
        rayleighwave=loginterp(dataSosEnattos.T[0],dataSosEnattos.T[1],frequency)
    elif ifo.Seismic.Site=='Terziet':
        rayleighwave=loginterp(dataTerziet.T[0],dataTerziet.T[1],frequency)
        
    return rayleighwave
    

def freq_to_time(ifo):
    '''
    This function generates a time series with a given PSD.
    
    Inputs:
        - f: Vector of frequencies [Hz]
        - sqrtPSD: Vector of seismic noise [m/sqrt(Hz)]
        
    Outputs:
        - t: vector conatining the times [s]
        - xt: vector with the time somain samples
    '''
    # Lenght of the time domain series
    frequency = np.linspace(1e-2,100,4096)
    M = 2*len(frequency)
    # Frequency step
    df = frequency[1]-frequency[0]
    # Time step
    dt = 1./(M*df)
    # PSD
    PSD = spectrum_rayleigh_horizontal(ifo)
    # Coefficients of each frequency
    Ak = np.sqrt(2*PSD*df)
    # Time vector
    t = np.arange(0,(M-1)*dt,dt)
    # Intialize the vector of time domain data
    xt = np.zeros(len(t))

    # Main loop
    for i in range(len(frequency)):
        # random phase between 0 and 2pi
        phi = np.random.uniform(0,2*np.pi)
        # Add the contribution of each frequency
        xt += Ak[i]*np.sin(2*np.pi*frequency[i]*t+phi)
    return t,xt

    
def Upconversion(ifo):
    frequency = np.linspace(1e-2,100,4096)
    lam=ifo.Laser.Wavelength
    t, Xt = freq_to_time(ifo)
    Ns = 512/2

    # Frequency step
    df = 1/(t[1]-t[0])

    # Auximilary variables
    arg1 = np.sin(4*np.pi/lam*Xt)
    arg2 = np.cos(4*np.pi/lam*Xt)

    # Compute the PSD of these auxiliary variables
    fp, psd1 = welch(arg1, df,nperseg=Ns,average='mean',scaling='density')
    fp, psd2 = welch(arg2, df,nperseg=Ns,average='mean',scaling='density')

    # Obtain the upcovnverted noise
    Xup = lam/(4*np.pi)*np.sqrt(psd1+psd2)

    # Ignore frequencies below the minimum
    Xup = Xup[fp>np.min(frequency)]
    fp  = fp[fp>np.min(frequency)]
    
    
    return fp,Xup


def Backscattering(ifo):
    '''
    This function computes the backscattering noise assuming a simple analytical model.
    
    Inputs:
        - f: Vector of frequencies [Hz]
        - X: Vector of seismic noise [m/sqrt(Hz)]
        - Pcirc: Circulating power of  the FP cavity [W]
        - k: BDRF parameter of the mirror when assumed of the form BRDF = k*theta^{-n} [str^{-1}]
        - n: BDRF parameter of the mirror when assumed of the form BRDF = k*theta^{-n}
        - m: Mass of the mirrors [kg]
        - L: Length of the FP cavity [m]
        - R: Radius of the tube [m]
        - Gain: Gain of the cavity formet bu the ITM and SRM
        - BRDFbaff: BRDF of the baffles [str^{-1}]
        - lam: Wavelength of the laser [m]
        - l0: Position of the first baffle [m]
        
    Outputs:
        - h: vector conatining the backscattered light noise for eadch f [1/sqrt(Hz)]
    '''
    frequency = np.linspace(1e-2,100,4096)
    Gain = 15.7
    L=ifo.Infrastructure.Length
    k=ifo.Materials.Substrate.k
    n=ifo.Materials.Substrate.nn
    m=ifo.Suspension.Stage[0].Mass
    R=ifo.Infrastructure.PipR
    BRDFbaff=ifo.Infrastructure.BRDFbaff
    l0=ifo.Infrastructure.l0
    lam=ifo.Laser.Wavelength
    Pcirc=ifo_power(ifo).parm
    
    f,X=Upconversion(ifo)
        
    Prad = 8*Gain*Pcirc/(const.c*m*np.pi*f**2)
    TF = lam**2+Prad**2

    tmin = R/(L-l0)
    tmax = R/(L/2)
  
    
    if n == 2:
        integral = np.log(tmax/tmin)
    else:
        integral = (tmax**(4-2*n)-tmin**(4-2*n))/(4-2*n)
        
    sumK = 2*np.pi*k**2/R**2*integral
    Disp = np.sqrt(TF*BRDFbaff)*X*np.sqrt(2*sumK) #displacement
    return f, Disp
    
def BackResult(f,ifo):
    ff, Disp =Backscattering(ifo)
    noise=loginterp(ff, Disp,f)
    return noise
