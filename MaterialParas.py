import numpy as np
import gwinc.noise.seismic
from scipy import constants
from scipy.interpolate import interp1d

# Interpolate spectrum data in log-space
def loginterp(x, y, f):
    return np.power(10, np.interp(np.log10(f),np.log10(x),np.log10(y)))
    #10**interp1d(np.log10(x), np.log10(y), kind="slinear", fill_value="extrapolate")(np.log10(f))

################################################################################
## Spectrum data
################################################################################
def Si_thermcond(T): #www.ioffe.ru/SVA/NSM/Semicond/Si/thermal.html
    Si_thermcond_T = np.array([6.46,7.37,10.5,13.7,19.42,20.64,24.19,23.07,26.96,31.81,41.68,53.16,79.96,110.16,164.02,219.17,328.6])
    Si_thermcond_K = np.array([2.13,4.46,8.27,16.5,23.3,27.26,26.14,32.74,36.53,34.3,35.39,31.2,17.45,10.4,5.15, 2.61,1.59])
    
    return loginterp(Si_thermcond_T,Si_thermcond_K,T)*100
    
def Si_specheat(T): #www.ioffe.ru/SVA/NSM/Semicond/Si/thermal.html
    Si_specheat_T = np.array([2.47,3.94, 6.037062194459862, 9.24,13.67, 19.53, 26.95,37.84,44.97, 58.26,83.22,131.1, 210.06,325.19])
    Si_specheat_Cp = np.array([0.00405,0.0161,0.0551,0.228,0.9247,3.52,13.97,44.98,72.7,120,211,348.3,540,723.5]) #J/g/K
    
    return loginterp(Si_specheat_T,Si_specheat_Cp,T) #J/kg/K
    
    
def Si_loss(T): #R Nawrodt et al 2008 J. Phys.: Conf. Ser. 122 012008
    Si_loss_T = np.array([4.81,7.83, 7.84, 9.66,14.84, 23.38, 34.65,56.60, 68.79,76.41, 77.62, 77.61, 79.13, 86.15, 92.56, 104.46, 114.83, 133.42,153.53,176.40, 198.66,219.39, 230.97, 250.79, 274.26, 291.64, 300.48])
    Si_loss_phi = np.array([2.18e-9,2.64e-9,4.39e-9,6e-9,6.97e-9,7.33e-9,8.96e-9,9.427e-9,1.07e-8,1.54e-8,2.34e-8,3.3e-8,4.19e-8,3.04e-8,1.94e-8,1.41e-8, 1.29e-8, 1.67e-8, 2.19e-8, 1.83e-8, 1.91e-8, 2.34e-8, 2.97e-8, 2.56e-8, 3.32e-8, 3.1e-8, 3.57e-8]) #@14kHz
    
    return loginterp(Si_loss_T,Si_loss_phi,T)
    
def Si_RefractiveIndex(T): #arxiv.org/pdf/physics/0606168.pdf
    Si_RefractiveIndex_T = np.array([30,40,50,60,70,80,90,100,150,200,250,295])
    Si_RefractiveIndex_n = np.array([3.45309,3.45319,3.45340,3.45373,3.45417,3.45471,3.45535,3.45609,3.46100,3.46760,3.47540,3.48314])
    
    return loginterp(Si_RefractiveIndex_T,Si_RefractiveIndex_n,T)
    
    
def Si_thermoptic(T): #doi.org/10.1063/1.4738989
    Si_thermoptic_T = np.array([5.207, 9.980, 15.17, 19.49, 25.07,29.80, 39.67, 52.95, 71.80,92.36, 118.9,147.6, 175.4, 205.4,234.0,263.6, 289.7,300])
    Si_thermoptic_dndT = np.array([9.693e-9,8.220e-8,4.426e-7,0.000001473,0.000002874,0.000006941,0.00001589,0.00002571,0.00004387,0.00006550,0.00008788,0.0001060,0.0001422,0.0001540,0.0001714,0.0001625,0.0001907,0.0001808])
    
    return loginterp( Si_thermoptic_T,Si_thermoptic_dndT,T)
    
        
def Si_thermalexpension(T): ##doi.org/10.1103/PhysRevB.92.174113
    Si_thermalexpension_T = np.linspace(8.15,293.15,58)
    Si_thermalexpension_K = np.array([0,-0.1,-2,-13,-40.9,-86.2,-144,-207.9,-271.9,-331.7,-383.7,-425.5,-455.5,-472.6,-476.2,-466.5,-443.8,-409,-363,-306.9,-242.1,-169.6,-90.8,-6.8,81.4,172.8,266.4,361.7,457.7,554,650.1,745.5,839.9,932.9,1024.4,1114.1,1201.9,1287.6,1371.2,1452.6,1531.7,1608.6,1683.3,1755.7,1825.9,1893.9,1959.7,2023.5,2085.1,2144.8,2202.4,2258.2,2312.2,2364.4,2414.7,2463.4,2510.5,2556.1])
    
    return np.interp(T,Si_thermalexpension_T,Si_thermalexpension_K)*1e-9


