# -*- coding: utf-8 -*-
"""
Created on Mon Jan 07 10:28:13 2019
@author: Marie Gruenbein, Max-Planck-Institute for Medical Research
         marie.gruenbein@mr.mpg.de

This programme calculates the absorbed number of photons per molecule inside a
light-absorbing sample based on the given sample and illumination properties.
It assumes equal probabilities for single-and multi-photon absorption as a
strong over-simplification.

Incident pump laser fluence is assumed to be the same for all area illuminated.
This corresponds to the situation where the probe beam is very small, probing a
cross sectional area of the excited sample over which the incident laser
power did not change significantly. E.g., this is the case in typical SFX
experiments where the optical pump pulse is focused to ~100 um diameter and the
X-ray probe beam is focused to ~1 um diameter.
As default, it is assumed that the probe beam is centred on the pump beam, thus
sampling molecules that were excited with the peak intensity of the pump pulse. 
The offset parameter can be used to shift the probing region off the peak
pump intensity.
Generally, these calculations assume that the pump pulse intensity has no
temporal variation and that the spatial intensity distribution is Gausssian.

The average number of absorbed molecules is evaluated as a function of depth 
into the sample (along the pump beam propagation axis) by subdividing the sample
into thin 'layers', calculating the absorption per layer according to Lambert-
Beer and subtracting the number of absorbed photons in the preceding k-1 layers
to calculate the remaining pump intensity and resulting number of absorbed
photons in layer k.


Input parameters:

    use_param:  string, either 'conc' or 'xtal_param'. If 'conc', use the
                concentration of absorbing molecules as input. If 'xtal_param' 
                calculate the concentration of absorbing molecules from the 
                input cell parameters of the crystal.
                
    a,b,c:      non-negative float, required if use_param='xtal_param'.
                Lengths of the principal axes of the crystal unit cell [m].
                
    alpha,beta,gamma:   non-negative float, required if use_param='xtal_param'.
                        Angles between principal axes of the unit cell [rad].
                        
    N_mol_UC:   non-negative int, required if use_param='xtal_param'.
                Number of absorbing molecules per unit cell.
    
    conc_mol:   non-negative float, required if use_param='conc'.
                Concentration of absorbing molecules in sample [mol/l].
                
    epsilon:    non-negative float. Molar absorption coefficient at given
                pump laser wavelength [l/(mol cm)].
    
    L_sample:   non-negative float. Sample thickness [m].
    
    moltype:    string. Description of sample. Used for title of Figures and 
                output files.
    
    E_pulse:    non-negative float. Pulse energy of the pump laser [J].
    
    d_beam:     non-negative float. Beam diameter of the pump laser [m].
    d_type:     string, either '1/e2' or 'FWHM'. Defines whether the input beam
                laser diameter d_beam describes the full width half
                maximum (FWHM) or the width where the intensity has decreased
                to 1/e^2 of the peak intensity.
                
    T_pulse:    non-negative float. Pump laser pulse duration (FWHM) [s]. To
                calculate the laser power density at the interaction point (IP).
                
    offset:     non-negative float. Distance between the interaction point (IP)
                (i.e. the focus of the probe pulse) and the centre of the pump
                laser beam [m]. If offset>0, the fluence and power density at
                the IP are lower than the peak fluence and peak power density.
                The decrease is a function of the ratio offset/d_beam.
                
    wavel:      non-negative float. Pump laser wavelength [m].
    
    savefig:    boolean. If True: Save figures as type specified in savefigas
                with name specified in savebasename.
                If False: Figures are not saved.
                
    savedata:   boolean. If True: Save input and output data as .txt file with
                name specified in savebasename. If False: Do not save data.
                
    savebasename: name figures shall be saved as.
    
    savefigas:  file type for saving figure (e.g., png, pdf, eps)
    

Output:
    
Output printed to screen:
    
    -   Fluence and power density at the interaction point (IP). The temporal 
        shape of the pump laser pulse is not considered such that
        power density = fluence / pulse duration.
    -   Average number of absorbed photons/molecule (averaged over all
        molecules in the sample of thickness L_sample).
    -   Fraction of molecules absorbing 0,1,2,3,... photons in the oversimpli-
        fied assumption of equal absorption probabilities to absorb N=1 or N>1
        photon. Fraction of molecules absorbing k photons is approximated by
        the ratio of the layer thickness absorbing k-0.5 <= k < k+0.5 photons
        to the total sample thickness.
        
Figures:
    
    1.  Number of photons absorbed per molecule as a function of depth into the
        sample (measured along the pump beam propagation axis).
        
    2.  Photon regimes, plotting the fraction of molecules absorbing 0,1,2,...
        photons.
        

Output saved to files:
    
    - Figures 1 and 2 (if savefig=True)
    
    - txt file with input parameters and results (if savedata=True)


"""

import numpy as np
import matplotlib.pyplot as plt

avogadro    = 6.02214086e23 # Avogadro's number [1/mol]
planck      = 6.62607004e-34 # Planck's constant [m^2 kg / s]
c_light     = 299792458 # speed of light in vacuum [m/s]


#%% required input parameters

# -------- Sample parameters --------

# Use crystal parameters or concentration of molecules as input
use_param   = 'xtal_param' # either 'xtal_param' or 'conc'

# Unit cell axes lengths [m] (required if use_param='xtal_param')
a           = 62e-10
b           = 62e-10
c           = 111e-10

# Unit cell angles [rad] (required if use_param='xtal_param')
alpha       = np.deg2rad(90)
beta        = np.deg2rad(90)
gamma       = np.deg2rad(120)

# Number of absorbing molecules per unit cell (required if use_param='xtal_param')
N_mol_UC    = 6

# Concentration of absorbing molecules [mol/l] (required if use_param='conc')
conc_mol    = 0.027

# Absorption coefficient [ l / (mol cm) ]
epsilon     = 45600

# Sample thickness [m]
L_sample    = 5e-6

# Name of system / protein / crystal
moltype     = 'MySample'



# -------- Pump laser paramters --------

# Pulse energy [J]
E_pulse      = 0.25e-6

# Beam waist diameter [m]
d_beam      = 100e-6
d_type      = "1/e2" # 1/e2 or FWHM

# Pulse length [s]
T_pulse     = 145e-15

# Laser offset [m]
offset      = 1e-6

# Wavelength [m]
wavel       = 532e-9



# -------- save output? --------
# If true, will save figures and output text file
savefig         = True     # True or False
savedata        = True     # True or False
savebasename    = 'mysample' # any string of your choice
savefigas       = 'png'     # pdf, png, eps


#%% Derived parameters

# ---- crystal

# Absorption coefficient unit conversion
epsilon         *= 100 # [ l / (mol cm) ] -> [ l / (mol m) ]


if use_param == 'xtal_param':

    # Unit cell volume [m^3]
    V_UC        = a*b*c*np.sqrt(1+2*np.cos(alpha)*np.cos(beta)*np.cos(gamma)-\
                            np.cos(alpha)**2-np.cos(beta)**2-np.cos(gamma)**2)
    V_UC_l      = V_UC*1000 # in [l]
    
    # Density of absorbing molecules 
    rho_mol     = N_mol_UC/V_UC # [molecules/m^3]
    rho_mol_l   = N_mol_UC/V_UC_l # [molecules/l]
    
    # Concentration of absorbing molecules [mol/l]
    conc_mol    = rho_mol_l/avogadro
    
elif use_param == 'conc':
    
    rho_mol_l   = conc_mol*avogadro # [molecules/l]
    rho_mol     = rho_mol_l*1000    # [molecules/m^3]

# penetration depth
d_pen = -1*np.log10(1/np.exp(1))/(conc_mol*epsilon)    

print("\nConcentration of absorbing molecules: {:.3f} mol/l".format(conc_mol))
print("Penetration depth of sample: {:.1e} m".format(d_pen))



# ---- laser

# if necesary convert FWHM beam diameter to 1/e^2 diameter
if d_type == "FWHM": d_beam *= 1.699

# beam radius [m]
r_beam          = 0.5*d_beam

# effective beam area (spatially gaussian pulse) [m^2]
A_beam          = 0.5*np.pi*r_beam**2

# reduction factor due to laser offset (assume spatially Gaussian intensity distribution)
offsetred       = np.exp(-2*offset**2/r_beam**2)

# peak fluence at centre of Gaussian beam and at IP
F_peak          = E_pulse/A_beam # [J/m^2]
F_IP            = F_peak*offsetred # [J/m^2]

# power density at centre of Gaussian beam and at IP
P_peak          = F_peak/T_pulse # [W/m^2]
P_IP            = P_peak*offsetred # [W/m^2]

# Photon energy [J]
E_photon        = planck*c_light/wavel

# Photon density at IP [photons/m^2]
rho_photon_IP   = F_IP/E_photon

print("\nFluence at IP: {} mJ/mm^2".format(F_IP/1000))
print("Power density at IP: {} GW/cm^2\n".format(P_IP*1e-9/(100*100)))

#%% Apply Lambert-Beer iteratively to layers at different sample depths


# Define Lambert-Beer law
def LambertBeer(eps,conc,L):
    return 10**(-1*eps*conc*L)


# average distance between molecules [m] used as layer thickness
d_layer     =  np.power(1/rho_mol,(1./3)) 
# Number of layers for sample of thickness L_sample
N_layers    = L_sample/d_layer
# Sampling points (depth [m]) for each layer of thickness d_layer
pts         = np.arange(0,L_sample,d_layer)

# fraction of photons absorbed per layer of thickness d_layer
f_LB        = 1-LambertBeer(epsilon,conc_mol,d_layer)

# area density of absorbing molecules within thickness d_layer [1/m^2]
N_mol_A     = rho_mol*d_layer
# Number of molecules in a cuboid of 1x1um^2 cross section x d_layer thickness
N_mol_1um   = N_mol_A*1e-12
# Number of photons incident onto a square with 1x1 um^2 cross section
N_phot_1um  = rho_photon_IP*1e-12


# Initialize array
pts_photons_per_mol_absd    = np.zeros(pts.shape)

# Iteratively calculate number of absorbed photons/molecule in different layers
for kk in np.arange(len(pts)):
    # Number of photons absorbed in layer kk
    pts_photons_per_mol_absd[kk] = f_LB*N_phot_1um/N_mol_1um
    # Number of photons remaining to propagate to subsequent layers
    N_phot_1um -= f_LB*N_phot_1um
del kk

# Print average number of photons per molecule
print('Sample thickness: {} um'.format(L_sample*1e6))
print('Average number of photons absorbed/molecule: {a:.1f}'.format(a=np.average(pts_photons_per_mol_absd)))


#%% Plot Figure 1
# Number of photons absorbed per molecule as function of depth into sample.

plt.figure()
plt.plot(pts*1e6,pts_photons_per_mol_absd,'.', label= 'Eps = {} l/(mol cm)'.format(epsilon/100))
plt.xlabel(r'Position in sample [$\rm\mu$m]')
plt.ylabel('Absorbed # photons/molecule')
plt.suptitle('Number of photons absorbed per molecule',y=1.0,fontsize=14)
plt.title(r'Laser: E={e:.2f}$\rm\mu$J, {a:.0f}$\rm\mu$m offset, w={w:.0f}$\rm\mu$m, $\rm\lambda$={l:.0f}nm. {s} parameters'.format(e=E_pulse*1e6,a=offset*1e6,\
          w=r_beam*1e6, l=wavel*1e9, s=moltype ))
plt.grid()
plt.legend(loc='best')
if savefig: plt.savefig(savebasename+'.'+savefigas, transparent=True)


#%% Calculate fraction of molecules absorbing 0,1,2,3,... photons

# Max. number of absorbed photons/molecule (round to next integer)
Nphot_max   = int(np.ceil(np.max(pts_photons_per_mol_absd)))

# Array: absorbing 0,1,2,...,Nphot_max photons.
Nphot_all   = np.arange(0,Nphot_max+1)
# Initialize array: fraction absorbing 0,1,2,... photons
Nphot_perc  = np.zeros(len(Nphot_all))

# Iterate over # absorbed photons & calculate correponding fraction of molecules
for Nphot in Nphot_all:#range(0,Nphot_max+1):
    
    # Find indices of layers where each molecule absorbs ~Nphot photons
    indices = (pts_photons_per_mol_absd>=Nphot-0.5) & (pts_photons_per_mol_absd<Nphot+0.5)
    
    if any(indices):
        
        # Fraction of molecules absorbing Nphot photons is approximated by the
        # ratio of the total layer thickness absorbing ~Nphot photons to the
        # total sample thickness
        begin               = min(pts[indices])
        end                 = max(pts[indices])
        frac                = (end-begin)/L_sample
        Nphot_perc[Nphot]   = frac
        print('{a:.3f} percent absorb {b:.0f} photons.'.format(a=frac*100,b=Nphot))
        
        del begin, end, frac
        
    else:
        
        print('{a:.3f} percent absorb {b:.0f} photons.'.format(a=0.000,b=Nphot))
        Nphot_perc[Nphot] = 0
        
    del indices
        
#%% Plot Figure 2
# Photon regimes, plotting the fraction of molecules absorbing 0,1,2,... photons.

plt.figure()
plt.plot(Nphot_all,Nphot_perc,'x-', label= 'Eps = {} l/(mol cm)'.format(epsilon/100))
plt.xlabel('Number of absorbed photons/molecule')
plt.ylabel('Fraction of molecules')
plt.suptitle('Dominant photon regimes',y=1.0,fontsize=14)
plt.title(r'Laser: E={e:.2f}$\rm\mu$J, {a:.0f}$\rm\mu$m offset, w={w:.0f}$\rm\mu$m, $\rm\lambda$={l:.0f}nm. {s} xtal parameters'.format(e=E_pulse*1e6,a=offset*1e6,\
          w=r_beam*1e6, l=wavel*1e9, s=moltype ))
plt.legend(loc='best')
plt.grid()
if savefig: plt.savefig(savebasename+'_photon-regimes.'+savefigas, transparent=True)


#%% save input and output to .txt file (if savedata=True)

if savedata:
    
    fo = open(savebasename+'.txt', 'w')
    
    fo.write('########## Result for {:s} ##########\n'.format(moltype))
             
    fo.write('\n***Sample properties***')
    fo.write('\nConcentration of absorbing molecule = {:.3e} mol/l'.format(conc_mol))
    fo.write('\nPenetration depth = {:.1e} m'.format(d_pen))
    
    fo.write('\n\n***Pump characteristics***')
    fo.write('\nPeak fluence = {:.2e} J/m^2 = {:.2e} mJ/cm^2 = {:.2e} uJ/mm^2'.format(F_peak, F_peak/10, F_peak))
    fo.write('\nPeak power density = {:.3e} W/m^2 = {:.3e} GW/cm^2'.format(P_peak, P_peak*1e-13))
    if offset != 0:
        fo.write('\nFluence and power density at the interaction point (IP) reduced to {:.2e} % of peak values due to offset between pump and probe.'.format(100.*offsetred))
        fo.write('\nFluence at IP = {:.2e} J/m^2 = {:.2e} mJ/cm^2 = {:.2e} uJ/mm^2'.format(F_IP, F_IP/10, F_IP))
        fo.write('\Power density at IP = {:.3e} W/m^2 = {:.3e} GW/cm^2'.format(P_IP, P_IP*1e-13))
    
    fo.write('\n\n***Photon absorption distribution***')
    fo.write('\nAverage number of photons absorbed/molecule = {:.2e}  '.format(np.average(pts_photons_per_mol_absd)))
    for Nphot in Nphot_all:
        fo.write('\n{a:.1f} % absorb {b:.0f} photons.'.format(a=Nphot_perc[Nphot]*100,b=Nphot))
        
    fo.write('\n\n\n\n########## Input ##########\n')
    fo.write('\n***Laser parameters***')
    fo.write('\nPulse energy = {:.2e} J'.format(E_pulse))
    fo.write('\nBeam diameter (FWHM) = {:.2e} m'.format(d_beam/1.699))
    fo.write('\nBeam diameter (1/e2) = {:.2e} m'.format(d_beam))
    fo.write('\nSpatial pump-probe offset = {:.2e} m'.format(offset))
    fo.write('\nPulse duration = {:.2e} s'.format(T_pulse))
    fo.write('\nWavelength = {:.2e} m'.format(wavel))
    
    fo.write('\n\n***Sample characteristics***')
    fo.write('\nAbsorption coeffcient = {:.0f} l/(mol cm)'.format(epsilon/100))
    fo.write('\nSample thickness = {:.2e}'.format(L_sample))
    if use_param == 'conc':
        fo.write('\nConcentration of absorbing molecule = {:.3e} mol/l'.format(conc_mol))
    else:
        fo.write('\nUnit cell axis length: a = {:.0f} Ang, b = {:.0f} Ang, c = {:.0f} Ang'.format(a*1e10,b*1e10,c*1e10))
        fo.write('\nUnit cell angles: alpha = {:.0f} deg, beta = {:.0f} deg, gamma = {:.0f} deg'.format(np.rad2deg(alpha),np.rad2deg(beta),np.rad2deg(gamma)))
        fo.write('\n{:.0f} absorbing molecules/unit cell'.format(N_mol_UC))    
    
    fo.close()
    