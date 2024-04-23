"""
Created Oct 7 2020
author: Hans Mohrmann (jkcm@uw.edu)

some utility functions for CSET and MERRA projects
"""

#Standard library
import collections
import datetime as dt
from functools import lru_cache
import math

#Specials
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytz
from scipy.special import erf
from scipy import stats 
from functools import lru_cache



############################################################################
#
#                general use functions and variables
#
############################################################################


def get_lon_prime(lat, lon, lon0=-140, lat0=30):
        lonp = lon0 + 0.8*(lon-lon0) + 0.4*(lat-lat0)
        return lonp
    
    
def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1
    else:
        return idx
    
    
lev_map = {'1': 0.0100, '2': 0.0200, '3': 0.0327, '4': 0.0476,
           '5': 0.0660, '6': 0.0893, '7': 0.1197, '8': 0.1595,
           '9': 0.2113, '10': 0.2785, '11': 0.3650, '12': 0.4758,
           '13': 0.6168, '14': 0.7951, '15': 1.0194, '16': 1.3005,
           '17': 1.6508, '18': 2.0850, '19': 2.6202, '20': 3.2764,
           '21': 4.0766, '22': 5.0468, '23': 6.2168, '24': 7.6198,
           '25': 9.2929, '26': 11.2769, '27': 13.6434, '28': 16.4571,
           '29': 19.7916, '30': 23.7304, '31': 28.3678, '32': 33.8100,
           '33': 40.1754, '34': 47.6439, '35': 56.3879, '36': 66.6034,
           '37': 78.5123, '38': 92.3657, '39': 108.6630, '40': 127.8370,
           '41': 150.3930, '42': 176.9300, '43': 208.1520, '44': 244.8750,
           '45': 288.0830, '46': 337.5000, '47': 375.0000, '48': 412.5000,
           '49': 450.0000, '50': 487.5000, '51': 525.0000, '52': 562.5000,
           '53': 600.0000, '54': 637.5000, '55': 675.0000, '56': 700.0000,
           '57': 725.0000, '58': 750.0000, '59': 775.0000, '60': 800.0000,
           '61': 820.0000, '62': 835.0000, '63': 850.0000, '64': 865.0000,
           '65': 880.0000, '66': 895.0000, '67': 910.0000, '68': 925.0000,
           '69': 940.0000, '70': 955.0000, '71': 970.0000, '72': 985.0000}


pres_map = {}
for k, v in lev_map.items():
    pres_map[v] = int(k)

def get_MERRA_level(pressure):
    a, b = zip(*[(float(k), v) for k, v in lev_map.items()])
    levels = sorted(a)
    pressures = sorted(b)

    return(interp1d(pressures, levels)(pressure))

def MERRA_lev(lev, invert=False, lev_map=lev_map):
    if invert:
        pres_map = {}
        for k, v in lev_map.items():
            pres_map[str(v)] = int(k)
        lev_map = pres_map
    if isinstance(lev, collections.Iterable):
        pres = [lev_map[str(int(i))] for i in lev]
    else:
        pres = lev_map[int(float(str(lev)))]
    return pres

merra_species_dict_colarco = {

      # NOTE: SST Rm is a guess, based on Re and lognormal
          # NOTE: Dust Rm is a guess, based on Re and lognormal, not great

        # NOTE: SST lower not included
    'OCPHILIC': dict(dist='trunc_lognormal', density=1800, geometric_std_dev=2.20, mode_radius=0.0212, upper=0.3),
    'OCPHOBIC': dict(dist='trunc_lognormal', density=1800, geometric_std_dev=2.20, mode_radius=0.0212, upper=0.3),
    'BCPHILIC': dict(dist='trunc_lognormal', density=1800, geometric_std_dev=2.00, mode_radius=0.0118, upper=0.3),
    'BCPHOBIC': dict(dist='trunc_lognormal', density=1800, geometric_std_dev=2.00, mode_radius=0.0118, upper=0.3),
    'SO4': dict(dist='trunc_lognormal', density=1700, geometric_std_dev=2.03, mode_radius=0.0695, upper=0.3),
    'DU001': dict(dist='power_special', effective_radius=0.73, density=2500, geometric_std_dev=2.00, mode_radius=0.220),  # weird bin
    'DU002': dict(dist='power', effective_radius=1.4, density=2650, geometric_std_dev=2.00, mode_radius=0.421, upper=1.8, lower=1.0),
    'DU003': dict(dist='power', effective_radius=2.4, density=2650, geometric_std_dev=2.00, mode_radius=0.7220, upper=3.0, lower=1.8),
    'DU004': dict(dist='power', effective_radius=4.5, density=2650, geometric_std_dev=2.00, mode_radius=1.3540, upper=6.0, lower=3.0),
    'DU005': dict(dist='power', effective_radius=8.0, density=2650, geometric_std_dev=2.00, mode_radius=2.4068, upper=10.0, lower=6.0),
    'SS001': dict(dist='trunc_MG', density=2200, geometric_std_dev=2.03, mode_radius=0.023, upper=0.1, lower=0.03),
    'SS002': dict(dist='trunc_MG', density=2200, geometric_std_dev=2.03, mode_radius=0.090, upper=0.5, lower=0.1),
    'SS003': dict(dist='trunc_MG', density=2200, geometric_std_dev=2.03, mode_radius=0.090, upper=1.5, lower=0.5),
    'SS004': dict(dist='trunc_MG', density=2200, geometric_std_dev=2.03, mode_radius=0.805, upper=5.0, lower=1.5),
    'SS005': dict(dist='trunc_MG', density=2200, geometric_std_dev=2.03, mode_radius=2.219, upper=10.0, lower=5.0)
    }




############################################################################
#
#                functions for CAMS conversions
#
############################################################################


era_name_map = {'aermr01': 'SS001', 'aermr02': 'SS002', 'aermr03': 'SS003', 'aermr04': 'DU001', 'aermr05': 'DU002', 'aermr06': 'DU003', 
                'aermr07': 'OCPHILIC', 'aermr08': 'OCPHOBIC', 'aermr09': 'BCPHILIC', 'aermr10': 'BCPHOBIC', 'aermr11': 'SO4'}


Muskatel_CAMS_lookup = {
'aermr01': [9.5491368106e-18, 9.5491368106e-18, 9.5491368106e-18, 9.5491368106e-18, 2.8508222886e-17, 3.5552048177e-17, 4.4141680140e-17, 5.5592627265e-17, 7.5022872211e-17, 9.2401351367e-17, 1.2548908641e-16, 2.2733094176e-16],
'aermr02': [4.2343711220e-16, 4.2343711220e-16, 4.2343711220e-16, 4.2343711220e-16, 1.2643975497e-15, 1.5768574457e-15, 1.9580637788e-15, 2.4653863535e-15, 3.3268891197e-15, 4.0977089657e-15, 5.5658632784e-15, 1.0083678901e-14],
'aermr03': [1.1890425978e-15, 1.1890425978e-15, 1.1890425978e-15, 1.1890425978e-15, 3.5504471910e-15, 4.4278271135e-15, 5.4982035375e-15, 6.9229124851e-15, 9.3420826653e-15, 1.1506545624e-14, 1.5628979584e-14, 2.8314815079e-14],
'aermr04': [2.8694347373e-16, 2.8694347373e-16, 2.8694347373e-16, 2.8694347373e-16, 2.8694347373e-16, 2.8694347373e-16, 2.8694347373e-16, 2.8694347373e-16, 2.8694347373e-16, 2.8694347373e-16, 2.8694347373e-16, 2.8694347373e-16],
'aermr05': [4.7291076037e-16, 4.7291076037e-16, 4.7291076037e-16, 4.7291076037e-16, 4.7291076037e-16, 4.7291076037e-16, 4.7291076037e-16, 4.7291076037e-16, 4.7291076037e-16, 4.7291076037e-16, 4.7291076037e-16, 4.7291076037e-16],
'aermr06': [1.5569834525e-15, 1.5569834525e-15, 1.5569834525e-15, 1.5569834525e-15, 1.5569834525e-15, 1.5569834525e-15, 1.5569834525e-15, 1.5569834525e-15, 1.5569834525e-15, 1.5569834525e-15, 1.5569834525e-15, 1.5569834525e-15],
'aermr07': [1.3411058245e-18, 1.5524975799e-18, 1.7367708373e-18, 1.9869071216e-18, 2.2599763045e-18, 2.5569842145e-18, 2.8789366803e-18, 3.2268395306e-18, 4.0045196995e-18, 4.8980713505e-18, 5.9155411131e-18, 8.9112054687e-18],
'aermr08': [4.9135101849e-19, 4.9135101849e-19, 4.9135101849e-19, 4.9135101849e-19, 4.9135101849e-19, 4.9135101849e-19, 4.9135101849e-19, 4.9135101849e-19, 4.9135101849e-19, 4.9135101849e-19, 4.9135101849e-19, 4.9135101849e-19],
'aermr09': [5.9773588135e-20, 6.9195324740e-20, 7.7408384085e-20, 8.8556952333e-20, 1.0072763231e-19, 1.1396525317e-19, 1.2831464408e-19, 1.4382063418e-19, 1.7848172858e-19, 2.1830716962e-19, 2.6365559054e-19, 3.9717161997e-19],
'aermr10': [5.9773588135e-20, 5.9773588135e-20, 5.9773588135e-20, 5.9773588135e-20, 5.9773588135e-20, 5.9773588135e-20, 5.9773588135e-20, 5.9773588135e-20, 5.9773588135e-20, 5.9773588135e-20, 5.9773588135e-20, 5.9773588135e-20],
'aermr11': [2.8658111470e-18, 2.8658111470e-18, 2.8658111470e-18, 2.8658111470e-18, 4.5781613964e-18, 5.2038774089e-18, 6.0382576467e-18, 7.2566325298e-18, 9.3848411035e-18, 1.1325129380e-17, 1.4889880922e-17, 2.5975608980e-17]
}


era_species_dict_reddy = {
    
    # getting these vals from https://atmosphere.copernicus.eu/sites/default/files/FileRepository/Resources/Documentation/Radiative_Forcing/CAMS74_2016SC1_D74.1-1_201612_Documentation_v1.pdf.
    # everything points to this paper:
    #     Reddy M.S., Boucher O., Bellouin N., Schulz M., Balkanski Y., Dufresne J.-L., and Pham M., 2005,
    # Estimates of global multicomponent aerosol optical depth and direct radiative perturbation in the Laboratoire de Meteorologie Dynamique general circulation model. J. Geophys. Res.-Atmospheres, 110,
    # D10S16, doi:10.1029/2004JD004757
    #also here: https://gmd.copernicus.org/preprints/gmd-2019-149/gmd-2019-149.pdf
    #only confusion is over sea salt - is the underlying distribution lognormal or bimodally lognomal with fixed concs of 70/cc and 3/cc??
    
    'OCPHILIC': dict(dist='lognormal', density=1760, geometric_std_dev=2.0, mode_radius=0.0355),
    'OCPHOBIC': dict(dist='lognormal', density=1760, geometric_std_dev=2.0, mode_radius=0.0355),
    'BCPHILIC': dict(dist='lognormal', density=1800, geometric_std_dev=2.0, mode_radius=0.0118),
    'BCPHOBIC': dict(dist='lognormal', density=1800, geometric_std_dev=2.0, mode_radius=0.0118),
    'SO4': dict(dist='lognormal', density=1840, geometric_std_dev=2.0, mode_radius=0.0355),
    'DU001': dict(dist='trunc_lognormal', mode_radius=0.29, median_radius=0.135, density=2160, geometric_std_dev=2.0, lower=0.1, upper=0.55), # lower is actually 0.03.
    'DU002': dict(dist='trunc_lognormal', mode_radius=0.29, median_radius=0.704, density=2160, geometric_std_dev=2.0, lower=0.55, upper=0.9),
    'DU003': dict(dist='trunc_lognormal', mode_radius=0.29, median_radius=4.4, density=2160, geometric_std_dev=2.0, lower=0.9, upper=20),
    'SS001': dict(dist='CAMS_SS', mode_radius=[0.1992, 1.992], density=1183, geometric_std_dev=[1.9,2.0], lower=0.03, upper=0.5), #lower is actually 0.03
    'SS002': dict(dist='CAMS_SS', mode_radius=[0.1992, 1.992], density=1183, geometric_std_dev=[1.9,2.0], lower=0.5, upper=5.0),
    'SS003': dict(dist='CAMS_SS', mode_radius=[0.1992, 1.992], density=1183, geometric_std_dev=[1.9,2.0], lower=5.0, upper=20),
    }


def mass_to_number_CAMS(mass, RH, air_density, name):
    vol_mean_mass = Muskatel_CAMS_lookup[era_name_map_rev[name]]

    RH_t = [0, 10, 20, 30, 40, 50, 60, 70, 80, 85, 90, 95]

    vfunc = np.vectorize(lambda x: vol_mean_mass[find_nearest(RH_t, x)])

    mean_masses = vfunc(RH)
    # mean_masses = np.array([vol_mean_mass[find_nearest(RH_t, i)] for i in RH])
    n = mass*air_density/mean_masses*1e-6
    return n



def mass_to_number_CAMS_SS(mass, particle_density, air_density, mode_radius, geo_std_dev, lower_lim, upper_lim):


    mu_g = mode_radius
    sigma_g = geo_std_dev

    mu = np.log(mu_g)
    sigma = np.log(sigma_g)

    mu_vol = mu + 3*sigma**2
    mu_vol_g = np.exp(mu_vol)

    novo = 3/(4*np.pi*np.exp(3*mu + 4.5*sigma**2))
    
    if upper_lim:
        lower_lim, upper_lim, sigma, mu_g, mu_vol_g, 

        x=np.linspace(lower_lim,upper_lim,1000)
        pdf = stats.lognorm.pdf(x, s=sigma, loc=0, scale=mu_g) # original distribution
        pdf_vol = stats.lognorm.pdf(x, s=sigma, loc=0, scale=mu_vol_g) # original distribution

        idx = np.logical_and(x<upper_lim, x>=lower_lim)
        frac = np.trapz(pdf[idx], x[idx])/np.trapz(pdf_vol[idx], x[idx])


        novo = novo*frac
    
    n_0 = novo*mass/particle_density*1e18 # per um3 to per m3
    num_cm3 = air_density*n_0*1e-6
    
    return num_cm3



############################################################################
#
#                functions for converting mass to number using erf
#
############################################################################

@lru_cache(maxsize=500)
def get_bounded_lognormal_frac(r_max, r_min, std_dev, mode_radius):
    """
    some notes: mode radius = same as 1/2 of median diameter, then using Zender eq 27
    so mode radius is Zender's D-tilde_n/2
    std dev here is NOT the same as std dev in the Grainger text. std dev should be somewhere around 2. 
    this will give n(rmin,rmax)/n0, or else v(rmin,rmax)/n0, if the volume median radius is provided instead of the median radius.

    """
    n_over_n0 = 0.5*(erf(np.log(r_max/mode_radius)/(np.sqrt(2)*np.log(std_dev))) - erf(np.log(r_min/mode_radius)/(np.sqrt(2)*np.log(std_dev))))
    return n_over_n0

@lru_cache(maxsize=500)
def get_v0_over_n0(r_max, r_min, std_dev, mode_radius):
    vol_mode_radius = mode_radius*np.exp(3*(np.log(std_dev)**2)) # this is D-tilde_v/2 from Zender.
    print(mode_radius)
    print(vol_mode_radius)
    v0_over_n0 = get_bounded_lognormal_frac(r_max, r_min, std_dev, vol_mode_radius)
    return v0_over_n0

def get_n0(mass, density, r_max, r_min, std_dev, mode_radius):
    v0_over_n0 = get_v0_over_n0(r_max, r_min, std_dev, mode_radius)
    n0 = mass*1e18/(density*v0_over_n0)
    return n0

def get_n_subset(n0, r_min, r_max, std_dev, mode_radius):
    """this is for subsetting for a bounded range after n0 has been calculated 
    """
    n_over_n0 = get_bounded_lognormal_frac(r_max, r_min, std_dev, mode_radius)
    n = n_over_n0 * n0
    return n

def get_m_subset(density, n0, r_min, r_max, std_dev, mode_radius):
    vol_mode_radius = mode_radius*np.exp(3*np.log(std_dev)**2) # this is D-tilde_v/2 from Zender.
    v_over_n0 = get_bounded_lognormal_frac(r_max, r_min, std_dev, vol_mode_radius)
    v = n0*v_over_n0/1e18
    m = v*density
    return m



############################################################################
#
#           functions for converting mass to number using integration
#
############################################################################

def mass_to_number(mass, air_density, shape_params):
    if shape_params['dist'] == 'trunc_lognormal':
        return mass_to_number_trunc_lognormal(mass=mass, particle_density=shape_params['density'], mode_radius=shape_params['mode_radius'], 
                                               geo_std_dev=shape_params['geometric_std_dev'], air_density=air_density, 
                                              upper_lim=shape_params['upper'])
    elif shape_params['dist'] == 'lognormal':
        return mass_to_number_lognormal(mass=mass, particle_density=shape_params['density'], mode_radius=shape_params['mode_radius'], 
                                               geo_std_dev=shape_params['geometric_std_dev'], air_density=air_density)
    elif shape_params['dist'] == 'power':
        return mass_to_number_trunc_power(mass=mass, particle_density=shape_params['density'], air_density=air_density, upper_lim=shape_params['upper'], lower_lim=shape_params['lower'])
    elif shape_params['dist'] == 'power_special':
        return mass_to_number_trunc_power_dust_smallest(mass=mass, particle_density=shape_params['density'], air_density=air_density)
    elif shape_params['dist'] == 'trunc_MG':
        return mass_to_number_trunc_MG(mass=mass, particle_density=shape_params['density'], 
                                       air_density=air_density, upper_lim=shape_params['upper'], lower_lim=shape_params['lower'])
    else:
        raise ValueError('shape params dist type not recognized')

        
def mass_to_number_lognormal(mass, particle_density, mode_radius, geo_std_dev, air_density):
    """ Calculates the number concentration for a lognormal mode, given mass and mode parameters
    
    Args:
        mass::float
            mode total aerosol mass, in kg/kg (sanity check: something like 1e-10, 1e-11)
        particle_density::float
            particle density in kg/m3 (sanity check: something like 1e3)
        modal_radius::float
            mode radius r_m, in m (meters) (sanity check: something like 1e-6, around a micron)
        geo_std_dev::float
            mode geometric standard deviation (sanity check: something around 2e-6, ie 2 microns)
            
    
    TODO doctstring plz
    """
    
    @lru_cache(maxsize=500)
    def num_per_mass(geo_std_dev, particle_density, mode_radius):
        exp = np.exp(-4.5*np.log(geo_std_dev)**2)
        num_per_mass = 3*exp/(particle_density*4*np.pi*mode_radius**3)*1e18 # per um3 to per m3
        return num_per_mass
    
    num_kg = mass*num_per_mass(geo_std_dev, particle_density, mode_radius)
    num_cm3 = air_density*num_kg*1e-6
    return num_cm3

def mass_to_number_trunc_MG(mass, particle_density, air_density, upper_lim, lower_lim):
    # using eqn 2 from here: https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2003GB002079
    
  
    @lru_cache(maxsize=500)
    def n0_per_v0(lower_lim, upper_lim):
        r=np.linspace(lower_lim,upper_lim,100)
        Theta = 30
        A = 4.7*(1+Theta*r)**(-0.017*r**-1.44)
        B = (0.433 - np.log10(r))/0.433
        dfdr = r**(-A)*(1+0.057*r**3.45)*10**(1.607*np.exp(-B**2))
        dfdn = dfdr*(r**3)
        n0_per_v0 = (3/(4*np.pi))*(np.trapz(dfdr, x=r)/np.trapz(dfdn, x=r))    
        return n0_per_v0
        
    n_0 = n0_per_v0(lower_lim, upper_lim)*mass/particle_density*1e18 # per um3 to per m3
    num_cm3 = air_density*n_0*1e-6
    
    return num_cm3
    



def mass_to_number_trunc_lognormal(mass, particle_density, mode_radius, geo_std_dev, air_density, upper_lim, lower_lim=0.1):
    #get the mass to number the old way, by integrating.
    @lru_cache(maxsize=500)
    def n0_per_v0(mode_radius, geo_std_dev, upper_lim, lower_lim):
        mu_g = mode_radius
        sigma_g = geo_std_dev

        mu = np.log(mu_g)
        sigma = np.log(sigma_g)

        mu_vol = mu + 3*sigma**2
        mu_vol_g = np.exp(mu_vol)

        x=np.linspace(0,upper_lim*2,1000)

        n0_per_v0_full = 3/(4*np.pi*np.exp(3*mu + 4.5*sigma**2))

        pdf = stats.lognorm.pdf(x, s=sigma, loc=0, scale=mu_g) # original distribution
        pdf_vol = stats.lognorm.pdf(x, s=sigma, loc=0, scale=mu_vol_g) # original distribution

        idx = np.logical_and(x<upper_lim, x>=lower_lim)
        n0_per_v0 = n0_per_v0_full*np.trapz(pdf[idx], x[idx])/np.trapz(pdf_vol[idx], x[idx])
        return n0_per_v0

    n_0 = n0_per_v0(mode_radius, geo_std_dev, upper_lim, lower_lim)*mass/particle_density*1e18 # per um3 to per m3
    num_cm3 = air_density*n_0*1e-6
    
    return num_cm3

# def mass_to_number_trunc_lognormal_bimodal(mass, particle_density, mode_radius, geo_std_dev, mode_radius_2, geo_std_dev_2, air_density, upper_lim=None, lower_lim=0.1):
#     mu_g = mode_radius
#     sigma_g = geo_std_dev

#     mu = np.log(mu_g)
#     sigma = np.log(sigma_g)

#     mu_vol = mu + 3*sigma**2
#     mu_vol_g = np.exp(mu_vol)

#     x=np.linspace(0,upper_lim*2,1000)

#     novo = 3/(4*np.pi*np.exp(3*mu + 4.5*sigma**2))
    
#     if upper_lim:
#         pdf = stats.lognorm.pdf(x, s=sigma, loc=0, scale=mu_g) # original distribution
#         pdf_vol = stats.lognorm.pdf(x, s=sigma, loc=0, scale=mu_vol_g) # original distribution

#         idx = np.logical_and(x<upper_lim, x>=lower_lim)
#         novo = novo*np.trapz(pdf[idx], x[idx])/np.trapz(pdf_vol[idx], x[idx])
    
#     n_0 = novo*mass/particle_density*1e18 # per um3 to per m3
#     num_cm3 = air_density*n_0*1e-6
    
#     return num_cm3




def mass_to_number_trunc_power(mass, particle_density, upper_lim, lower_lim, air_density):
    
    @lru_cache(maxsize=500)
    def n0_per_v0(particle_density, lower_lim, upper_lim):
        novo = 9/particle_density*(lower_lim**-3 - upper_lim**-3)/(4*np.pi*np.log(upper_lim/lower_lim))*1e18 # per um3 to per m3
        return novo
    
    novo = n0_per_v0(particle_density, lower_lim, upper_lim)
    n_0 = mass*novo
    num_cm3 = air_density*n_0*1e-6
    return num_cm3


    
def mass_to_number_trunc_power_dust_smallest(mass, particle_density, air_density):

    mass_bins = [mass*0.009, mass*0.081, mass*0.234, mass*0.676]
    lowers = [0.1, 0.18, 0.3, 0.6]
    uppers = [0.18, 0.3, 0.6, 1.0]
    
    numbers = [mass_to_number_trunc_power(mass=mass_i, particle_density=particle_density, upper_lim=upper_i, lower_lim=lower_i, air_density=air_density) 
               for (mass_i, lower_i, upper_i) in zip(mass_bins, lowers, uppers)]
    
    return sum(numbers)

    
    
    
@lru_cache(maxsize=500)
def sum_trunc_lognormal(Dg, sigma, Dmax=None, Dmin=None):
    """code from Matt Wyant, following Zender (https://patarnott.com/pdf/SizeDistributions.pdf)
    """
    sqrt2 = np.sqrt(2)
    from scipy.special import erf
    s = 0.
    if Dmax is not None:
        s += erf(np.log(Dmax/Dg)/(sqrt2*np.log(sigma)))       
    else:
        s += 1.
        
    if Dmin is not None:
        s -= erf(np.log(Dmin/Dg)/(sqrt2*np.log(sigma)))
    else:
        s += 1.    
    
    return s/2
