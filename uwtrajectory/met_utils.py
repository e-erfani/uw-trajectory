# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 16:42:46 2016

@author: jkcm
"""
import numpy as np
np.warnings.filterwarnings('ignore')

# import warnings
from scipy import integrate
# warnings.simplefilter("ignore")
p0 = 1000.  # reference pressure, hPa
Rdry = 287.  # gas const for dry air, J/K/kg
Rvap = 461.  # gas const for water vapor, J/K/kg
eps = Rvap/Rdry - 1
cp = 1004.  # cp_dry, specific heat of dry air at const pressure, J/K/kg
g = 9.81   # grav acceleration at sea level, m/s2
lv = 2.5*10**6  # latent heat of vaporization at 0C, J/kg






def smooth(x, window_len=11, window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal. Courtesy of scipy-Cookbook

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett',
            blackman', flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
    """

    if isinstance(x, list):
        x = np.array(x)
    if window_len % 2 == 0:
        raise ValueError("please use odd-numbered window_len only.")
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")
    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
    if window_len < 3:
        return x

    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is one of 'flat', 'hanning', "
                         "hamming', 'bartlett', 'blackman'")

    s = np.r_[x[window_len-1:0:-1], x, x[-1:-window_len:-1]]
    if window == 'flat':  # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')

    y = np.convolve(w/w.sum(), s, mode='valid')
    return y[int(window_len/2):-int(window_len/2)]

def qvs_from_p_T(p, T):
    """p in Pa, T in K. return is in kg/kg
    """
    es = 611.2*np.exp(17.67*(T-273.15)/(T-29.65))
    qvs = 0.622*es/(p-0.378*es)
    return qvs
    
def qv_from_p_T_RH(p, T, RH):
    """p in Pa, T in K, Rh in pct. return is in kg/kg
    """
#     es = 611.2*np.exp(17.67*(T-273.15)/(T-29.65))
    qvs = qvs_from_p_T(p, T)
    rvs = qvs/(1-qvs)
    rv = RH/100. * rvs
    qv = rv/(1+rv)
    return qv

def tvir_from_T_w(T, w):
    """T in L, w in kg/kg"""
    t_vir = T*(1+0.61*w)
    return t_vir


def theta_from_p_T(p, T, p0=1000):
    theta = T * (p0/p)**(Rdry/cp)
    return theta

def get_liquid_water_theta(temp, theta, q_l):
    """temp = air temp (K) theta = pot temp, q_l = liquid water MR"""
    theta_l = theta - (theta*lv*q_l/(temp*cp*1000))
    return theta_l

def density_from_p_Tv(p, Tv):
    return p/(Rdry*Tv)
    
#def thetae_from_theta_w_T(theta, w, T):
#    """theta = pot temp in K, w = mr in kg/kg, T = temp in K"""
#    returnb theta*np.exp(lv*)
    
def thetae_from_t_tdew_mr_p(t, tdew, mr, p):
    """From Bolton, 1980
    t, tdew in K, mr in kg/kg, p in Pa
    """
    t_lcl = 56 + 1/((1/(tdew-56))+(np.log(t/tdew)/800))
    e = p*mr/(mr + 0.622)
    K = 0.2854  # Rdry/cp
    theta_lcl = t*(100000/(p-e))**K*(t/t_lcl)**(0.28*mr)
    theta_e = theta_lcl*np.exp((3036/t_lcl - 1.78)*mr*(1+0.488*mr))
    return theta_e

def get_LCL(t, t_dew, z):
    raise NotImplementedError('use calculate_LCL instead')
#     if np.any(t_dew > t):
#                 t_dew = np.minimum(t, t_dew)
# #         raise ValueError('dew point temp above temp, that\'s bananas')
#     return z + 125*(t - t_dew)

def get_virtual_dry_static_energy(T, q, z):
    return cp*T*(1+eps*q) + g*z


def get_moist_adiabatic_lapse_rate(T, p): #K and hPa
    
    es = 611.2*np.exp(17.67*(T-273.15)/(T-29.65)) # Bolton formula, es in Pa
    qs = 0.622*es/(p*100-0.378*es)
    num = 1 + lv*qs/(Rdry*T)
    denom = 1 + lv**2*qs/(cp*Rvap*T**2)
    gamma = g/cp*(1-num/denom)
    return gamma

def get_moist_adiabat(t, p, p_arr):
    pass
    

def get_Ri_profile(u, v, q, T, z, T0=None, z0=None, q0=None, filt=False):
    if filt:
        u = smooth(u, window_len=filt)
        v = smooth(v, window_len=filt)
        q = smooth(q, window_len=filt)
        T = smooth(T, window_len=filt)
        z = smooth(z, window_len=filt)

    if T0 is None:
        T0 = T[0]
    if z0 is None:
        z0 = z[0]
    if q0 is None:
        q0 = q[0]

    del_U_sq = u**2 + v**2
    sv_0 = get_virtual_dry_static_energy(T0, q0, z0)
    sv_hbl = get_virtual_dry_static_energy(T, q, z)
    Ri_b = z*(2*g*(sv_hbl - sv_0))/(del_U_sq*(sv_hbl + sv_0 - g*z0 - g*z))
    return Ri_b


def Ri_pbl_ht(u, v, q, T, z, T0=None, z0=None, q0=None, smooth=False):
    indx = np.flatnonzero(z < 40)[-1]  # avg all values below this for sfc values
    if T0 is None:
        T0 = np.nanmean(T[:indx])
    if z0 is None:
        z0 = np.nanmean(z[:indx])
    if q0 is None:
        q0 = np.nanmean(q[:indx])

    Ri = get_Ri_profile(u, v, q, T, z, T0, z0, q0, smooth)
    try:
        indx = np.flatnonzero(np.array(Ri) > 0.25)[0]
        z_pbl = z[indx]
        if z_pbl > 4000:
            raise IndexError
    except IndexError:
        return 0, float('nan')
    return indx, z_pbl


def RH_fancy_pblht_1d( z, RH):
    """
    at least 80% of the time, one could 
    identify a 'RH inversion base' as the altitude of max RH for which RH(zi 
    + 300 m) - RH(zi) < -0.3.  If such a layer does not exist below 4 km or 
    the top of the sounding, we say an inversion is not present.
    """
    if z.shape != RH.shape:
        raise ValueError('z and RH must have the same shape')
    if len(z.shape) != 1:  # height axis
        raise ValueError('data has an invalid number of dimensions')
    if not (z < 4000).any():
        raise ValueError('must have data below 4000m')
        
    z = smooth(z, window_len=5)
    RH = smooth(RH, window_len=5)
#    
#    flipped = False
##    if np.all(z[100:-100] != sorted(z[100:-100])):  # not in ascending order
##        if np.all(z[100:-100:-1] == sorted(z[100:-100])):  # in descenting order
#    if z[0] > z[-1]: # starts off higher
#        if True:
#            z = z[::-1]
#            theta = theta[::-1]
#            flipped = True
#        else:
#            raise ValueError("data not in ascending or descending order")
#   
# %%
    z_p300 = z + 300
#    i_arr = np.empty_like(z)
    RH_diff  = np.empty_like(z)
    for n, zpn in enumerate(z_p300):
        i_arr = np.abs(z - zpn).argmin()
        RH_diff[n] = RH[n] - RH[i_arr]
#        print(RH_diff)
    inv_cands = np.where(RH_diff > 30)[0]
    RH_cands = RH[inv_cands]
    if len(RH_cands) == 0:
        if np.all(np.isnan(RH_diff)):
            return {'z': np.nan, 'RH': np.nan, 'i': np.nan, 'inversion': False}
        biggest_drop = np.nanargmax(RH_diff)
        z_drop=z[biggest_drop]
        RH_drop = RH[biggest_drop]
        return {'z': z_drop, 'RH': RH_drop, 'i': biggest_drop, 'inversion': False}
    best = np.argmax(RH_cands)
    best_index = inv_cands[best]
    RH_bot = RH[best_index]
    z_bot = z[best_index]
    return {'z': z_bot, 'RH': RH_bot, 'i': best_index, 'inversion': True}
    
# %%
#    inv_cands = RH
    
    
    
    
#    for n,i in enumerate(i_arr):
#        print(z[n] - z[i])
    

def RH_50_pblht_1d(z, RH):
    """
    Given z and RH, return height where RH drops below 50
    """
    if z.shape != RH.shape:
        raise ValueError('z and RH must have the same shape')
    if len(z.shape) != 1:  # height axis
        raise ValueError('data has an invalid number of dimensions')
    if not (z < 4000).any():
        raise ValueError('must have data below 4000m')
        
    z = smooth(z, window_len=5)
    RH = smooth(RH, window_len=5)
    
    nanfrac = sum(np.isnan(RH)/len(RH))
    
    if np.nanmin(RH) > 50 or np.all(np.isnan(RH)) or np.nanmax(RH) < 50:
#        print('no inversion')
        return {"z": np.nan, "i": np.nan, 'inversion': False}
    
    i_min = np.where(z == z[RH < 50].min())[0][0]
    z_i = z[i_min]
#    print(z_i)
#    RH_i = RH[i_min]
    return {"z": z_i, "i":i_min, 'inversion': True}
    
    pass


def Peter_inv(z, rh, theta, polyfit_range=1):
    
    idx = z<5000
    
    THETA = theta[idx]
    RH = rh[idx]
    z = z[idx]
    
    dz = np.diff(z)#, prepend=(2*z[0]-z[1]))
    d_THETA = np.diff(THETA)
    d_RH = np.diff(RH)
#     d_THETA = THETA[:-1] - THETA[1:]
#     d_RH = RH[:-1] - RH[1:]
#     z  = zz[i,:]  if len(zz.shape) > 1 else zz
#     dz = dzz[i,:] if len(zz.shape) > 1 else dzz
    dTHETA_dz = d_THETA / dz
    dRH_dz = d_RH / dz
    dTHETA_dz[dTHETA_dz < 0] = 0
    dRH_dz[dRH_dz > 0] = 0    
    func = dTHETA_dz * dRH_dz
    indx = np.argmin(func)                
    # approximate func as a parabola around the minimum and find the height where that parabola is minimized. 
    # This will allow the inversion height to vary continuously as the input profiles change.
    # inversion_test is defined at midpoints of grid
    zavg = 0.5 * (z[:-1] + z[1:])
    rnge = range(indx-polyfit_range, min(indx+1+polyfit_range, len(func))) # edited to add min() to avoid overrun
    # we define the parabola, converting from m to km.
    try:
        pp = np.polyfit(1e-3 * zavg[rnge], func[rnge], 2)
        # take the derivative of the parabola in coeffient space.
        pp_prime = np.array([2 * pp[0], pp[1]]) # this is its derivative
        # find the zero-crossing of the derivative. This is the inversion height in meters
        z_inv = -1e3 * pp_prime[1] / pp_prime[0]         
    except np.linalg.LinAlgError as e:
        print(z[indx])
        print(zavg[rnge])
        print(func[rnge])
        print(dTHETA_dz[rnge])
        print(dRH_dz[rnge])
        print(z[rnge])
        print(rh[rnge])
        print(theta[rnge])

        return z[indx]
    return z_inv

def Peter2_inv(z, rh, theta):
    
    idx = z<3000
    rh = rh[idx]
    theta = theta[idx]
    z = z[idx]
    
    
    grad_rh = np.gradient(rh, z)
    grad2_rh = np.gradient(grad_rh, z)
    grad2_rh[np.where(grad2_rh>0)] = 0 # looking for a strong decrease in rh grad
    grad2_rh[np.where(grad_rh>0)] = 0 # must be negative grad in rh
    grad_theta = np.gradient(theta, z)
    grad2_theta = np.gradient(grad_theta, z)
    grad2_theta[np.where(grad2_theta<0)] = 0 #looking for a strong increase in theta grad
    grad2_theta[np.where(grad_theta<0)] = 0 #must be positive grad in theta
    grad2_prod = grad2_rh*grad2_theta
    
    return(z[np.argmin(grad2_prod)])
    
    
def moist_static_energy(t, z, q):
    return cp*t + g*z + lv*q


def get_inversion_layer_2d(z, t, p, axis=0, handle_nans=False):
    res_dict = {key: np.empty(z.shape[axis]) for key in ["z_top", "z_mid", "z_bot",  "i_top", "i_mid", "i_bot",
                                                        "t_above_inv", "t_below_inv", "d_t_inv"]}

    for i,(z_i,t_i,p_i) in enumerate(zip(z,t,p)):
        try:
            res = quick_inversion(z_i,t_i,p_i)
        except ValueError as e:
            if handle_nans:
                res = {"z_top": np.nan, "z_mid": np.nan, "z_bot": np.nan, 
                       "i_top": np.nan, "i_mid": np.nan, "i_bot": np.nan}
            else:
                import matplotlib.pyplot as plt
                plt.plot(t_i[z_i<4000], z_i[z_i<4000])
                raise e
        for key, value in res.items():
            res_dict[key][i] = value
    return res_dict
    
def quick_inversion(z, t, p, smooth_t=False): # z in meters, t in K, p in hPa

    #getting layers
    gamma_moist = get_moist_adiabatic_lapse_rate(T=t, p=p)*1000
    if smooth_t:
        gamma = -np.gradient(smooth(t, window_len=31), z)*1000
    else: 
        gamma = -np.gradient(t, z)*1000
    gamma[np.gradient(z)>-1] = np.nan
    gamma[z<330] = np.nan
#     gamma[z>3000] = np.nan    
    gamma[np.abs(gamma)>100] = np.nan
    gamma_diff = (gamma-gamma_moist)/1000

    return_dict = {"z_top": np.nan, "z_mid": np.nan, "z_bot": np.nan, 
                   "i_top": np.nan, "i_mid": np.nan, "i_bot": np.nan,
                   "t_above_inv": np.nan, "t_below_inv": np.nan, "d_t_inv": np.nan}
    
    #inversion center
    #quick hack to only look below 3km, but allow the inversion top to be above 3km:
    gamma_lower = gamma.copy()
    gamma_lower[z>3000] = np.nan
    i_mid = np.nanargmin(gamma_lower)  # middle of inversion is where the lapse rate is the strongest    
    if np.isnan(i_mid):
        print('no i_mid')
        return buncha_nans
    z_mid = z[i_mid]
    return_dict['i_mid'] = i_mid
    return_dict['z_mid'] = z_mid
    
    #inversion base
    max_gap = gamma[i_mid] - gamma_moist[i_mid]
    try: # first way to get the inversion base: where the lapse rate is sufficiently close to the moist adiabat again
        z_bot = np.max(z[np.logical_and(z<z[i_mid], gamma-gamma_moist>max_gap/4)])
    except ValueError as v: # no crossing of the max_gap/4 line go for smallest gap below zmid
        cands = z<z[i_mid] # second way to get the inversion base: wherever it gets closest.
        if not np.any(cands):
            raise ValueError("no values below inversion middle!")            
        z_bot = z[cands][np.argmin(np.abs(gamma[cands]-gamma_moist[cands]))]
    i_bot = np.argwhere(z==z_bot)[0][0]
    return_dict['i_bot'] = i_bot
    return_dict['z_bot'] = z_bot

    #inversion top
    top_candidates = np.logical_and(z>z[i_mid], gamma-gamma_moist>max_gap/4)
    if np.any(top_candidates):
        z_top = np.min(z[top_candidates]) # first way to get inversion top: where the lapse rate is sufficiently close to the moist adiabat again
        i_top = np.argwhere(z==z_top)[0][0]
    else: #second way to get inversion top: wherever it gets closest
        cands = z>z[i_mid]
        if not np.any(cands):
            raise ValueError("no values above inversion middle!")
        z_top = z[cands][np.argmin(np.abs(gamma[cands]-gamma_moist[cands]))]
        i_top = np.argwhere(z==z_top)[0][0]
    return_dict['i_top'] = i_top
    return_dict['z_top'] = z_top
    
    t_below_inv = t[i_bot]
    i_inv = np.logical_and(z>z_bot, z<z_top)
    d_t_inv = integrate.trapz(gamma_diff[i_inv], z[i_inv])
    t_above_inv = t_below_inv + d_t_inv
    return_dict['t_above_inv'] = t_above_inv
    return_dict['t_below_inv'] = t_below_inv
    return_dict['d_t_inv'] = d_t_inv
    
    
    return return_dict
    

    
def calc_decoupling_and_inversion_from_sounding(sounding_dict, usetheta=False, get_jumps=True, smooth_t=True):
    #Setting up variables
    z = sounding_dict['GGALT']
    theta = sounding_dict['THETA']
    theta_e = sounding_dict['THETAE']
    qv = sounding_dict['QV']
    t = sounding_dict['ATX']
    if 'PSXC' in sounding_dict.keys():
        p = sounding_dict['PSXC']
    else:
        p = sounding_dict['PSX']
    if not usetheta:
        theta_l = sounding_dict['THETAL']
        ql = sounding_dict['QL']
        if np.all(np.isnan(ql)):
            qt = qv
        else:
            qt = qv + ql
    else:
        theta_l = sounding_dict['THETA']
        qt = qv

    #failing quietly
    buncha_nans = {"d_qt": np.nan, "d_theta_e": np.nan, "d_theta_l": np.nan,
                "alpha_thetal": np.nan, "alpha_qt":np.nan, "alpha_thetae": np.nan,
                "d_q_inv": np.nan, "d_t_inv": np.nan,
                "t_below_inv": np.nan, "t_above_inv": np.nan, "q_below_inv": np.nan, "q_above_inv": np.nan,
                "z_top": np.nan, "z_mid": np.nan, "z_bot": np.nan, "i_top": np.nan, "i_mid": np.nan, "i_bot": np.nan}
        
        
    buncha_nans['lat'] = np.nanmean(sounding_dict['GGLAT'])
    buncha_nans['lon'] = np.nanmean(sounding_dict['GGLON'])
    buncha_nans['lon_p'] =-140 + 0.8*(buncha_nans['lon']+140) + 0.4*(buncha_nans['lat']-30)   

    buncha_nans['time'] = sounding_dict['TIME'][0]
        
    inv_levs = quick_inversion(z, t, p, smooth_t=smooth_t)
    buncha_nans.update(inv_levs)
    z_top, z_mid, z_bot = inv_levs['z_top'], inv_levs['z_mid'], inv_levs['z_bot']
    i_top, i_mid, i_bot = inv_levs['i_top'], inv_levs['i_mid'], inv_levs['i_bot']
    
#     for key, value in inv_levs.items():
#             buncha_nans[key] = value # better with dict.update()?
    
    #jumps in q, t    
#     i_upper = np.logical_and(z<=z_top, z>=z_mid) #this is the upper half of the inversion
#     if np.sum(i_upper) == 0:
#         print("error: no upper inv layer: z_top: {}  z_mid: {}".format(z_top, z_mid))
#         return buncha_nans
#     i_lower = np.logical_and(z>z_bot, z<z_mid) #this is the lower half of the inversion

    q_above_inv = qt[i_top]
    q_below_inv = qt[i_bot]
    d_q_inv = q_above_inv - q_below_inv
    buncha_nans['q_above_inv'] = q_above_inv
    buncha_nans['q_below_inv'] = q_below_inv
    buncha_nans['d_q_inv'] = d_q_inv
    


    #decoupling ests
    upper_25 = z_bot - (z_bot - min(z))/4. #top quarter of the MBL
    u_i = np.logical_and(z > upper_25, z < z_bot)
    lower_25 = min(z) + (z_bot - min(z))/4. #bottom quarter of the MBL
    l_i = np.logical_and(z < lower_25, z > min(z))
    
    ft_base = z_top
    ft_top = ft_base + 500
    l_ft = np.logical_and(z < ft_top, z > ft_base) #lower_free tropospheric values

    if z_bot - min(z) < 300 or np.sum(l_ft) == 0:
        return buncha_nans # can't calculate decouplng values if there is not enough MBL vertical or free-tropospheric 
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        theta_e_sml = np.nanmean(theta_e[l_i])
        theta_e_bzi = np.nanmean(theta_e[u_i])
        theta_e_uzi = np.nanmean(theta_e[l_ft])

        theta_l_sml = np.nanmean(theta_l[l_i])
        theta_l_bzi = np.nanmean(theta_l[u_i])
        theta_l_uzi = np.nanmean(theta_l[l_ft])
    
        qt_sml = np.nanmean(qt[l_i])
        qt_bzi = np.nanmean(qt[u_i])
        qt_uzi = np.nanmean(qt[l_ft])

    d_theta_e = theta_e_bzi - theta_e_sml
    d_theta_l = theta_l_bzi - theta_l_sml
    d_qt = qt_bzi - qt_sml
    buncha_nans['d_qt'] = d_qt
    buncha_nans['d_theta_l'] = d_theta_l
    buncha_nans['d_theta_e'] = d_theta_e
    
    alpha_thetal = (theta_l_bzi - theta_l_sml)/(theta_l_uzi - theta_l_sml)
    alpha_qt = (qt_bzi - qt_sml)/(qt_uzi - qt_sml)
    alpha_thetae = (theta_e_bzi - theta_e_sml)/(theta_e_uzi - theta_e_sml)
    buncha_nans['alpha_thetal'] = alpha_thetal
    buncha_nans['alpha_qt'] = alpha_qt
    buncha_nans['alpha_thetae'] = alpha_thetae
    
        
        
    return buncha_nans # hopefully no longer nans
    
    
    
def calc_decoupling_and_zi_from_flight_data(flight_data, usetheta=False):
    
    var_list = ['GGLAT', 'GGLON', 'GGALT', 'RHUM', 'ATX', 'MR', 'THETAE', 'THETA', 'PSX', 'DPXC', 'PLWCC']    
    
    
    sounding_dict = {}
    sounding_dict['TIME'] = flight_data.time.values
    for i in var_list:
        sounding_dict[i] = flight_data[i].values
    if 'ATX' in var_list:
        sounding_dict['ATX'] = sounding_dict['ATX'] + 273.15
    sounding_dict['DENS'] = density_from_p_Tv(flight_data['PSX'].values*100, flight_data['TVIR'].values+273.15)  
    sounding_dict['QL'] = flight_data['PLWCC'].values/sounding_dict['DENS']
    sounding_dict['THETAL'] = get_liquid_water_theta(
        sounding_dict['ATX'], sounding_dict['THETA'], sounding_dict['QL'])
    sounding_dict['QV'] = flight_data['MR'].values/(1+flight_data['MR'].values/1000)
    
    decoupling_dict = calc_decoupling_and_inversion_from_sounding(sounding_dict, usetheta=usetheta)
#     zi_dict = calc_zi_from_sounding(sounding_dict)
    return {**decoupling_dict}

def calculate_LTS(t_700, t_1000):
    """calculate lower tropospheric stability
    t_700: 700 hPa temperature in Kelvin
    t_1000: 1000 hPa temperature in Kelvin
    
    returns: lower tropospheric stability in Kelvin
    """
    theta_700 = theta_from_p_t(p=700.0, t=t_700)
    lts = theta_700-t_1000
    return lts
    
    
def calculate_moist_adiabatic_lapse_rate(t, p): 
    """calculate moist adiabatic lapse rate from pressure, temperature
    p: pressure in hPa
    t: temperature in Kelvin
    
    returns: moist adiabatic lapse rate in Kelvin/m
    """
    es = 611.2*np.exp(17.67*(t-273.15)/(t-29.65)) # Bolton formula, es in Pa
    qs = 0.622*es/(p*100-0.378*es)
    num = 1 + lv*qs/(Rdry*t)
    denom = 1 + lv**2*qs/(cp*Rvap*t**2)
    gamma = g/cp*(1-num/denom)
    return gamma
    
def theta_from_p_t(p, t, p0=1000.0):
    """calculate potential temperature from pressure, temperature
    p: pressure in hPa
    t: temperature in Kelvin
    
    returns: potential temperature in Kelvin
    """
    theta = t * (p0/p)**(Rdry/cp)
    return theta


def calculate_LCL(t, t_dew, z=0.0):
    """calculate lifting condensation level from temperature, dew point, and altitude
    t: temperature in Kelvin
    t_dew: dew point temperature in Kelvin
    z: geopotential height in meters. defaults to 0
    
    returns: lifting condensation level in meters
    
    raises: ValueError if any dew points are above temperatures (supersaturation)
    """
    if np.any(t_dew > t):
        t_dew = np.minimum(t, t_dew)
#         raise ValueError('dew point temp above temp, that\'s bananas')
    return z + 125*(t - t_dew)

def calculate_EIS(t_1000, t_850, t_700, z_1000, z_700, r_1000):
    """calculate estimated inversion strength from temperatures, heights, relative humidities
    t_1000, t_850, t_700: temperature in Kelvin at 1000, 850, and 700 hPa
    z_1000, z_700: geopotential height in meters at 1000 and 700 hPa
    r_1000: relative humidity in % at 1000 hPa
    
    returns: estimated inversion strength (EIS) in Kelvin
    """
    if hasattr(r_1000, '__iter__'):
        r_1000[r_1000>100] = 100  # ignoring supersaturation for lcl calculation
    t_dew = t_1000-(100-r_1000)/5.0
    lcl = calculate_LCL(t=t_1000, t_dew=t_dew, z=z_1000)
    lts = calculate_LTS(t_700=t_700, t_1000=t_1000)
    gamma_850 = calculate_moist_adiabatic_lapse_rate(t=t_850, p=850)
    eis = lts - gamma_850*(z_700-lcl)
    return eis





# def DEC_inv_layer_from_sounding(sounding):
#     rh = sounding['RHUM']
#     z = sounding['GGALT']
    
    
#     i_above_inv = np.where(rh<60)[0]
#     z_above_inv = z[i_above_inv]
#     if np.any(i_above_inv):
#         z_mid = np.min(z_above_inv)
#     else:
#         z_mid = np.nan
#     return {'z_mid': z_mid}
    
    
    
    
# def DEC_calc_decoupling_from_sounding(sounding_dict, usetheta=False, get_jumps=True, smooth_t=True):
#     z = sounding_dict['GGALT']
#     theta = sounding_dict['THETA']
#     theta_e = sounding_dict['THETAE']
#     qv = sounding_dict['QV']
#     t = sounding_dict['ATX']
#     if 'PSXC' in sounding_dict.keys():
#         p = sounding_dict['PSXC']
#     else:
#         p = sounding_dict['PSX']

#     if not usetheta:
#         theta_l = sounding_dict['THETAL']
#         ql = sounding_dict['QL']
#         if np.all(np.isnan(ql)):
#             qt = qv
#         else:
#             qt = qv + ql
#     else:
#         theta_l = sounding_dict['THETA']
#         qt = qv

        
        
#     zi = heffter_pblht_1D(z, theta)

    
#     upper_25 = zi['z_bot'] - (zi['z_bot'] - min(z))/4.
#     u_i = np.logical_and(z > upper_25, z < zi['z_bot'])
#     lower_25 = min(z) + (zi['z_bot'] - min(z))/4.
#     l_i = np.logical_and(z < lower_25, z > min(z))
    
#     ft_base = zi['z_bot']+500
#     ft_top = ft_base + 500
#     l_ft = np.logical_and(z < ft_top, z > ft_base)
    
#     buncha_nans = {"d_qt": np.nan, "d_theta_e": np.nan, "d_theta_l": np.nan,
#                 "alpha_thetal": np.nan, "alpha_qt":np.nan, "alpha_thetae": np.nan,
#                 "d_q_inv": np.nan, "d_t_inv": np.nan,
#                 "t_below_inv": np.nan, "t_above_inv": np.nan, "q_below_inv": np.nan, "q_above_inv": np.nan,
#                 "z_top": np.nan, "z_mid": np.nan, "z_bot": np.nan, "i_top": np.nan, "i_mid": np.nan, "i_bot": np.nan}
#     if zi['z_bot'] - min(z) < 300 or np.sum(l_ft) == 0:
#         return buncha_nans
#     with warnings.catch_warnings():
#         warnings.simplefilter("ignore", category=RuntimeWarning)
#         theta_e_sml = np.nanmean(theta_e[l_i])
#         theta_e_bzi = np.nanmean(theta_e[u_i])
#         theta_e_uzi = np.nanmean(theta_e[l_ft])

#         theta_l_sml = np.nanmean(theta_l[l_i])
#         theta_l_bzi = np.nanmean(theta_l[u_i])
#         theta_l_uzi = np.nanmean(theta_l[l_ft])
    
#         qt_sml = np.nanmean(qt[l_i])
#         qt_bzi = np.nanmean(qt[u_i])
#         qt_uzi = np.nanmean(qt[l_ft])

#     d_theta_e = theta_e_bzi - theta_e_sml
#     d_theta_l = theta_l_bzi - theta_l_sml
#     d_qt = qt_bzi - qt_sml
    
#     alpha_thetal = (theta_l_bzi - theta_l_sml)/(theta_l_uzi - theta_l_sml)
#     alpha_qt = (qt_bzi - qt_sml)/(qt_uzi - qt_sml)
#     alpha_thetae = (theta_e_bzi - theta_e_sml)/(theta_e_uzi - theta_e_sml)
    
#     # getting jumps across the inversion, old-fashioned way. bad for fuzzy inversions
# #     z_inv = inv_layer_from_sounding(sounding_dict)['z_mid']

# #     i_below_inv = np.logical_and(z > z_inv-200, z < z_inv)
# #     i_above_inv = np.logical_and(z > z_inv, z < z_inv+200)
# #     q_below_inv = np.nanmax(qt[i_below_inv])
# #     q_above_inv = np.nanmin(qt[i_above_inv])
# #     d_q_inv = q_above_inv - q_below_inv
    
# #     t_below_inv = np.nanmin(theta[i_below_inv])
# #     t_above_inv = np.nanmax(theta[i_above_inv])
# #     d_t_inv = t_above_inv - t_below_inv
    
    
#     if get_jumps:
#         ### moist adiabatic way
#         gamma_moist = get_moist_adiabatic_lapse_rate(T=t, p=p)*1000
#         if smooth_t:
#             gamma = -np.gradient(smooth(t, window_len=31), z)*1000
#         else: 
#             gamma = -np.gradient(t, z)*1000
#         gamma[np.gradient(z)>-1] = np.nan
#         gamma[z<330] = np.nan
#         gamma[z>3000] = np.nan
#         gamma[np.abs(gamma)>100] = np.nan
#         gamma_diff = (gamma-gamma_moist)/1000
# #         import matplotlib.pyplot as plt
# #         plt.plot(gamma, z)
# #         plt.ylim([0, 3000])
# #         raise ValueError('hahahah')
#         i_mid = np.nanargmin(gamma)
#         if np.isnan(i_mid):
#             print('no i_mid')
#             return buncha_nans
#         z_mid = z[i_mid]
#         max_gap = gamma[i_mid] - gamma_moist[i_mid]

#     #     z_bot = np.max(z[np.logical_and(z<z[i_mid], gamma>gamma_moist)])
#         try:
#             z_bot = np.max(z[np.logical_and(z<z[i_mid], gamma-gamma_moist>max_gap/4)])
#         except ValueError as v: # no crossing of the max_gap/4 line go for smallest gap below zmid
#             cands = z<z[i_mid]
#             if not np.any(cands):
#                 raise ValueError("no values below inversion middle!")            
#             z_bot = z[cands][np.argmin(np.abs(gamma[cands]-gamma_moist[cands]))]
#             i_bot = np.argwhere(z==z_bot)[0][0]
            

#         i_bot = np.argwhere(z==z_bot)[0][0]

#         top_candidates = np.logical_and(z>z[i_mid], gamma-gamma_moist>max_gap/4)
#         if np.any(top_candidates):
#             z_top = np.min(z[top_candidates])
#             i_top = np.argwhere(z==z_top)[0][0]
#         else: 
#             cands = z>z[i_mid]
#             if not np.any(cands):
#                 import matplotlib.pyplot as plt
#                 plt.plot(theta, z)
#                 plt.figure()
#                 plt.plot(gamma, z)
#                 raise ValueError("no values above inversion middle!")

#             z_top = z[cands][np.argmin(np.abs(gamma[cands]-gamma_moist[cands]))]
#             i_top = np.argwhere(z==z_top)[0][0]

#         i_upper = np.logical_and(z<=z_top, z>=z_mid)
#         if np.sum(i_upper) == 0:
#             print("error: no upper inv layer: z_top: {}  z_mid: {}".format(z_top, z_mid))
#             return buncha_nans
#         i_lower = np.logical_and(z>z_bot, z<z_mid)

#         q_above_inv = qt[i_top]
#         q_below_inv = qt[i_bot]
#         d_q_inv = q_above_inv - q_below_inv

#         t_below_inv = t[i_bot]
#         i_inv = np.logical_and(z>z_bot, z<z_top)

#         d_t_inv = integrate.trapz(gamma_diff[i_inv], z[i_inv])
#         t_above_inv = t_below_inv + d_t_inv
#     else:
#         i_bot, i_mid, i_top, z_bot, z_mid, z_top, q_above_inv, q_below_inv, t_above_inv, t_below_inv, d_q_inv, d_t_inv = np.nan, \
#             np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        
        
#     return{"d_qt": d_qt, "d_theta_e": d_theta_e, 'd_theta_l': d_theta_l,
#            "alpha_thetal": alpha_thetal, "alpha_qt":alpha_qt, "alpha_thetae": alpha_thetae,
#            "d_q_inv": d_q_inv, "d_t_inv": d_t_inv, 
#            "t_below_inv": t_below_inv, "t_above_inv": t_above_inv, "q_below_inv": q_below_inv, "q_above_inv": q_above_inv,
#            "z_top": z_top, "z_mid": z_mid, "z_bot": z_bot, "i_top": i_top, "i_mid": i_mid, "i_bot": i_bot }

# def DEC_calc_zi_from_sounding(sounding_dict):
#     z = sounding_dict['GGALT']
#     theta = sounding_dict['THETA']
#     RH = sounding_dict['RHUM']
#     T = sounding_dict['ATX']
#     zi_dict = {}
# #    zi_dict['Rich'] = mu.Ri_pbl_ht(u, v, q, T, z, smooth=True)
#     zi_dict['RH50'] = RH_50_pblht_1d(z, RH)
#     zi_dict['RHCB'] = RH_fancy_pblht_1d(z, RH)
#     zi_dict['Heff'] = heffter_pblht_1D(z, theta)
#     zi_dict['Heff']['T_bot'] = T[zi_dict['Heff']['i_bot']]
#     zi_dict['Heff']['T_top'] = T[zi_dict['Heff']['i_top']]
#     zi_dict['lat'] = np.nanmean(sounding_dict['GGLAT'])
#     zi_dict['lon'] = np.nanmean(sounding_dict['GGLON'])
#     zi_dict['time'] = sounding_dict['TIME'][0]
#     zi_dict['lon_p'] = -140 + 0.8*(zi_dict['lon']+140) + 0.4*(zi_dict['lat']-30)   
#     return zi_dict


    
def DEC_heffter_pblht_1D(z, theta, find_top=False):

    def moving_average(a, n=3):
        ret = np.cumsum(a, dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        return ret[n - 1:] / n
    """
    [courtesy of jmcgibbon]
    [made a little better by jkcm]
    Given height and theta returns
    the planetary boundary layer height from the Heffter criteria and the
    index of that height in the z-array. Assumes the data is 1-D with the
    axis being height. Assumes height in meters, and theta in K.
    """
    if z.shape != theta.shape:
        raise ValueError('z and theta must have the same shape')
    if len(z.shape) != 1:  # height axis
        raise ValueError('data has an invalid number of dimensions')
    if not (z < 4000).any():
        raise ValueError('must have data below 4000m')
        
    
#    z = moving_average(z, n=5)
#    theta = moving_average(theta, n=3)
    z = smooth(z, window_len=15)
    theta = smooth(theta, window_len=15)
    flipped = False
#    if np.all(z[100:-100] != sorted(z[100:-100])):  # not in ascending order
#        if np.all(z[100:-100:-1] == sorted(z[100:-100])):  # in descending order
    if z[0] > z[-1]: # starts off higher
        if True:
            z = z[::-1]
            theta = theta[::-1]
            flipped = True
        else:
            raise ValueError("data not in ascending or descending order")
    
    
    dtheta = np.diff(theta)
    dz = np.diff(z)
    dtheta_dz = np.zeros_like(dtheta)
    valid = dz != 0
    dtheta_dz[valid] = dtheta[valid]/dz[valid]
    del valid
    in_inversion = False
    found_inversion = False
    found_top = False
    
    theta_bot = np.nan
    z_bot = np.nan
    i_bot = np.nan
    theta_top = np.nan
    i_top = np.nan
    z_top = np.nan
    for i in range(z.shape[0]-1):  # exclude top where dtheta_dz isn't defined
        if z[i] > 4000.:
            # not allowed to have inversion height above 4km
            break
        if in_inversion:
            # check if we're at PBL top
            if theta[i] - theta_bot > 2:
                found_inversion = True
                theta_top = theta[i]
                i_top = i
                z_top = z[i]
                if not find_top:
                    break
                else:
                    break
                    #keep going up until we break the 
            # check if we're still in an inversion
#            layer_dtheta_dz = (theta[i] - theta_bot)/(z[i]-z_bot)
#            if layer_dtheta_dz > 0.005:
            if dtheta_dz[i] > 0.005:  # criterion for being in inversion
                pass  # still in inversion, keep going
            else:
                in_inversion = False
                theta_bot = np.nan
                z_bot = np.nan
                i_bot = np.nan
        else:
            if dtheta_dz[i] > 0.005:  # just entered inversion
                theta_bot = theta[i]
                i_bot = i
                z_bot = z[i]
                in_inversion = True
            else:
                # still not in inversion, keep going
                pass
    if found_inversion:
        if flipped:
            i_top = len(z)-i_top-1
            i_bot = len(z)-i_bot-1
        return {"z_top": z_top, "theta_top": theta_top, "i_top": i_top,
                "z_bot": z_bot, "theta_bot": theta_bot, "i_bot": i_bot,
                "inversion": True}
    else:
        # we didn't find a boundary layer height
        # return height of highest dtheta_dz below 4000m
        i_max = np.where(dtheta_dz == dtheta_dz[z[:-1] < 4000].max())[0][0]
        z_max = z[i_max]
        theta_max = theta[i_max]
        if flipped:
            i_max = len(z)-i_max-1
        return {"z_top": z_max, "theta_top": theta_max, "i_top": i_max,
                "z_bot": z_max, "theta_bot": theta_max, "i_bot": i_max,
                "inversion": False}
def DEC_heffter_pblht_2d(z, theta, axis=0, handle_nans=False):
    dummy = heffter_pblht_1D(np.arange(100), np.arange(100))
    res_dict = {key: np.empty(z.shape[axis]) for key in dummy.keys()}

    result = np.empty(z.shape[axis])
    for i,(z_i,theta_i) in enumerate(zip(z, theta)):
        try:
            res = heffter_pblht_1D(z_i,theta_i)
        except ValueError as e:
            if handle_nans:
                res = {"z_top": float('nan'), "theta_top": float('nan'), "i_top": float('nan'),
                "z_bot": float('nan'), "theta_bot": float('nan'), "i_bot": float('nan'),
                "inversion": False}
            else:
                raise e
        for key, value in res.items():
            res_dict[key][i] = value
    return res_dict