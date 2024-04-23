"""
Created April 6 2021
author: Hans Mohrmann (jkcm@uw.edu)
author: Ehsan Erfani (most functions, where marked)
"""

#standard library
import os
import pandas as pd


#Specials
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scipy.optimize import fsolve, fmin


#constants
Cp = 1005.7 # specific heat of dry air at constant pressure: J/kg/K
Lv = 2.5e6  # latent heat of vaporization: J/kg
Rd_cp = 0.286
g = 9.81


###############################################################
### Calculate inversion height (min of dtheta/dz * dRH/dz)
def inv_h(THETA, RH, dzz, zz):
    #from Ehsan's notebook, minor adjustments
#    d_THETA = delta(THETA)
#    d_RH = delta(RH)
    d_THETA = THETA[:-1] - THETA[1:]
    d_RH = RH[:-1] - RH[1:]
    z  = zz[i,:]  if len(zz.shape) > 1 else zz
    dz = dzz[i,:] if len(zz.shape) > 1 else dzz
    dTHETA_dz = d_THETA / dz[1:]
    dRH_dz = d_RH / dz[1:]
    dTHETA_dz[dTHETA_dz < 0] = 0
    dRH_dz[dRH_dz > 0] = 0    
    func = dTHETA_dz * dRH_dz
    indx = np.min(np.where(func == np.nanmin(func))[0])                    
    # approximate func as a parabola around the minimum and find the height where that parabola is minimized. 
    # This will allow the inversion height to vary continuously as the input profiles change.
    # inversion_test is defined at midpoints of grid
    zavg = 0.5 * (z[:-1] + z[1:])
    rnge = range(indx-1, indx+2)
    # we define the parabola, converting from m to km.
    pp = np.polyfit(1e-3 * zavg[rnge], func[rnge], 2)
    # take the derivative of the parabola in coeffient space.
    pp_prime = np.array([2 * pp[0], pp[1]]) # this is its derivative
    # find the zero-crossing of the derivative. This is the inversion height in meters
    z_inv = -1e3 * pp_prime[1] / pp_prime[0]         
    return z_inv

# Function to calculate saturation mixing ratio from es and P, and es from T based on Bulk equation
# inputs: P in hPa, and Tk in K.
def sat_q(Tk, P):
    #Tc    = Tk - 273.15 
    Tc    = max(-85, min(Tk - 273.15, 70))
    #es    = 6.1121 * np.exp( (18.678 - Tc / 234.5) * (Tc / (257.14 + Tc)) ) # Bulk Eq.
    pcoef_esw = [-0.3704404e-15,   0.2564861e-13,   0.6936113e-10, \
              0.2031998e-7,    0.2995057e-5,    0.2641412e-3,  \
              0.1430341e-1,    0.4440316,       6.105851]
    es = np.polyval(pcoef_esw, Tc)
#    qsat = 621.97 * es / (P - es)
    qsat = 621.97 * es / max( es, (P - es) )

    return abs(qsat)

# Function to compute the temperature of moist air that is induced by the condensation/evaporation of liquid water
# (NO ICE PERMITTED), resulting in a mixture that is at equilibrium. 
# The equilibrium mixture will either have no liquid water or its water vapor mixing ratio will be equal to the saturation value.
def sat_adj(Z, TL, QT, P):
    QV    = QT.copy()
    QL    = QT.copy()
    T0    = QT.copy()
    T     = QT.copy()
    QV[:] = np.nan
    QL[:] = np.nan
    T0[:] = np.nan
    T[:]  = np.nan
    for i in range(QT.shape[0]):
        T0[i]  = TL[i] - g * Z[i] / Cp# + (Lv / Cp) * QT[i] / 1000
        qsat = sat_q(T0[i], P[i])
        if qsat >= QT[i]:
            QV[i] = QT[i]
            QL[i] = 0
            T[i]  = T0[i]
        #elif qsat < QT[i]:
        else:
            T0[i]  = TL[i] - g * Z[i] / Cp + (Lv / Cp) * QT[i] / 1000
            # Solve for T, based on the notion that moist static energy is unchanged by condensation/evaporation: 
            # Build a function that takes in a guess for T along with fixed values of TL, QT, and Z. 
            # Search for a value of T that returns f(T) as 0. 
            # This basically says that our guess for T should give the same moist static energy as that implied by 
            # the input profiles TL(z) and QT(z).
            f = lambda x: x + g * Z[i] / Cp + (Lv / Cp) * min( QT[i], sat_q(x, P[i]) ) / 1000 - TL[i] - (Lv / Cp) * QT[i] / 1000
            T[i] = fsolve(f, T0[i])
            QV[i] = sat_q(T[i], P[i])
            QL[i] = QT[i] - QV[i]

    return QL, T

## find delta after calculating mid-point values
def delta(THETA):
    THETA = THETA[::-1]
    nz = THETA.shape[0]
    THETAi = np.concatenate((THETA, np.array([0.])))
    THETAi[:] = np.nan    
    THETAi[0] = 0
    THETAi[1:nz] = 0.5 * (THETA[:nz-1] + THETA[1:nz])
    THETAi[nz] = 1.5 * THETA[nz-1] - 0.5 * THETA[nz-2]
    d_THETA = THETAi[1:] - THETAi[:-1]    
    return d_THETA[::-1]

# Function to calculate rho-weighted vertical integral
def z_integral(var, ERA_z, ERA_RHO, upbound):
#     idx = np.nanmin(np.where( abs(ERA_z - upbound) == np.nanmin(abs(ERA_z - upbound)) ))
    idx = np.nanargmin(np.abs(ERA_z - upbound))
    delta_ERA_z = delta(ERA_z)
    try:
        var_integ = np.nansum(var[idx:] * delta_ERA_z[idx:] * ERA_RHO[idx:])
    except IndexError as e:
        print(idx)
        print(var)
        print(delta_ERA_z)
        print(ERA_RHO)
    return var_integ






## Function to produce modified initial ERA5 profile
def produce_ERA5_profile(ERA_Z, ERA_TL, ERA_QT, CTH, TL_ML, QT_ML, L_FT, TL_diff, QT_diff):

    #Ehsan's code, some variables added by Hans
    
    f_FT    = 3   # factor to determine the top of the region where we're fitting the line (zi + f_FT * L_FT ). 

    
    Zi = CTH
    idx_zi    = np.nanmax( np.where( abs(ERA_Z - Zi) == np.nanmin(abs(ERA_Z - Zi)) ) )
    # Divide variables in two parts: FT and BL
    ERA_Z_BL  = ERA_Z[idx_zi:]
    ERA_Z_FT  = ERA_Z[:idx_zi]
    ERA_TL_BL = ERA_TL[idx_zi:]
    ERA_TL_FT = ERA_TL[:idx_zi]
    ERA_QT_BL = ERA_QT[idx_zi:]
    ERA_QT_FT = ERA_QT[:idx_zi]
    
    # Define outputs:
    TL_BL_o = ERA_TL_BL.copy()
    TL_FT_o = ERA_TL_FT.copy()
    QT_BL_o = ERA_QT_BL.copy()
    QT_FT_o = ERA_QT_FT.copy()

    TL_BL_o[:] = np.nan
    TL_FT_o[:] = np.nan
    QT_BL_o[:] = np.nan
    QT_FT_o[:] = np.nan
    
    ## FT:
    idx_L_FT    = np.nanmax( np.where( abs(ERA_Z_FT - (Zi + L_FT) ) == np.nanmin(abs(ERA_Z_FT - (Zi + L_FT) )) ) )
    TL_FT_o[:idx_L_FT+1] = ERA_TL_FT[:idx_L_FT+1]
    QT_FT_o[:idx_L_FT+1] = ERA_QT_FT[:idx_L_FT+1]
    
    # fit a line to the ERA TL & QT profiles away from the inversion, and extrapolate down to the inversion. 
    zind  = np.where( (ERA_Z_FT >= (Zi + L_FT) ) & (ERA_Z_FT <= (Zi + f_FT * L_FT)) )[0]
    zind2 = np.where( (ERA_Z_FT >=  Zi )         & (ERA_Z_FT <= (Zi + L_FT)) )[0]
    # Include a few points above the zi+L_FT to avoid jump in the FT profiles.
    zind2 = np.concatenate( (np.array([np.nanmin(zind2) - 4, np.nanmin(zind2) - 3, np.nanmin(zind2) - 2,\
                                   np.nanmin(zind2) - 1]), zind2) )
    TL_FT_o[zind2] = np.polyval(np.polyfit(ERA_Z_FT[zind], ERA_TL_FT[zind], 1), ERA_Z_FT[zind2])
    QT_FT_o[zind2] = np.polyval(np.polyfit(ERA_Z_FT[zind], ERA_QT_FT[zind], 1), ERA_Z_FT[zind2])
    
    # Make sure the profile has strong inversion and well-mixed a few hundreds meters below the inversion
    #TL_FT_o = np.array([max(TL_FT_o[i], TL_ML + TL_diff) for i in range(len(TL_FT_o))])   #commented out by Ehsan, 4/12
    QT_FT_o = np.array([min(QT_FT_o[i], QT_ML - abs(QT_diff)) for i in range(len(QT_FT_o))])

    # Make sure that dq/dz <= 0:   #Added by Ehsan, 4/12
    QT_FT_o = QT_FT_o[::-1]
    for i in range(len(QT_FT_o) -1):
        if QT_FT_o[i+1] > QT_FT_o[i]:
            QT_FT_o[i+1] = QT_FT_o[i]
    QT_FT_o = QT_FT_o[::-1]

    
    ## BL:
    TL_BL_o = np.array([min(ERA_TL_BL[i], TL_ML) for i in range(len(ERA_TL_BL))])
    QT_BL_o = np.array([max(ERA_QT_BL[i], QT_ML) for i in range(len(ERA_QT_BL))])
    
    ## Concatenate BL and FT
    TL_o = np.concatenate((TL_FT_o, TL_BL_o))
    QT_o = np.concatenate((QT_FT_o, QT_BL_o))
    
    return TL_o, QT_o

def optimize_ERA5(ERA_Z, ERA_TL, ERA_QT, ERA_p, ERA_t, ERA_q, ERA_RHO, ERA_lwc, LWP_target, CTH_target, CTH, TL_ML, QT_ML, ff, TL_diff, QT_diff, L_FT):
    # Ehsan's code, some variables added by Hans
    
    upbound = 4000 # upper bound (in meter) for the vertical integral
    
    TL_o, QT_o = produce_ERA5_profile(ERA_Z, ERA_TL, ERA_QT, CTH, TL_ML, QT_ML,L_FT, TL_diff, QT_diff)    
    
    # compute LWC and Tadj by using saturation adjustment at each height
    QL_o_adj, T_o_adj = sat_adj(ERA_Z, TL_o, QT_o, ERA_p)

    # Calculate Trho, which is the density temperature. It's the same as the virtual temperature except for the inclusion of 
    # the density change due to liquid water loading. This will keep the vertically-integrated density from drifting too far 
    # from that of ERA5, which is a more physical constraint than the one on Tl.
    #T_o_adj  = TL_o - g * ERA_Z / Cp + (Lv / Cp) * QL_o_adj / 1000
    T_o_rho    = T_o_adj * (1 + 0.61 * (QT_o - QL_o_adj) / 1000 - QL_o_adj / 1000)
    ERA_T_rho  = ERA_t   * (1 + 0.61 * ERA_q / 1000 - ERA_lwc / 1000)

    # Compute the LWP of the modified profiles and other vertical integrals
    LWP_adj       = z_integral(QL_o_adj, ERA_Z, ERA_RHO, upbound)
    TWP_o         = z_integral(QT_o,     ERA_Z, ERA_RHO, upbound)
    TWP_ERA       = z_integral(ERA_QT,   ERA_Z, ERA_RHO, upbound)
    TL_o_zint     = z_integral(TL_o,     ERA_Z, ERA_RHO, upbound)
    TL_ERA_zint   = z_integral(ERA_TL,   ERA_Z, ERA_RHO, upbound)
    To_rho_zint   = z_integral(T_o_rho,  ERA_Z, ERA_RHO, upbound)
    T_rho_ERA_zint= z_integral(ERA_T_rho,  ERA_Z, ERA_RHO, upbound)    
    ones          = ERA_Z.copy()
    ones[:]       = 1
    RHO_zint      = z_integral(ones,   ERA_Z, ERA_RHO, upbound)
            
    # Output a (positive) number that tells how well the resulting profile matches LWP_target 
    # we have a few targets that we want to match: LWP_adj, vertically-integrated TL_o(z), vertically-integrated QT_o(z) & CTH 
    
    # First, factors; which have a general formula of [1 / (variable uncertainty) ] ^ 2
    # uncertainties are small values for the deviation of the integrated variables 
    # (with each scaled by int_z rho*dz and then squared) 
    # This would essentially force the vertical integrals to be conserved in the new profiles.
    
    # For fac2, using the factor of Lv/Cp forces some equivalence between moisture and temperature errors.  
    # The default ff value is equal to 0.3 kg/m2 because it's probably a couple percent of the TWP.
    # Modifying it has minimal effect on the solution.  
    fac1 = 1 / ( 10.  ** 2 )
    fac2 = (1 / ( ff ** 2 ) ) * ( Cp / Lv ) ** 2
    fac3 = 1 / ( (ff * 1000) ** 2 )
    fac4 = 1 / ( 20.  ** 2 )

    output = fac1 * (LWP_target - LWP_adj) ** 2 + fac2 * (To_rho_zint - T_rho_ERA_zint) ** 2 + \
             fac3 * (TWP_o      - TWP_ERA) ** 2 + fac4 * (CTH_target - CTH) ** 2
    #+ fac2 * (TL_o_zint - TL_ERA_zint) ** 2 + \
    
    return output

def get_adjusted_ERA_profile(profile, zi_init, target_lwp, plot=False):
    #Hans' adaptation of Ehsan's truncated code
    
    ERA_z = profile['z']
    ERA_Tl = profile['theta_l']
    ERA_qt = profile['qt']
    ERA_p = profile['p']
    ERA_lwc = profile['lwc']
    ERA_RHO = profile['rho']
    ERA_t = profile['t']
    ERA_q = profile['qv']
    CTH_target = zi_init
    ff0    = 0.3                          #WHAT IS THIS
    TL_diff= 10                           #WHAT IS THIS
    QT_diff= 2.53                         #WHAT IS THIS
    CTH0 = zi_init                        #INITIAL guess for inversion top
    L_FT    = 500 # height above the inversion where ERA profiles don't "feel" the BL anymore
    
    idx_zi = np.nanargmin(np.abs(ERA_z-zi_init))
    TL_ML0 = ERA_Tl[idx_zi]               #initial guess for theta_l of 
    QT_ML0 = ERA_qt[idx_zi]               #initial guess for 
    
    #HANS: added ERA_p, ERA_t, ERA_q, ERA_RHO, ERA_lwc and L_FT to optimize_ERA5 args
    #HANS: changed "SSMI_lwp" to "target_lwp", which is set to ERA lwp if no SSMI in args to ERA5_args
    #HANS: cosmetic: changed fmin output to include some flags (not necessary change)
    #HANS: moved upbound into optimize_ERA5 (could also declare global)
    #HANS: moved f_FT into produce_ERA_profile (could also declare global)
    
    func_ERA = lambda x: optimize_ERA5(ERA_z, ERA_Tl, ERA_qt, ERA_p, ERA_t, ERA_q, ERA_RHO, ERA_lwc,
                                       target_lwp, CTH_target, x[0], x[1], x[2], x[3], TL_diff, x[4], L_FT)    
    opt_vals, fopt, iters, _, flag = fmin(func=func_ERA, x0=[CTH0, TL_ML0, QT_ML0, ff0, QT_diff], disp=False, full_output=True)
    flag_str = '' if flag==0 else f'                         FAILED ({flag})'
    
    print('Initial values: ', CTH0, TL_ML0, QT_ML0)
    print('Final values: ', opt_vals, flag_str)
    
    TL_o_f2, QT_o_f2 = produce_ERA5_profile(ERA_z, ERA_Tl, ERA_qt, opt_vals[0], opt_vals[1], opt_vals[2],L_FT, TL_diff, QT_diff)

    QL_o_f, T_o_f = sat_adj(ERA_z, TL_o_f2, QT_o_f2, ERA_p)
    LWP_0  = z_integral(ERA_lwc, ERA_z, ERA_RHO, 10000)
    LWP_f  = z_integral(QL_o_f, ERA_z, ERA_RHO, 10000)
    
    if plot:
        QL_o_f, T_o_f = sat_adj(ERA_z, TL_o_f2, QT_o_f2, ERA_p)
        LWP_0  = z_integral(ERA_lwc, ERA_z, ERA_RHO, 10000)
        LWP_f  = z_integral(QL_o_f, ERA_z, ERA_RHO, 10000)

        str_CTH0  = '{:.0f}'.format(CTH0)
        str_TL_ML0 = '{:.1f}'.format(TL_ML0)
        str_QT_ML0 = '{:.2f}'.format(QT_ML0)
        str_LWP_0  = '{:.1f}'.format(LWP_0)

        str_CTHf  = '{:.0f}'.format(opt_vals[0])
        str_TL_MLf = '{:.1f}'.format(opt_vals[1])
        str_QT_MLf = '{:.2f}'.format(opt_vals[2])
        str_LWP_f  = '{:.1f}'.format(LWP_f)

        ### range of the x-axis:
        xx1 = 265 
        xx2 = 315
        ###

        fig  = plt.figure(figsize=(11,8))
        axis = fig.add_subplot(111)
        cc5  = axis.plot(ERA_Tl, ERA_z, c='k', label='ERA5 ctrl $T_l$', linewidth= 2)
        #cc50 = axis.plot(TL_o_f, Z_o_f, c='tab:red', label='ERA5 adj1 $T_l$', linewidth= 1)
        cc500= axis.plot(TL_o_f2, ERA_z, '-.', c='tab:blue', label='ERA5 adj $T_l$', linewidth= 2)
        #cc60 = axis.hlines(CTH_f, xx1, xx2, colors='tab:red', linestyles='dashed', linewidth=1, label='ERA5 adj1 $Z_{target}$')
        cc600= axis.hlines(opt_vals[0], xx1, xx2, colors='tab:cyan', linestyles='solid', linewidth=1, label='ERA5 adj $Z_{f}$')#, *, data=None, **kwargs)[source]
        axis.set_ylabel('Height (m)', fontsize = '18')
        axis.set_xlabel('$T_l$ ($K$)', fontsize = '18', horizontalalignment='right', x=.85)
        axis.set_title('$Z_0=$'+str_CTH0+', $T_{l_0}$='+str_TL_ML0+', $q_{t_0}$='+str_QT_ML0+', $LWP_0$='+str_LWP_0+'\n'+\
                       '$Z_f=$'+str_CTHf+', $T_{l_f}$='+str_TL_MLf+', $q_{t_f}$='+str_QT_MLf+', $LWP_f$='+str_LWP_f,\
                       fontsize = '18', horizontalalignment='right', x=.85)
        axis.set_ylim([0, 4000])
        axis.set_xlim([xx1, xx2])
        axis.set_xticks(range(290,316,5))
        axis.legend(fontsize = '16', loc='lower right')
        axis.grid(linestyle=':', axis='y', linewidth=1)
        axis.grid(linestyle=':', axis='x', linewidth=1)
        axis.tick_params(axis='both', which='major', labelsize=14)

        ax = axis.twiny()
        cc7  = ax.plot(ERA_qt, ERA_z, c='k', label='ERA5 ctrl $q_t$', linewidth= 2)
        #cc70 = ax.plot(QT_o_f, Z_o_f, c='tab:red', label='ERA5 adj1 $q_t$', linewidth= 1)
        cc700= ax.plot(QT_o_f2, ERA_z, '-.', c='tab:blue', label='ERA5 adj $q_t$', linewidth= 2)
        ax.set_xlabel('$q_t$ ($g$ $kg^{-1}$)', fontsize = '18', horizontalalignment='right', x=.3)
        ax.set_xlim([0, 22])
        ax.set_xticks(range(4,14,2))
        ax.legend(fontsize = '16', loc='lower left')
        ax.grid(linestyle=':', axis='x', linewidth=1)
        ax.tick_params(axis='both', which='major', labelsize=14)


    
    else:
        fig, ax, axis = None, None, None
    

    return {'zi': opt_vals[0], 'theta_l': TL_o_f2, 'q_t': QT_o_f2, 'z': ERA_z, 
            'opt_flag': flag, 'fig': fig, 'axl': ax, 'axr': axis}

#     return TL_o_f2, QT_o_f2, QL_o_f, T_o_f