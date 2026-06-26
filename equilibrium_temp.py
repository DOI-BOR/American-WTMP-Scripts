import datetime as dt
import math

import tz_offset
reload(tz_offset)

# TODO: hardcoded lat/lon sould be variable

zero_K = 273.15
Tk = zero_K
stefan_boltz_const = 5.67e-8
lw_emissivity_water = 0.97
reference_density = 1000.
wind_a = 2.8e-9
wind_b = 1.2e-9
Magnus_a = 6.1078
Magnus_b = 17.27
Magnus_c = 237.7
Cb = 0.61

def bisection_solve(fTe_min, a, b, H1, uw, Ta_in, td, sat_vp_ta, tol=1e-2, max_iter=1000):
    
    # Evaluate the function at the lower and upper interval bounds
    fa = fTe_min(a, H1, uw, Ta_in, td, sat_vp_ta)
    fb = fTe_min(b, H1, uw, Ta_in, td, sat_vp_ta)

    # Initialize iteration counter for convergence protection
    iter_count = 0
    
    while (b - a) / 2.0 > tol:
        
        # Compute midpoint of the current search interval
        c = (a + b) / 2.0
        
        # Protect against divide-by-zero during iteration
        if (c-td) == 0.0:
            
            c += tol
            
        # Evaluate function at midpoint    
        fc = fTe_min(c, H1, uw, Ta_in, td, sat_vp_ta)

        # Stop if exact root found or interval is within tolerance
        if fc == 0 or (b - a) / 2.0 < tol:
            return c  # Root found
            
        # Track iterations and prevent infinite loops
        iter_count += 1
        
        if iter_count > max_iter:
            
            raise ValueError("bisection_solve: Maximum iterations exceeded")

        # Select subinterval that contains the root
        if fa * fc < 0.0:
            
            b = c
            fb = fc
            
        else:
            
            a = c
            fa = fc

    # Return midpoint of final interval as best estimate
    return (a + b) / 2.0  


def newton_raphson_solve(f, df, x0, H1, uw, Ta_in, Td, sat_vp_ta, tol=1e-2, max_iter=1000):
    # Initialize iteration counter and starting estimate
    iter_count = 0
    x_old = x0
    
    # Iterate until convergence or maximum iterations reached
    while iter_count < max_iter:
        
        # Evaluate function and derivative at current estimate
        fx_old = f(x_old, H1, uw, Ta_in, Td, sat_vp_ta)
        dfx_old = df(x_old, H1, uw, Ta_in, Td, sat_vp_ta)

        # Newton-Raphson cannot proceed with near-zero derivative
        if abs(dfx_old) < tol:
            raise ValueError("Derivative near zero")

        # Compute next approximation
        x_new = x_old - fx_old / dfx_old
        
        # Convergence achieved when successive estimates are close
        if abs(x_new - x_old) < tol:
            return x_new  # Root found

        # Update state for next iteration
        x_old = x_new
        iter_count += 1

    # Raise error if convergence was not achieved
    raise ValueError("newton_raphson_solve: Maximum iterations exceeded")


def get_decimal_day_of_year(tin):
    
    # Define the first moment of the current year
    foy = dt.datetime(tin.year, 1, 1, 0, 0, 0)
    
    # Convert elapsed time since Jan 1 into fractional day-of-year
    ddoy = (tin - foy).total_seconds() / 86400. + 1. + tz_offset.days
    
    # Return decimal day-of-year value
    return ddoy


def solar_alt_angle(tin):
    
    # Site-specific geographic coordinates and time-zone reference longitude
    latitude = 38.1
    longitude = 121.8
    time_zone_longitude = 120.

    # Calculate longitude correction from standard time meridian
    delta = (longitude - time_zone_longitude) / 15.
    doy = get_decimal_day_of_year(tin)

    # Estimate solar declination angle from day of year
    solar_decl = 0.4092 * math.cos(2. * math.pi / 365. * (172. - doy))
    
    # Components used in solar altitude angle calculation
    t1 = math.sin(solar_decl) * math.sin(latitude * math.pi / 180.)
    t2 = math.cos(solar_decl) * math.cos(latitude * math.pi / 180.)
    
    # Determine decimal hour and corresponding solar hour angle
    hr = (doy - math.floor(doy)) * 24.
    hr_angle = math.pi * (hr - delta - 12.) / 12.
    
    # Compute sine of solar altitude angle
    sin_solar_alt = t1 + t2 * math.cos(hr_angle)
    
    # Convert result from radians to degrees
    solar_alt_angle = math.asin(sin_solar_alt) * 180. / math.pi
    
    # Return solar altitude above the horizon
    return solar_alt_angle


def sw_water_reflectance(tin, cloud_in):
    
    # Determine solar altitude angle for the specified time
    alt_angle = solar_alt_angle(tin)
    
    # At night, assume all incoming shortwave radiation is reflected
    if alt_angle <= 0.:
        
        return 1.
        
    else:
        
        # Select empirical reflectance coefficients based on cloud cover
        if cloud_in > 0.95:
            
            acoef = 0.33
            bcoef = -0.45
            
        elif cloud_in > 0.55:
            
            acoef = 0.95
            bcoef = -0.75
            
        elif cloud_in > 0.05:
            
            acoef = 2.20
            bcoef = -0.97
            
        else:
            
            acoef = 1.18
            bcoef = -0.77
            
        # Compute water surface reflectance from solar altitude
        water_sfc_reflectance = acoef * alt_angle ** bcoef
        
        # Limit reflectance to a maximum physical value of 1.0
        water_sfc_reflectance = min(water_sfc_reflectance, 1.)
        
        # Return estimated shortwave reflectance fraction
        return water_sfc_reflectance


def heat_flux_surface_longwave_down(air_temp_in, cloud_in):
    
    # Constants controlling surface reflection and cloud effects
    sfc_reflect_fract = 0.03
    cloud_cover_fract_coef = 0.17
    
    # Empirical coefficient used for atmospheric emissivity
    con1 = 0.937e-5
    Tkel = air_temp_in + zero_K
    
    # Estimate effective atmospheric emissivity
    emissivity_air = con1 * (1. + cloud_cover_fract_coef * cloud_in ** 2) * Tkel ** 2
    
    # Calculate downward longwave radiation flux at the surface
    hlwd = emissivity_air * stefan_boltz_const * (1. - sfc_reflect_fract) * Tkel ** 4
    
    # Return downward longwave heat flux [W/m2]
    return hlwd


def latent_heat_vaporization(temp_in):
    
    # Temperature-dependent latent heat of vaporization [J/kg]
    return 1000. * (2499. - 2.36 * temp_in)


def sat_water_vapor_pres(temp_in):
    
    # Magnus-Tetens empirical formula for saturation vapor pressure
    return Magnus_a * math.exp(Magnus_b * temp_in / (temp_in + Magnus_c))


def fTe(Te, H1, uw, Ta_in, Td, sat_vp_ta):
    '''
    Calculation of fTe with only Te-dependent terms
    '''
    
    # Compute latent heat of vaporization at the equilibrium temperature
    Lw = latent_heat_vaporization(Te)
    
    # Compute saturation vapor pressure at the equilibrium temperature
    sat_vp_te = sat_water_vapor_pres(Te)
    
    # Estimate vapor pressure gradient between water surface and dewpoint
    beta = (sat_vp_te - sat_vp_ta) / (Te - Td)
    
    # Combine radiative and evaporative heat transfer terms
    return H1 - lw_emissivity_water * stefan_boltz_const * (Tk ** 4 + 4. * Tk ** 3 * Te + 6. * Tk ** 2 * Te ** 2) + \
        reference_density * Lw * (wind_a + wind_b * uw) * ((Cb * Ta_in + beta * Td) - (Cb + beta) * Te)


def dfTe_dTe(Te, H1, uw, Ta_in, Td, sat_vp_ta):
    '''
    Calculation of fTe with only Te-dependent terms
    '''

    # Recompute latent heat and saturation vapor pressure at Te
    Lw = 1000. * (2499. - 2.36 * Te)
    sat_vp_te = Magnus_a * math.exp(Magnus_b * Te / (Te + Magnus_c))

    # Derivative of latent heat with respect to temperature
    dLw_dTe = -1000. * 2.36
    
    # Derivative of saturation vapor pressure from Magnus-Tetens equation
    dsat_vp_te_dTe = (Magnus_a * Magnus_b * math.exp(Magnus_b * Te / (Te + Magnus_c))) / (Te + Magnus_c) - \
                     (Magnus_a * Magnus_b * Te * math.exp(Magnus_b * Te / (Te + Magnus_c))) / (Te + Magnus_c) ** 2
    
    # Derivative of beta using the quotient rule
    dbeta_dTe = (dsat_vp_te_dTe * (Te - Td) - (sat_vp_te - sat_vp_ta)) / (Te - Td) ** 2

    # Derivative of the longwave radiation component
    dTe_term = - 4.0 * lw_emissivity_water * stefan_boltz_const * Tk ** 3 - 12.0 * lw_emissivity_water * stefan_boltz_const * Tk ** 2 * Te
    
    # Contribution from temperature dependence of latent heat transfer
    Lw_term = reference_density * (wind_a + wind_b * uw) * ((Cb * Ta_in + (sat_vp_te - sat_vp_ta) / (Te - Td) * Td) - (
            Cb + (sat_vp_te - sat_vp_ta) / (Te - Td)) * Te)
    
    # Contribution from the explicit beta term
    beta_term = reference_density * Lw * (wind_a + wind_b * uw) * (-(Cb + (sat_vp_te - sat_vp_ta) / (Te - Td)))

    # Return the complete analytical derivative of fTe
    return dTe_term + Lw_term + beta_term + reference_density * dLw_dTe * (wind_a + wind_b * uw) * (
            (Cb * Ta_in + (sat_vp_te - sat_vp_ta) / (Te - Td) * Td) - (
            Cb + (sat_vp_te - sat_vp_ta) / (Te - Td)) * Te) + reference_density * Lw * (
            wind_a + wind_b * uw) * dbeta_dTe * Td


def fTe_min(Te, H1, uw, Ta_in, Td, sat_vp_ta):
    
    # Reformulate root-finding problem for bisection solver
    return Te - fTe(Te, H1, uw, Ta_in, Td, sat_vp_ta)

def fTe_abs(Te, H1, uw, Ta_in, Td, sat_vp_ta):
    
    # Absolute objective function used by scipy minimization
    return abs(fTe(Te, H1, uw, Ta_in, Td, sat_vp_ta))


def equilibrium_temp(dtt1, at, cl, sr, ws, td, te_guess, type='nr'):
    
    # Calculate fraction of shortwave radiation reflected by the water surface
    reflectance = sw_water_reflectance(dtt1, cl)
    
    # Net absorbed shortwave radiation
    Hsw = sr * (1. - reflectance)
    
    # Downward atmospheric longwave radiation
    HH = heat_flux_surface_longwave_down(at, cl)
    
    # Combined radiative heat input term
    H1 = Hsw + HH
    
    # Saturation vapor pressure at air temperature
    sat_vp_ta = sat_water_vapor_pres(at)
    
    # Solve using Newton-Raphson method (fastest when convergent)
    if type == 'nr':
        return newton_raphson_solve(fTe, dfTe_dTe, te_guess, H1, ws, at, td, sat_vp_ta, tol=1.0e-2, max_iter=1000)
    
    # Solve using bisection method (more robust but slower)
    elif type == 'bs':
        return bisection_solve(fTe_min, -30.0, 55.0, H1, ws, at, td, sat_vp_ta, tol=1.0e-2, max_iter=1000)
    
    # Legacy scipy implementation retained for reference
    elif type == 'scipy':
       
       # Original equilibrium temp calc using scipy.optimize/minimize.  
       # Can't use scipy in jython
       args = (H1, ws, at, td, sat_vp_ta)
       
       # Minimize absolute residual to estimate equilibrium temperature
       res = minimize(fTe_abs, x0=te_guess, args=args)
       return res.x
      
    # Reject unsupported solution methods  
    else:
        raise ValueError('equilibrium_temp() - calculation type not known:', type)


def calc_equilibrium_temp(dtt, at, cl, sr, td, ws):
    '''
    adapted from numpy/scipy routine by Steve, Ben, Scott, which was generated from:
    "Stratification and heat transfer in lakes and reservoirs" chapter in "Hydrodynamics and Transport
    for Water Quality Modeling" by Martin and McCutcheon (pg. 373)

    dtt = list of datetime objects
    at = list of airtemps [C]
    cl = list of cloud cover fractions [0-1]
    sr = shortwave irradiance [W/m2]
    dp = dewpoint [C]
    ws = wind speed [m/s]
    '''
    
    # Determine the number of timesteps and initialize output list
    nt = len(dtt)
    te = []

    # Process each timestep independently
    for j in range(nt):

        # Use air temperature as the initial guess for the first timestep
        if j == 0:
            x0 = at[j]
        
        else:
            
            # Use the previous equilibrium temperature as the next initial guess
            x0 = te[j - 1]

        # The Newton-Raphson solver sometimes fails (when the derivitive goes big/non-linear) or produces a really bad
        # value.  We calulcate a secondary solution (which is sometimes inaccurate be over a degree, but usually is
        # within 0.2 degrees of the solution using scipy.optimize.minimize) and use it in those cases.
        
        # Initialize output slot with a placeholder value
        te.append(-999)
        
        
        # Always compute a bisection-based solution as a fallback
        te_bs = equilibrium_temp(dtt[j], at[j], cl[j], sr[j], ws[j], td[j], x0, type='bs')
        
        try:
            # Attempt the faster Newton-Raphson solution
            te[j] = equilibrium_temp(dtt[j], at[j], cl[j], sr[j], ws[j], td[j], x0, type='nr')
        
        except:
            
            # If Newton-Raphson fails to converge, use the bisection result
            te[j] = te_bs
        
        # Detect unrealistic Newton-Raphson results by comparison to bisection
        if abs(te[j] - te_bs) > 1.0:
            
            # Replace suspect Newton-Raphson values with the more robust solution
            te[j] = te_bs

    # Return equilibrium temperature series for all timesteps
    return te


