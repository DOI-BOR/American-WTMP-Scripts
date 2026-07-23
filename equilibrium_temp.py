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
    """
    Finds the root of fTe_min within a bracketed interval [a, b] using the
    bisection method.

    The bisection method repeatedly halves the search interval, selecting the
    subinterval that contains a sign change (and therefore a root). It is more
    robust than Newton-Raphson but converges more slowly. A small perturbation
    is applied when the midpoint coincides with the dewpoint temperature to avoid
    a divide-by-zero in the residual function.

    Inputs:
      fTe_min   -- callable; the residual function to solve, with signature
                   fTe_min(Te, H1, uw, Ta_in, td, sat_vp_ta)
      a         -- float; lower bound of the search interval (deg C)
      b         -- float; upper bound of the search interval (deg C)
      H1        -- float; combined radiative heat input term (W/m2)
      uw        -- float; wind speed (m/s)
      Ta_in     -- float; air temperature (deg C)
      td        -- float; dewpoint temperature (deg C)
      sat_vp_ta -- float; saturation vapor pressure at air temperature (mbar)
      tol       -- float; convergence tolerance; iteration stops when the interval
                   width is less than tol (default 1e-2)
      max_iter  -- int; maximum number of bisection iterations before raising an
                   error (default 1000)

    Output:
      Returns the estimated root (equilibrium temperature, deg C) as a float.
      Raises ValueError if the maximum number of iterations is exceeded.
    """
    
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
    """
    Finds the root of function f using the Newton-Raphson iterative method,
    starting from an initial estimate x0.

    Each iteration refines the estimate using the ratio of the function value
    to its derivative. Converges quadratically near the root but may fail when
    the derivative is near zero or when the function is highly nonlinear.
    Use bisection_solve as a fallback when this method does not converge.

    Inputs:
      f         -- callable; the residual function to solve, with signature
                   f(Te, H1, uw, Ta_in, Td, sat_vp_ta)
      df        -- callable; the analytical derivative of f with respect to Te,
                   with the same signature as f
      x0        -- float; initial estimate of the root (deg C)
      H1        -- float; combined radiative heat input term (W/m2)
      uw        -- float; wind speed (m/s)
      Ta_in     -- float; air temperature (deg C)
      Td        -- float; dewpoint temperature (deg C)
      sat_vp_ta -- float; saturation vapor pressure at air temperature (mbar)
      tol       -- float; convergence tolerance; iteration stops when the change
                   between successive estimates is less than tol (default 1e-2)
      max_iter  -- int; maximum number of Newton-Raphson iterations before raising
                   an error (default 1000)

    Output:
      Returns the estimated root (equilibrium temperature, deg C) as a float.
      Raises ValueError if the derivative is near zero or if the maximum number
      of iterations is exceeded.
    """
    
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
    """
    Computes the decimal day-of-year for a given Python datetime object,
    including a timezone offset correction from the tz_offset module.

    Day 1.0 corresponds to midnight on January 1st. Fractional days represent
    the elapsed portion of the day. The tz_offset.days correction accounts for
    the difference between local time and the reference time zone.

    Inputs:
      tin -- Python datetime object representing the local date and time

    Output:
      Returns a float representing the decimal day-of-year
      (e.g., 1.5 = noon on Jan 1, 32.0 = midnight on Feb 1),
      adjusted by tz_offset.days.
    """
    
    # Define the first moment of the current year
    foy = dt.datetime(tin.year, 1, 1, 0, 0, 0)
    
    # Convert elapsed time since Jan 1 into fractional day-of-year
    ddoy = (tin - foy).total_seconds() / 86400. + 1. + tz_offset.days
    
    # Return decimal day-of-year value
    return ddoy


def solar_alt_angle(tin):
    """
    Computes the solar altitude angle (degrees above the horizon) for a hardcoded
    site location (Natoma/Folsom area: 38.1 N, 121.8 W) at a given local datetime.

    Uses a simplified astronomical formula based on solar declination, latitude,
    and the solar hour angle. The time zone reference longitude is set to 120 W
    (Pacific Standard Time meridian).

    NOTE: Latitude, longitude, and time_zone_longitude are currently hardcoded.
    A TODO exists in the module header to make these configurable.

    Inputs:
      tin -- Python datetime object representing the local date and time

    Output:
      Returns a float representing the solar altitude angle in degrees.
      Negative values indicate the sun is below the horizon (nighttime).
    """
    
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
    """
    Estimates the shortwave (solar) radiation reflectance fraction of the water
    surface as a function of solar altitude angle and cloud cover fraction.

    Uses empirical power-law coefficients that vary by cloud cover category.
    Returns 1.0 (full reflection) at night when the solar altitude is at or
    below the horizon. The result is clamped to a maximum of 1.0.

    Inputs:
      tin      -- Python datetime object representing the local date and time,
                  used to compute the solar altitude angle
      cloud_in -- float; fractional cloud cover (0.0 = clear sky, 1.0 = fully overcast)

    Output:
      Returns a float in the range [0, 1.0] representing the fraction of incident
      shortwave radiation reflected by the water surface.
      Returns 1.0 at night (solar altitude <= 0).
    """
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
    """
    Estimates the downward longwave (thermal infrared) radiation flux at the
    water surface using an empirical atmospheric emissivity model.

    Atmospheric emissivity is estimated from air temperature and cloud cover
    using the formula:
        emissivity = 0.937e-5 * (1 + 0.17 * cloud^2) * Tair_K^2

    The net downward flux accounts for a small surface reflection fraction (3%).

    Inputs:
      air_temp_in -- float; air temperature (deg C)
      cloud_in    -- float; fractional cloud cover (0.0 = clear sky, 1.0 = overcast)

    Output:
      Returns a float representing the net downward longwave heat flux at the
      water surface (W/m2).
    """
    
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
    """
    Computes the latent heat of vaporization of water as a linear function of
    temperature using an empirical approximation.

    The approximation is:
        Lw = 1000 * (2499 - 2.36 * T)    [J/kg]

    Inputs:
      temp_in -- float; water or air temperature (deg C)

    Output:
      Returns a float representing the latent heat of vaporization (J/kg).
    """
    
    # Temperature-dependent latent heat of vaporization [J/kg]
    return 1000. * (2499. - 2.36 * temp_in)


def sat_water_vapor_pres(temp_in):
    """
    Computes the saturation vapor pressure of water at a given temperature
    using the Magnus-Tetens empirical formula:

        e_s = Magnus_a * exp(Magnus_b * T / (T + Magnus_c))

    where Magnus_a = 6.1078, Magnus_b = 17.27, Magnus_c = 237.7 (module constants).

    Inputs:
      temp_in -- float; temperature (deg C)

    Output:
      Returns a float representing the saturation vapor pressure (mbar).
    """
    
    # Magnus-Tetens empirical formula for saturation vapor pressure
    return Magnus_a * math.exp(Magnus_b * temp_in / (temp_in + Magnus_c))


def fTe(Te, H1, uw, Ta_in, Td, sat_vp_ta):
    """
    Computes the net heat flux residual at a candidate equilibrium temperature Te.

    This function implements the surface energy balance equation for a water body.
    The equilibrium temperature is the value of Te that makes fTe = 0, balancing
    incoming radiation (H1) against outgoing longwave radiation and evaporative
    heat loss. Used as the target function for the Newton-Raphson solver.

    The energy balance form is:
        fTe = H1 - H_lw_up(Te) + H_evap(Te)

    where H_lw_up is the upward longwave emission and H_evap is the evaporative
    heat flux (latent + sensible, via Bowen ratio approximation).

    Inputs:
      Te        -- float; candidate equilibrium temperature (deg C)
      H1        -- float; combined radiative heat input (shortwave + downward longwave) (W/m2)
      uw        -- float; wind speed (m/s)
      Ta_in     -- float; air temperature (deg C)
      Td        -- float; dewpoint temperature (deg C)
      sat_vp_ta -- float; saturation vapor pressure at air temperature (mbar)

    Output:
      Returns a float representing the net heat flux residual (W/m2).
      A value of zero indicates that Te is the equilibrium temperature.
    """
    
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
    """
    Computes the analytical derivative of fTe with respect to Te (dF/dTe),
    used by the Newton-Raphson solver to refine the equilibrium temperature estimate.

    Differentiates all Te-dependent terms in the surface energy balance:
      - Upward longwave radiation (Stefan-Boltzmann, linearized)
      - Latent heat flux (via temperature-dependent Lw and beta)
      - Sensible heat flux (via Bowen ratio and beta gradient)

    Inputs:
      Te        -- float; candidate equilibrium temperature (deg C)
      H1        -- float; combined radiative heat input (W/m2); not directly used
                   in the derivative but required to match the solver's call signature
      uw        -- float; wind speed (m/s)
      Ta_in     -- float; air temperature (deg C)
      Td        -- float; dewpoint temperature (deg C)
      sat_vp_ta -- float; saturation vapor pressure at air temperature (mbar)

    Output:
      Returns a float representing the analytical derivative of fTe with respect
      to Te (W/m2/degC). Used by newton_raphson_solve to compute the Newton step.
    """

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
    """
    Reformulates the surface energy balance as a root-finding problem suitable
    for the bisection solver by returning Te - fTe(Te, ...).

    The bisection solver (bisection_solve) requires a function whose sign changes
    at the root. This wrapper converts the energy balance residual into the form
    needed by that solver.

    Inputs:
      Te        -- float; candidate equilibrium temperature (deg C)
      H1        -- float; combined radiative heat input (W/m2)
      uw        -- float; wind speed (m/s)
      Ta_in     -- float; air temperature (deg C)
      Td        -- float; dewpoint temperature (deg C)
      sat_vp_ta -- float; saturation vapor pressure at air temperature (mbar)

    Output:
      Returns a float equal to Te - fTe(Te, H1, uw, Ta_in, Td, sat_vp_ta).
      A value of zero indicates that Te is the equilibrium temperature.
    """
    
    # Reformulate root-finding problem for bisection solver
    return Te - fTe(Te, H1, uw, Ta_in, Td, sat_vp_ta)

def fTe_abs(Te, H1, uw, Ta_in, Td, sat_vp_ta):
    """
    Returns the absolute value of the surface energy balance residual fTe,
    used as a scalar objective function for scipy-based minimization.

    NOTE: scipy.optimize is not available in Jython. This function is retained
    for reference and potential future use in a CPython environment only.

    Inputs:
      Te        -- float; candidate equilibrium temperature (deg C)
      H1        -- float; combined radiative heat input (W/m2)
      uw        -- float; wind speed (m/s)
      Ta_in     -- float; air temperature (deg C)
      Td        -- float; dewpoint temperature (deg C)
      sat_vp_ta -- float; saturation vapor pressure at air temperature (mbar)

    Output:
      Returns a non-negative float equal to |fTe(Te, H1, uw, Ta_in, Td, sat_vp_ta)|.
      Minimizing this value finds the equilibrium temperature.
    """
    
    # Absolute objective function used by scipy minimization
    return abs(fTe(Te, H1, uw, Ta_in, Td, sat_vp_ta))


def equilibrium_temp(dtt1, at, cl, sr, ws, td, te_guess, type='nr'):
    """
    Computes the equilibrium (surface) temperature of a water body at a single
    time step by solving the surface energy balance equation.

    Assembles the net radiative heat input from shortwave and downward longwave
    fluxes, then solves for the temperature Te at which the total heat flux
    (radiation + evaporation + sensible heat) is zero.

    Three solution methods are supported:
      'nr'    -- Newton-Raphson (fastest, may fail on nonlinear regions)
      'bs'    -- Bisection over [-30, 55] deg C (slower but robust)
      'scipy' -- scipy.optimize.minimize on |fTe| (CPython only; not available in Jython)

    Inputs:
      dtt1     -- Python datetime object for the current time step (used for solar angle)
      at       -- float; air temperature (deg C)
      cl       -- float; fractional cloud cover (0.0 = clear, 1.0 = overcast)
      sr       -- float; incoming shortwave irradiance (W/m2)
      ws       -- float; wind speed (m/s)
      td       -- float; dewpoint temperature (deg C)
      te_guess -- float; initial guess for the equilibrium temperature (deg C)
      type     -- string; solution method: 'nr', 'bs', or 'scipy' (default 'nr')

    Output:
      Returns a float representing the equilibrium water surface temperature (deg C).
      Raises ValueError if the solution method is not recognized.
    """
    
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
    """
    Computes the equilibrium water surface temperature for a time series of
    meteorological inputs by iterating over each time step and applying the
    surface energy balance model.

    At each time step, the Newton-Raphson solver is attempted first for speed.
    A bisection solution is always computed as a fallback. If Newton-Raphson
    fails to converge or produces a result more than 1.0 deg C from the bisection
    estimate, the bisection result is used instead.

    The initial guess for the first time step is the air temperature. Subsequent
    time steps use the previous equilibrium temperature as the initial guess.

    Adapted from a NumPy/SciPy routine by Steve, Ben, and Scott, based on:
    "Stratification and heat transfer in lakes and reservoirs" in
    "Hydrodynamics and Transport for Water Quality Modeling"
    by Martin and McCutcheon (pg. 373).

    Inputs:
      dtt -- list of Python datetime objects; one per time step
      at  -- list of floats; air temperature (deg C), one per time step
      cl  -- list of floats; fractional cloud cover (0.0-1.0), one per time step
      sr  -- list of floats; shortwave irradiance (W/m2), one per time step
      td  -- list of floats; dewpoint temperature (deg C), one per time step
      ws  -- list of floats; wind speed (m/s), one per time step

    Output:
      Returns a list of floats containing the estimated equilibrium water surface
      temperature (deg C) for each time step. Length equals len(dtt).
    """
    
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


