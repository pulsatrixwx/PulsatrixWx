"""
derived.py
Purpose:    Utility package of functions for calculating common derived parameters
Started:    09 September 2010 by David John Gagne (djgagne#@ou.edu)
Completed:  [not yet]
Modified:   [not yet]
"""
import numpy as np

#from callCount import FunctionCallCount

eo = 611.0 # Pa
g = 9.806 # m/s^2
Lv = 2.5e6 # J/kg
Rv = 461.5 # J/(K*kg)
Rd = 287.047 # J/(K*kg)
cp = 1005. # J/(K*kg)
Po = 100000. # Pa
epsilon = Rd / Rv
sigma = 5.67e-8
_static_units = { 
    'potentialTemperature':{'return':'K',    'args':['K', 'Pa'] },
    'equivalentPotentialTemperatureQV':{'return':'K', 'args':['K','Pa','kg kg-1']},
    'equivalentPotentialTemperatureRH':{'return':'K', 'args':['K','Pa','']},
    'dewpointQV':          {'return':'K',     'args':['K', 'kg kg-1', 'Pa'] },
    'dewpointRH':        {'return':'K',     'args':['K', ''] },
    'dewpointE':         {'return':'K',     'args':['Pa'] },
    'temperature':       {'return':'K',     'args':['K', 'Pa'] }, 
    'mslp':              {'return':'Pa',    'args':['Pa', 'm', 'K', 'kg kg-1'] }, 
    'virtualTemperature':{'return':'C',     'args':['K', 'kg kg-1'] }, 
    'relativeHumidity':  {'return':'%',     'args':['K', 'kg kg-1', 'Pa'] },
    'lclHeightRH':       {'return':'m',     'args':['K', ''] },
    'lclHeightDP':       {'return':'m',     'args':['K', 'K'] },
    'pwtrToMM':          {'return':'m',     'args':['kg m-2'] },

    'g':                 {'return':'m s-2', 'args':[] },
    'radToTemp':         {'return':'K',     'args':['W m-2'] },
    'tempToRad':         {'return':'W m-2',     'args':['K'] }
}

#@FunctionCallCount
def dewpointQV(temperature_k,qv,P):
    """
    dewpoint()
    Purpose:    Calculate the dewpoint temperature
    Parameters: temperature_k - Temperature in Kelvin
                qv - Water vapor mixing ratio in kg/kg
    Returns:    Dewpoint in degrees Kelvins 
    """
    return np.power(1/273.0 - Rv/Lv * np.log(vaporPressure(qv,P) / eo), -1)

def dewpointRH(temperature, rel_hum):
    vapor_pres = rel_hum * satVaporPressure(temperature)
    return np.minimum(dewpointE(vapor_pres), temperature)

def dewpointE(vapor_pressure):
    return 1 / (1 / 273.15 - Rv / Lv * np.log(vapor_pressure / eo))

#@FunctionCallCount
def temperature(theta, pressure):
    """
    temperature()
    Purpose:    Calculate temperature from potential temperature
    Parameters: theta - potiential temperature (in K)
                pressure - air pressure (in Pa)
    Returns:    Air temperature in Kelvins
    """
    return theta * np.power(pressure / Po, Rd / cp)

#@FunctionCallCount
def mslp(station_pressure,station_height,temperature_k,qv):
    """
    mslp()
    Purpose:    Calculate the mean sea level pressure
    Parameters: station_pressure - surface pressure in Pa
                station_height - Height above sea level in meters
                temperature_k - Temperature in Kelvin
                qv - Water vapor mixing ratio in kg/kg
    Returns:    Mean Sea Level Pressure in Pa
    """
    return station_pressure * np.exp( ( g * station_height) / ( Rd * virtualTemperature(temperature_k,qv)))

#@FunctionCallCount
def virtualTemperature(temperature_k,qv):
    """
    virtualTemperature()
    Purpose:    Calculate the virtual temperature
    Parameters: temperature_k - Temperature in Kelvin
                qv - Water vapor mixing ratio in kg/kg
    Returns:    Virtual Temperature in degrees Celsius
    """
    virtual_temp_const = 1.0/(epsilon) - 1
    return (1 + virtual_temp_const * specificHumidity(qv)) * temperature_k

def potentialTemperature(temperature,pressure):
    """
    potentialTemperature()
    Purpose:  Calculate the potential temperature from a given temperature and pressure
    Parameters:     temperature - Temperature in Kelvin
                    pressure  - Pressure in Pa
    Returns:  Potential temperature in the units of temperature
    """
    r_cp = Rd / cp
    return temperature * (100000.0/pressure) ** r_cp

def equivalentPotentialTemperatureQV(temperature,pressure,qv):
    """
    potentialTemperature()
    Purpose:  Calculate the potential temperature from a given temperature and pressure
    Parameters:     temperature - Temperature in Kelvin
                    pressure  - Pressure in Pa
                    qv - Mixing ratio in kg/kg
    Returns:  Potential temperature in the units of temperature
    """
    

    thtm = temperature * (100000.0/pressure) ** (2.0/7.0 * (1 - (.28 * qv)))
    lclt = lclTemperature(temperature,pressure,qv)
    thetae = thtm * np.exp((3.376 / lclt - .00254) * (qv/0.001 * (1 + 0.81 * qv)))

    return np.where(qv == 0, potentialTemperature(temperature, pressure), thetae)

def equivalentPotentialTemperatureRH(temperature, pressure, rh):
    qv = mixingRatio(rh * satVaporPressure(temperature), pressure)
    return equivalentPotentialTemperatureQV(temperature, pressure, qv)

def lclTemperature(temperature,pressure,qv):
    """Using GEMPAK formula"""
    dewp = dewpointQV(temperature,qv,pressure)
    return (1.0 / (1.0 / (dewp - 56.0) + np.log(temperature/dewp) / 800.0)) + 56.0

def lclHeightRH(temperature, rel_hum):
    dewp = dewpointRH(temperature, rel_hum)
    return lclHeightDP(temperature, dewp)

def lclHeightDP(temperature, dewp):
    return np.maximum(125 * (temperature - dewp), np.zeros(dewp.shape))

#@FunctionCallCount
def specificHumidity(qv):
    """
    specificHumidity()
    Purpose:    Calculate the specific humidity from mixing ratio
    Parameters: qv - Water vapor mixing ratio in kg/kg
    Returns:    Specific humidity in kg/kg
    """    
    return qv/(1 + qv)

#@FunctionCallCount
def vaporPressure(qv,P):
    """
    vaporPressure()
    Purpose:    Calculate the vapor pressure (e) from mixing ratio and pressure
    Parameters: qv - Water vapor mixing ratio in kg/kg, P - pressure in Pa
    Returns:    vapor pressure in Pa
    """
    return (P * qv) / ( qv + epsilon)

def mixingRatio(vapor_pres, pres):
    return (epsilon * vapor_pres) / (pres - vapor_pres)

#@FunctionCallCount
def satVaporPressure(temperature_k):
    """
    satVaporPressure()
    Purpose:    Calculate the saturated vapor pressure (es) from temperature
    Parameters: temperature_k - Temperature in Kelvin
    Returns:    saturated vapor pressure in Pa
    """
    return eo * np.exp(Lv/Rv * (1 / 273.0 - 1/temperature_k))

#@FunctionCallCount
def relativeHumidity(temperature_k,qv,P):
    """
    relativeHumidity()
    Purpose:    Calculate the relative humidity from model variables
    Parameters: temperature_k - Temperature in Kelvin
                qv - Water vapor mixing ratio in kg/kg
                P - pressure in Pascals
    Returns:    relative humidity as a percentage
    """
    return vaporPressure(qv,P) / satVaporPressure(temperature_k) * 100

def pwtrToMM(precip_water):
    return precip_water / 1000

def vectorDecompose(mag, dir):
    """
    vectorDecompose()
    Purpose:    Decompse vector from magnitude and direction to u and v components
    Parameters: mag - Magnitude of the vector
                dir - Bearing of the vector in degrees (0 degrees is north, increasing clockwise).
    Returns:    A tuple containing the u and v components of the given vector
    """
    return (mag * np.cos(np.pi * (90 - dir + 180) / 180.0), mag * np.sin(np.pi * (90 - dir + 180) / 180.0))

#@FunctionCallCount
def vectorMagnitude(u, v):
    """
    vectorMagnitude()
    Purpose:    Return the magnitude of the given vector
    Parameters: u - The u component of the vector
                v - The v component of the vector
    Returns:    The magnitude of the vector
    """
    return np.sqrt(u ** 2 + v ** 2)

def layerAverage(grid):
    return grid.mean(axis=1)

def layerMin(grid):
    return grid.min(axis=1)

def layerMax(grid):
    return grid.max(axis=1)

def layerStd(grid):
    return grid.std(axis=1)

def layerDifference(grid):
    """
    layerDifference()
    Purpose:    Calculate the difference in layers in a 4D gridded file
    Parameters: grid - A 4d array where the order is assumed to be t, z, y, x
    Returns:    The array that contains the layer difference
    """
    dif_grid = grid[:,0] - grid[:,-1]
    return dif_grid

def maxCloudArea(*args):
    all_clouds = np.vstack(args)
    return all_clouds.max(axis=0)

def radToTemp(radiance):
    return (radiance / sigma ) ** 0.25

def tempToRad(temperature):
    return sigma * temperature ** 4

#@FunctionCallCount
def _units(func, *args):
    """
    _units()
    Purpose:    Returns information about what the above ("HootPy") functions take and return for units.
    Parameters: args [type=*list]
                    Units of the arguments that will be given to the function.  Needed only for things like vectorMagnitude(), where
                        the units of the output depend on the units of the input.
    Returns:    A dictionary giving the units returned ('return' key) and expected ('args' key).
    """
    try:
        return _static_units[func]
    except KeyError:
        _units = {
              'vectorMagnitude':   {'return':args[0], 'args':[ args[0] for a in args ] },
              'vectorDecompose':   {'return':args[0], 'args':[ args[0] for a in args ] },
              'maxCloudArea':      {'return':'',      'args':[ '' for a in args ] },
              'layerAverage':     {'return':args[0], 'args':[ args[0] for a in args ] },
              'layerMin':     {'return':args[0], 'args':[ args[0] for a in args ] },
              'layerMax':     {'return':args[0], 'args':[ args[0] for a in args ] },
              'layerDifference':     {'return':args[0], 'args':[ args[0] for a in args ] },
              'layerStd':     {'return':args[0], 'args':[ args[0] for a in args ] },
        }

        return _units[func]

def main():
    print "Testing"
    import Nio as nio
    data = nio.open_file('/home/wrfuser/hootpy/data/wrfout_d01_PLEV.nc')
#   print mslp(data.variables['PSFC'][:],data.variables['HGT'][:],data.variables['T2'][:],data.variables['Q2'][:])
    dewp = dewpoint(data.variables['T'][0,:],data.variables['QVAPOR'][0,:],data.variables['P'][:]*100)
    print dewp.max()
    print dewp.min()
    print dewp.mean()

if __name__ == "__main__":
    main()
