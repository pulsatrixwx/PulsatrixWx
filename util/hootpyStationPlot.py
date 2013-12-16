import re

from units import Units
from util import fatalError
import weatherSymbol as wsym

class HootPyStationPlot(object):
    """
    HootPyStationPlot
    Purpose:    Base class handling the plotting of meteorological station plots on surface and upper air maps
    Started:    16 June 2010 by Tim Supinie (tsupinie@ou.edu)
    Completed:  [not yet]
    Modified:   [not yet]
    """

    _internal_units = {'temp':'C', 'dewp':'C', 'pres':'mb', 'wnds':'kts', 'wndg':'kts', 'hght':'m'}

    _layout_keys = ['T', 'TR', 'R', 'BR', 'B', 'BL', 'L', 'TL']

    def __init__(self, location_info, valid_attrs, **kwargs):
        """
        __init__()
        Purpose:    Constructor for the HootPyStationPlot class.
        Parameters:    location_info [type=dict]
                        Dictionary containing the station's location information
                    Key-value pairs (station variable to value, e.g. tmpf=50)
        """

        self._valid_attrs = valid_attrs

        # Set station location information
        self._station_id = location_info['ID']
        self._latitude = location_info['Latitude']
        self._longitude = location_info['Longitude']
        self._elevation = location_info['Elevation']

        # Set type and station variables
        self.set(**kwargs)
        return

    def set(self, **kwargs):
        """
        set() [public]
        Purpose:    Sets station plot variables (such as temperature, dewpoint, etc.)
        Parameters:    Key-value pairs (station variable to value, e.g. tmpf=50.
        Returns:    [nothing]
        """

        re_value_string = re.compile("^(\\-?[\\d]+\\.?[\\d]*)[\\s]*([\\w\\d\\- ]+)$")

        for attr, value in kwargs.iteritems():
            if attr in self._valid_attrs:
                re_value_match = re_value_string.match(value)
                if re_value_match is not None:
                    obs_value, obs_units = re_value_match.groups()
                    obs_value = Units.convert(float(obs_value), obs_units, HootPyStationPlot._internal_units[attr])
                    setattr(self, "_%s" % attr, obs_value)
                else:
                    setattr(self, "_%s" % attr, value)
            else:
                fatalError("Attribute %s is not a valid attribute for this type of station plot." % attr)
        return

    def plot(self, layout):
        """
        plot() [public]
        Purpose:    Plots the station on an image.
        Parameters: layout [type=dict]
                        Dictionary that specifies where to put each value in the station plot
        Returns:    [nothing]
        """
        print "Plotting station plot for %s ..." % self._station_id
        for location in hootpyStationPlot._layout_keys:
            try:
                value = getattr(self, "_%s" % layout[location])
            except AttributeError:
                continue

            
        return

#   def _convert(self, attr, value):
#       """
#       _convert() [private]
#       Purpose:    Convert incoming values to the correct type (float) and internal units.
#       Parameters:    attr [type=str]
#                       What the value is (e.g. tmpf = 'Temperature in Fahrenheit')
#                   value [type=int,float,str]
#                       The value to be converted
#       Returns:    A tuple in which the first element is the converted attribute type and the 
#                       second element is the converted value.
#       """

#       return HootPyStationPlot._convert_attrs[attr], HootPyStationPlot._convert_fns[attr](value)

    def _dump(self, file):
        """
        _dump() [private]
        Purpose:    Debug function that dumps the meterological attributes to a file.
        Parameters:    file [type=file]
                        File object to write the output to.
        Returns:    [nothing]
        """

        print "Station %s: %f,%f, Elev.: %d" % (self._station_id, self._latitude, self._longitude, self._elevation)
        for attr in self._valid_attrs:
            if hasattr(self, "_%s" % attr):
                file.write("%s: %s\n" % (attr, str(getattr(self, '_' + attr))))

if __name__ == "__main__":
    import sys

    sr = StationReader('data/sfc_stations.csv')
    hpsp = HootPyStationPlot(sr.seachByID('KOUN'), ['temp', 'dewp', 'pres', 'wnds', 'wndg'], temp="50 F", dewp="32 F")
    hpsp._dump(sys.stdout)
