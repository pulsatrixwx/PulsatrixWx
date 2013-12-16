from util.units import *

from hootpyStationPlot import HootPyStationPlot

class UpperAirStationPlot(HootPyStationPlot):
    """
    UpperAirStationPlot
    Purpose:    Handles the plotting of meteorological station plots on upper air maps
    Started:    17 June 2010 by Tim Supinie (tsupinie@ou.edu)
    Completed:  [not yet]
    Modified:   [not yet]
    """

    # Static private list of valid attributes for the station plot
    _valid_attrs = ['tmpc', 'tmpf', 'tmpk', 'dwpc', 'dwpf', 'dwpk', 'hgtd', 'hgtm', 'wdir', 'wsdk', 'wsdm', 'wsdi']

    def __init__(self, location_info, **kwargs):
        """
        __init__()
        Purpose:    Constructor for the UpperAirStationPlot class.
        Parameters: location_info [type=dict]
                        Dictionary containing the location information for the station.
                    Key-value pairs (station variable to value, e.g. tmpf=50)
        """
        
        # Call superclass's constructor
        super(UpperAirStationPlot, self).__init__(location_info, UpperAirStationPlot._valid_attrs, **kwargs)
        return

    def parseFromSounding(self, sounding_string):
        """
        parseFromSounding() [public]
        Purpose:    Initialize the station plot by parsing a sounding and extracting the 
                        relevant variables.
        Parameters: sounding_string [type=str]
                        A string containing the sounding to be parsed.
        Returns:    [nothing]
        """
        return

if __name__ == "__main__":
    import sys
    from data.StationReader import StationReader

    sr = StationReader('data/ua_stations.csv')
    uasp = UpperAirStationPlot(sr.searchByID('KOUN'), tmpf=50, dwpf=32, wdir=180, wsdk=10, hgtm=5700)
    uasp._dump(sys.stdout)
