from hootpyStationPlot import HootPyStationPlot
from metar import METAR

class SurfaceStationPlot(HootPyStationPlot):
    """
    SurfaceStationPlot
    Purpose:    Handles the plotting of meteorological station plots on surface maps
    Started:    17 June 2010 by Tim Supinie (tsupinie@ou.edu)
    Completed:  25 June 2010
    Modified:   [not yet]
    """

    # Static private list of valid attributes for the station plot
    _valid_attrs = ['timh', 'timm', 'tmpc', 'tmpf', 'tmpk', 'dwpc', 'dwpf', 'dwpk', 'pslm', 'psli', 'wdir', 'wsdk', 'wsdm', 'wsdi', 'skyc', 'wsym']

    def __init__(self, location_info, **kwargs):
        """
        __init__()
        Purpose:    Constructor for the SurfaceStationPlot class.
        Parameters: location_info [type=dict]
                        Dictionary containing the location information for the station
                    Key-value pairs (station variable to value, e.g. tmpf=50)
        """

        # Call superclass's constructor
        super(SurfaceStationPlot, self).__init__(location_info, SurfaceStationPlot._valid_attrs, **kwargs)
        return

    def parseFromMETAR(self, metar_string):
        """
        parseFromMETAR() [public]
        Purpose:    Initialize the station plot by parsing a METAR and extracting the 
                        relevant variables.
        Parameters: metar_string [type=str]
                        A string containing the METAR to be parsed.
        Returns:    [nothing]
        """

        metar = METAR(metar_string)
        if (metar.get('station_id') != self._station_id):
            raise ValueError("METAR is for " + metar.get('station_id') + ", not " + self._station_id)
        self.set(timh=metar.get('time_hours'),
                 timm=metar.get('time_minutes'),
                 wdir=metar.get('wind_direction'),
                 wsdk=metar.get('wind_speed'),
                 skyc=metar.get('sky_conditions'),
                 wsym=metar.get('current_weather'),
                 tmpc=metar.get('temperature'),
                 dwpc=metar.get('dewpoint'),
                 psli=metar.get('pressure'))

        return

if __name__ == "__main__":
    import sys
    from data.StationReader import StationReader

    sr = StationReader('data/sfc_stations.csv')
    ssp = SurfaceStationPlot(sr.searchByID('KOUN'), tmpf=50, dwpf=32, wdir=180, wsdk=10, psli=29.92)
    ssp._dump(sys.stdout)
    ssp.parseFromMETAR("METAR KOUN 240212Z AUTO 16008KT 10SM CLR 32/21 A2998 RMK A02=")
    ssp._dump(sys.stdout)
    ssp = SurfaceStationPlot(sr.searchByID('KMDW'))
    ssp.parseFromMETAR("METAR KMDW 232256Z 30028G39KT 1SM R31C/P6000FT -TSRA SQ FEW028 BKN036CB OVC042 26/20 A2985 RMK AO2 PK WND 30039/2255 RAB55 PRESRR FRQ LTGICCG OHD TS OHD MOV E-SE P0003=")
    ssp._dump(sys.stdout)
    ssp = SurfaceStationPlot(sr.searchByID('KOMA'))
    ssp.parseFromMETAR("METAR KOMA 231253Z VRB04KT -TSSN BR OVC001 M02/M02 A2985 RMK A02=")
    ssp._dump(sys.stdout)
