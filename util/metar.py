import re

class METAR:
    """
    METAR
    Purpose:    Parses and keeps track of variables related to METAR observations
    Started:    25 June 2010 by Tim Supinie (tsupinie@ou.edu)
    Completed:  26 June 2010
    Modified:   [not yet]
    """

    # Private static list giving the priority for reporting the sky conditions (higher priorities are near
    #    the front of the list).
    _skyc_priority = ['OVC', 'BKN', 'SCT', 'FEW', 'CLR', 'SKC']

    # Private static list giving the priority for reporting precipitation types (higher priorities are near
    #    the front of the list).
    _wsym_priority = ['PL', 'SN', 'RA', 'SG', 'DZ', 'IC', 'BR', 'FG', 'VA', 'DU', 'HZ']

    # Private static list giving the priority for reporting precipitation descriptors (higher priorities
    #    are near the front of the list).
    _wmod_priority = ['TS', 'FZ', 'BL', 'DR', 'SH', 'MI', 'PR', 'BC']

    # Privat static list giving the priority for reporting precipitation intensity (higher priorities are
    #    near the front of the list).
    _wint_priority = ['+', '', '-']

    def __init__(self, metar_string):
        """
        __init__()
        Purpose:    Constructor for the METAR class.  Calls the _parse() method on the METAR.
        Parameters:    metar_string
                        A string containing the METAR to be parsed.
        """

        # Parse the METAR
        self._parse(metar_string)
        return

    def get(self, attr):
        """
        get() [public]
        Purpose:    Get a meteorological value from the METAR.
        Parameters:    attr [type=str]
                        Which attribute to return.
        Returns:    The value of the attribute for this METAR.
        """
        if hasattr(self, attr):
            return getattr(self, attr)
        else:
            return None

    def _parse(self, metar_string):
        """
        _parse() [private]
        Purpose:    Parse the metar string and extract relevant variables.
        Parameters:    metar_string [type=str]
                        A string containing the METAR to be parsed.
        Returns:    [nothing]
        """

        # Split the observation from the remarks (they're handle differently.  And the remarks are weird.)
        observation, remarks = metar_string.split('RMK')

        # Note: before you try to understand this code, you may want to read up on regular expressions (regex).
        #    Regex information and its Python API can be found here: http://docs.python.org/library/re.html

        # Match the METAR string at the very beginning (or SPECI, indicating a special observation).
        #    The (?=[\\s]) looks for spaces ahead of the string, but doesn't include it in the match (necessary later).
        metr_match = re.search('(METAR|SPECI)(?=[\\s])', observation)

        # Match the ICAO identifier (K followed by 3 characters).  The (?<=[\\s]) looks for spaces behind
        #    the string, but doesn't include it in the match (again, necessary later).
        icao_match = re.search('(?<=[\\s])(K[\\w]{3})(?=[\\s])', observation)

        # Match the time stamp (three sets of two digits followed by a Z.  String multiplication is done on the 
        #    two-digit sets because the regex engine needs three sets of parenthesis to capture three values.)
        time_match = re.search('(?<=[\\s])' + '([\\d]{2})' * 3 + 'Z(?=[\\s])', observation)

        # Match the wind observation (a set of three digits or "VRB" followed by a set of two digits, followed by an optional
        #    "gust" cluster, which is a G and two digits, followed by a KT).
        wind_match = re.search('(?<=[\\s])([\\d]{3}|VRB)([\\d]{2})(?:G([\\d]{2}))?KT(?=[\\s])', observation)

        # Match the visibility (not strictly necessary for the station plot, and untested.  I may have broken it.)
        visb_match = re.search('(?<=[\\s])([\\d][\\s]?[\\d]?/?[\\d]?)SM(?=[\\s])', observation)

        # Match the sky conditions (any of the strings in _skyc_priority followed by three digits for the ceiling
        #    height, followed by an optional CB cluster).
        skyc_matches = re.findall('(?<=[\\s])(' + "|".join(METAR._skyc_priority) + 
            ')([\\d]{3})?(?:CB)?(?=[\\s])', observation)

        # Match the weather observation (any two of the strings in _wmod_priority, both of which are optional,
        #    followed by any two of the strings in _wsym_priority, one of which is optional).
        wsym_matches = re.findall('(?<=[\\s])(\\+|\\-)?' + ('(' + "|".join(METAR._wmod_priority) + ')?') * 2 +
            ('(' + "|".join(METAR._wsym_priority) + ')') * 2 + '?(?=[\\s])', observation)

        # Match the temperature and dewpoint observations (two two digit numbers, each of which may be preceded
        #    by an M, separated by a slash).
        temp_match = re.search('(?<=[\\s])(M?[\\d]{2})/(M?[\\d]{2})(?=[\\s])', observation)

        # Match the pressure observation (an A followed by four digits).
        pres_match = re.search('(?<=[\\s])A([\\d]{4})(?=[\\s])', observation)

        if not metr_match:
            # METAR doesn't start with the METAR or SPECI, so it's not a valid METAR
            raise ValueError(metar_string + ' is not a valid METAR (no METAR or SPECI header).')
        if icao_match:
            setattr(self, 'station_id', icao_match.group(1))
        else:
            # METAR doesn't have an ICAO station identifier, so it's not a valid METAR
            raise ValueError(metar_string + ' is not a valid METAR (no ICAO station identifier).')
        if time_match:
            # Set the observation time
            setattr(self, 'time_hours', int(time_match.group(2)))
            setattr(self, 'time_minutes', int(time_match.group(3)))
        if wind_match:
            # Set the wind observations
            wdir = wind_match.group(1)
            if (wdir != 'VRB'): wdir = int(wdir)
            setattr(self, 'wind_direction', wdir)
            setattr(self, 'wind_speed', int(wind_match.group(2)))
            if wind_match.lastindex == 3:
                setattr(self, 'wind_gust', int(wind_match.group(3)))
        if visb_match:
            visb_string = visb_match.group(1)
            visb = 0
            if ' ' in visb_string:
                whole, partial = visb_string.split(' ')
                numerator, denominator = partial.split('/')
                visb = int(whole) + float(numerator) / float(denominator)
            elif '/' in visb_string:
                numerator, denominator = visb_string.split('/')
                visb = float(numerator) / float(denominator)
            else:
                visb = int(visb_string)
                
            setattr(self, 'visibility', visb)
            pass
        if skyc_matches:
            # Parse out the sky cover and pick the most important one.  The most important has the lowest
            #    index in the _skyc_priority list.
            skyc = zip(*skyc_matches)[0]
            for cover in METAR._skyc_priority:
                if cover in skyc:
                    setattr(self, 'sky_conditions', cover)
                    break
        if wsym_matches:
            # Parse out the weather symbols and pick the most important one.  The most important is scored
            #    by index in the _wsym_priority list plus a factor for intensity.  The factor is -1 for a +
            #    intensity, 0 for no intensity and +1 for a - intensity (yeah, backwards, I know ...).  The
            #    most important weather symbol has the lowest score.
            scores = ([], [])
            for wsym in wsym_matches:
                intensity = wsym[0]
                for idx in range(len(METAR._wsym_priority)):
                    weather = METAR._wsym_priority[idx]
                    if weather in wsym:
                        scores[0].append(idx + METAR._wint_priority.index(intensity) - 1)
                        scores[1].append("".join(wsym))
                        break

            setattr(self, 'current_weather', scores[1][scores[0].index(min(scores[0]))])
        if temp_match:
            # Parse the temperature and dewpoint and set them.
            def parseTemp(temp):
                """
                parseTemp() [local]
                Purpose:    Parse METAR temperature to an actual value
                Parameters:    temp [type=str]
                                Temperature (in degrees Celsius, with negative values denoted by a preceding M)
                Returns:    Temperature of type float in degrees Celsius.
                """

                if temp[0] == 'M':
                    return -1 * int(temp[1:])
                else:
                    return int(temp)

            setattr(self, 'temperature', parseTemp(temp_match.group(1)))
            setattr(self, 'dewpoint', parseTemp(temp_match.group(2)))
        if pres_match:
            # Convert the pressure to in Hg and set it.
            setattr(self, 'pressure', float(pres_match.group(1)) / 100)
        return

if __name__ == "__main__":
#    m = METAR("METAR KOUN 240212Z AUTO 16008G17KT 10SM CLR 32/21 A2998 RMK A02=")
    m = METAR("METAR KMDW 232256Z 30028G39KT 1SM R31C/P6000FT -TSRA SQ FEW028 BKN036CB OVC042 26/20 A2985 RMK AO2 PK WND 30039/2255 RAB55 PRESRR FRQ LTGICCG OHD TS OHD MOV E-SE P0003=")
    print "Station:", m.get('station_id')
    print "Time: %d:%d UTC" % (m.get('time_hours'), m.get('time_minutes'))
    print "Wind Direction:", m.get('wind_direction')
    print "Wind Speed:", m.get('wind_speed')
    print "Wind Gust:", m.get('wind_gust')
    print "Visibility:", m.get('visibility')
    print "Sky Conditions:", m.get('sky_conditions')
    print "Current Weather:", m.get('current_weather')
    print "Temperature:", m.get('temperature')
    print "Dewpoint:", m.get('dewpoint')
    print "Pressure:", m.get('pressure')
