from copy import copy

import numpy as np

class StationReader:
    """
    StationReader
    Purpose:    Parse the CSV file containing the station location information for all the surface
                    and upper air stations.
    Started:    19 June 2010 by Tim Supinie (tsupinie@ou.edu)
    Completed:  [not yet]
    Modified:   [not yet]
    """

    # List of keys (fields) in the CSV files
    _keys = ['State', 'Name', 'ID', 'Latitude', 'Longitude', 'Elevation']

    # Static private dictionary of keys pointing to their python data types
    _data_types = {'State':str, 'Name':str, 'ID':str, 'Latitude':float, 'Longitude':float, 'Elevation':int}

    def __init__(self, file_name):
        """
        __init__() [public]
        Purpose:    Constructor for the StationReader class.  Parses the CSV file containing the
                        station location information
        Parameters: file_name [type=str]
                        Name of the file to parse
        """

        self._stations = []
        convert_fns = [ StationReader._data_types[k] for k in StationReader._keys ]

        for line in open(file_name, 'rb'):
            # Loop through the file.
            values = line.strip().split(',')

            # Convert each value to its respective python data type.
            values = [ convert_fns[idx](values[idx]) for idx in range(len(values)) ]

            # Create a dictionary containing the values for each station and append it to the list
            self._stations.append(dict(zip(StationReader._keys, values)))

        return

    def searchByID(self, value):
        """
        searchByID() [public]
        Purpose:    Search for station information by ID
        Parameters: value [type=str]
                        What to search for
        Returns:    A dictionary containing the station's location information
        """

        return self._search(ID=value)

    def searchByLatLon(self, bounds):
        """
        searchByLatLon() [public]
        Purpose:    Search for station information by location
        Parameters: bounds [type=dict]
                        A dictionary pointing to the latitude and longitude information for the search.  The information can either be a point or a tuple for bounds.
        Returns:    Either the station at the given latitude and longitude or the stations within the latitude and longitude bounds.
        """
        return self._search(Latitude=bounds['lat'], Longitude=bounds['lon'])

    def searchByDomain(self, map):
        """
        searchByDomain() [public]
        Purpose:    Search for station information in a given domain. 
        Parameters: map [type=Basemap]
                        A Basemap object corresponding to a domain in which to search for stations.
        Returns:    The station information for stations inside the given domain.
        """
        stn_lats = np.array([ stn['Latitude'] for stn in self._stations ])
        stn_lons = np.array([ stn['Longitude'] for stn in self._stations ])

        ll_corner_x, ll_corner_y = map(map.llcrnrlon, map.llcrnrlat)
        ur_corner_x, ur_corner_y = map(map.urcrnrlon, map.urcrnrlat)

        stn_xs, stn_ys = map(stn_lons, stn_lats)

        idxs = np.where((ll_corner_x < stn_xs) & (stn_xs < ur_corner_x) & (ll_corner_y < stn_ys) & (stn_ys < ur_corner_y))
        return [ self._stations[idx] for idx in idxs[0] ]

    def _search(self, **kwargs):
        """
        _search() [private]
        Purpose:    Search for a station by the key given
        Parameters: key [type=str]
                        Which data field to search
                    value [type=str,tuple]
                        What to search for.  If the type is string, it will match the value exactly.  If the type is tuple, it will match all objects between the elements of the tuple.
        Returns:    A dictionary containing the station's location information
        """

        # We want to search for things that have bounds first, then search for individual items
        search_order = [ key for key, val in kwargs.iteritems() if type(val) in [ list, tuple ] ]
        search_order.extend([ key for key, val in kwargs.iteritems() if type(val) not in [ list, tuple ] ])

        sub_list = copy(self._stations)

        # Create a list of values that go with the key (field) and find the index of the value in the list.
        for key in search_order:
            if type(kwargs[key]) in [ list, tuple ]:
                sub_list = [ stn for stn in sub_list if stn[key] >= kwargs[key][0] and stn[key] <= kwargs[key][1] ]
            else:
                sub_list = [ stn for stn in sub_list if stn[key] == kwargs[key] ]

        if len(sub_list) == 0:
            return
        elif len(sub_list) == 1:
            return sub_list[0]
        return sub_list

if __name__ == "__main__":
    sr = StationReader("sfc_stations.csv")
    print sr.searchByID("KOUN")
    print sr.searchByLatLon({'lat':(34, 36), 'lon':(-99, -96)})
