from datetime import datetime, timedelta
import calendar
from copy import copy

from math import sqrt, ceil

import numpy as np

import matplotlib
from matplotlib.transforms import Bbox
matplotlib.use('agg')
import pylab
from mpl_toolkits.basemap import Basemap
import pyproj

from util.config import Config
import util.derived as derived
from util.units import Units
from util.dataIO import DataIO
from util.util import any, uniquify, timeParse, warning, fatalError
from util.timezone import utc
from util.weatherSymbol import drawSkyCover, parseAndPlot
#from util.callCount import FunctionCallCount
from data.StationReader import StationReader
from hootpy import HootPy

class MapPlotter(HootPy):
    """
    MapPlotter
    Purpose:    Handles the plotting of data on cartesian maps (such as surface/upper air observations, model output,
                    radar data, and satellite data).
    Started:    14 June 2010 by Tim Supinie (tsupinie#@ou.edu)
    Completed:  [not yet]
    Modified:   [not yet]
    """

    _plot_order = [ 'sat_fn', 'grid_fn_fill', 'grid_fn_contour', 'grid_fn_vector', 'obs_fn' ]

    #@FunctionCallCount
    def __init__(self, config):
        """
        __init__()
        Purpose:    Constructor for the MapPlotter class.
        Parameters:    config [type=dictionary]
                        Dictionary containing the configuration parameters for this run.
        """

        super(MapPlotter, self).__init__(config)
        return

    #@FunctionCallCount
    def loadData(self):
        """
        loadData() [public]
        Purpose:    Handles the loading in of data.
        Parameters:    [none]
        Returns:    [nothing]
        """

        # Initialize dictionary of data sources to open netCDF files.
        self._netcdf = {}

        # Set up the basemap instance (used for specifying the bounds for loading data in)
        map_projection_conversion = {'ll_corner_longitude':'llcrnrlon',
                                     'll_corner_latitude':'llcrnrlat',
                                     'ur_corner_latitude':'urcrnrlat',
                                     'ur_corner_longitude':'urcrnrlon',
                                     'center_latitude':'lat_0',
                                     'center_longitude':'lon_0',
                                     'standard_latitude_1':'lat_1',
                                     'standard_latitude_2':'lat_2',
                                     }
        map_projection = {}
        for m,k in self._map_projection.iteritems():
            if m in map_projection_conversion.keys():
                map_projection[map_projection_conversion[m]] = k
            else:
                map_projection[m] = k
        if 'area_thresh' not in map_projection.keys():
            map_projection['area_thresh'] = 2000
        if 'resolution' not in map_projection.keys():
            map_projection['resolution'] = 'i'
        self._basemap = Basemap(**map_projection)

        # Loop over all the plot variables
        for plot in MapPlotter._plot_order:
            # Check to make sure the user specified this plot variable
            if not hasattr(self, "_%s" % plot):
                continue

            if plot == "obs_fn":
                self._clearCache()

            plot_info = getattr(self, "_%s" % plot)

            # Put plain dictionaries into a list so we can loop over it regardless of whether it was in a list in the 
            #    input file.
            if type(plot_info) == dict:
                setattr(self, "_%s" % plot, [ plot_info ])
                plot_info = getattr(self, "_%s" % plot)

            # Split the function lists in plots into different plots for each function.
#           MapPlotter._splitPlots(plot_info)

            # Loop over all the dictionaries (may be just a single iteration)
            for idx, plot_dict in enumerate(plot_info):

                try:
                    data_src = self._findAttribute(plot_dict, 'data_src')
                except AttributeError:
                    fatalError("No 'data_src' attribute found for plot %s, index %d." % (plot, idx))

                try:
                    data_scheme = self._findAttribute(plot_dict, 'data_scheme')
                except AttributeError:
                    fatalError("No 'data_scheme' attribute found for plot %s, index %d." % (plot, idx))

                try:
                    vertical_level = self._findAttribute(plot_dict, 'vertical_level')
                except AttributeError:
                    fatalError("No 'vertical_level' attribute found for plot %s, index %d." % (plot, idx))

                try:
                    plot_times = self._findAttribute(plot_dict, 'plot_times')
                except AttributeError:
                    fatalError("No 'plot_times' attribute found for plot %s, index %d." % (plot, idx))

                try:
                    data_src_map = self._data_src_map[data_scheme]
                except KeyError:
                    fatalError("Error: Data scheme '%s' unknown." % data_scheme)

                # Don't open a file we've already opened
                if data_src not in self._netcdf:
                    print "Opening file %s ..." % data_src
                    self._netcdf[data_src] = DataIO(data_src)

                t_coord, plot_dict['z_coord'] = self._loadTemporalVerticalCoords(self._netcdf[data_src], plot_dict)

                if data_src_map['internal_time'] == 'absolute':
                    if data_src_map['time_format'] == "integer":
                        seconds = Units.convert(t_coord, data_src_map['time_unit'], 's')
                        plot_dict['t_coord'] = np.array([ datetime.utcfromtimestamp(t).replace(tzinfo=utc) for t in seconds ])
                    else:
                        plot_dict['t_coord'] = np.array([ timeParse("".join(t), data_src_map['time_format']) for t in t_coord ])

                    valid_time = plot_dict['t_coord'][0]

                elif data_src_map['internal_time'] == 'relative':
#                   t_coord = t_coord - t_coord[0] + data_src_map['time_offset']
                    valid_time = timeParse(self._base_time, data_src_map['time_format']).replace(tzinfo=utc)
#                   print self._base_time
                    plot_dict['t_coord'] = [] 
                    for data_time in t_coord:
                        offset_seconds = Units.convert(data_time - t_coord[0], data_src_map['time_unit'], 's')
                        plot_dict['t_coord'].append(valid_time + timedelta(seconds=offset_seconds))

                    plot_dict['t_coord'] = np.array(plot_dict['t_coord'])

                self._valid_time = valid_time

                if not hasattr(self, '_master_plot'):
                    self._master_plot = plot_dict

                # Load element configuration files (how the lines, vectors, etc. are drawn)
                try:
                    plot_dict['element_config'] = Config(self._config_base + plot_dict['element_config'])
                except KeyError:
                    fatalError("element_config attribute not found for plot %s (%d)" % (plot, idx))

                try:
                    # If the units in the config file are in a dictionary, extract the units associated with the function.
                    plot_dict['element_config']['units'] = plot_dict['element_config']['units'][plot_dict['function']]
                except TypeError, IndexError:
                    pass

                if type(plot_dict['plot_times']) not in [ list, tuple ]: 
                    plot_dict['plot_times'] = [ plot_dict['plot_times'] ]

                if type(plot_dict['plot_times'][0]) == int:
                    plot_dict['plot_times'] = [ valid_time + timedelta(seconds=Units.convert(pt, self._plot_time_units, 's')) for pt in plot_dict['plot_times'] ]

                if data_src_map['grid_type'] is not None:
                    dindex = self._getDataIndexGrid(self._netcdf[data_src], plot_dict)
                else:
                    # grid_type is set to None, so assume it's an observation
                    dindex = self._getDataIndexObs(self._netcdf[data_src], plot_dict, data_src_map)

                # Load data to memory and evaluate the function
                self._loadFunctionString(self._netcdf[data_src], plot_dict, data_src_map, data_index=dindex)

        # Clear the data cache after everything has been loaded and derived variables have been computed
        self._clearCache()
        return

    #@FunctionCallCount
    def plot(self):
        """
        plot() [public]
        Purpose:    Handles the plotting of data.
        Parameters: [none]
        Returns:    [nothing]
        """

        self._initializeMap()

        # Loop to find which times for which we actually have data.
        for plot_time in self._master_plot['plot_times']:
            # Reset product for this image.
            self._resetProduct()

            n_plots_done = 0

            # Loop over all types of plots.
            for plot_list in MapPlotter._plot_order:
                # If we don't have the plot in the input file, don't try to plot it.
                if not hasattr(self, "_%s" % plot_list):
                    continue

                plots = getattr(self, "_%s" % plot_list)

                # Plot all grids for a given plot type.
                for idx, plot in enumerate(plots):
                    if plot_time in plot['t_coord']:
 #                      print plot['t_coord']
                        print "Plotting %s at time %s ..." % (plot['plot_name'], plot_time.strftime(self._product_time_format))
                        file_index = np.where(plot['t_coord'] == plot_time)[0][0]
                        n_plots_done += 1

                        # Plot contours ...
                        if plot_list == "grid_fn_contour":
                            self._plotContour(file_index, plot)

                        # Plot vectors ...
                        elif plot_list == "grid_fn_vector":
                           self._plotVector(file_index, plot)

                        # Plot fills ...
                        elif plot_list == "grid_fn_fill":
                            self._plotFill(file_index, plot)

                        # Plot conventional observations
                        elif plot_list == "obs_fn":
                            self._plotObservations(file_index, plot)
                        
                        # Plot Satellite image
                        elif plot_list == 'sat_fn':
                            self._plotSatellite(file_index, plot)


            if n_plots_done > 0:
                if hasattr(self, "_obs_fn") or hasattr(self, "_sat_fn") or hasattr(self, "_radar_fn"):
                    # If we're plotting any observed data, assume we want it to be an observed product
                    self._finalizeMap(plot_time, False)
                else:
                    self._finalizeMap(plot_time, True)
        return

    def _sanitize(self):
#       self._sanitizeCheck('map_n_lat', [int, float], ( lambda l: l > -90 and l < 90, "Parameter 'map_n_lat' must be between -90 and 90 (was given as '%(map_n_lat)s')." ))
#       self._sanitizeCheck('map_c_lat', [int, float], ( lambda l: l > -90 and l < 90, "Parameter 'map_c_lat' must be between -90 and 90 (was given as '%(map_c_lat)s')." ))
#       self._sanitizeCheck('map_s_lat', [int, float], ( lambda l: l > -90 and l < 90, "Parameter 'map_s_lat' must be between -90 and 90 (was given as '%(map_s_lat)s')." ))
#       self._sanitizeCheck('map_e_lon', [int, float], ( lambda l: l > -180 and l < 180, "Parameter 'map_e_lon' must be between -180 and 180 (was given as '%(map_e_lon)s')." ))
#       self._sanitizeCheck('map_c_lon', [int, float], ( lambda l: l > -180 and l < 180, "Parameter 'map_c_lon' must be between -180 and 180 (was given as '%(map_c_lon)s')." ))
#       self._sanitizeCheck('map_w_lon', [int, float], ( lambda l: l > -180 and l < 180, "Parameter 'map_w_lon' must be between -180 and 180 (was given as '%(map_w_lon)s')." ))
#       self._sanitizeCheck(['map_n_lat', 'map_s_lat'], [int, float], (lambda n, s: n > s, "Parameter 'map_n_lat' (given as '%(map_n_lat)s') must be greater than 'map_s_lat' (given as '%(map_s_lat)s')." ))

        self._sanitizeCheck('plot_times', [int, list, tuple], ( lambda t: type(t) in [list, tuple] and all([ type(x) == int for x in t]) or type(t) == int, "Parameter 'plot_times' must be either an integer or a list of integers." ))
        self._sanitizeCheck('vertical_level', str)

        if self._sanitizeCheck('grid_fn_contour', [list, tuple, dict], required=False):
            if type(self._grid_fn_contour) not in [ list, tuple ]:
                self._grid_fn_contour = [ self._grid_fn_contour ]

            for idx in range(len(self._grid_fn_contour)):
                self._sanitizeDict(self._grid_fn_contour[idx], 'plot_names', [str, dict], "the dictionary 'grid_fn_contour'", private=False, required=False)
                self._sanitizeDict(self._grid_fn_contour[idx], 'function', [str, list, tuple], "the dictionary 'grid_fn_contour'", private=False)
                self._sanitizeDict(self._grid_fn_contour[idx], 'element_config', str, "the dictionary 'grid_fn_contour'", private=False)

        if self._sanitizeCheck('grid_fn_fill', [list, tuple, dict], required=False):
            if type(self._grid_fn_fill) not in [ list, tuple ]:
                self._grid_fn_fill = [ self._grid_fn_fill ]

            for idx in range(len(self._grid_fn_fill)):
                self._sanitizeDict(self._grid_fn_fill[idx], 'plot_names', [str, dict], "the dictionary 'grid_fn_fill'", private=False, required=False)
                self._sanitizeDict(self._grid_fn_fill[idx], 'function', [str, list, tuple], "the dictionary 'grid_fn_fill'", private=False)
                self._sanitizeDict(self._grid_fn_fill[idx], 'element_config', str, "the dictionary 'grid_fn_fill'", private=False)

        if self._sanitizeCheck('grid_fn_vector', [list, tuple, dict], required=False):
            if type(self._grid_fn_vector) not in [ list, tuple ]:
                self._grid_fn_vector = [ self._grid_fn_vector ]

            for idx in range(len(self._grid_fn_vector)):
                self._sanitizeDict(self._grid_fn_vector[idx], 'plot_names', [str, dict], "the dictionary 'grid_fn_vector'", private=False, required=False)
                self._sanitizeDict(self._grid_fn_vector[idx], 'function', [str, list, tuple], "the dictionary 'grid_fn_vector'", private=False)
                self._sanitizeDict(self._grid_fn_vector[idx], 'element_config', str, "the dictionary 'grid_fn_vector'", private=False)

        if self._sanitizeCheck('obs_fn', [list, tuple, dict], required=False):
            if type(self._obs_fn) not in [ list, tuple ]:
                self._obs_fn = [ self._obs_fn ]

            for idx in range(len(self._obs_fn)):
                self._sanitizeDict(self._obs_fn[idx], 'plot_names', [str, dict], "the dictionary 'obs_fn'", private=False, required=False)
                self._sanitizeDict(self._obs_fn[idx], 'function', [str, list, tuple], "the dictionary 'obs_fn'", private=False)
                self._sanitizeDict(self._obs_fn[idx], 'element_config', str, "the dictionary 'obs_fn'", private=False)

        if self._sanitizeCheck('sat_fn', [list, tuple, dict], required=False):
            if type(self._sat_fn) not in [ list, tuple ]:
                self._sat_fn = [ self._sat_fn ]

            for idx in range(len(self._sat_fn)):
                self._sanitizeDict(self._sat_fn[idx], 'plot_names', [str, dict], "the dictionary 'sat_fn'", private=False, required=False)
                self._sanitizeDict(self._sat_fn[idx], 'function', [str, list, tuple], "the dictionary 'sat_fn'", private=False)
                self._sanitizeDict(self._sat_fn[idx], 'element_config', str, "the dictionary 'sat_fn'", private=False)

        if self._sanitizeCheck('radar_fn', [list, tuple, dict], required=False):
            if type(self._radar_fn) not in [ list, tuple ]:
                self._radar_fn = [ self._radar_fn ]

            for idx in range(len(self._radar_fn)):
                self._sanitizeDict(self._radar_fn[idx], 'plot_names', [str, dict], "the dictionary 'radar_fn'", private=False, required=False)
                self._sanitizeDict(self._radar_fn[idx], 'function', [str, list, tuple], "the dictionary 'radar_fn'", private=False)
                self._sanitizeDict(self._radar_fn[idx], 'element_config', str, "the dictionary 'radar_fn'", private=False)

        self._sanitizeCheck('image_rotation', int, (lambda r: r >= 0, "Parameter 'image_rotation' must be positve or zero (was given as '%(image_rotation)s')." ), required=False)

        return

    def _loadTemporalVerticalCoords(self, nc_file, plot_dict):
        """
        _loadTemporalVerticalCoords() [protected]
        Purpose:    Load the time and vertical coordinates for a given plot dictionary.
        Parameters: nc_file [ type=DataIO ]
                        A DataIO object for reading the data in.
                    plot_dict [ type=dict ]
                        A dictionary containing attributes for this plot.
        Returns:    A tuple contianing the t- and z-coordinates (in that order).  If the data does not have z-coordinates, return None in its place.
                        If the data does not have a t-coordinate, then find all the times in the file and uniquify them.
        """

        def allElementsEqual(L):
            if len(L) > 0:
                return L.count(L[0]) == len(L)
            else:
                return None

        def loadCoord(coord_names, name):
            # Verify the coordinates for the all the function dependencies are the same
            coords_equal = allElementsEqual(coord_names)
            if coords_equal is not None:
                if not coords_equal:
                    # They're not the same.  Uh-oh ...
                    fatalError("Loaded %s-coordinate names for function '%s' differ (are %s)." % (name, plot_dict['function'], ", ".join(coord_names)))
                else:
                    # They're the same load the coordinate ...
                    variable_list = nc_file.get_variable_list()
                    if coord_names[0] in variable_list:
                        # The dimension has a corresponding variable.  Just load it ...
                        coords = nc_file.get_variable(coord_names[0])
                    elif name in ['t', 'z']:
                        # The time dimension does not have a corresponding variable.  Look at the time coordinates we were given and find the one with the same dimensions as our dependencies
                        if name == 't':
                            coord_variables = data_src_map['time_var_name']
                        elif name == 'z':
                            coord_variables = [ ]
                            for key, val in data_src_map.iteritems():
                                if 'coord_name' in key:
                                    if type(val) in [ list, tuple ]:
                                        coord_variables.extend(val)
                                    else:
                                        coord_variables.append(val)

                        if type(coord_variables) not in [ list, tuple ]: coord_variables = [ coord_variables ]

                        file_coord_names = [ nc_file.get_variable_dimensions(c)[0] for c in coord_variables ]
                        coord_index = file_coord_names.index(coord_names[0])
                        coords = nc_file.get_variable(coord_variables[coord_index])
            else:
                coords = None
            return coords

        data_src_map = self._data_src_map[plot_dict['data_scheme']]

        if type(plot_dict['function']) not in [ list, tuple ]:
            plot_dict['function'] = [ plot_dict['function'] ]

        file_data_names = uniquify([ c for f in plot_dict['function'] for c in self._parseFunctionConstituents(f, parse_consts=False, data_scheme=data_src_map) ])
        file_data_dims = [ nc_file.get_variable_dimensions(v) for v in file_data_names ]

        z_coord_names = [ c[-3] for c in file_data_dims if len(c) == 4 ]
        z_coords = loadCoord(z_coord_names, 'z')

        if data_src_map['time_coordinate']:
            t_coord_names = [ c[0] for c in file_data_dims if len(c) > 2 ]
            t_coords = loadCoord(t_coord_names, 't')
            print t_coord_names,t_coords
        else:
            t_coords = np.unique(nc_file.get_variable(data_src_map['time_var_name']))

        return t_coords, z_coords

    def _getDataIndexGrid(self, nc_file, plot_dict):
        # Construct grid
        data_src_map = self._data_src_map[plot_dict['data_scheme']]
        nc_lats, nc_lons = self._constructGrid(nc_file, data_src_map)
        lats_shape = nc_lats.shape
        lons_shape = nc_lons.shape

        # Find bounds of the grid
        #ll_corner_x, ll_corner_y = self._basemap(self._map_projection['ll_corner_longitude'], self._map_projection['ll_corner_latitude'])
        #ur_corner_x, ur_corner_y = self._basemap(self._map_projection['ur_corner_longitude'], self._map_projection['ur_corner_latitude'])
        ll_corner_x, ll_corner_y = (self._basemap.llcrnrx, self._basemap.llcrnry)
        ur_corner_x, ur_corner_y = (self._basemap.urcrnrx, self._basemap.urcrnry)

        nc_x, nc_y = self._basemap(nc_lons, nc_lats)

        idxs = np.where((ll_corner_x < nc_x) & (nc_x <= ur_corner_x) & (ll_corner_y < nc_y) & (nc_y <= ur_corner_y))

        lat_idx_start = max(0,                idxs[1].min() - 1)
        lat_idx_end   = min(nc_lats.shape[1], idxs[1].max() + 1)
        lon_idx_start = max(0,                idxs[2].min() - 1)
        lon_idx_end   = min(nc_lons.shape[2], idxs[2].max() + 1)

        plot_dict['y_data_bound'] = (lat_idx_start, lat_idx_end)
        plot_dict['x_data_bound'] = (lon_idx_start, lon_idx_end)

        n_time_steps = plot_dict['t_coord'].shape[0]
        t_slice = slice(n_time_steps)
        x_slice = slice(*plot_dict['x_data_bound'])
        y_slice = slice(*plot_dict['y_data_bound'])

        plot_dict['grid_lats'] = nc_lats[0, y_slice, x_slice]
        plot_dict['grid_lons'] = nc_lons[0, y_slice, x_slice]
        plot_dict['grid_xs'] = nc_x[0, y_slice, x_slice]
        plot_dict['grid_ys'] = nc_y[0, y_slice, x_slice]

        if plot_dict['z_coord'] is None:
            dindex = (t_slice, y_slice, x_slice)
        else:
            # Load vertical coordinate levels
            level, coord = plot_dict['vertical_level'].strip().split(" ", 1)
            if level.find("-") > -1:
                levels = level.split("-")
            else:
                levels = (level, level)

            # Convert the level (or level bounds) to the same units as the data in the file.
            if Units.quantity(coord) == "pressure":
                levels = tuple([ Units.convert(int(l), coord, data_src_map['p_default_units']) for l in levels ])
            elif Units.quantity(coord) == "length":
                levels = tuple([ Units.convert(int(l), coord, data_src_map['gp_default_units']) for l in levels ])

            z_data_bound = tuple([ np.where(plot_dict['z_coord'] == l)[0][0] for l in levels ])
            z_slice = slice(z_data_bound[0], z_data_bound[1] + 1)

            if np.all((plot_dict['z_coord'] < levels[0]) | (plot_dict['z_coord'] > levels[1])):
                fatalError("Vertical level '%s' out of bounds." % plot_dict['vertical_level'])

            z_coord_levels = plot_dict['z_coord'][z_slice]
            z_shape = list(z_coord_levels.shape)
            for dim in range(len(z_shape), 3):
                z_shape.append(1)
            z_coord_levels.reshape(z_shape)

            if Units.quantity(coord) == "pressure" and not self._inCache('P'):
                self._updateCache('P', z_coord_levels, coord)
            elif Units.quantity(coord) == "length" and not self._inCache('GP'):
                self._updateCache('GP', z_coord_levels, coord)

            # Generate data indexes
            dindex = (t_slice, z_slice, y_slice, x_slice)
        return dindex

    def _getDataIndexObs(self, nc_file, plot_dict, data_src_map):
        if plot_dict['vertical_level'] in ['surface', 'sfc']:
            obs_meta_filename = "data/sfc_stations.csv"
            upper_air = False
        else:
            obs_meta_filename = "data/ua_stations.csv"
            upper_air = True                        

        stn_reader = StationReader(obs_meta_filename)

        # Get the information for all the stations within the bounds of the map
        # Get time information from the file and find all the observations at the time we want to plot
        data_times = nc_file.get_variable(data_src_map['time_var_name'])
        obs_at_plot_time = (data_times == calendar.timegm(plot_dict['plot_times'][0].timetuple()))
        plot_obs_index = np.where(obs_at_plot_time)[0]

        # Search for stations by latitude and longitude bounds of the domain
        map_stations = stn_reader.searchByDomain(self._basemap)

        # Retrieve the stations at the time we want to plot from the file
        data_stations = nc_file.get_variable(data_src_map['stn_var_name'], data_index=plot_obs_index)
        data_stations = [ "".join([ c for c in data_stations[idx] if type(c) != np.ma.core.MaskedConstant ]) for idx in range(data_stations.shape[0]) ]

        missing = []

        stn_indexes = np.zeros((0,), dtype=np.int32)
        stn_lats = np.zeros((0,))
        stn_lons = np.zeros((0,))

        if upper_air:
            pres_indexes = np.zeros((0,), dtype=np.int32)

        for stn in map_stations:
            # Search for every station in the domain and do "filtering" (make sure we're not plotting stations close together).
            stn_id = stn['ID']
            if data_src_map['stn_id_type'] == "IATA": stn_id = stn_id[1:]

            # Create a boolean array that says whether each station has the same ID as this station
            try:
                stn_idx = data_stations.index(stn_id)
            except ValueError:
                stn_idx = -1

            if len(stn_lats) > 0:
                # If we already have stations, compute the minimum distance between the new stations and all the others.
                min_dist = np.min(np.sqrt((stn['Latitude'] - stn_lats) ** 2 + (stn['Longitude'] - stn_lons) ** 2))
            else:
                # If not, set it to an arbitrarily large value
                min_dist = 999999

            # Find the index station in the data file and keep it only if it's not "too close" to what we've already got.
            if stn_idx >= 0 and min_dist > 0.40:

                if upper_air:
                    level, coord = plot_dict['vertical_level'].strip().split(" ", 1)

                    pressure_levels = nc_file.get_variable(data_src_map['p_var_name'])[plot_obs_index[stn_idx]]
#                   print np.where(pressure_levels == Units.convert(int(level), coord, data_src_map['p_default_units']))
                    pres_idx = np.where(pressure_levels == Units.convert(int(level), coord, data_src_map['p_default_units']))[0]

                    if len(pres_idx) > 0:
                        stn_indexes = np.append(stn_indexes, stn_idx)
                        pres_indexes = np.append(pres_indexes, pres_idx[0])

                        stn_lats = np.append(stn_lats, stn['Latitude'])
                        stn_lons = np.append(stn_lons, stn['Longitude'])
                    else:
                        missing.append(stn['ID'])
                else:
                    stn_indexes = np.append(stn_indexes, stn_idx)

                    stn_lats = np.append(stn_lats, stn['Latitude'])
                    stn_lons = np.append(stn_lons, stn['Longitude'])
            else:
                missing.append(stn['ID'])

        # Create the data index array from the obs_at_plot_time array (take the indexes from that array that we found in the domain)
        if upper_air:
            dindex = (plot_obs_index[stn_indexes], pres_indexes)
        else:
            dindex = (plot_obs_index[stn_indexes], )

        # Keep around the list of stations IDs and their latittudes and longitudes
        plot_dict['stns'] = [ s['ID'] for s in map_stations if s['ID'] not in missing ]
        plot_dict['stn_lats'] = stn_lats
        plot_dict['stn_lons'] = stn_lons

        plot_dict['stn_xs'], plot_dict['stn_ys'] = self._basemap(stn_lons, stn_lats)

        return dindex

    #@FunctionCallCount
    def _constructGrid(self, nc_file, data_scheme):
        """
        _constructGrid() [protected]
        Purpose:    Extract grid information from the data file and return an array of latitudes and longitudes describing the grid.
        Parameters: nc_file [type=DataIO]
                        The DataIO instance used to load in the file
                    data_scheme [type=dictionary]
                        A dictionary describing the mappings from internal variable names to those in the file.
        Returns:    A 2D array of latitudes and a 2D array of longitudes.
        """

        def findVariable(variable_name):
            # Find a variable in the file, taking into account the "var:attr" syntax for a variable attribute.
            if ':' in data_scheme[variable_name]:
                var_name, attr_name = data_scheme[variable_name].split(':')
                return nc_file.get_variable_attribute(var_name, attr_name)[0]
            else:
                return nc_file.get_variable(data_scheme[variable_name])[0]

        if data_scheme['grid_type'] == 'explicit':
            # The grid is explicitly in the file, so extract it directly
            lats = nc_file.get_variable(data_scheme['lt_var_name'])
            lons = nc_file.get_variable(data_scheme['ln_var_name'])
        elif data_scheme['grid_type'] == 'explicit_1d':
            lats_1d = nc_file.get_variable(data_scheme['lt_var_name'])
            lons_1d = nc_file.get_variable(data_scheme['ln_var_name'])

            lons, lats = np.meshgrid(lons_1d, lats_1d)

            self._data_map = Basemap(projection='cyl', resolution=None,
                llcrnrlat=lats.min(), llcrnrlon=lons.min(), urcrnrlat=lats.max(), urcrnrlon=lons.max()
            )
        elif data_scheme['grid_type'] == 'lcc':
            # The parameters to describe a Lambert Conic Conformal grid are in the file, so extract them and reconstruct the grid.
            proj_grid = { 'projection':'lcc', 'resolution':None }
            proj_conf= {}
            for param in [ 'std_lat1', 'std_lat2', 'std_lon', 'grid_dx', 'grid_dy', 'ngrid_x', 'ngrid_y' ]:
                proj_conf[param] = findVariable(param)

            width = proj_conf['ngrid_x'] * proj_conf['grid_dx']
            height = proj_conf['ngrid_y'] * proj_conf['grid_dy']

            p = pyproj.Proj(proj='lcc', lon_0=proj_conf['std_lon'], lat_1=proj_conf['std_lat1'], lat_2=proj_conf['std_lat2'], x_0=width, y_0=height)

            if 'ctr_lat' in data_scheme and 'ctr_lon' in data_scheme:
                ctrx, ctry = p(findVariable('ctr_lon'), findVariable('ctr_lat'))
                llcrnrx = ctrx - width / 2.
                llcrnry = ctry - height / 2.
            elif 'll_lat' in data_scheme and 'll_lon' in data_scheme:
                llcrnrx, llcrnry = p(findVariable('ll_lon'), findVariable('ll_lat'))
            else:
                fatalError("Data scheme must contain both of either ll_lat and ll_lon or ctr_lat and ctr_lon for the Lambert data grid")

            x = llcrnrx + proj_conf['grid_dx'] * np.arange(proj_conf['ngrid_x'])
            y = llcrnry + proj_conf['grid_dy'] * np.arange(proj_conf['ngrid_y'])

            x, y = np.meshgrid(x, y)

            lons,lats = p(x, y, inverse=True)

            self._data_map = Basemap(projection='lcc', resolution=None, 
                llcrnrlat=lats[0, 0], llcrnrlon=lons[0, 0], urcrnrlat=lats[-1, -1], urcrnrlon=lons[-1, -1],
                lon_0=proj_conf['std_lon'], lat_1=proj_conf['std_lat1'], lat_2=proj_conf['std_lat2'])

#           grid_x, grid_y = np.meshgrid(np.arange(proj_conf['ngrid_x']), np.arange(proj_conf['ngrid_y']))
#           lons, lats = self._data_map(grid_x, grid_y, inverse=True)

        elif data_scheme['grid_type'] == "Stereographic":
            x_values = nc_file.get_variable('x')
            y_values = nc_file.get_variable('y')
            proj_center = nc_file.get_variable_attribute('Stereographic', 'latitude_of_projection_origin')[0]
            std_lon = nc_file.get_variable_attribute('Stereographic', 'longitude_of_projection_origin')[0]
            scale_fac = nc_file.get_variable_attribute('Stereographic', 'scale_factor_at_projection_origin')[0]

            width = x_values[-1] - x_values[0]
            height = y_values[0] - y_values[-1]

            map = Basemap(projection='stere', resolution=None,
                width=width, height=height, lat_0=proj_center, lon_0=std_lon, lat_ts=60, boundinglat=0)
            self._data_map = map

            map_ur_x, map_ur_y = map(map.urcrnrlon, map.urcrnrlat)

            delta_x = map_ur_x - width
            delta_y = map_ur_y - height

            grid_x, grid_y = np.meshgrid(x_values - x_values[0] + delta_x / 2., y_values - y_values[-1] + delta_y / 2.)
            lons, lats = map(grid_x, grid_y, inverse=True)

        dims = len(lats.shape)
        extra_dims = []
        for idx in xrange(dims, 3): 
            extra_dims.append(1)
        extra_dims.extend(lats.shape)
        return lats.reshape(extra_dims), lons.reshape(extra_dims)

    #@FunctionCallCount
    def _initializeMap(self):
        """
        _initializeMap() [protected]
        Purpose:    Set up the plot grids for the map and initialize the Basemap instance
        Parameters: [none]
        Returns:    [nothing]
        """
        self._initializeProduct()

        # Find our master data source and the plot scheme for our master data source
#       nc = self._netcdf[self._master_plot['data_src']]
#       data_src_map = self._data_src_map[self._master_plot['data_scheme']]
#       if data_src_map['grid_type'] is not None:
#           # We have a grid, so get the x and y data bounds for subsetting the grid and then subset the grid.
#           data_lats = self._master_plot['grid_lats']
#           data_lons = self._master_plot['grid_lons']
#       else:
#           # No grid, so just get the list of station latitudes and longitudes.
#           data_lats = self._master_plot['stn_lats']
#           data_lons = self._master_plot['stn_lons']
#
#       # This sets the standard grid point structure at full resolution
#       self._x_grid, self._y_grid = self._basemap(data_lons, data_lats)
#
#       if hasattr(self, "_grid_fn_vector"):
#           stride = min([ plot['stride'] for plot in self._grid_fn_vector ]) 
#           # This sets a thinned out grid point structure for plotting
#           # wind barbs at the interval specified in "stride"
#
#           self._x_wind_grid = self._x_grid[::stride, ::stride]
#           self._y_wind_grid = self._y_grid[::stride, ::stride]

        return

    #@FunctionCallCount
    def _finalizeMap(self, plot_time, is_forecast):
        """
        _finalizeMap() [protected]
        Purpose:    Add final things to the map, such as the state/country borders, color bar (if it's a fill plot), 
                        title, valid time, and image border, and then save the image.
        Parameters:    forecast_hour [type=int]
                        Forecast hour for model products (pass in None for an observed product).
        Returns:    [nothing]
        """

        # Draw the lines (should add parameters to the file to say which ones to draw and how).
        self._basemap.drawcoastlines()
        self._basemap.drawcountries()
        self._basemap.drawstates()

        plots = []
        plot_names = []
        # Assemble a list of plots
        for plot in MapPlotter._plot_order:
            try:
                plots.append(getattr(self, "_%s" % plot))
            except AttributeError:
                pass

        # Assemble a list of plot names (to pass to the finalize product function)
        for grid_fn in plots:
            for plot_dict in grid_fn:
                if type(plot_dict['plot_name']) != dict:
                    if type(plot_dict['element_config']['units']) != dict and plot_dict['element_config']['units'] is not None:
                        units = Units(plot_dict['element_config']['units'])
                        plot_names.append("%s (%s)" % (plot_dict['plot_name'], units.toLaTexString()))
                    else:
                        plot_names.append("%s" % plot_dict['plot_name'])
                else:
                    if type(plot_dict['element_config']['units']) != dict:
                        pn = plot_dict['plot_name'].values()

                        if len(pn) > 1:
                            pn[-1] = "and " + pn[-1]

                        if plot_dict['element_config']['units'] is not None:
                            units = Units(plot_dict['element_config']['units'])
                            plot_names.append("%s (%s)" % (", ".join(pn), units.toLaTexString()))
                        else:
                            plot_names.append("%s" % ", ".join(pn))
                    else:
                        for key in plot_dict['plot_name'].keys():
                            if plot_dict['element_config']['units'][key] is not None:
                                units = Units(plot_dict['element_config']['units'][key])
                                plot_names.append("%s (%s)" % (plot_dict['plot_name'][key], units.toLaTexString()))
                            else:
                                plot_names.append("%s" % plot_dict['plot_name'][key])

        # Finalize the product.  Should be last.
        self._finalizeProduct(plot_time, is_forecast, plot_names=plot_names)
        return

    #@FunctionCallCount
    def _plotContour(self, nc_fh_index, plot_info):
        """
        _plotContour() [private]
        Purpose:    Plot contours on the map.
        Parameters:    nc_fh_index [type=int]
                        Temporal index into the NetCDF data array.
                    plot_info [type=dictionary]
                        Information about the plot (as specified in the input file).
        Returns:    [nothing]
        """

        plevels = []
        pcolors = []
        plstyles = []
        plwidths = []

        if len(plot_info['data'].shape) == 4:
            plot_data = plot_info['data'][nc_fh_index, 0, ...]
        elif len(plot_info['data'].shape) == 3:
            plot_data = plot_info['data'][nc_fh_index, ...]

        for line in plot_info['element_config']['lines']:
             
            if line['lower_bound'] is None:
                line_set_lb = ceil(plot_data.min() / line['contour_interval']) * line['contour_interval']
            else:
                line_set_lb = line['lower_bound']

            if line['upper_bound'] is None:
                line_set_ub = ceil(plot_data.max() / line['contour_interval']) * line['contour_interval']
            else:
                line_set_ub = line['upper_bound']

            for level in np.arange(line_set_lb, line_set_ub, line['contour_interval']):
                plevels.append(level)
                pcolors.append(line['line_color'])
                plstyles.append(line['line_style'])
                plwidths.append(line['line_weight'])

        if len(plevels) == 0:
            warning("No contours will be plotted with the current configuration.")
        else:
            CS = pylab.contour(plot_info['grid_xs'], plot_info['grid_ys'], plot_data, levels=plevels, colors=pcolors, linestyles=plstyles, linewidths=plwidths)
            pylab.clabel(CS, fontsize=8, inline=1, fmt='%d', inline_spacing=1)
        return

    #@FunctionCallCount
    def _plotFill(self, nc_fh_index, plot_info):
        """
        _plotFill() [private]
        Purpose:    Plot filled areas on the map.
        Parameters:    nc_fh_index [type=int]
                        Temporal index into the NetCDF data array.
                    plot_info [type=dictionary]
                        Information about the plot (as specified in the input file).
        Returns:    [nothing]
        """

        if len(plot_info['data'].shape) == 4:
            plot_data = plot_info['data'][nc_fh_index, 0, ...]
        elif len(plot_info['data'].shape) == 3:
            plot_data = plot_info['data'][nc_fh_index, ...]

        plevels = [ ]
        pcolors = [ ]

        extend_lower = False
        extend_upper = False
        for x in plot_info['element_config']['fill']:
            add_me = True
            if x['lower_bound'] is None:
                add_me = False 
                extend_lower = True

            if x['upper_bound'] is None:
                add_me = False 
                extend_upper = True

            if add_me:
                plevels.extend([x['lower_bound'], x['upper_bound']])

        plevels = uniquify(plevels)
        plevels.sort()

        for idx in xrange(len(plevels) - 1):
            for fill in plot_info['element_config']['fill']:
                if fill['lower_bound'] is not None and fill['upper_bound'] is not None and fill['lower_bound'] == plevels[idx] and fill['upper_bound'] == plevels[idx + 1]:
                    pcolors.append(fill['color'])
                    break

        if extend_upper:
            for fill in plot_info['element_config']['fill']:
                if fill['upper_bound'] is None:
                    pcolors.append(fill['color'])
                    break

        if len(pcolors) == 0:
            warning("No fills will be plotted with the current configuration.")
        else:
            extend = 'neither'
            if extend_lower and extend_upper: extend = 'both'
            elif extend_lower:                extend = 'min'
            elif extend_upper:                extend = 'max'

            CS = pylab.contourf(plot_info['grid_xs'], plot_info['grid_ys'], plot_data, levels=plevels, colors=pcolors, extend=extend)

            if extend_lower:
                for fill in plot_info['element_config']['fill']:
                    if fill['lower_bound'] is None:
                        CS.cmap.set_under(fill['color'])
                        break

            if extend_upper:
                for fill in plot_info['element_config']['fill']:
                    if fill['upper_bound'] is None:
                        CS.cmap.set_over(fill['color'])
                        break
            plevels_arr = np.array(plevels,dtype=float)
            mean_diff = np.mean(plevels_arr[1:] - plevels_arr[0:-1])
            if mean_diff < 1:
                fmt = '%0.2f'
            else:
                fmt = '%d'
            cbar = pylab.colorbar(orientation='horizontal', extend=extend, format=fmt, fraction=.06, aspect=65, shrink=.9, pad=0, ticks=plevels)
            pylab.setp(pylab.getp(cbar.ax.axes,'xticklabels'),fontsize=8)
            
                
        return

    #@FunctionCallCount
    def _plotVector(self, nc_fh_index, plot_info):
        """
        _plotVector() [private]
        Purpose:    Plot vectors on the map.
        Parameters: nc_fh_index [type=int]
                        Temporal index into the NetCDF data array.
                    plot_info [type=dictionary]
                        Information about the plot (as specified in the input file).
        Returns:    [nothing]
        """

        u, v = plot_info['data']
        stride = plot_info['stride']
        if len(u.shape) == 4:
            u = u[nc_fh_index, 0, ::stride, ::stride]
        elif len(u.shape) == 3:
            u = u[nc_fh_index, ::stride, ::stride]

        if len(v.shape) == 4:
            v = v[nc_fh_index, 0, ::stride, ::stride]
        elif len(v.shape) == 3:
            v = v[nc_fh_index, ::stride, ::stride]
        
        wind_lats = plot_info['grid_lats'][::stride, ::stride]
        wind_lons = plot_info['grid_lons'][::stride, ::stride]
        wind_xs = plot_info['grid_xs'][::stride, ::stride]
        wind_ys = plot_info['grid_ys'][::stride, ::stride]

        # This is a hack.  Fix me for real, which will involve finding a way to rotate vectors between arbitrary map projections.
        data_src_map = self._data_src_map[plot_info['data_scheme']]
#       if data_src_map['grid_type'] not in [ 'explicit' ]:
##          u, v = self._basemap.rotate_vector(u, v, plot_info['grid_lons'][::stride, ::stride], plot_info['grid_lats'][::stride, ::stride])
##          vec_start_x, vec_start_y = self._data_map(plot_info['grid_lons'][::stride, ::stride], plot_info['grid_lats'][::stride, ::stride])

#           if data_src_map['grid_type'] in [ 'explicit_1d' ]:
#                vec_end_x_data = wind_lons + u * 0.0001 # Position of the parcels 1 second later
#                vec_end_y_data = wind_lats + v * 0.0001
#           else:
#               print "Rotating wind vectors ..."
#               vec_end_x_data = self._x_wind_grid + u # Position of the parcels 1 second later
#               vec_end_y_data = self._y_wind_grid + v

#           vec_end_lon, vec_end_lat = self._data_map(vec_end_x_data, vec_end_y_data, inverse=True)

#           vec_end_x_plot, vec_end_y_plot = self._basemap(vec_end_lon, vec_end_lat)

#           orig_speed = np.hypot(u, v)
#           u_rot = vec_end_x_plot - self._x_wind_grid
#           v_rot = vec_end_y_plot - self._y_wind_grid

#           u = u_rot / np.hypot(u_rot, v_rot) * orig_speed
#           v = v_rot / np.hypot(u_rot, v_rot) * orig_speed
#       else:
#           def computeGridRot(grid1_x, grid1_y, grid2_x, grid2_y):
#               grid1_diff_x = grid1_x[1:, :] - grid1_x[:-1, :]
#               grid1_diff_y = grid1_y[1:, :] - grid1_y[:-1, :]

#               grid1_diff_vec = np.array((grid1_diff_x, grid1_diff_y))

#               grid2_diff_x = grid2_x[1:, :] - grid2_x[:-1, :]
#               grid2_diff_y = grid2_y[1:, :] - grid2_y[:-1, :]

#               grid2_diff_vec = np.array((grid2_diff_x, grid2_diff_y))

#               dotted_diffs = np.empty(grid2_diff_x.shape)
#               for idx in np.ndindex(dotted_diffs.shape):
#                   new_idx = [ slice(None) ]
#                   new_idx.extend(idx)
#                   dotted_diffs[idx] = np.dot(grid1_diff_vec[tuple(new_idx)], grid2_diff_vec[tuple(new_idx)])

#               rot_angle = np.arccos(dotted_diffs / (np.hypot(*grid1_diff_vec) * np.hypot(*grid2_diff_vec)))
#               return rot_angle

#           print wind_lons.shape

#           rot_left = computeGridRot(wind_lons[:-1, :], wind_lats[:-1, :], self._x_wind_grid[:-1, :], self._y_wind_grid[:-1, :])
#           rot_right = computeGridRot(wind_lons[1:, :], wind_lats[1:, :], self._x_wind_grid[1:, :], self._y_wind_grid[1:, :])

#           print rot_left.shape

#           rot_avg = np.empty(wind_lons.shape)
#           rot_avg[0] = rot_left[0]
#           rot_avg[-1] = rot_right[-1]
#           rot_avg[1:-1] = np.mean([rot_left, rot_right], axis=0)

#           new_angle = np.arctan2(v, u) + rot_avg

#           u_rot = np.hypot(u, v) * np.cos(new_angle)
#           v_rot = np.hypot(u, v) * np.sin(new_angle)

#           u = u_rot
#           v = v_rot

        if plot_info['element_config']['type'] == "barb":
            pylab.barbs(wind_xs, wind_ys, u, v, length = 5.5, linewidth = .5)
        elif plot_info['element_config']['type'] == "arrow":
            pylab.quiver(wind_xs, wind_ys, u, v, linewidth = 1)
        else:
            fatalError("Vector type %s not recognized" % plot_info['element_config']['type'])

        return

    #@FunctionCallCount
    def _plotObservations(self, nc_index, plot_info):
        """
        _plotObservations() [private]
        Purpose:    Plot observations on the map.
        Parameters:
        Returns:    [nothing]
        """

        # A few things to setup at the beginning.  Maybe make global parameters?
        offset = 1e4
        clip_box = Bbox([[0, 0], [1, 1]])

        for function, loc in plot_info['element_config']['location'].iteritems():
            # Loop over the location dictionary to find where to plot everything

            # Set plot points in the x-direction
            if loc in ['TopLeft', 'Left', 'BottomLeft']:
                plot_xs = plot_info['stn_xs'] - offset
                horiz_align = 'right'
            elif loc in ['TopRight', 'Right', 'BottomRight']:
                plot_xs = plot_info['stn_xs'] + offset
                horiz_align = 'left'
            else:
                plot_xs = plot_info['stn_xs'] 
                horiz_align = 'center'

            # Set plot points in the y-direction
            if loc in ['TopLeft', 'Top', 'TopRight']:
                plot_ys = plot_info['stn_ys'] + offset
                vert_align = 'bottom'
            elif loc in ['BottomLeft', 'Bottom', 'BottomRight']:
                plot_ys = plot_info['stn_ys'] - offset
                vert_align = 'top'
            else:
                plot_ys = plot_info['stn_ys']
                vert_align = 'center'

            if function in plot_info['function']:
                # The function to be plotted is one of the functions that we computed earlier
                if type(plot_info['data'][function]) in [ list, tuple ]:
                    # Plot vectors (which have been converted from magnitude and direction to u and v earlier).
                    u, v = plot_info['data'][function]

                    # Rotate the vector field to take into account the map projection.
                    rot_u, rot_v = self._basemap.rotate_vector(u, v, plot_info['stn_lons'], plot_info['stn_lats'])

                    # Do the plot
                    pylab.barbs(plot_xs, plot_ys, rot_u, rot_v, length=6.0, linewidth=.5)
                else:
                    if 'SC' in function:
                        for x, y, skc in zip(plot_xs, plot_ys, plot_info['data'][function]):
                            if skc < 0 : skc = 0
                            oktas = int(round(skc * 8))
                            drawSkyCover((x, y), 0.75 * offset, oktas=oktas)
                    elif 'WX' in function:
                        for x, y, wx in zip(plot_xs, plot_ys, plot_info['data'][function]):
                            if vert_align == "bottom": y += 0.75 * offset
                            elif vert_align == "top": y -= 0.75 * offset

                            if horiz_align == "left": x += 0.75 * offset
                            elif horiz_align == "right": x -= 0.75 * offset

                            parseAndPlot("".join(wx), (x, y), 0.75 * offset, color=plot_info['element_config']['colors'][function])
                    elif function == 'P':
                        kwargs = { 'size':'x-small', 'color':plot_info['element_config']['colors'][function], 'va':vert_align, 'ha':horiz_align, 
                                   'clip_box':clip_box, 'clip_on':True }

                        for x, y, pres in zip(plot_xs, plot_ys, plot_info['data'][function]):
                            if not np.isnan(pres):
                                if plot_info['element_config']['abbreviate_pressure']:
                                    pres = Units.convert(pres, plot_info['element_config']['units'][function], 'daPa')
                                    pres = str(int(round(pres)))[-3:]
                                else:
                                    pres = str(int(round(pres)))

                                pylab.text(x, y, pres, **kwargs)
                    else:
                        # After some more elif blocks, assume it's just text we want to put on there.
                        kwargs = { 'size':'x-small', 'color':plot_info['element_config']['colors'][function], 'va':vert_align, 'ha':horiz_align, 
                                   'clip_box':clip_box, 'clip_on':True }


                        for x, y, text in zip(plot_xs, plot_ys, plot_info['data'][function]):
                            if not np.isnan(text):
                                if type(text) in [ float, np.float32, np.float64 ]:
                                    text = int(round(text))
                                pylab.text(x, y, str(text), **kwargs)

            elif function == 'STN':
                # This means we want to plot the station ID's, so plot them.
                kwargs = { 'size':'x-small', 'color':plot_info['element_config']['colors'][function], 'va':vert_align, 'ha':horiz_align, 
                           'clip_box':clip_box, 'clip_on':True }

                for x, y, text in zip(plot_xs, plot_ys, plot_info['stns']):
                    # Loop over all station IDs and do the plot
                    pylab.text(x, y, text, **kwargs)

        return

    def _plotSatellite(self, file_idx, plot_info):
        print plot_info.keys()
#       img = plot_info['img']
#       xs = plot_info['']
#       img = 
#       pylab.
        return

if __name__ == "__main__":
    
    cfg = {
        'forecast_hours':[0, 1, 2, 3, 4, 5, 6, 9, 12],
        'product_title':"RUC 500mb Forecast",
        'image_file_name':"ruc_500mb_f%02d.png"
    }

    hpm = MapPlotter(cfg)
    hpm.loadData()
    hpm.plot()
