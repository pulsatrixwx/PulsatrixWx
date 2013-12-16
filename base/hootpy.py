from datetime import datetime, timedelta
from copy import copy
import re
import os

import pylab 

from util.config import Config
from util.units import Units
import util.derived as derived
from util.util import uniquify, any, fatalError, warning
#from util.callCount import FunctionCallCount

class HootPy(object):
    """
    HootPy
    Purpose:    Base class for HootPy.  Keeps track of product valid times, adding the finishing touches,
                    and saving the images.    
    Started:    14 June 2010 by Tim Supinie (tsupinie@ou.edu)
    Completed:  [not yet]
    Modified:   [not yet]
    """

    

#   @FunctionCallCount
    def __init__(self, config):
        """
        __init__()
        Purpose:    Constructor for the HootPy class.
        Parameters:    config [type=dictionary]
                        Dictionary containing the configuration parameters for this run.
        """

        self._valid_time = datetime.utcnow()

        for key, value in config.iteritems():
            setattr(self, "_%s" % key, value)

        try:
            meta_filename = config['meta_filename']
        except KeyError:
            meta_filename = "default/meta.hp"

        meta = Config(meta_filename)

        for key, value in meta.iteritems():
            setattr(self, "_%s" % key, value)

        self._var_name_cache = []

        self._sanitizeHootPy()
        self._sanitize()

        return

    def _sanitize(self):
        """
        sanitize() [protected, abstract]
        Purpose:    Abstract sanitize method.  Implement in a subclass to handle the sanitizing (quality control) of the configuration.
        Parameters: [none]
        Returns:    [nothing]
        """
        self._abstract()

#   @FunctionCallCount
    def loadData(self):
        """
        loadData() [public, abstract]
        Purpose:    Abstract loadData method.  Implement in a subclass to handle the loading in of data.
        Parameters: [none]
        Returns:    [nothing]
        """
        self._abstract()

#   @FunctionCallCount
    def plot(self):
        """
        plot() [public, abstract]
        Purpose:    Abstract plot method.  Implement in a subclass to handle the plotting of data.
        Parameters: [none]
        Returns:    [nothing]
        """
        self._abstract()

    def _sanitizeHootPy(self):
        """
        _sanitizeHootPy() [private]
        Purpose:    Sanitizes the plotting configuration for variables that are common to all subclasses.
        Parameters: [none]
        Returns:    [nothing]
        """

        self._sanitizeCheck('product_title', str)
        #self._sanitizeCheck('data_src', str, ( lambda n: os.path.exists(n), "Data file '%(data_src)s' does not exist." ))
        self._sanitizeCheck('data_scheme', str)
        self._sanitizeCheck('image_file_name', str, ( lambda n: os.path.exists(os.path.dirname(n)), "Path to file '%(image_file_name)s' does not exist." ))
        self._sanitizeCheck('image_size_x', int, ( lambda x: x > 0, "Parameter 'image_size_x' must be positive (was given as '%(image_size_x)d')." ))
        self._sanitizeCheck('image_size_y', int, ( lambda y: y > 0, "Parameter 'image_size_y' must be positive (was given as '%(image_size_y)d')." ))
        return

    def _sanitizeCheck(self, variables, dtype, constraint=None, required=True, private=True):
        return self._sanitizeDict(self.__dict__, variables, dtype, "the configuration file", constraint, required, private)

    def _sanitizeDict(self, dictionary, variables, dtype, src_name, constraint=None, required=True, private=True):
        if constraint is not None:
            func, error = constraint

        if type(variables) not in [ list, tuple ]:
            variables = [ variables ]

        if type(dtype) not in [ list, tuple ]:
            dtype = [ dtype ]

        values = []
        for var in variables:
            try:
                if private: key = "_%s" % var
                else:       key = var

                value = dictionary[key]
            except KeyError:
                if required:
                    fatalError("Parameter '%s' must be specified in %s." % (var, src_name))
                else:
                    return False

            if type(value) not in dtype:
                dtype_strings = [ str(t)[1:-1] for t in dtype ]
                if len(dtype) > 1:
                    dtype_strings[-1] = "or %s" % dtype_strings[-1]

                if len(dtype) > 2:
                    dtype_string = ", ".join(dtype_strings)
                else:
                    dtype_string = " ".join(dtype_strings)

                fatalError("Parameter '%s' in %s must have %s" % (var, src_name, dtype_string))

            values.append(value)

        if constraint is not None:
            if not func(*values):
                fatalError(error % dict(zip(variables, values)))
        return True

#   @FunctionCallCount
    def _initializeProduct(self):
        """
        _initializeProduct() [protected]
        Purpose:    Create the initial product and set it up to fill as much of the figure as possible.
        Parameters: [none]
        Returns:    [nothing]
        """

        dpi = 80 * 1.25
        pylab.figure(figsize=(float(self._image_size_x) / dpi, float(self._image_size_y) / dpi), dpi=dpi)
        pylab.axes((0, 0, 1, 1))
        return

    def _resetProduct(self):
        """
        _resetProduct() [protected]
        Purpose:    Reset the product so we don't have contours and fills, etc, bleeding over on time steps
        Parameters: [none]
        Returns:    [nothing]
        """
        dpi = 80 * 1.25

        pylab.clf()
        pylab.axes((0, 0, 1, 1))
        pylab.gcf().set_size_inches(float(self._image_size_x) / dpi, float(self._image_size_y) / dpi)
        return

#   @FunctionCallCount
    def _finalizeProduct(self, plot_time, is_forecast, plot_names=[]):
        """
        _finalizeProduct() [protected]
        Purpose:    Add final things to the product, such as the title, valid time, and border, and then save.
        Parameters: forecast_hour [type=int]
                        Forecast hour for model products (pass in None for an observed product).
        Returns:    [nothing]
        """

        plot_names = uniquify(plot_names)

        # Modify the last plot name for joining for the title string
        if len(plot_names) > 1:
            plot_names[-1] = "and " + plot_names[-1]

        # Create the forecast hour string according to whether or not we're passed a forecast hour.
        plot_time_delta = plot_time - self._valid_time
        hour = Units.convert(plot_time_delta.microseconds, 'us', 'hr') + Units.convert(plot_time_delta.seconds, 's', 'hr') + Units.convert(plot_time_delta.days, 'dy', 'hr')
        file_name = self._image_file_name % { 'plot_time':hour }

        if is_forecast:
            fh_string = " (F%03d)" % hour
        else:
            fh_string = ""

        if self._vertical_level in ['surface', 'sfc', "None"]:
            vert_level_str = ""
        else:
            vert_level_str = " %s" % self._vertical_level

        # Create the valid time string and assemble the title string
        valid_time_string = plot_time.strftime(self._product_time_format)
        title_string = "%s%s %s Valid: %s%s" % (self._product_title, vert_level_str, ", ".join(plot_names), valid_time_string, fh_string)

        # Put the title on the image
        pylab.title(title_string, weight="bold", size="x-small", bbox=dict(facecolor="#ffffff", alpha=0.7),x=0.5,y=0.95)

        # Save the figure
        try:
            pylab.savefig(file_name)
        except IOError:
            fatalError("Couldn't save image to %s" % file_name)

        print "Saved product '%s', valid at %s%s, to file %s" % (self._product_title, 
            valid_time_string, fh_string, file_name)

        pylab.close()
        return

#   @FunctionCallCount
    def _findAttribute(self, plot, attr_name):
        """
        _findAttribute() [protected]
        Purpose:    Find an attribute in the plot dictionary.  If the attribute isn't in the dictionary, look in the member variables of the class. Failing that,
                        raise an error.
        Parameters: plot [type=dictionary]
                        Dictionary of plot attributes to their values that is searched for an attribute.
                    attr_name [type=string]
                        The name of the attrbute to seach for.
        Returns:    The value of the attribute given by attr_name.
        """
        attribute = ""
        # Find the proper data source
        try:
            attribute = plot[attr_name]
        except KeyError:
            # The user didn't specify a data_src attribute in the plot dictionary, so look for a global attribute.
            try:
                attribute = getattr(self, "_%s" % attr_name)
                plot[attr_name] = attribute
            except:
                # The user didn't specify a global data_src attribute, either.  Uh-oh ...
                raise AttributeError("Attribute %s not found")
        return attribute
   
    @classmethod
    def _splitPlots(klass, plot_dicts):
        """
        """
        # Split plot_dicts with lists of functions into different entries in plot_dict
        indexes = [ idx for idx in xrange(len(plot_dicts)) if type(plot_dicts[idx]['function']) in [ list, tuple ] ]
        indexes.sort(reverse=True)

        for idx in indexes:
            plot = plot_dicts.pop(idx)
            for func in plot['function']:
                new_plot = copy(plot)
                new_plot['function'] = func
                plot_dicts.insert(idx, new_plot)
        return

    def _parseFunctionConstituents(self, function, parse_consts=True, data_scheme=None):
        """
        _parseFunctionConstituents() [protected]
        Purpose:    Take a function string and parse out the variables and constants that the function needs for computations.  If the data_scheme 
                        variable is given, convert the list to a file-based variable names.
        Parameters: function [ type=str ]
                        String containing the function.
                    parse_consts [ type=bool ] [ optional ]
                        Boolean value specifying whether or not to parse out the constants in the file name.  Default is True.
                    data_scheme [ type=dict ] [ optional ]
                        A dictionary containing the mapping of internal variables names to the variable names in the data file.
        Returns:    A list of internal constant names and/or a list of internal variable names or file variable names, depending on whether or not 
                        parse_consts was set and data_scheme was given.
        """
        # Put the variable list into a regexp-like format
        nc_variable_list = "(?:^|(?<=[\\W]))(?:" + "|".join(self._var_map.keys()) + ")(?:(?=[\\W])|$)"

        # Find all NetCDF variables in the function, removing duplicates
        var_list = uniquify(re.findall(nc_variable_list, function))

        if parse_consts:
            hp_const_list = "|".join([ const for const in dir(derived) if const[0] != "_" ])

            # Find all the constants in the function, removing duplicates
            const_list = uniquify(re.findall(hp_const_list, function))

        if data_scheme is not None:
            var_list = [ data_scheme[self._var_map[v]] for v in var_list ]

        if parse_consts:
            return var_list, const_list
        else:
            return var_list

#   @FunctionCallCount
    def _loadFunctionString(self, nc, plot_dict, data_scheme, data_index=None, scratch=False):
        """
        _loadFunctionString() [protected]
        Purpose:    Parses a function string from the input file, loads the data from the file, and converts the data to the proper units for plotting.
        Parameters: nc [type=DataIO]
                        DataIO object that loads in the data. 
                    plot_dict [type=dictionary]
                        Dictionary containing the plot attributes and their values.
                    data_scheme [type=dictionary]
                        Dictionary containing the mapping of internal variable names to the variable names in the data file.
                    data_index [type=np.array]
                        An array containing the indexes into the data array to return [not yet implemented].
        Returns:    [nothing]
        """

        if type(plot_dict['function']) not in [ list, tuple ]:
            plot_dict['function'] = [ plot_dict['function'] ]

        plot_dict['data'] = {}        

        for function in plot_dict['function']:
            parse_function = False
            try:
                default_units = data_scheme[self._unit_map[function]]
            except KeyError:
                parse_function = True

            if parse_function or default_units is not None:
                units_function = function

                var_list, const_list = self._parseFunctionConstituents(function)

                parsed_function = function
                for const in const_list:
                    # Replace each constant in the function and units function with the proper source code to get the value
                    parsed_function = re.sub("(?:^|(?<=[\\W]))(%s)(?:(?=[\\W])|$)" % const, "derived.\\1", parsed_function)

                for nc_variable in var_list:
                    if scratch:
                        plot_dict['scratch'] = nc.get_variable(data_scheme[self._var_map[nc_variable]], data_index)
                        return

                    if self._inCache(nc_variable):
                        # Check the cache to make sure the units in the cache are what we think they are (they might have been converted before putting them in).
                        file_units = data_scheme[self._unit_map[nc_variable]]
                        cache_units = self._getFromCache(nc_variable, 'units')
                        if file_units != cache_units:
                            nc_data = self._getFromCache(nc_variable)
                            self._updateCache(nc_variable, Units.convert(nc_data, cache_units, file_units), file_units)
                    else:
                        # Put data in the global namespace for easy access (will be deleted later)
                        self._updateCache(nc_variable, nc.get_variable(data_scheme[self._var_map[nc_variable]], data_index),
                            data_scheme[self._unit_map[nc_variable]])

                for const in const_list:
                    # Find each constant and HootPy function in the string
                    match = re.search("(?:^|(?<=[\\W]))%s\\(([\\w\\, ]+)\\)?(?:(?=[\\W])|$)" % const, units_function)
                    # Parse out the arguments to each function.  If it's not a function (and really a constant, such as g) give it an empty list for arguments.
                    if match is not None and match.group(1) is not None:
                        args = re.split("\,[\s]*", match.group(1))
                    else:
                        args = []

                    # Determine what the units of the data for the arguments are.  If the argument's variable name is not in the cache, 
                    #   then that probably means it's the units being output from another HootPy function that's already been subbed 
                    #   into the string.  Put None in its place.
                    arg_units = [ self._getFromCache(a, 'units') if self._inCache(a) else None for a in args ]

                    # Determine what the function is expecting
                    func_units = derived._units(const, *arg_units)

                    # A bit of idiot-proofing on the arguments
                    if len(arg_units) != len(func_units['args']): fatalError("Incorrect number of arguments for function %s." % const)

                    for idx in xrange(len(args)):
                        # Convert all argument data to the units the function is expecting (only do it if we actually have units there, and they don't need to be converted.
                        if arg_units[idx] is not None and arg_units[idx] != func_units['args'][idx]:
                            self._updateCache(args[idx], Units.convert(self._getFromCache(args[idx], 'value'), arg_units[idx], func_units['args'][idx]), 
                                func_units['args'][idx])

                    # Substitute the units output from this function back into the units string
                    units_function = re.sub("(?:^|(?<=[\\W]))((?:%s)(?:\\([\w\\, ]+\\))?)(?:(?=[\\W])|$)" % const, 
                        func_units['return'], units_function)

                for nc_variable in var_list:
                    # Sub individual variables' units into the units string
                    if units_function.find(nc_variable) > -1:
                        units_function = re.sub("(?:^|(?<=[\\W]))(%s)(?:(?=[\\W])|$)" % nc_variable, 
                            data_scheme[self._unit_map[nc_variable]], units_function)
    
                plot_dict['default_units'] = Units.evaluateFunction(units_function)
            else:
                parsed_function = function
                self._updateCache(parsed_function, nc.get_variable(data_scheme[self._var_map[parsed_function]], data_index), None)
                plot_dict['default_units'] = None
    
            # Load data
            if len(plot_dict['function']) == 1:
                if not scratch:
                    if type(plot_dict['plot_name']) == dict:
                        print "Loading/computing data for %s ..." % plot_dict['plot_name'][function]
                    else:
                        print "Loading/computing data for %s ..." % plot_dict['plot_name']

                exec "plot_dict['data'] = %s " % parsed_function in globals(), locals()

                # Do units conversion
                if plot_dict['element_config']['units'] is not None:
                    if type(plot_dict['data']) in [ list, tuple ]:
                        plot_dict['data'] = tuple([ Units.convert(d, plot_dict['default_units'], plot_dict['element_config']['units']) for d in plot_dict['data'] ])
                    else:
                        plot_dict['data'] = Units.convert(plot_dict['data'], plot_dict['default_units'], plot_dict['element_config']['units'])
            else:
                if not scratch:
                    if type(plot_dict['plot_name']) == dict:
                        print "Loading/computing data for %s (%s) ..." % (plot_dict['plot_name'][function], function)
                    else:
                        print "Loading/computing data for %s (%s) ..." % (plot_dict['plot_name'], function)

                exec "plot_dict['data']['%s'] = %s " % (function, parsed_function) in globals(), locals()

                # Do units conversion
                if plot_dict['element_config']['units'][function] is not None:
                    if type(plot_dict['data'][function]) in [ list, tuple ]:
                        plot_dict['data'][function] = tuple([ Units.convert(d, plot_dict['default_units'], plot_dict['element_config']['units'][function]) for d in plot_dict['data'][function] ])
                    else:
                        plot_dict['data'][function] = Units.convert(plot_dict['data'][function], plot_dict['default_units'], plot_dict['element_config']['units'][function])

        if scratch:
            for nc_variable in var_list:
                self._clearCache(nc_variable)

        return 

#   @FunctionCallCount
    def _updateCache(self, key, value, units):
        """
        _updateCache() [protected]
        Purpose:    Put data into a global cache.  This has a couple of advantages.  First, data that is loaded once can be kept around so it isn't 
                        loaded again.  Second, putting a variable into the cache means it is accessible by its name only.  So putting, say, temperature
                        data into the cache under 'T' means it can be accesed through the T variable.  This provides a convenient way to tell python
                        about the data when it's evaluating a function, and doesn't involve as many messy regex substitutions.
        Parameters: key [type=string]
                        The name of the variable (maybe 'T' for temperature).
                    value [type=int,float,np.array]
                        The data to store in the cache.
                    units [type=str]
                        The units of the data (maybe K for temperature).
        Returns:    [nothing]
        """
        globals()[key] = value
        globals()["%s_units" % key] = units
        if key not in self._var_name_cache:
            self._var_name_cache.append(key)
        return

#   @FunctionCallCount
    def _getFromCache(self, key, item=None):
        """
        _getFromCache() [protected]
        Purpose:    Retrieve an item from the global cache.
        Parameters: key [type=str]
                        The name of the variable to retrieve from the cache.  Equivalent to just calling whatever the name of the variable is in 
                            plaintext.
                    item [type=str,None]
                        String specifying what to retrieve.  "value" or None will return the data, and "units" will return the units of the data.
        Returns:    Whatever is specified by the key and item arguments.
        """
        if item is None or item == "value":
            # Return data
            return globals()[key]
        elif item == "units":
            # Returns the units of the data
            return globals()["%s_units" % key]

#   @FunctionCallCount
    def _inCache(self, key):
        """
        _inCache() [protected]
        Purpose:    Return a boolean specifying whether or not a variable is in the global cache.
        Parameters: key [type=str]
                        The name of the variable we want to check for.
        Returns:    A boolean specifying whether or not the variable is in the cache.
        """
        return key in self._var_name_cache

#   @FunctionCallCount
    def _clearCache(self,variable=None):
        """
        _clearCache() [protected]
        Purpose:    Clear the global cache of raw data to free up the memory.
        Parameters: [none]
        Returns:    [nothing]
        """
        if variable is None:
            for nc_variable in self._var_name_cache:
                # Delete data in the global namespace
                del globals()[nc_variable] 
                del globals()["%s_units" % nc_variable]
            self._var_name_cache = [ ]
        else:
            try:
                del globals()[variable]
                del globals()['%s_units' % variable]
                self._var_name_cache.remove(variable)
            except KeyError:
                warning("Variable %s does not exist in the cache" % variable)
        return

    def _abstract(self):
        """
        _abstract() [protected]
        Purpose:    Emulate abstraction behavior of C++/Java.  Raises an exception in abstract methods.
        Parameters: [none]
        Returns:    [nothing]
        """
        raise NotImplementedError('Abstract method must be implemented in subclass.')

if __name__ == "__main__":
    cfg = { 
        'forecast_hours':[0, 3, 6, 9, 12],
    }
    hp = HootPy(cfg)
#   hp.loadData()
#   hp.plot()    
