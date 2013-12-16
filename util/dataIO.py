import numpy as np
from util import warning, fatalError

#from callCount import FunctionCallCount

try:
    # Look for PyNIO ...
    import Nio as reader
    print "Using PyNIO as data reader."
except ImportError:
    # PyNIO is not installed, so ...
    try:
        from netCDF4 import Dataset as reader
        warning("Using netCDF4-python as data reader.  NetCDF format supported only.")
    except ImportError:
        try:
            # Look for scipy's NetCDF reader ...
            import scipy.io.netcdf as reader
            warning("Using scipy.io as data reader.  NetCDF format supported only.")
        except ImportError:
            # scipy's NetCDF reader isn't installed either.  Uhoh ...
            fatalError("No scientific data reader found.  Exiting ...")

class DataIO:
    """
    DataIO
    Purpose:    Handles the reading and writing of data for HootPy using whichever reader is installed
    Started:    30 September 2010 by Tim Supinie (tsupinie#@ou.edu)
    Completed:    [not yet]
    Modified:    [not yet]
    """

    # Probably about time to write down a running list of quirks I've noticed with scipy's netcdf reader:
    # 1.    Can't have an attribute be an empty string (nc.variables['stuff'].stuff = "" will fail on file close)
    # 2.    Doesn't have a "read and write" mode, so if you want to do any sort of reading, you must open the file in read mode, close it, then open it in write mode.

    def __init__(self, file_name=None, mode=None):
        """
        __init__()
        Purpose:    Constructor for the DataIO class.
        Parameters:    file_name [type=string]
                        Name of the file to open.  If None, don't open a file.
                    mode [type=string]
                        How to open the file (e.g. 'r' for reading, 'rw' for reading and writing, etc.)
                        Passing None defaults to opening for reading only.
        """
        self._df = None
        self._rm_variables = []
        self._sp_variables = {}
        self._ov_variables = {}
        self._rm_attributes = { '__file__':[] }
        self._attempted_close = False

        if file_name is not None:
            self.open_file(file_name, mode)
        return

    #@FunctionCallCount
    def open_file(self, file_name, mode=None):
        """
        open_file() [public]
        Purpose:    Opens a file with the given mode.
        Parameters:    file_name [type=string]
                        Name of the file to open.
                    mode [type=string]
                        How to open the file (e.g. 'r' for reading, 'rw' for reading and writing, etc.)
                        Passing None defaults to opening for reading only.
        Returns:    [nothing]
        """
        # Set the mode default
        if mode is None: mode = 'r'

        # If the file is already open, close it.
        if self._df is not None:
            self.close()

        self._file_name = file_name

        # Figure out whether we're reading, writing, or both.
        self._read = (mode.find('r') > -1)
        self._write = (mode.find('w') > -1 or mode == 'r+')

#       if self._write:
#           warning("DataIO: Writing is currently not supported.  Opening as read-only.")
#           self._write = False

        if reader.__name__ == "Nio":
            # If the reader is PyNIO open it with PyNIO's open command
            self._df = reader.open_file(file_name, mode=mode)
        elif reader.__name__ == "scipy.io.netcdf":
            # If the reader is scipy, open it with scipy's open command
            if self._read:
                self._df = reader.netcdf_file(file_name, 'r')
            elif self._write:
                self._df = reader.netcdf_file(file_name, 'w')
        elif reader.__name__ == "Dataset":
            # If the reader is netCDF4, open it with that open command
            if self._read:
                self._df = reader(file_name, mode='r')
            elif self._write:
                self._df = reader(file_name, mode='a')
        return

    def create_dimension(self,dimName,value):
        """
        create_dimension() [public]
        Purpose:  Creates a dimension in the file.
        Parameters:     dimName [type=string]
                            Name of the dimension being created.
                        value [type=int]
                            Value of the dimension.
        Returns:  None
        """
        if dimName not in self._df.dimensions.keys():
            if reader.__name__=="Nio":
                self._df.create_dimension(dimName,value)
            else:
                self._df.createDimension(dimName,value)
        else:
            print 'Error:  Dimension already created!'
    def get_dimensions(self):
        """
        get_dimensions() [public]
        Purpose:    Get the dictionary showing all the dimensions and their values
        Parameters: None
        Returns:    Dictionary containing dimension information
        """
        return self._df.dimensions
    def create_variable(self,variable,type,dimensions=None):
        """
        create_variable() [public]
        Purpose:    Creates a variable in the file.
        Parameters: variable [type=string]
                        Name of NetCDF variable to be created.
                    type [type=string]
                        The variable primitive datatypes correspond to the dtype attribute of a numpy array. You can specify the 
                        datatype as a numpy dtype object, or anything that can be converted to a numpy dtype object. Valid 
                        datatype specifiers include: 'f4' (32-bit floating point), 'f8'  (64-bit floating point), 'i4' (32-bit 
                        signed integer), 'i2' (16-bit signed integer), 'i8' (64-bit singed integer), 'i1' (8-bit signed integer), 
                        'u1' (8-bit unsigned integer), 'u2' (16-bit unsigned integer), 'u4' (32-bit unsigned integer), 'u8' 
                        (64-bit unsigned integer), or 'S1' (single-character string). The old Numeric single-character typecodes 
                        ('f','d','h', 's','b','B','c','i','l'), corresponding to ('f4','f8','i2','i2','i1','i1','S1','i4','i4'),
                        will also work.
                    dimensions [type=tuple of strings]
                        The dimensions of the variable as a tuple of strings that name dimensions already created in the netCDF file.
                        If the variable is a scalar, only specify the variable and type.
        Returns:  An instance of the Variable class that points to the variable in the file.
        """
        if variable not in self._df.variables.keys():
            if reader.__name__ == "Nio":
                print 'Creating variable for Nio'
                if dimensions == None:
                    self._df.create_variable(variable,type,dimensions=())
                else:
                    self._df.create_variable(variable,type,dimensions)
            elif reader.__name__ in ['scipy.io.netcdf','Dataset']:
                print 'Creating variable for netcdf'
                if dimensions == None:
                    self._df.createVariable(variable,type)
                else:
                    self._df.createVariable(variable,type,dimensions)
        else:
            warning("DataIO: Variable '%s' already exists in file '%s'." % (variable, self._file_name))

    #@FunctionCallCount
    def get_variable(self, variable, data_index=None):
        """
        get_variable() [public]
        Purpose:    Returns a masked data array for a variable in the file.
        Parameters: variable [type=string,tuple,list]
                        NetCDF variable to return.  Could be multiple variables in a list.
                    data_index [type=tuple,np.array]
                        Either a tuple of slice objects (created by slice()) or a numpy array of indexes that specify which data
                        to return
        Returns:    A numpy masked array containing the variable's data.
        """
        if type(variable) in [list, tuple]:
            var_data = []
            for var in variable:
                try:
                    if data_index is None:
                        var_data.append(self._df.variables[var][:])
                    else:
                        var_data.append(self._df.variables[var][data_index])
                except KeyError:
                    fatalError("DataIO: File '%s' does not contain a variable '%s'." % (self._file_name, var))
        else:
            try:
            # Assign the data to a variable to play with later ...
                if data_index is None:
                    var_data = self._df.variables[variable][:]
                else:
                    if reader.__name__ == "Nio":
                        var_data = self._df.variables[variable].get_value()[data_index]
                    else:
                        var_data = self._df.variables[variable][data_index]
            except KeyError:
                fatalError("DataIO: File '%s' does not contain a variable '%s'." % (self._file_name, variable))

        if reader.__name__ == "Nio":
            # If we're using PyNIO, it already returns a masked array, so just 
            #    return the data we already have.
            return var_data
        elif reader.__name__ in ["scipy.io.netcdf","Dataset"]:
            # If we're using scipy, we have to mask the data ourselves ...
            if type(variable) in [list, tuple]:
                var_data_masked = []
                for idx in xrange(len(var_data)):
                    try:
                        missing_value = self.get_variable_attribute(variable[idx], 'missing_value', force_exception=True)
                        mask = (var_data[idx] == missing_value) 
                    except ValueError:
                        mask = np.zeros(var_data[idx].shape)
                    var_data_masked.append(np.ma.array(var_data[idx], mask=mask))
            else:
                try:
                    missing_value = self.get_variable_attribute(variable, 'missing_value', force_exception=True)
                    mask = (var_data == missing_value) 
                except ValueError:
                    mask = np.zeros(var_data.shape)

                var_data_masked = np.ma.array(var_data, mask=mask)
            return var_data_masked
            

    #@FunctionCallCount
    def get_file_attribute(self, attribute):
        """
        get_file_attribute() [public]
        Purpose:    Returns an attribute of the NetCDF file.
        Parameters:    attribute [type=string]
                        The attribute to return
        Returns:    A string containing the attribute in the file.
        """
        try:
            if reader.__name__ == "Dataset":
                return getattr(self._df, attribute)
            else:    
                return self._df.attributes[attribute]
        except KeyError:
            fatalError("DataIO: File '%s' does not contain a variable '%s'." % (self._file_name, attribute))

    #@FunctionCallCount
    def get_variable_attribute(self, variable, attribute, force_exception=False):
        """
        get_variable_attribute() [public]
        Purpose:    Returns an attribute of a variable in the NetCDF file.
        Parameters:    variable [type=string]
                        The NetCDF variable to find an attribute from.
                    attribute [type=string]
                        The attribute to return.
                    force_exception [type=boolean]
                        Force the function to raise an exception on error, instead of just exiting.
        Returns:    A string containing the attribute of the NetCDF variable.
        """
        try:
            return getattr(self._df.variables[variable], attribute)
        except AttributeError:
            if force_exception:
                raise ValueError("DataIO: Variable '%s' in file '%s' does not have an attribute '%s'." % (variable, self._file_name, attribute))
            else:
                fatalError("DataIO: Variable '%s' in file '%s' does not have an attribute '%s'." % (variable, self._file_name, attribute))
        except KeyError:
            if force_exception:
                raise ValueError("DataIO: Variable '%s' in file '%s' does not have an attribute '%s'." % (variable, self._file_name, attribute))
            else:
                fatalError("DataIO: File '%s' does not contain a variable '%s'." % (self._file_name, variable))

    def get_variable_dimensions(self, variable):
        """
        get_variable_dimensions() [public]
        Purpose:    Returns the dimensions tuple associated with a variable.
        Parameters: variable
                        Variable name.
        Returns:    A tuple containing the names of the dimensions used by the variable
        """
        return self._df.variables[variable].dimensions

    def get_variable_attribute_list(self, variable):
        """
        get_variable_dimensions() [public]
        Purpose:    Returns a list of all the attributes for the given variable.
        Parameters: variable
                        Variable name.
        Returns:    A list of all the public attributes for the variable
        """
        if reader.__name__ == "Dataset":
            return self._df.variables[variable].ncattrs()
        else:
            return [ a for a in dir(self._df.variables[variable]) if a[0] != '_' ]

    #@FunctionCallCount
    def get_variable_list(self):
        """
        get_variable_list() [public]
        Purpose:    Returns a list of the variable names in the NetCDF file.
        Parameters: [none]
        Returns:    A list of the variables names.
        """
        return self._df.variables.keys()

    #@FunctionCallCount
    def get_file_attribute_list(self):
        """
        get_file_attribute_list() [public]
        Purpose:    Returns a list of attribute values of the NetCDF file.
        Parameters: [none]
        Returns:    A list of the values of the attributes of the NetCDF file.
        """
        if reader.__name__ == "Dataset":
            return self._df.ncattrs()
        else:
            return [ a for a in dir(self._df) if a[0] != '_' and a not in ['variables', 'dimensions'] ]

    def get_dimension(self,dimName):
        """
        get_dimension() [public]
        Purpose:  Returns the value of a specified dimension
        Parameters:  dimName [string]
                        Name of the dimension
        Returns:  The value of the dimension
        """
        if reader.__name__ == "Nio":
            return self._df.dimensions[dimName]
        else:
            return len(self._df.dimensions[dimName])


    def set_variable(self, variable, value):
        """
        set_variable() [public]
        Purpose:    Put an array of data in the file, but only if write flag has been set.
        Parameters: variable [type=string]
                        Name of the variable to set.
                    value [type=np.array]
                        Array to put in the file.
        Returns:    [nothing]
        """
        if not self._write:
            warning("DataIO: Write flag has not been set on file '%s'.  No action taken." % self._file_name)
            return
        
        self._ov_variables[variable] = value
        return

    def set_file_attribute(self, attribute, value):
        """
        set_file_attribute() [public]
        Purpose:    Set a file attribute, but only if the write flag is set.
        Parameters: attribute [type=string]
                        The name of the attribute to set.
                    value [type=int,float,char,string]
                        The value to put in the file.
        Returns:    [nothing]
        """
        if not self._write:
            warning("DataIO: Write flag has not been set on file '%s'.  No action taken." % self._file_name)
            return
        setattr(self._df, attribute, value)
        return

    def set_variable_attribute(self, variable, attribute, value):
        """
        set_variable_attribute() [public]
        Purpose:    Set a variable attriute, but only if the write flag is set.
        Parameters: variable [type=string]
                        Name of the variable whose attribute to set.
                    attribute [type=string]
                        Name of the attribute to set.
                    value [type=int,float,char,string]
                        The value to put in the file.
        Returns:    [nothing]
        """
        if not self._write:
            warning("DataIO: Write flag has not been set on file '%s'.  No action taken." % self._file_name)
            return
        if reader.__name__ in ["Nio","scipy.io.netcdf"]:
            setattr(self._df.variables[variable], attribute, value)
        elif reader.__name__ == "Dataset":
            setattr(self._df.variables[variable],attribute,value)
        return

    def splice_variable(self, variable, new_dim_name, data_index):
        """
        splice_variable() [public]
        """

        unique_indexes = []
        for index_list in data_index:
            unique_indexes.append(np.unique(index_list))
        self._sp_variables[variable] = (new_dim_name, tuple(unique_indexes))
        return

    def remove_variable(self, variable):
        """
        remove_variable() [public]
        Purpose:    Set a file variable to the remove list, to be "removed" on close.
        Parameters: variable [type=string]
                        Name of the variable to remove.
        Returns:    [nothing]
        """
        if not self._write:
            warning("DataIO: Write flag has not been set on file '%s'.  No action taken." % self._file_name)
            return

        if variable not in self._rm_variables:
            self._rm_variables.append(variable)
        else:
            warning("DataIO: Variable '%s' has already been removed from file '%s'." % (variable, self._file_name))
        return

    def remove_file_attribute(self, attribute):
        """
        remove_file_attribute() [public]
        Purpose:    Set a file attribute to the remove list, to be "removed" on close.
        Parameters: attribute [type=string]
                        Name of the attribute to remove.
        Returns:    [nothing]
        """
        if not self._write:
            warning("DataIO: Write flag has not been set on file '%s'.  No action taken." % self._file_name)
            return

        if attribute not in self._rm_attributes['__file__']:
            self._rm_attributes['__file__'].append(attribute)
        else:
            warning("DataIO: Attribute '%s' has already been removed from file '%s'." % (attribute, self._file_name))
        return

    def remove_variable_attribute(self, variable, attribute):
        """
        remove_variable_attribute () [public]
        Purpose:    Set a variable attribute to the remove list, to be "removed" on close.
        Parameters: variable [type=string]
                        Name of the variable whose attribute to remove.
                    attribute [type=string]
                        Name of the attribute to remove.
        Returns:    [nothing]
        """
        if not self._write:
            warning("DataIO: Write flag has not been set on file '%s'.  No action taken." % self._file_name)
            return

        try:
            if attribute not in self._rm_attributes[variable]:
                self._rm_attributes[variable].append(attribute)
            else:
                warning("DataIO: Attribute '%s' has already been removed from variable '%s' in file '%s'." % (attribute, variable, self._file_name))
        except KeyError:
            self._rm_attributes[variable] = []
            self._rm_attributes[variable].append(attribute)
        return

    def close(self,out_file_name=''):
        """
        close() [public]
        Purpose:    Closes the file and resets the class to its initial state.  The file is set for writing, and any changes have occurred, create a 
                        new file and copy all the variables over, without anything that's been "removed."
        Parameters: [none]
        Returns:    [nothing]
        """
        self._attempted_close = True
        if self._write:
            # Copy the data to another file, excluding the removed variables ...
            tmp_file_name = "%s.tmp" % self._file_name
            tmp_file = None
            if reader.__name__ == "Nio":
                tmp_file = reader.open_file(tmp_file_name, format='nc', mode='w')
            elif reader.__name__ == "scipy.io.netcdf":
                tmp_file = reader.netcdf_file(tmp_file_name, mode='w')
            elif reader.__name__ == "Dataset":
                tmp_file = reader(tmp_file_name,'w')

            # Copy dimensions over
            for dimension, length in self._df.dimensions.iteritems():
                if reader.__name__ == "Nio":
                    tmp_file.create_dimension(dimension, length)
                elif reader.__name__ in ["scipy.io.netcdf", "Dataset"]:
                    tmp_file.createDimension(dimension, len(length))

            # Create new dimensions added by splicing
            for variable, splice_info in self._sp_variables.iteritems():
                dim_names, index_lists = splice_info
                for dim_name, indexes in zip(dim_names, index_lists):
                    if dim_name not in tmp_file.dimensions.keys():
                        if reader.__name__ == "Nio":
                            tmp_file.create_dimension(dim_name, len(np.unique(indexes)))
                        elif reader.__name__ in ["scipy.io.netcdf", "Dataset"]:
                            tmp_file.createDimension(dim_name, len(np.unique(indexes)))

            # Copy variable definitions over
            for variable, data in self._df.variables.iteritems():
                if variable not in self._rm_variables:
                    # If this variable is going to be spliced, find its new dimensions
                    if variable not in self._sp_variables:
                        dim_names = data.dimensions
                    else:
                        dim_names, index_lists = self._sp_variables[variable]

                    if reader.__name__ == "Nio":
                        tmp_file.create_variable(variable, data.typecode(), dim_names)
                    elif reader.__name__ == "scipy.io.netcdf":
                        tmp_file.createVariable(variable, data.dtype, dim_names)
                    elif reader.__name__ == "Dataset":
                        if hasattr(data, "_FillValue"):
                            tmp_file.createVariable(variable, data.dtype, dim_names, fill_value=data._FillValue)
                        else:
                            tmp_file.createVariable(variable, data.dtype, dim_names)

            # Copy file attributes over
            for attribute in dir(self._df):
                if not hasattr(tmp_file, attribute) and attribute not in self._rm_attributes['__file__']:
                    setattr(tmp_file, attribute, getattr(self._df, attribute))

            # Copy variable attributes over
            for variable, data in tmp_file.variables.iteritems():
                try:
                    for attribute in dir(self._df.variables[variable]):
                        if not hasattr(data, attribute) and attribute not in self._rm_attributes[variable]:
                            setattr(data, attribute, getattr(self._df.variables[variable], attribute))
                except KeyError:
                    for attribute in dir(self._df.variables[variable]):
                        if not hasattr(data, attribute):
                            setattr(data, attribute, getattr(self._df.variables[variable], attribute))

            # Copy variable data over (NetCDF4 Python appears to require this be done after setting all the variable attributes)
            for variable, data in self._df.variables.iteritems():
                if variable not in self._rm_variables:
                    tmp_variable = np.copy(self._df.variables[variable][:])

                    # Overwrite the data if we need to ...
                    if variable in self._ov_variables:
                        tmp_variable = np.copy(self._ov_variables[variable])

                    # Do splicing ...
                    if variable in self._sp_variables:
                        dim_names, index_lists = self._sp_variables[variable]
                        tmp_variable = np.copy(self._df.variables[variable][np.ix_(*index_lists)])

                    # Put it in the file ...
                    tmp_file.variables[variable][:] = tmp_variable

#           for dimension, length in tmp_file.dimensions.iteritems():
#               print dimension, length

        self._df.close()

        if self._write:
            # Move temporary file to original file's location
            import shutil
            if len(out_file_name) > 0:
                shutil.move(tmp_file_name,out_file_name)
            else:
                shutil.move(tmp_file_name, self._file_name) 

        self._df = None
        self._file_name = ""
        self._read = False
        self._write = False
        self._rm_variables = []
        self._sp_variables = {}
        self._rm_attributes = { '__file__':[] }
        return

    def __del__(self):
        if self._df is not None and not self._attempted_close:
            self.close()
        return

if __name__ == "__main__":
    import shutil
    shutil.copy("../data/wrfout_d01_PLEV.nc", "../data/backup_wrfout_d01_PLEV.nc") 
    nc = DataIO("../data/backup_wrfout_d01_PLEV.nc", mode='rw')
#   p = nc.get_variable('P')
    missing_value = nc.get_variable_attribute('T', 'missing_value')
    nc.set_variable_attribute('P', 'missing_value', missing_value)

    nc.close()
