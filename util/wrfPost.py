from dataIO import DataIO
from derived import *
from os.path import expanduser
class wrfPost(object):
    """
    wrfPost
    Purpose:  wrfPost ingests a WRF output netCDF file and splits it into three separate netCDF files for surface, upper air, and composite variables.  For surface data it saves the WRF surface grids and can calculate other derived variables from the basic ones.  For upper air data it interpolates the grid coordinates from model levels to pressure and isentropic levels.  For composite variables it does operations on each grid column to calculate variables such as CAPE and Precipitable Water.
    Started:    25 June 2010 by David John Gagne (djgagne@ou.edu)
    Completed:    [not yet]
    Modified:    [not yet]
    """
    
    def __init__(self,wrfFileName):
        """
        __init__
        Purpose:  Constructor for wrfPost.
        Parameters:
            wrfFileName:  Full path and filename for WRF output file to be processed.
        """
        self.filename = wrfFileName
        # load WRF data into class from file
        self.data = self.loadData()
        #self.processSurface()
        self.processUpperair()

    def loadData(self):
        """
        loadData()
        Purpose:  Load WRF data file into class using the Nio module
        """
        if self.filename[-3:] != '.nc':
            self.filename = self.filename + '.nc'
        data = DataIO(self.filename,mode='w')
        return data

    def processSurface(self):
        """
        processSurface()
        Purpose:  Separate surface grids from WRF data, calculate derived parameters, and write to self.filename_surface.nc netCDF file.
        """
        # Create output surface netcdf file and copy dimensions and global attributes
        outfilename = self.filename[:-3] + '_surface.nc'
        print outfilename
        
        # Calculate derived parameters and output to surface netCDF file    
        # Mean Sea Level Pressure
        #print self.data.get_variable_dimensions('PSFC')
        print self.data.get_variable_list()
        if 'MSLP' in self.data.get_variable_list():
            print 'MSLP in list already'
        else:    
            self.data.create_variable('MSLP','f4',dimensions=self.data.get_variable_dimensions('PSFC'))
        print 'Setting attributes'
        self.data.set_variable_attribute('MSLP','description','Mean Sea Level Pressure')
        print 'desc'
        self.data.set_variable_attribute('MSLP','MemoryOrder','XY')
        print 'XY'
        self.data.set_variable_attribute('MSLP','coordinates','XLONG XLAT')
        print 'XLONG XLAT'
        self.data.set_variable_attribute('MSLP','stagger','')
        print 'stagger'
        self.data.set_variable_attribute('MSLP','FieldType',np.array([104],dtype=np.int32))
        print 'FieldType'
        self.data.set_variable_attribute('MSLP','units','Pa')
        print 'Units'
        print 'Set Value'
        self.data.set_variable('MSLP',mslp(self.data.get_variable('PSFC'),self.data.get_variable('HGT'),self.data.get_variable('T2'),self.data.get_variable('Q2')))
        print 'Closing file'
        # Dewpoint
        #outfile.createVariable('DEWP','f',outfile.variables['T2'].dimensions)
        #outfile.variables['DEWP'].attributes = {'description':'Dewpoint temperature','MemoryOrder':'XY','coordinates': 'XLONG XLAT','stagger':'','FieldType':np.array([104],dtype=np.int32),'units':'K'}
        #outfile.variables['DEWP'].assignValue(dewpoint(outfile.variables['T2'][:],outfile.variables['Q2'][:],outfile.variables['PSFC'][:]))
        
        self.data.close(outfilename)

    def processUpperair(self,p_levels=[100000,92500,85000,70000,50000,30000,25000]):
        """
        processUpperair()
        Purpose:  Call pressure and isentropic interpolation functions, calculate upperair derived parameters, and write to self.filename_upperair.nc netCDF file.
        """
        outfilename = self.filename[:-3] + '_upperair.nc'
        if 'p_levels' not in self.data.get_dimensions().keys():
            self.data.create_dimension('p_levels',len(p_levels)) 
        self.data.create_variable('p_levs','f',('p_levels'))
        self.data.set_variable('p_levs',p_levels)
        pressure = self.data.get_variable('P') + self.data.get_variable('PB')
        print 'Processing Height'
        ph = self.data.get_variable('PHB') + self.data.get_variable('PH')
        ph = self.unstagger_grid(ph,1)
        self.data.create_variable('PH_plev','f',('Time','p_levels','south_north','west_east'))
        self.data.set_variable('PH_plev',self.pressure_interp(pressure,ph,p_levels))
        
        print 'Processing Temperature'
        tair = self.data.get_variable('T')
        self.data.create_variable('T_plev','f',('Time','p_levels','south_north','west_east'))
        self.data.set_variable('T_plev',self.pressure_interp(pressure,tair,p_levels))
        
        print 'Processing U'
        u = self.unstagger_grid(self.data.get_variable('U'),3)
        self.data.create_variable('U_plev','f',('Time','p_levels','south_north','west_east'))
        self.data.set_variable('U_plev',self.pressure_interp(pressure,u,p_levels))
        
        print 'Processing V'
        v = self.unstagger_grid(self.data.get_variable('V'),2)
        self.data.create_variable('V_plev','f',('Time','p_levels','south_north','west_east'))
        self.data.set_variable('V_plev',self.pressure_interp(pressure,v,p_levels))
        
        self.data.close(outfilename)

    def pressure_interp(self,pressure,var,p_levels):
        """
        pressure_interp()
        Purpose:  Interpolate raw 3d levels to pressure coordinates using a linear interpolation        scheme.  
        """
        var_interpolated = np.zeros((pressure.shape[0],len(p_levels),pressure.shape[2],pressure.shape[3]))
        for (t,y,x),v in np.ndenumerate(var_interpolated[:,0,:,:]):
            var_interpolated[t,:,y,x] = np.interp(p_levels,pressure[t,::-1,y,x],var[t,::-1,y,x])
        return var_interpolated

    def unstagger_grid(self,grid,staggered_dim):
        """
        unstagger_grid()
        Purpose:  Unstagger grid by averaging in direction of stagger"""
        if staggered_dim == 0:
            return 0.5 * (grid[0:-1] + grid[1:])
        elif staggered_dim == 1:
            return 0.5 * (grid[:,0:-1] + grid[:,1:])
        elif staggered_dim == 2:
            return 0.5 * (grid[:,:,0:-1] + grid[:,:,1:])
        else:
            return 0.5 * (grid[:,:,:,0:-1] + grid[:,:,:,1:])
def main():
    post = wrfPost(expanduser('~') + '/wrfout_d01')

if __name__ == "__main__":
    main()
