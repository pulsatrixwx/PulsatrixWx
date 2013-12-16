
class Config(dict):
    """
    Config
    Purpose:    Handles the loading in of configuration files (the meta data file and the configuration files for each run of HootPy).
    Started:    Fall 2010 by Tim Supinie (tsupinie@ou.edu)
    Completed:  14 January 2011 by Tim Supinie
    Modified:   [not yet]
    """

    def __init__(self, file_name):
        """
        __init__()
        Purpose:    Constructor for the Config class.
        Parameters: file_name [type=string]
                        Name of the file to read in.
        """
        config = {}
        print "Reading configuration from %s ..." % file_name

        # Treat the file as valid Python and execute it.
        execfile(file_name, {}, config)

        for key, value in config.iteritems():
            self.__setitem__(key, value)
        return

if __name__ == "__main__":
    c = Config("default.hp")
