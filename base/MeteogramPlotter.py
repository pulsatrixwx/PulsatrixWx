from datetime import datetime

from hootpy import HootPy

class MeteogramPlotter(HootPy):
    """
    MeteogramPlotter
    Purpose:    Handles the plotting of meteogram products (both observed and forecast).
    Started:    14 June 2010 by Tim Supinie (tsupinie@ou.edu)
    Completed:  [not yet]
    Modified:   [not yet]
    """

    def __init__(self, config):
        """
        __init__()
        Purpose:    Constructor for the MeteogramPlotter class.
        Parameters: config [type=dictionary]
                        Dictionary containing configuration parameters for the run.
        """

        super(MeteogramPlotter, self).__init__(config)
        return

    def loadData(self):
        """
        loadData() [public]
        Purpose:    Handles the loading in of data.
        Parameters: [none]
        Returns:    [nothing]
        """
        return

    def plot(self):
        """
        plot() [public]
        Purpose:    Plot meteogram products.
        Parameters: [none]
        Returns:    [nothing]
        """

        self._finalizeMeteogram()
        return

    def _finalizeMeteogram(self):
        """
        _finalizeMeteogram() [protected]
        Purpose:    Add final things to the profile, such as the background, 
                        title, valid time, and image border, and then save the file.
        Parameters: forecast_hour [type=int]
                        Forecast hour for model products (pass in None for an observed product).
        Returns:    [nothing]
        """

        # Finish creating the product.  Should be last.
        self._finalizeProduct(None)
        return

if __name__ == "__main__":
    cfg = {
        'forecast_hours':[0, 3, 6, 9, 12],
        'product_title':"NAM Forecast Meteogram for KOUN",
        'image_file_name':"nam_fcstmeteo_KOUN.png"
    }
    hpm = MeteogramPlotter(cfg)
    hpm.loadData()
    hpm.plot()
