from datetime import datetime

from hootpy import HootPy

class ProfilePlotter(HootPy):
    """
    ProfilePlotter
    Purpose:    Handles the plotting of products on Skew-T diagrams (such as observed and forecast soundings).
    Started:    14 June 2010 by Tim Supinie (tsupinie@ou.edu)
    Completed:  [not yet]
    Modified:   [not yet]
    """

    def __init__(self, config):
        """
        __init__()
        Purpose:    Constructor for the ProfilePlotter class.
        Parameters: forecast_hours [type=dictionary]
                        Dictionary containing the configuration parameters for the run.
        """

        super(ProfilePlotter, self).__init__(config)
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
        Purpose:    Plot Skew-T products.  For model products, plots products for all forecast hours.
        Parameters: [none]
        Returns:    [nothing]
        """

        if self._forecast_hours is None:
            # Plot profile here ...

            self._finalizeProfile(None)
        else:
            for fh in self._forecast_hours:
                # Plot the profiles here ...

                self._finalizeProfile(fh)
        return

    def _finalizeProfile(self, forecast_hour):
        """
        _finalizeProfile() [protected]
        Purpose:    Add final things to the profile, such as the Skew-T background, 
                        title, valid time, and image border, and then save the image.
        Parameters: forecast_hour [type=int]
                        Forecast hour for model products (pass in None for an observed product).
        Returns:    [nothing]
        """

        # Finish creating the product.  Should be last.
        self._finalizeProduct(forecast_hour)
        return

if __name__ == "__main__":
    cfg = {
        'forecast_hours':[0, 3, 6, 9, 12],
        'product_title':"NAM Forecast Sounding for KOUN",
        'image_file_name':"nam_fsound_KOUN_f%02d.png"
    }

    hpp = ProfilePlotter(cfg)
    hpp.loadData()
    hpp.plot()
