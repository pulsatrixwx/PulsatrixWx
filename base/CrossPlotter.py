from datetime import datetime

from hootpy import HootPy

class CrossPlotter(HootPy):
    """
    CrossPlotter
    Purpose:    Handles the plotting of cross section products.
    Started:    14 June 2010 by Tim Supinie (tsupinie@ou.edu)
    Completed:  [not yet]
    Modified:   [not yet]
    """

    def __init__(self, config):
        """
        __init__()
        Purpose:    Constructor for the CrossPlotter class.
        Parameters: config [type=dictionary]
                        Dictionary containing configuration parameters for the run.
        """

        super(CrossPlotter, self).__init__(config)
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
        Purpose:    Plot cross section products.  For model products, plots products for all forecast hours.
        Parameters: [none]
        Returns:    [nothing]
        """

        if self._forecast_hours is None:
            # Plot cross section here ...

            self._finalizeCrossSection(None)
        else:
            for fh in self._forecast_hours:
                # Plot the cross section here ...

                self._finalizeCrossSection(fh)
        return

    def _finalizeCrossSection(self, forecast_hour):
        """
        _finalizeCrossSection() [protected]
        Purpose:    Add final things to the profile, such as the background, 
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
        'product_title':"NAM Forecast Cross Section KDRT-KGRB",
        'image_file_name':"nam_fcross_KDRT-KGRB_f%02d.png"
    }

    hpc = CrossPlotter(cfg)
    hpc.loadData()
    hpc.plot()
