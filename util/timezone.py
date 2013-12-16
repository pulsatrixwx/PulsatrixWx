from datetime import timedelta, tzinfo

ZERO = timedelta(0)

class UTC(tzinfo):
    """
    UTC
    Purpose:    Class to handle the Python datetime module's notion of the UTC time zone.
    Started:    16 March 2011 by Tim Supinie (tsupinie@ou.edu)
    Completed:  16 March 2011
    Modified:   [not yet]
    """

    def utcoffset(self, dst):
        """
        utcoffset() [public]
        Purpose:    Returns the UTC offset of this time zone
        Parameters: dst [type=boolean]
                        Specifies whether or not we're in daylight savings time (required by the datetime module, but unused, as UTC doesn't observe DST).
        Returns:    The UTC offset of UTC (... zero).
        """
        return ZERO

    def tzname(self, dst):
        """
        tzname() [public]
        Purpose:    Returns the human-readable name of this time zone
        Parameters: dst [type=boolean]
                        Specifies whether or not we're in daylight savings time (required by the datetime module, but unused, as UTC doesn't observe DST).
        Returns:    The name for the UTC time zone (will always be "UTC")
        """
        return "UTC"

    def dst(self, dst):
        """
        dst() [public]
        Purpose:    Returns the daylight savings time offset of this time zone
        Parameters: dst [type=boolean]
                        Specifies whether or not we're in daylight savings time (required by the datetime module, but unused, as UTC doesn't observe DST).
        Returns:    The daylight savings time offset for UTC (will always be a timedelta object set to 0).
        """
        return ZERO

# Create an instance of the UTC class to use
utc = UTC()
