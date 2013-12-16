#!/usr/local/epd/bin/python

import sys
sys.path.append("../")

from util.dataIO import DataIO
import numpy as np

df = DataIO("wrfout_d01_PLEV.nc", mode="rw")
missing_value = df.get_variable_attribute('T', 'missing_value')
df.set_variable_attribute('P', 'missing_value', missing_value)

df.close()
