
import Nio as nio
import numpy as np

import os
import urllib, urllib2
import re

def _stitchVariable(new_var_name, var_A, var_B, nc_stitch, stitch_dir):
    new_dimensions = list(var_A.dimensions[:-2])
    new_dimensions.extend(['lat', 'lon'])
    nc_stitch.create_variable(new_var_name, var_A.typecode(), tuple(new_dimensions))

    if stitch_dir == "lat":
        lat_size_A = var_A.shape[-2]
        lat_size_B = var_B.shape[-2]

        nc_stitch.variables[new_var_name][..., :lat_size_A, :] = var_A
        nc_stitch.variables[new_var_name][..., -lat_size_B:, :] = var_B
    elif stitch_dir == "lon":
        lon_size_A = var_A.shape[-1]
        lon_size_B = var_B.shape[-1]

        nc_stitch.variables[new_var_name][..., :, :lon_size_B] = var_B
        nc_stitch.variables[new_var_name][..., :, -lon_size_A:] = var_A

    return

def findUrlsECMWF(target_time):
    base_url = "http://motherlode.ucar.edu"
    path = ['Data', 'IDD+Data', 'native', 'grid', 'ECMWF', 'Global_2p5']

    base_page = urllib2.urlopen("%s/repository" % base_url).read()
    candidates = re.findall(r"(?<=folderClick\().*?(?=\))", base_page)

    urls = []

    for dir in path:
        next_url = ""
        dir = re.sub(r"\+", r"\+", dir)
        for candidate in candidates:
            url = candidate.split(',')[1][1:-1]
            if re.search(r"\/%s\?entryid" % dir, url):
                next_url = url

        if next_url == "":
            raise ValueError("Directory '%s' not found (in _find_ECMWF_urls)." % dir)

        xml = urllib2.urlopen("%s%s" % (base_url, next_url)).read()
        candidates = re.findall(r"(?<=folderClick\().*?(?=\))", xml)

    for candidate in candidates:
        url = candidate.split(',')[1][1:-1]
        if re.search(target_time, url):
            urls.append(url)

    for idx in xrange(len(urls)):
        url = urls[idx]

        print url
        xml = urllib2.urlopen("%s%s" % (base_url, url)).read()
        urls[idx] = urllib.unquote("%s%s" % (base_url, re.findall(r"(?<=a href=\")/repository.*?(?=\")", xml)[0]))

    return urls

def stitchECMWF(base_path, time):
    variables = ['HGT', 'PRMSL', 'TMP', 'U_GRD', 'V_GRD', 'R_H']
    coords = ['lv_ISBL3', 'lv_ISBL4', 'lv_ISBL5', 'forecast_time0']

    grib_A = nio.open_file("%s/ECMWF_Global_2p5_A_%s.grib1" % (base_path, time), mode='r')
    grib_B = nio.open_file("%s/ECMWF_Global_2p5_B_%s.grib1" % (base_path, time), mode='r')

    nc_file_name = "%s/ECMWF_Global_2p5_%s.nc" % (base_path, time)

    if os.path.exists(nc_file_name):
        os.system("rm %s" % nc_file_name)
    nc_stitch = nio.open_file(nc_file_name, mode='w')

    stitch_dir = ""
    if np.all(grib_A.variables['lat_1'][:] == grib_B.variables['lat_2'][:]):
        nc_stitch.create_dimension('lat', grib_A.dimensions['lat_1'])
        new_len = grib_A.dimensions['lon_1'] + grib_B.dimensions['lon_2'] - 1
        nc_stitch.create_dimension('lon', new_len)


        stitch_dir = 'lon'
#       print "Lats are the same"

    elif np.all(grib_A.variables['lon_1'][:] == grib_B.variables['lon_2'][:]):
        new_len = grib_A.dimensions['lat_1'] + grib_B.dimensions['lat_2'] - 1
        nc_stitch.create_dimension('lat', new_len)
        nc_stitch.create_dimension('lon', grib_A.dimensions['lon_1'])

        stitch_dir = 'lat'
#       print "Lons are the same"

    nc_stitch.create_variable('lat', grib_A.variables['lat_1'].typecode(), ('lat',))
    nc_stitch.create_variable('lon', grib_A.variables['lon_1'].typecode(), ('lon',))

    if stitch_dir == 'lat':
        lats_A = grib_A.variables['lat_1'][:]
        lats_B = grib_B.variables['lat_2'][:]
        lons = grib_A.variables['lon_1'][:]

        lons = np.where(lons < -180, lons + 360, lons)
        lons = np.where(lons >= 180, lons - 360, lons)

        nc_stitch.variables['lat'][:(lats_A.shape[0])] = lats_A 
        nc_stitch.variables['lat'][-(lats_B.shape[0]):] = lats_B 
        nc_stitch.variables['lon'][:] = lons
    elif stitch_dir == 'lon':
        lats = grib_A.variables['lat_1'][:]
        lons_A = grib_A.variables['lon_1'][:]
        lons_B = grib_B.variables['lon_2'][:]

        lons_A = np.where(lons_A < -180, lons_A + 360, lons_A)
        lons_A = np.where(lons_A >= 180, lons_A - 360, lons_A)
        lons_B = np.where(lons_B < -180, lons_B + 360, lons_B)
        lons_B = np.where(lons_B >= 180, lons_B - 360, lons_B)

        nc_stitch.variables['lat'][:] = lats
        nc_stitch.variables['lon'][:(lons_B.shape[0])] = lons_B
        nc_stitch.variables['lon'][-(lons_A.shape[0]):] = lons_A

    for coord_name in ['lat', 'lon']:
        coord_name_A = "%s_1" % coord_name
        for attr, value in grib_A.variables[coord_name_A].attributes.iteritems():
            setattr(nc_stitch.variables[coord_name], attr, value)

    for coord_name in coords:
        nc_stitch.create_dimension(coord_name, grib_A.dimensions[coord_name])
        nc_stitch.create_variable(coord_name, grib_A.variables[coord_name].typecode(), (coord_name,))
        nc_stitch.variables[coord_name][:] = grib_A.variables[coord_name][:]

        for attr, value in grib_A.variables[coord_name].attributes.iteritems():
            setattr(nc_stitch.variables[coord_name], attr, value)

#   print nc_stitch.variables['lon'][:]
#   print nc_stitch.variables.keys()

    for var_name in variables:
        if var_name == "PRMSL":
            A_var_name = "%s_1_SFC" % var_name
            B_var_name = "%s_2_SFC" % var_name
        else:
            A_var_name = "%s_1_ISBL" % var_name
            B_var_name = "%s_2_ISBL" % var_name

        _stitchVariable(var_name, grib_A.variables[A_var_name], grib_B.variables[B_var_name], nc_stitch, stitch_dir)

        for attr, value in grib_A.variables[A_var_name].attributes.iteritems():
            setattr(nc_stitch.variables[var_name], attr, value)

    nc_stitch.close()
    grib_A.close()
    grib_B.close()
    return

def main():
    print findUrlsECMWF("20120709_1200")
    stitchECMWF("/home/operations/hootpy/data", "20120709_1200")
    return

if __name__ == "__main__":
    main()
