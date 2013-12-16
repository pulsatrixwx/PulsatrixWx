from currentModelDate import getCurrentModelDate
from time import sleep
import argparse
import subprocess
import os
from glob import glob

from dataLibECMWF import findUrlsECMWF, stitchECMWF

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--model',default='NAM',
        help="Model being downloaded from Motherlode")
    parser.add_argument('--dir', default = os.environ["HOME"] + '/hootpy/data',
        help="Directory where model file is to be downloaded")
    parser.add_argument('--max',default=4,type=int,
        help="Number of model files to keep in data directory.  Oldest will be deleted if the specified number is exceeded")
    args = parser.parse_args()
    downloadModel(args.model,args.dir)
    removeOldModelRuns(args.model,args.dir,args.max)

def removeOldModelRuns(model,data_dir,max_files):
    """removeOldModelRuns()
        Purpose:  delete old model runs so hard drive is not backed up with extra data.
        Paramters:  
            model (str) - Name of model
            data_dir (str) - Location of model files
            max_files (int) - Maximum number of model files to be kept.  The oldest files will be deleted first.
    """
    model = model.upper()
    if model == "NAM":
        model_files = sorted(glob(data_dir + '/NAM*.grib2'))
    elif model == "RUC":
        model_files = sorted(glob(data_dir + '/RUC*.grib2'))
    elif model == "SREF":
        model_files = sorted(glob(data_dir + '/SREF*.grib2'))
    elif model == "RAP":
        model_files = sorted(glob(data_dir + '/RR*.grib2'))
    elif model == "GFS":
        model_files = sorted(glob(data_dir + '/GFS*.grib1'))
    elif model == "ECMWF":
        model_files = sorted(glob(data_dir + '/ECMWF*.nc'))
    else:
        model_files = []
    print model_files
    if len(model_files) > max_files:
        num_deleted = len(model_files) - max_files
        for i in xrange(num_deleted):
            subprocess.call(["rm",model_files[i]])
    return

def downloadModel(model,data_dir=os.environ["HOME"] + "/hootpy/data"):
    """downloadModel()
        Purpose:  Given a model type, download the latest grib file from 
        motherlode.ucar.edu
        Parameters:
            model (str) - Name of model (NAM,GFS,RUC)
            data_dir (str) - directory where data files are downloaded
    """
    model = model.upper()
    print data_dir
    model_date = getCurrentModelDate(model,"%Y%m%d_%H%M")
    if model == "NAM":
        curr_files = glob('%s/NAM_CONUS_40km_conduit_%s.grib2' % (data_dir,model_date))
        if len(curr_files) > 0:
            print 'Deleting',curr_files[0]
            subprocess.call(['rm',curr_files[0]])
        link = "http://thredds.ucar.edu/thredds/fileServer/grib/NCEP/NAM/CONUS_40km/conduit/files/NAM_CONUS_40km_conduit_%s.grib2" % model_date 
        #link = "http://motherlode.ucar.edu/thredds/fileServer/fmrc/NCEP/NAM/CONUS_40km/conduit/files/NAM_CONUS_40km_conduit_%s.grib2" % model_date
        cmd = ["wget","-P",data_dir,link]
        subprocess.call(cmd)
    elif model == "RAP":
        curr_files = glob('%s/RR_CONUS_20km_%s.grib2' % (data_dir,model_date))
        if len(curr_files) > 0:
            print 'Deleting',curr_files[0]
            subprocess.call(['rm',curr_files[0]])
        #link = "http://motherlode.ucar.edu/thredds/fileServer/fmrc/NCEP/RAP/CONUS_20km/files/RR_CONUS_20km_%s.grib2" % model_date
        link = "http://thredds.ucar.edu/thredds/fileServer/grib/NCEP/RAP/CONUS_20km/files/RR_CONUS_20km_%s.grib2" % model_date
        cmd = ["wget","-P",data_dir,link]
        subprocess.call(cmd)
    elif model == "SREF":
        curr_files = glob('%s/SREF_CONUS_40km_ensprod_%s.grib2' % (data_dir,model_date))
        if len(curr_files) > 0:
            print 'Deleting',curr_files[0]
            subprocess.call(['rm',curr_files[0]])
        link = "http://thredds.ucar.edu/thredds/fileServer/grib/NCEP/SREF/CONUS_40km/ensprod/files/SREF_CONUS_40km_ensprod_%s.grib2" % model_date
        cmd = ["wget","-P",data_dir,link]
        subprocess.call(cmd)
    elif model == "GFS":
        curr_files = glob('%s/GFS_CONUS_80km_conduit_%s.grib2' % (data_dir,model_date))
        if len(curr_files) > 0:
            print 'Deleting',curr_files[0]
            subprocess.call(['rm',curr_files[0]])
        link ="http://thredds.ucar.edu/thredds/fileServer/grib/NCEP/GFS/CONUS_80km/files/GFS_CONUS_80km_%s.grib1" % model_date
        cmd = ["wget","-P",data_dir,link]
        subprocess.call(cmd)
    elif model == "ECMWF":
        curr_files = glob('%s/ECMWF_Global_2p5_%s.nc' % (data_dir, model_date))
        if len(curr_files) > 0:
            print 'Deleting',curr_files[0]
            subprocess.call(['rm',curr_files[0]])

        pathECMWF = '/data1/native/grid/ECMWF/Global_2p5'
        print(model_date)
        print(data_dir)

        os.system('cp %s/ECMWF_Global_2p5_A_%s.grib1 %s' %(pathECMWF, model_date, data_dir))
        os.system('cp %s/ECMWF_Global_2p5_B_%s.grib1 %s' %(pathECMWF, model_date, data_dir))

        stitchECMWF(data_dir, model_date)

        for file in glob("%s/ECMWF_Global_2p5_?_%s.grib1" % (data_dir, model_date)):
            cmd = ["rm", file]
            subprocess.call(cmd)
    else:
        print "Model " + model + " is not available"
        exit()

if __name__ == "__main__":
    main()
