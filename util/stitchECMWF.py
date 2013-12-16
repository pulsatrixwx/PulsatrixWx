from currentModelDate import getCurrentModelDate
from time import sleep
import argparse 
import subprocess
import os
from glob import glob
from dataLibECMWF import findUrlsECMWF, stitchECMWF


print('Stitching')

data_dir=os.path.abspath('.') + '/data'
model_date = "20131112_0000"
#model_date = getCurrentModelDate(model,"%Y%m%d_%H%M")


stitchECMWF(data_dir, model_date)
print('Done')
