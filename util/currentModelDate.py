from datetime import datetime,timedelta
import numpy as np
import os

def getCurrentModelDate(model,time_format='%Y%m%d_%H:%M:%S'):
    """getCurrentModelDate()
       Purpose:  Given a model, find the most current operational run
       Parameters:  model:  (String) name of model
                    time_format:  (String) datetime strftime format 
    """
    model = model.upper()
    model_date_file = open(os.environ['HOME'] + '/hootpy/util/modelDates.hp')
    exec(model_date_file.read())
    curr_time = datetime.utcnow() - timedelta(minutes=modelDateInfo[model]['offset'])
    model_hours = np.arange(modelDateInfo[model]['start'],24,24/modelDateInfo[model]['runs'])
    model_diff = model_hours - curr_time.hour
    model_time = datetime(curr_time.year,curr_time.month,curr_time.day,model_hours[np.argmin(np.where(model_diff > 0,1000,np.abs(model_diff)))],0,0)
    return model_time.strftime(time_format)

def getCurrentUpperairDate():
    """Gets current upperair time""" 
   
def getSeason():
    month = datetime.utcnow().month
    if month < 3 or month >= 12:
        season = "wi"
    elif month >= 3 and month < 6:
        season = "sp"
    elif month >= 6 and month < 9:
        season = "su"
    elif month >= 9 and month < 12:
        season = "fa"

    return season

def main():
    model = 'NAM'   
    time_format = '%Y%m%d_%H:%M:%S'
    print getCurrentModelDate(model,time_format)
if __name__ == "__main__":
    main()
