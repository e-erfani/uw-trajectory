import os
import numpy as np

## Specify the time range that data is available
years  = [2018, 2019, 2020, 2021]
months = [6, 7, 8, 9]
days   = [30, 31, 31, 30]  # number of days in each month

######
for year in years:
    for mm in range(len(months)):
        for dd in range(1, days[mm]+1):
            cmd2     = 'nohup ipython Direct_Download_Pressure.py ' + str(year) + ' ' + str(months[mm]) + ' ' + str(dd) + '  &'
            cmd_run2 = os.popen(cmd2)
            data_get = cmd_run2.read()         
