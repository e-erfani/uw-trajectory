This directory contains the python code for reading SSMI and SSMIS V7 data.
These files are needed: bytemaps.py, ssmi_daily_v7.py, ssmi_averaged_v7.py 
ssmis_daily_v7.py, ssmis_averaged_v7.py, ssmis_uncertainty_v7.py and example_usage.py

In order to test the programs, you need to 
1) download the appropriate verify.txt file located in the SSMI/ssmi_support/verify/ directory
2) download the 4 files required to test:
For ssmi:
Time              Directory                        File
daily passes      F10/bmaps_v07/Y1995/M01        F10_19950120v7.gz
3-day mean        F10/bmaps_v07/Y1995/M01        F10_19950120v7_d3d.gz
weekly mean       F10/bmaps_V07/WEEKS            F10_19950121v7.gz
monthly mean      F10/bmaps_v07/Y1995/M01        F10_199501v7.gz

Each of these files contains F10 SSM/I data for:
the day January 20th, 1995 (daily)
the days January 19,20 and 21st, 1995 (3-day mean)
January 15th (Sunday) to January 21st (Saturday), 1995 (weekly mean)
or the month of January 1995 (monthly mean).


For SSMIS
Time              Directory                        File
daily passes      F17/bmaps_v07/Y2009/M01        F17_20090120v7.gz
3-day mean        F17/bmaps_v07/Y2009/M01        F17_20090120v7_d3d.gz
weekly mean       F17/bmaps_v07/WEEKS            F17_20090124v7.gz
monthly mean      F17/bmaps_v07/Y2009/M01        F17_200901v7.gz

Each of these files contains F17 SSMIS data for:
the day January 20th, 2009 (daily)
the days January 19,20 and 21st, 2009 (3-day mean)
January 18th (Sunday) to January 24th (Saturday), 2009 (weekly mean)
or the month of January 2009 (monthly mean). 


For SSMIS F17 uncertainty:
Time              Directory                        	File
daily passes      F17/bmaps_v07/Y2014/M02        	F17_20140221v7.gz
uncertainty file  F17/bmaps_v07_uncertainty/Y2014/M02 F17_20140221v7.unc.gz

Each of these files contains F17 SSMIS data for:
the day February 21, 2014 (daily and associated uncertainty file)



3) place the files from step 1 and 2 in the same directory as the programs   
   
First run the daily and averaged routines to be sure they execute correctly.  You will get a 'verification failed' message if there is a problem.  If they work correctly, the message 'all tests completed successfully' will be displayed.

After confirming the routines work, use the example_usage.py routine as your base program and adapt to your needs.  This code shows you how to call the needed subroutines and diplays an example image.  Once you change the program, make sure to run it on the test files and check that the results match those listed in the verify_ssmi_v7.txt file or verify_ssmis_v7.txt file or verify_ssmis_v7_uncertainty.txt file.

If you have any questions regarding these programs 
or the RSS binary data files, contact RSS support:
http://www.remss.com/support