This directory contains the python code for reading AMSR-2 (F34) data.
Four files are needed: bytemaps.py, amsr2_daily.py, amsr2_averaged.py and example_usage.py

In order to test the programs, you need to 
1) download the appropriate verify.txt file AND test files located in the AMSR/support/verify_v*

2) place these files from step 1 and 2 in the same directory as the programs   
   
First run the daily and averaged routines to be sure they execute correctly.  You will get a 'verification failed' message if there is a problem.  If they work correctly, the message 'all tests completed successfully' will be displayed.

After confirming the routines work, use the example_usage.py routine as your base program and adapt to your needs.  This code shows you how to call the needed subroutines and diplays an example image.  Once you change the program, make sure to run it on the test files and check that the results match those listed in the verify_amsr2_v*txt file.

If you have any questions regarding these programs 
or the RSS binary data files, contact RSS support:
http://www.remss.com/support