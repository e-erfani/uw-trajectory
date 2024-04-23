# -*- coding: utf-8 -*-
"""Provides LoopTimer class for printing useful output on a for loop iteration
"""
import datetime as dt
import sys

class LoopTimer:
    """Generic loop timer.
    Given a start time, current index(i.e. where in the loop we are), and
    total_index, (i.e. total length of loop) and optional offset (how many
    loop iterations we skipped), this prints our where we are in the loop,
    percentage completion, and estimated completion time assuming each loop
    iteration takes the same amount of time.

    Usage:
    lt = LoopTimer(len(x))
    for i,elem in enumerate(x):
        lt.update(i)
    """
    def __init__(self, loop_length, loop_offset=0):
        self.start_time = dt.datetime.now()
        self.loop_length = loop_length
        self.loop_offset = loop_offset
        self.init_frac = loop_offset/float(loop_length)
        self.current_frac = self.init_frac
        self.frac_done = 0
        self.time_left = -1
        self.current_index=loop_offset

    def update(self, new_index=None, overwrite=True, it=' '):
        """
        current_index should refer to which iteration you are through your
        for loop; using enumerate([ITERITEM]) is the quickest way to get this
        overwrite sets whether the output string overwrites the previous line
        on stdout. So if you are printing things in your for loop, use
        overwrite=False (and prepare for a lot of output on stdout); if you
        are looping over a large iterable and don't want many lines on stdout,
        overwrite=True will keep the output string on the same line of stdout.
        """
        if new_index is None:
            new_index = self.current_index + 1
        if new_index > self.loop_offset:
            current_time = dt.datetime.now()
            time_passed = current_time - self.start_time
            self.current_frac = (new_index/float(self.loop_length))
            self.frac_done = self.current_frac - self.init_frac
            frac_todo = 1.0-(self.frac_done+self.init_frac)
            rate = self.frac_done/time_passed.total_seconds()
            self.time_left = frac_todo/rate
        self.current_index = new_index

        pct_str = '{:.2f}%      '.format(self.current_frac*100)
        mins_left = self.time_left/60.
        hours_left = mins_left/60.
        if hours_left > 2:
            lstr = 'time left: {:.1f} hours        '.format(hours_left)
        elif mins_left > 3:
            lstr = 'time left: {:.1f} minutes      '.format(mins_left)
        else:
            lstr = 'time left: {:.0f} seconds      '.format(self.time_left)
        eta = dt.datetime.now() + dt.timedelta(seconds=self.time_left)

        if overwrite:
            print('\r' + pct_str + 'ETA {:%H:%M:%S}      '.format(eta) + lstr + it, end=' ')
        else:
            print(pct_str + 'ETA {:%H:%M:%S}      '.format(eta) + lstr + it)
        sys.stdout.flush()
                
        
if __name__ == "__main__":
    import time
    import numpy as np
    some_iterable = np.arange(80)
    lt = LoopTimer(len(some_iterable))
    for i in some_iterable:
        lt.update()
        time.sleep(0.1)
        
        


__author__ = "Johannes Mohrmann"
__copyright__ = "Copyright 2016 Johannes Mohrmann"
__license__ = "MIT"
