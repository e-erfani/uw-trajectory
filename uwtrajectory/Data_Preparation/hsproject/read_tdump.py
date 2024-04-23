# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 14:39:53 2015

@author: Hans
read_tdump code courtesy of Jayson Stemmler
"""


def read_tdump(tdump):
    """ Read a tdump file as output by the HYSPLIT Trajectory Model

        Returns a pandas DataFrame object.
    """

    import pandas as pd
    from datetime import datetime as dt

    def parseFunc(s):
        return dt.strptime('-'.join([i.zfill(2) for i in s.split()]),
                           '%y-%m-%d-%H-%M')

    def parseFunc_icept(y, m, d, H, M):
        return dt(int('20'+y), int(m), int(d), int(H), int(M))

    columns = ['tnum', 'gnum', 'y', 'm', 'd', 'H', 'M',
               'fhour', 'age', 'lat', 'lon', 'height', 'pres']

    tmp = pd.read_table(tdump, nrows=100, header=None)
    l = [len(i[0]) for i in tmp.values]
    skiprows = l.index(max(l))

    D = pd.read_table(tdump, names=columns,
                      skiprows=skiprows,
                      engine='python',
                      sep=r'\s*',
                      parse_dates={'dtime': ['y', 'm', 'd', 'H', 'M']},
                      date_parser=parseFunc_icept,
                      index_col='dtime')

    return D
