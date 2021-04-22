#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd


def hash_df(in_df, column_name):

    import collections
    hashed_col = collections.defaultdict(list)

    for column_val in in_df[column_name].iteritems():
        hashed_col[column_val[1]].append(column_val[0])
    # for index, column_val in enumerate(in_df[column_name].iteritems()):
    #     hashed_col[column_val].append(index)
    
    return hashed_col
