#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd


def hash_df(in_df, column_name):
    """
    Creates a dictionary of de-duplicated column values as keys, whose value is a list of 
    all of the indices that the column exists, providing a O(1) lookup complexity for finding column rows
    """
    import collections
    hashed_col = collections.defaultdict(list)

    for column_val in in_df[column_name].iteritems():
        hashed_col[column_val[1]].append(column_val[0])
    
    return hashed_col
