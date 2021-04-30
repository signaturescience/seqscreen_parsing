#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import collections


def hash_df(in_df_column):
    """
    Creates a dictionary of de-duplicated column values as keys, whose value is a list of 
    all of the row indices that the column value exists. Provides a O(1) lookup complexity 
    for finding all the rows that contain a specificed value in the provided column name.
    """
    import collections
    hashed_col = collections.defaultdict(list)

    for column_val in in_df_column.iteritems():
        hashed_col[column_val[1]].append(column_val[0])
    
    return hashed_col


def intersection(list1, list2):
    """
    Find the interesection of two lists using set()
    """
    # return list(set(list1) & set(list2))
    return set(list1).intersection(list2)


def slice_dict(main_dict, rows):
    """
    Returns a slice of dictionary setup as a "column, row" format, typically converted from pandas
    """
    if not isinstance(rows, list):
        rows = [rows]

    current_slice = collections.defaultdict(list)
    for col, row_dict in main_dict.items():
        for row_index in rows:
            current_slice[col].append(row_dict[row_index])
    return current_slice
