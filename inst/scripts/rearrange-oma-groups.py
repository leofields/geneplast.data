# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 03:03:16 2020

@author: Leonardo
"""
import gzip

with open('../../working_directory/oma-groups-rearranged.txt', 'w') as res:
    with gzip.open('../../working_directory/oma-groups.txt.gz') as f:
        lines = f.readlines()
        for l in lines:
            cols = l.split('\t')
            for index, col in enumerate(cols):
                if index >= 2:
                    res.write(''.join([cols[0],'\t',cols[1],'\t',col,'\n']))
            
