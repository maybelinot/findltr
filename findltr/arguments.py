#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: maybelinot
# @Email: edik.trott@yandex.ru
# @Date:   2015-09-12 13:31:59
# @Last Modified by:   maybelinot
# @Last Modified time: 2015-09-12 14:34:07

import sys, getopt

def arguments(argv):
    tmp = {'input_file'     :'',
           'algorithm'      :'original',
           'output_file'    :'',
           'min_pattern_len':40,
           'min_distance'   :1000,
           'max_distance'   :20000,
           'max_ltr_len'    :1000,
           'min_ltr_len'    :100}
    try:
        opts, args = getopt.getopt(argv,"a:i:o:",['min_pattern_len=', 'min_distance=', 'max_distance=', 'max_ltr_len=', 'min_ltr_len=' ])
    except getopt.GetoptError:
        print ('run.py -i <inputfile> -o <outputfile> -a <algorithm> <parametrs>')
        sys.exit(2)
    for opt, arg in opts: 
        if opt == '-o':
            tmp['output_file'] = arg
        elif opt =='-i':
            tmp['input_file'] = arg
        elif opt == '-a':
            tmp['algorithm'] = arg
        elif opt == '--min_pattern_len':
            tmp['min_pattern_len'] = int(arg)
        elif opt == '--min_distance':
            tmp['min_distance'] = int(arg)
        elif opt == '--max_distance':
            tmp['max_distance'] = int(arg)
        elif opt == '--max_ltr_len':
            tmp['max_ltr_len'] = int(arg)
        elif opt == '--min_ltr_len':
            tmp['min_ltr_len'] = int(arg)
    if tmp['input_file']=='':
        print('Type full path to .fasta input file or .\input_file.fa if it is located in current directory')
        sys.exit(2)
    return tmp