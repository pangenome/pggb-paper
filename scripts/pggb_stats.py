#!/bin/python

import sys

dataset = ''
runtime_s = ''
memory_kb = ''
tool = ''
suffix = ''

for line in sys.stdin:
    line_list = line.strip().split()
    #print(line_list)
    if 'max memory' in line:
        runtime_s = line_list[6]
        memory_kb = line_list[8]
    else:
        if 'wfmash' in line_list[0]:
            dataset = line.split('--tmp-base')[-1].split()[0]
            if '-m' in line:
                suffix='-mapping'
            else:
                suffix='-alignment'
        elif 'odgi' in line_list[0]:
            suffix = '-' + line_list[1]
        else:
            suffix = ''
        tool = line_list[0].split(':')[-1].split('/')[-1].split('-')[0]

    if dataset and runtime_s and memory_kb and tool:
        # dataset, tool, runtime_s, memory_kb,
        print(dataset, tool+suffix, runtime_s.rstrip('s'), memory_kb.rstrip('Kb'), sep = '\t')
        runtime_s = ''
        memory_kb = ''
        tool = ''
