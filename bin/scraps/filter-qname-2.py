#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 16:57:03 2022

@author: kalavattam
"""

import gzip
import sys

with gzip.open(sys.argv[1], 'r') as indexfile:
    ids = set(l.rstrip('\r\n') for l in indexfile)

for line in sys.stdin:
    qname, _ = line.split('\t', 1)
    if qname in ids:
        sys.stdout.write(line)
