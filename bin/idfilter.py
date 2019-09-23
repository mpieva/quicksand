#!/usr/bin/env python3

import sys

with open(sys.argv[1], 'r') as indexfile:
    ids = set(l.rstrip('\r\n') for l in indexfile)

for line in sys.stdin:
    rec_id, _ = line.split('\t', 1)
    if rec_id in ids:
        sys.stdout.write(line)
