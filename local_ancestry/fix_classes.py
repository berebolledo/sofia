#! /bin/env python

import sys

newvalues = []
for line in sys.stdin:
    values = line.strip().split()
    for value in values:
        newvalues.append(str(value))
        newvalues.append(str(value))
newvalues.append('\n')
sys.stdout.write(' '.join(newvalues))
