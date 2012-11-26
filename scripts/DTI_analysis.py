# pmx  Copyright Notice
# ============================
#
# The pmx source code is copyrighted, but you can freely use and
# copy it as long as you don't change or remove any of the copyright
# notices.
#
# ----------------------------------------------------------------------
# pmx is Copyright (C) 2006-2013 by Daniel Seeliger
#
#                        All Rights Reserved
#
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name of Daniel Seeliger not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
#
# DANIEL SEELIGER DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
# SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS.  IN NO EVENT SHALL DANIEL SEELIGER BE LIABLE FOR ANY
# SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
# RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
# CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------
import sys, os
from glob import glob
from numpy import *


def read_data(fn):
    l = open(fn).readlines()
    data = []
    for line in l:
        if line[0] not in ['@','#']:
            entr = line.split()
            try:
                data.append( float(entr[1] ) )
            except:
                pass
    print >>sys.stderr, 'Read file:', fn, ' with %d data points' % len(data)
    return data

def datapoint_from_time(time):
    return time*500

def block_aver( data, block_size = 1000, offset = 0):
    # make blocks of 1 ns length
    total_time = len(data) / 500.
    next_time = block_size
    results = []
    while next_time < total_time:
        beg = datapoint_from_time(offset)
        end = datapoint_from_time(next_time)
        res = average( data[beg:end] )
        results.append( (offset+block_size*.5, res ) )
        offset = next_time
        next_time += block_size
    return results

def convergence( data, block_size = 1000, offset = 0):
    total_time = len(data) / 500.
    next_time = block_size
    results = []
    while next_time < total_time:
        beg = datapoint_from_time(offset)
        end = datapoint_from_time(next_time)
        res = average( data[beg:end] )
        results.append( (next_time, res ) )
        next_time += block_size
    return results
    


run_dirs = sys.argv[1:]
for d in run_dirs:
    os.chdir( d )
    print 'Processing -> %s' % d,
    data = read_data( 'dhdl.xvg' )
    fp = open('fe_sum.txt','w')
    print >>fp, average(data),  len(data)
    fp.close()
    res = block_aver( data, 100 )
    fp = open('fe_block_100.txt','w')
    for t, r in res:
        print >>fp,  t, r
    fp.close()
    res = block_aver( data, 200 )
    fp = open('fe_block_200.txt','w')
    for t, r in res:
        print >>fp, t, r
    fp.close()

    res = convergence( data, 100 )
    fp = open('fe_convergence.txt','w')
    for t, r in res:
        print >>fp, t, r
    fp.close()
    print '.............done'
    os.chdir('..')

