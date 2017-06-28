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

def block_aver( data, block_size = 1000, offset = 200):
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

def convergence( data, block_size = 1000, offset = 200):
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
    data = read_data( 'dgdl.xvg' )
    fp = open('fe_sum.txt','w')
    print >>fp, average(data),  len(data)
    fp.close()
    res = block_aver( data, 1000 )
    fp = open('fe_block_1000.txt','w')
    for t, r in res:
        print >>fp,  t, r
    fp.close()
    res = block_aver( data, 2000 )
    fp = open('fe_block_2000.txt','w')
    for t, r in res:
        print >>fp, t, r
    fp.close()

    res = convergence( data, 1000 )
    fp = open('fe_convergence.txt','w')
    for t, r in res:
        print >>fp, t, r
    fp.close()
    print '.............done'
    os.chdir('..')
