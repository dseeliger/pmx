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
from pmx import *

def read_data(fn, b = 0, e = -1):
    if e == -1: e = 9999999999
    l = open(fn).readlines()
    data = []
    for line in l:
        if line[0] not in ['@','#']:
            entr = line.split()
            try:
                time = float(entr[0])
                if time > b and time < e:
                    data.append( float(entr[1] ) )
            except:
                pass
#    print >>sys.stderr, 'Read file:', fn, ' with %d data points' % len(data)
    return data

def datapoint_from_time(time):
    return time*500

def block_aver( data, block_size = 1000):

    total_time = len(data) / 500.
    next_time = block_size
    results = []
    offset = 0
    while next_time < total_time:
        beg = datapoint_from_time(offset)
        end = datapoint_from_time(next_time)
        res = average( data[beg:end] )
        results.append( (str(offset)+'-'+str(next_time), res ) )
        offset = next_time
        next_time += block_size
    return results

def convergence( data, block_size = 1000):
    total_time = len(data) / 500.
    next_time = block_size
    results = []
    offset = 0
    while next_time < total_time:
        beg = datapoint_from_time(offset)
        end = datapoint_from_time(next_time)
        res = average( data[beg:end] )
        results.append( (next_time, res ) )
        next_time += block_size
    return results
    


help_text = ('Calculate delta G from multiple DTI runs',)

options = [
        Option( "-b", "real", 500, "Start time [ps]"),
        Option( "-e", "real", -1, "End time[ps]"),
        Option( "-block1", "int", 100, "Time[ps] for block average"),
        Option( "-block2", "int", 500, "Time[ps] for block average"),
#        Option( "-r2", "rvec", [1,2,3], "some vector that does wonderful things and returns always segfaults")
        ]

files = [
    FileOption("-dgdl", "r",["xvg"],"run", "Input file with dH/dl values"),
    FileOption("-o", "w",["txt"],"results.txt", "Results"),
    FileOption("-oc", "w",["txt"],"convergence.txt", "text file with mutations to insert"),
    FileOption("-ob", "w",["txt"],"block.txt", "files with block averages"),
    
]


cmdl = Commandline( sys.argv, options = options,
                    fileoptions = files,
                    program_desc = help_text,
                    check_for_existing_files = False )

dgdl_file = cmdl['-dgdl']
start_time = cmdl['-b']
end_time = cmdl['-e']

print 'DTI_analysis__> Reading: ', dgdl_file
print 'DTI_analysis__> Start time = ', start_time, ' End time = ', end_time
data = read_data( dgdl_file, b = start_time, e = end_time )
av = average(data)
st = std(data)
size = len(data)
print 'DTI_analysis__> <dH/dl> = %8.4f'% av, ' | #data points = ', size
fp = open(cmdl['-o'],'w')
print >>fp, av, st, size
fp.close()

block1 = cmdl['-block1']
fn =os.path.splitext(cmdl['-ob'])[0]+str(block1)+os.path.splitext(cmdl['-ob'])[1]
print 'DTI_analysis__> Block averaging 1: ',  block1
res = block_aver( data, block1 )
fp = open(fn,'w')
for a, b in res:
    print >>fp, a, b
fp.close()
block2 = cmdl['-block2']
fn =os.path.splitext(cmdl['-ob'])[0]+str(block2)+os.path.splitext(cmdl['-ob'])[1]
print 'DTI_analysis__> Block averaging 2: ',  block2
res = block_aver( data, block2 )
fp = open(fn,'w')
for a, b in res:
    print >>fp, a, b
fp.close()

res = convergence( data, 100 )
fp = open(cmdl['-oc'],'w')
for t, r in res:
    print >>fp, t, r
fp.close()

