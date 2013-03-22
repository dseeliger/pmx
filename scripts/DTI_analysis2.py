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
#import random
from pylab import *
from pmx import *
from pmx.odict import OrderedDict

def simple_xy_plot(x,y,xlab, ylab, name, result, err, ylimit = False):
    figure(figsize=(8,8))
    plot(x,y,'rd-', lw=2)
    xlabel(xlab, fontsize=20)
    ylabel(ylab, fontsize=20)
    title(r"$\int_0^1\langle\frac{dH}{d\lambda}\rangle_{\lambda}d\lambda$ = %.2f +/- %.2f kJ/mol" % (result, err))
    grid(lw=2)
    xx = gca()
    if ylimit:
        ylim( ylimit)
    for x in xx.spines.values():
        x.set_lw(2)
    savefig(name, dpi=600)


def plot_random_blocks( outf, run_dic, block_file,avx, avy, ylimit = False):
    figure(figsize=(8,5))
    xlab = r'$\lambda$'
    ylab = r"dH/d$\lambda$ [kJ/mol]"
    block_res = []
    std_lda = []
    for lda, p in run_dic.items():
        r = read_convergence_file( os.path.join(p,'convergence.txt') )
        std_lda.append( std(map( lambda a: a[1], r)) )
        blck = read_blocks(p, block_file)
        block_res.append( (lda, blck ))
    block_res.sort(lambda a, b: cmp(a[0],b[0]) )
    xx = []
    for i in range(1000):
        dG, x, y = random_value_from_block( block_res, with_data = True )
        plot(x,y,'k-', lw=.5, alpha = .01)
    errorbar(avx,avy,yerr=std_lda, fmt ='kd', markersize=5)
    plot(avx,avy,'r-', lw=2)
    if ylimit:
        ylim( ylimit)
    xlabel(xlab, fontsize=22)
    ylabel(ylab, fontsize=22)
    xticks(fontsize=20)
    yticks(fontsize=20)
    subplots_adjust(bottom=0.25, left=0.15)
#    title("Curves from random blocks")
    grid(lw=2)
    xx = gca()
    for x in xx.spines.values():
        x.set_lw(2)
    savefig(outf, dpi=600)
    

def check( lst ):
    
    files = ['convergence.txt',
             'results.txt',
             'block100.txt',
             'block500.txt']
    missing = []
    for d in lst:
        for f in files:
            p = os.path.join(d, f)
            if not os.path.isfile( p ):
                if d not in missing:
                    missing.append(d)
#                    print 'file %s missing in directory %s -> Exiting' % (f,d)
#                missing = True
    if missing:
        for d in missing:
            print 'Run %s KIA or MIA!!' % (d)
        sys.exit(1)
        
def read_result(p):
    f = os.path.join(p,'results.txt')
    r = float(open(f).read().split()[0])
    return r

def read_blocks(p, block_file):
    f = os.path.join( p, block_file)
    l = open(f).readlines()
    r = []
    for line in l:
        r.append( float(line.split()[1] ))
    return r

def do_dgdl( results ):
    x = map(lambda a: a[0], results)
    y = map(lambda a: a[1], results)
    dG =  trapz(y,x)
    return dG, x, y

def random_value_from_block( block_res, with_data = False ):
    data = []
    n_blocks = len(block_res)
    for i in range(n_blocks):
        size = len(block_res[i][1])
        ri = randint(0,size-1)
        data.append( (block_res[i][0], block_res[i][1][ri]) )
    dg, x, y = do_dgdl( data )
    if with_data:
        return dg, x, y
    else:
        return dg

def error_from_block_aver(run_dic, block_file):
    block_res = []
    for lda, p in run_dic.items():
        blck = read_blocks(p, block_file)
        block_res.append( (lda, blck ))
    block_res.sort(lambda a, b: cmp(a[0],b[0]) )
    xx = []
    for i in range(1000):
        dG = random_value_from_block( block_res )
        xx.append( random_value_from_block( block_res ) )
    return std(xx)

def read_dGdl( run_dic) :
    results_all = []
    for lda, p in run_dic.items():
        r = read_result(p)
        print 'DTI_analysis2__> lambda = %4.3f'% lda, 'dH/dl = %8.4f' %  r
        results_all.append( (lda, r ))
        
    results_all.sort(lambda a, b: cmp(a[0],b[0]) )
    return results_all

def read_convergence_file( f ):
    l = open(f).readlines()
    r = []
    for line in l:
        entr = line.split()
        r.append( (float(entr[0]), float(entr[1]) ) )
    return r
                
def dG_over_time( run_dic ):
    fp = open('dG_vs_time.txt','w')
    dG_vs_time = []
    for lda, p in run_dic.items():
        r = read_convergence_file( os.path.join(p,'convergence.txt') )
        dG_vs_time.append( (lda, r) )
    dG_vs_time.sort(lambda a, b: cmp(a[0],b[0]) )
    min_size = 100
    for lda, lst in dG_vs_time:
        if len(lst) < min_size:
            min_size = len(lst)
    for i in range(min_size):
        time = dG_vs_time[0][1][i][0]
        lda_vals = map(lambda a: a[0], dG_vs_time)
        dgdl_vals = map(lambda a: a[1][i][1], dG_vs_time)
        dG = trapz( dgdl_vals, lda_vals)
        print >>fp, time, dG 

def get_max_number_of_blocks(run_dic, block_file):
    block_res = []
    for lda, p in run_dic.items():
        blck = read_blocks(p, block_file)
        block_res.append( (lda, blck ))
    block_res.sort(lambda a, b: cmp(a[0],b[0]) )
    min_size = min( map(lambda a: len(a[1]), block_res ) )
    return min_size

def dG_from_last_blocks( run_dic, block_file, nblocks ):
    block_res = []
    for lda, p in run_dic.items():
        blck = read_blocks(p, block_file)
        block_res.append( (lda, blck ))
    block_res.sort(lambda a, b: cmp(a[0],b[0]) )
    min_size = min( map(lambda a: len(a[1]), block_res ) )
    first_block = min_size - nblocks
#    first_block = 15
    res = []
    for b in block_res:
        lda = b[0]
        dGdl = average( b[1][first_block:] )
        res.append( (lda, dGdl) )
    dG = trapz( map(lambda a:a[1], res), map(lambda a:a[0], res) )
    return dG, nblocks
                        
###########################################################################################################


help_text = ('Calculate delta G from multiple DTI runs. Input files are assumed to have default names',)

options = [
#        Option( "-b", "real", 0, "Start time [ps]"),
#        Option( "-e", "real", -1, "End time[ps]"),
#        Option( "-b", "bool", True, "bool"),
#        Option( "-r2", "rvec", [1,2,3], "some vector that does wonderful things and returns always segfaults")
        ]

files = [
    FileOption("-d", "r/m",["dir"],"run", "Run directories of DTI runs"),
    FileOption("-o", "w",["txt"],"dti_result.txt", "Result with dG"),
    FileOption("-dGdl", "w",["png"],"lda_vs_dGdl.png", "plot with dG/dl values"),
    FileOption("-rblocks", "w",["png"],"random_blocks.png", "plot with dG/dl values and error estimation"),
]


cmdl = Commandline( sys.argv, options = options,
                    fileoptions = files,
                    program_desc = help_text,
                    check_for_existing_files = False )


run_dirs = cmdl['-d']



#exit()

    

#run_name = sys.argv[1]

#lst = glob('%s_*.*' % sys.argv[1])
#print lst
check(run_dirs)

run_dic = OrderedDict()
tmp_list = []
for d in run_dirs:
    tmp_list.append(d)
tmp_list = sorted(tmp_list, key=lambda f: float(f.split('_')[-1]) )
for d in tmp_list:
    lda = float(d.split('_')[-1])
    run_dic[lda] = d
    
results_all = read_dGdl( run_dic )
try:
    err = error_from_block_aver( run_dic, 'block100.txt')
    err2 = error_from_block_aver( run_dic, 'block500.txt')
except:
    err = 0
    err2 = 0
    
dgdl_all, x, y = do_dgdl(results_all)

outf = 'lda_vs_dGdl.txt'
fp = open(outf, 'w')
for i in range(len(x)):
    print >>fp, x[i], y[i]

outf = cmdl['-dGdl']
#try:
#    yl = ( float(sys.argv[1]), float(sys.argv[2]) )
#except:
#    yl = False
yl = False           
simple_xy_plot(x,y,r'$\lambda$', r"dG/d$\lambda$",outf, dgdl_all, err)#, yl)
if err != 0:
    outf1 = os.path.splitext(cmdl['-rblocks'])[0]+str(100)+'.png'
    outf2 = os.path.splitext(cmdl['-rblocks'])[0]+str(500)+'.png'
    print 'DTI_analysis2__> Plotting random block figures'
    plot_random_blocks(outf1, run_dic, 'block100.txt', x, y, yl)
    plot_random_blocks(outf2, run_dic, 'block500.txt', x, y, yl)

final_result = round(dgdl_all,2)
err100 =  round(err,2)
err500 =  round(err2,2)
print 'DTI_analysis2__> dG = %8.2f' % final_result, '+/- %8.4f' % err500
fp = open(cmdl['-o'],'w')
print >>fp, 'Result dG = ', round(dgdl_all,2), 'kJ/mol', 'Err(100) = ', round(err,2), 'kJ/mol', 'Err(500) = ', round(err2,2), 'kJ/mol'
dG_over_time( run_dic )

max_block = get_max_number_of_blocks( run_dic, 'block100.txt')
for i in range(1, max_block+1):
    dG, ndata = dG_from_last_blocks( run_dic, 'block100.txt', i)
    print >>fp, 'dG (block) = ', round(dG,2), 'kJ/mol', 'from last ', ndata, 'blocks'

