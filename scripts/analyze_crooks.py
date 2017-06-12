#!/usr/bin/env python
# pmx  Copyright Notice
# ============================
#
# The pmx source code is copyrighted, but you can freely use and
# copy it as long as you don't change or remove any of the copyright
# notices.
#
# ----------------------------------------------------------------------
# pmx is Copyright (C) 2006-2011 by Daniel Seeliger
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

import sys
import os
import time
from copy import deepcopy
from pmx import *  # TODO: use explicit imports: hard to find out what modules are being used
from pmx.parser import *
from pylab import *
from scipy.integrate import simps
from scipy.optimize import fmin
from scipy.special import erf
from random import gauss, randint, choice

# TODO:
# - consolidate functions: e.g. BAR* and Jarz* functions
# - define internal-use functions
# - add docstrings
# - cleanup main:
#     - move messages to stderr (tee) inside the different functions?
#     - make more readable by using less and more self-contained funcs, e.g.:
#           data = read_dgdl()
#           if args.do_jarz:
#               dg_jarz = jarz(data)
#           etc etc
#
#  The above will make it easier to use this both as command line tool as well
#  as within other scripts. E.g. one could get a list of the dgdl.xvg files
#  then estimate the free energy along with its error just like BAR(dgdl_list)
#  and do whatever they want with the numbers, or something like that
#
# - use logger for info and debug

debug = True  # unused: remove?

# move this inside plot functions?
params = {'legend.fontsize': 12}
rcParams.update(params)

# Constants
kb = 0.00831447215


def tee(fp, s):
    print >>fp, s
    print s


def cgi_error_from_mean(nruns, mu1, sig1, n1, mu2, sig2, n2):
    iseq = []

    for k in range(nruns):
        g1 = []
        g2 = []
        for i in range(n1):
            g1.append(gauss(mu1, sig1))
        for i in range(n2):
            g2.append(gauss(mu2, sig2))
        m1 = average(g1)
        s1 = std(g1)
        m2 = average(g2)
        s2 = std(g2)
        p1 = 1./(s1*sqrt(2*pi))  # unused variable: remove?
        p2 = 1./(s2*sqrt(2*pi))  # unused variable: remove?
        iq = (m1+m2)/2.
        iseq.append(iq)
    mean = average(iseq)  # unused variable: remove?
    err = std(iseq)
    return err


def cgi_error(nruns, mu1, sig1, n1, mu2, sig2, n2):
    iseq = []
    for k in range(nruns):
        g1 = []
        g2 = []
        for i in range(n1):
            g1.append(gauss(mu1, sig1))
        for i in range(n2):
            g2.append(gauss(mu2, sig2))
        m1 = average(g1)
        s1 = std(g1)
        m2 = average(g2)
        s2 = std(g2)
        p1 = 1./(s1*sqrt(2*pi))
        p2 = 1./(s2*sqrt(2*pi))
        iq = gauss_intersection([p1, m1, s1], [p2, m2, s2])
        iseq.append(iq)
    mean = average(iseq)  # unused variable: remove?
    err = std(iseq)
    return err


def sort_file_list(lst):

    # we assume that the directory is numbered
    # guess directory base name first
    dir_name = lst[0].split('/')[-2]
    base_name = ''
    for i, x in enumerate(dir_name):
        if x.isdigit():
            check = True
            for k in range(i, len(dir_name)):
                if not dir_name[k].isdigit():
                    check = False
            if check:
                base_name = dir_name[:i]
                break
    if base_name:
        def get_num(s): return int(s.split('/')[-2].split(base_name)[1])
        lst.sort(lambda a, b: cmp(get_num(a), get_num(b)))
        return lst
    else:
        return lst


def process_dgdl(fn, ndata=-1, lambda0=0):
    sys.stdout.write('\r------>  %s' % fn)
    sys.stdout.flush()
    l = open(fn).readlines()
    if not l:
        return None, None
    r = []
    for line in l:
        if line[0] not in '#@&':
            try:
                r.append([float(x) for x in line.split()])
            except:
                print ' !! Skipping %s ' % (fn)
                return None, None

    if ndata != -1 and len(r) != ndata:
        try:
            print(' !! Skipping %s ( read %d data points, should be %d )'
                  % (fn, len(r), ndata))
        except:
            print ' !! Skipping %s ' % (fn)
        return None, None
    # convert time to lambda
    ndata = len(r)
    dlambda = 1./float(ndata)
    if lambda0 == 1:
        dlambda *= -1
    data = []

    for i, (ttime, dgdl) in enumerate(r):
        data.append([lambda0+i*dlambda, dgdl])
    x = map(lambda a: a[0], data)
    y = map(lambda a: a[1], data)

    if lambda0 == 1:
        x.reverse()
        y.reverse()
    return simps(y, x), ndata


def check_first_dgdl(fn, lambda0):
    l = open(fn).readlines()
    if not l:
        return None
    r = []
    for line in l:
        if line[0] not in '#@&':
            r.append([float(x) for x in line.split()])
    ndata = len(r)
    dlambda = 1./float(ndata)
    if lambda0 == 1:
        dlambda *= -1
    print '---------------------------------------------'
    print '\t\t Checking simulation data.....'
    print '\t\t File: %s' % fn
    print '\t\t # data points: %d' % ndata
    print '\t\t Length of trajectory: %8.3f ps' % r[-1][0]
    print '\t\t Delta lambda: %8.5f' % dlambda
    print '---------------------------------------------'


def work_from_crooks(lst, lambda0, reverse=False):
    print '\nProcessing simulation data......'
    output_data = []
    check_first_dgdl(lst[0], lambda0)
    first_res, ndata = process_dgdl(lst[0], lambda0=lambda0)
    output_data.append([lst[0], first_res])
    results = [first_res]
    for f in lst[1:]:
        res, tmp = process_dgdl(f, ndata=ndata, lambda0=lambda0)
        if res is not None:
            results.append(res)
            output_data.append([f, res])

    if reverse is True:
        results = [x*(-1) for x in results]
        output_data = [[x, y*(-1)] for x, y in output_data]

    return results, output_data


def data_to_gauss(data):
    m = average(data)
    dev = std(data)
    A = 1./(dev*sqrt(2*pi))
    return m, dev, A


def gauss_intersection(g1, g2):
    A1, m1, s1 = g1
    A2, m2, s2 = g2
    p1 = m1/s1**2-m2/s2**2
    p2 = sqrt(1/(s1**2*s2**2)*(m1-m2)**2+2*(1/s1**2-1/s2**2)*log(s2/s1))
    p3 = 1/s1**2-1/s2**2
    x1 = (p1+p2)/p3
    x2 = (p1-p2)/p3
    # determine which solution to take
    if x1 > m1 and x1 < m2 or x1 > m2 and x1 < m1:
        return x1
    elif x2 > m1 and x2 < m2 or x2 > m2 and x2 < m1:
        return x2
    else:
        return False  # we do not take the intersection


def ksref():
    f = 1
    potent = 10000
    lamb = arange(0.25, 2.5, 0.001)
    q = array(zeros(len(lamb), float))
    res = []
    for k in range(-potent, potent):
        q = q + f*exp(-2.0*(k**2)*(lamb**2))
        f = -f
    for i in range(len(lamb)):
        res.append((lamb[i], q[i]))
    return res


def ksfunc(lamb):
    f = 1
    potent = 10000
    q = 0
    for k in range(-potent, potent):
        q = q + f*exp(-2.0*(k**2)*(lamb**2))
        f *= -1
    return q


def ks(data, alpha=.05, refks=None):
    N = len(data)
    nd, ed = edf(data)
    cd = cdf(data)
    siglev = 1-alpha
    dval = []
    for i, val in enumerate(ed):
        d = abs(val-cd[i])
        dval.append(d)
        if i:
            d = abs(ed[i-1]-cd[i])
            dval.append(d)
    dmax = max(dval)
    check = math.sqrt(N)*dmax
    if not refks:
        refks = ksref()
    lst = filter(lambda x: x[1] > siglev, refks)
    lam0 = lst[0][0]
    if check >= lam0:
        bOk = False
    else:
        bOk = True

    q = ksfunc(check)
    return (1-q), lam0, check, bOk


def edf(dg_data):
    edf_ = []
    ndata = []
    data = deepcopy(dg_data)
    data.sort()
    N = float(len(data))
    cnt = 0
    for item in data:
        cnt += 1
        edf_.append(cnt/N)
        ndata.append(item)
    ndata = array(ndata)
    edf_ = array(edf_)
    return ndata, edf_


def cdf(dg_data):
    data = deepcopy(dg_data)
    data.sort()
    mean = average(data)
    sig = std(data)
    cdf = 0.5*(1+erf((data-mean)/float(sig*sqrt(2))))
    return cdf


def data_from_file(fn):
    data = read_and_format(fn, 'sf')
    return map(lambda a: a[1], data)


def dump_integ_file(fn, data):
    fp = open(fn, 'w')
    for fn, w in data:
        print >>fp, fn, w
    fp.close()


# =================
# Estimator Classes
# =================
# TODO: docstrings

class Jarz(object):
    '''Jarzynski estimator.

    Description...

    Parameters
    ----------
    wf : array_like
        array of forward work values.
    wr : array_like
        array of reverse work values.
    T : float or int
    nboots : int
        number of bootstrap samples to use for error estimation.

    Examples
    --------


    Attributes
    ----------

    '''

    def __init__(self, wf, wr, T, nboots=0):
        self.wf = np.array(wf)
        self.wr = np.array(wr)
        self.T = float(T)
        self.nboots = nboots

        # Calculate all Jarz properties available
        self.dg_for = self.calc_dg(w=self.wf, c=1.0, T=self.T)
        self.dg_rev = -1.0 * self.calc_dg(w=self.wr, c=-1.0, T=self.T)
        self.dg_mean = (self.dg_for + self.dg_rev) * 0.5

        if nboots > 0:
            self.err_boot_for = self.calc_err_boot(w=self.wf, T=self.T,
                                                   c=1.0, nboots=nboots)
            self.err_boot_rev = self.calc_err_boot(w=self.wr, T=self.T,
                                                   c=-1.0, nboots=nboots)

    @staticmethod
    def calc_dg(w, T, c):
        beta = 1./(kb*T)
        n = float(len(w))
        mexp = 0.0
        m = 0.0
        m2 = 0.0
        for i in w:
            mexp = mexp + exp(-beta*c*i)
            m = m + c*i
            m2 = m2 + i*i
        mexp = mexp/n
        m = m/n
        m2 = m2/n
        var = (m2-m*m)*(n/(n-1))
        # Jarzynski estimator
        dg = -kb*T*log(mexp)
        # Fluctuation-Dissipation estimator
        # FIXME: unused atm, remove or return?
        dg2 = m - beta*var/2.0

        return dg

    @staticmethod
    def calc_err_boot(w, T, c, nboots):
        out = []
        n = int(len(w))
        for k in range(nboots):
            sys.stdout.write('\rJarzynski error bootstrap: iteration %s/%s'
                             % (k+1, nboots))
            sys.stdout.flush()
            for i in range(n):
                val = [choice(w) for _ in xrange(n)]
            foo = -1.0 * Jarz.calc_dg(val, T, c)
            out.append(foo)
        sys.stdout.write('\n')
        err = std(out)

        return err


class BAR(object):
    '''Bennett acceptance ratio (BAR).

    Description...

    Parameters
    ----------

    Examples
    --------
    '''

    def __init__(self, wf, wr, T, nboots=0):
        self.wf = array(wf)  # is array() needed? they are both lists
        self.wr = array(wr)
        self.T = float(T)
        self.nboots = nboots

        self.nf = len(wf)
        self.nr = len(wr)
        self.beta = 1./(kb*self.T)
        self.M = kb * self.T * log(float(self.nf) / float(self.nr))

        # Calculate all BAR properties available
        self.dg = self.calc_dg(self.wf, self.wr, self.T)
        self.err = self.calc_err(self.dg, self.wf, self.wr, self.T)
        if nboots > 0:
            self.err_boot = self.calc_err_boot(self.wf, self.wr, nboots,
                                               self.T)
        self.conv = self.calc_conv(self.dg, self.wf, self.wr, self.T)
        if nboots > 0:
            self.conv_err_boot = self.calc_conv_err_boot(self.dg, self.wf,
                                                         self.wr, nboots,
                                                         self.T)

    @staticmethod
    def calc_dg(wf, wr, T):
        '''Estimates and returns the free energy difference.

        Parameters
        ----------

        Returns
        ----------
        '''

        nf = float(len(wf))
        nr = float(len(wr))
        beta = 1./(kb*T)
        M = kb * T * log(nf/nr)

        def func(x, wf, wr):
            sf = 0
            for v in wf:
                sf += 1./(1+exp(beta*(M+v - x)))
                sr = 0
                for v in wr:
                    sr += 1./(1+exp(-beta*(M+v - x)))
            r = sf-sr
            return r**2

        avA = average(wf)
        avB = average(wr)
        x0 = (avA+avB)/2.
        dg = fmin(func, x0=x0, args=(wf, wr), disp=0)

        return float(dg)

    @staticmethod
    def calc_err(dg, wf, wr, T):
        '''Calculates the analytical error estimate.'''

        nf = float(len(wf))
        nr = float(len(wr))
        beta = 1./(kb*T)
        M = kb * T * log(nf/nr)

        err = 0
        for v in wf:
            err += 1./(2+2*cosh(beta * (M+v-dg)))
        for v in wr:
            err += 1./(2+2*cosh(beta * (M+v-dg)))
        N = nf + nr
        err /= float(N)
        tot = 1/(beta**2*N)*(1./err-(N/nf + N/nr))

        err = float(sqrt(tot))
        return err

    @staticmethod
    def calc_err_boot(wf, wr, nruns, T):
        '''Calculates the error by bootstrapping.'''

        nf = len(wf)
        nr = len(wr)
        res = []
        for k in range(nruns):
            sys.stdout.write('\rBAR error bootstrap: iteration %s/%s'
                             % (k+1, nruns))
            sys.stdout.flush()
            for i in range(nf):
                valA = [choice(wf) for _ in xrange(nf)]
            for i in range(nr):
                valB = [choice(wf) for _ in xrange(nr)]
            foo = BAR.calc_dg(valA, valB, T)
            res.append(foo)
        sys.stdout.write('\n')
        err_boot = std(res)

        return err_boot

    @staticmethod
    def calc_conv(dg, wf, wr, T):
        '''Evaluates BAR convergence'''

        wf = np.array(wf)
        wr = np.array(wr)

        beta = 1./(kb*T)
        nf = len(wf)
        nr = len(wr)
        N = float(nf + nr)

        ratio_alpha = float(nf)/N
        ratio_beta = float(nr)/N
        bf = 1.0/(ratio_beta + ratio_alpha * np.exp(beta*(wf-dg)))  # where is numpy coming from? probably from pmx import *
        tf = 1.0/(ratio_alpha + ratio_beta * np.exp(beta*(-wr+dg)))
        Ua = (np.mean(tf) + np.mean(bf))/2.0
        Ua2 = (ratio_alpha * np.mean(np.power(tf, 2)) +
               ratio_beta * np.mean(np.power(bf, 2)))
        conv = (Ua-Ua2)/Ua
        return conv

    @staticmethod
    def calc_conv_err_boot(dg, wf, wr, nruns, T):
        nf = len(wf)
        nr = len(wr)
        res = []
        for k in range(nruns):
            sys.stdout.write('\rConvergence error bootstrap: '
                             'iteration %s/%s' % (k+1, nruns))
            sys.stdout.flush()
            for i in range(nf):
                valA = [choice(wf) for _ in xrange(nf)]
            for i in range(nr):
                # FIXME: is "choice(wf)" correct? seems like it should be "choice(wr)"
                valB = [choice(wf) for _ in xrange(nr)]
            foo = BAR.calc_conv(dg, valA, valB, T)
            res.append(foo)
        sys.stdout.write('\n')
        err = std(res)
        return(err)


def gauss_func(A, mean, dev, x):
    x = array(x)
    y = A*exp(-(((x-mean)**2)/(2.0*(dev**2))))
    return y


def make_plot(fname, data1, data2, result, err, nbins, dpi):
    figure(figsize=(8, 6))
    mf, devf, Af = data_to_gauss(data1)
    mb, devb, Ab = data_to_gauss(data2)

    maxi = max(data1+data2)
    mini = min(data1+data2)
    n1, bins1, patches1 = hist(data1, range=(mini, maxi), bins=nbins,
                               facecolor='blue', alpha=0.75, normed=True,
                               label='0->1')
    n2, bins2, patches2 = hist(data2, range=(mini, maxi), bins=nbins,
                               facecolor='red', alpha=0.75, normed=True,
                               label='1->0')
    xlabel('W [kJ/mol]', fontsize=20)
    ylabel('Probability', fontsize=20)
    title(r'Work Distribution $\lambda$ 0->1 (blue) $\lambda$ 1->0 (red)')
    grid(lw=2)
    loc, lab = yticks()
    ll = []
    for i in range(len(lab)):
        ll.append("")
    yticks(loc, ll)
    x = arange(mini, maxi, .5)
    y1 = gauss_func(Af, mf, devf, x)
    y2 = gauss_func(Ab, mb, devb, x)

    plot(x, y1, 'b--', linewidth=2)
    plot(x, y2, 'r--', linewidth=2)

    size = max([max(y1), max(y2)])
    res_x = [result, result]
    res_y = [0, size*1.2]
    plot(res_x, res_y, 'k--', linewidth=2,
         label=r'$\Delta$G = %.2f $\pm$ %.2f kJ/mol' % (result, err))
    legend(shadow=True, fancybox=True)
    ylim(0, size*1.2)
    xl = gca()
    for val in xl.spines.values():
        val.set_lw(2)
    savefig(fname, dpi=dpi)


def make_W_over_time_plot(fname, data1, data2, result, err, nbins, dpi):
    figure(figsize=(8, 6))
    x1 = range(len(data1))
    x2 = range(len(data2))
    if x1 > x2:
        x = x1
    else:
        x = x2
    mf, devf, Af = data_to_gauss(data1)
    mb, devb, Ab = data_to_gauss(data2)

    maxi = max(data1+data2)
    mini = min(data1+data2)

    sm1 = smooth(array(data1))
    sm2 = smooth(array(data2))
    subplot(1, 2, 1)
    plot(x1, data1, 'g-', linewidth=2, label="Forward (0->1)", alpha=.3)
    plot(x1, sm1, 'g-', linewidth=3)
    plot(x2, data2, 'b-', linewidth=2, label="Backward (1->0)", alpha=.3)
    plot(x2, sm2, 'b-', linewidth=3)
    legend(shadow=True, fancybox=True, loc='upper center')
    ylabel(r'W [kJ/mol]', fontsize=20)
    xlabel(r'# Snapshot', fontsize=20)
    grid(lw=2)
    xlim(0, x[-1]+1)
    xl = gca()
    for val in xl.spines.values():
        val.set_lw(2)
    subplot(1, 2, 2)
    hist(data1, bins=nbins, orientation='horizontal', facecolor='green',
         alpha=.75, normed=True)
    hist(data2, bins=nbins, orientation='horizontal', facecolor='blue',
         alpha=.75, normed=True)

    x = arange(mini, maxi, .5)

    y1 = gauss_func(Af, mf, devf, x)
    y2 = gauss_func(Ab, mb, devb, x)

    plot(y1, x, 'g--', linewidth=2)
    plot(y2, x, 'b--', linewidth=2)
    size = max([max(y1), max(y2)])
    res_x = [result, result]
    res_y = [0, size*1.2]
    plot(res_y, res_x, 'k--', linewidth=2,
         label=r'$\Delta$G = %.2f $\pm$ %.2f kJ/mol' % (result, err))
    legend(shadow=True, fancybox=True, loc='upper center')
    xticks([])
    yticks([])
    xl = gca()
    for val in xl.spines.values():
        val.set_lw(2)
    subplots_adjust(wspace=0.0, hspace=0.1)
    savefig(fname, dpi=dpi)


def smooth(x, window_len=11, window='hanning'):

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")
    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
    if window_len < 3:
        return x
    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', "
                         "'bartlett', 'blackman'")
    s = r_[2*x[0]-x[window_len:1:-1], x, 2*x[-1]-x[-1:-window_len:-1]]
    # moving average
    if window == 'flat':
        w = ones(window_len, 'd')
    else:
        w = eval(window+'(window_len)')
    y = convolve(w/w.sum(), s, mode='same')
    return y[window_len-1:-window_len+1]


def select_random_subset(lst, n):
    ret = []
    idx = []
    while len(ret) < n:
        rn = randint(0, len(lst)-1)
        if rn not in idx:
            idx.append(rn)
            ret.append(lst[rn])
    print idx
    return ret


def parse_options(argv):

    # should we keep versioning just for the whole pmx?
    version = "1.2"

    options = [
        Option("-nbins", "int", 10, "number of histograms bins for plot"),
        Option("-T", "real", 298, "Temperature for BAR calculation"),
        Option("-dpi", "int", 300, "plot resolution"),
        Option("-reverseB", "bool", False, "reverse state B"),
        Option("-firstA", "int", 0,
               "first trajectory to analyze"
               " (by default all values are taken)"),
        Option("-lastA", "int", 100,
               "last trajectory to analyze"
               " (by default all values are taken)"),
        Option("-firstB", "int", 0,
               "first trajectory to analyze"
               " (by default all values are taken)"),
        Option("-lastB", "int", 100,
               "last trajectory to analyze (by default all values are taken)"),
        Option("-rand", "int", 50, "take a random subset of trajectories"),
        Option("-integ_only", "bool", False,
               "Do integration only. Skip analysis."),
        Option("-KS", "bool", True, "Do Kolmogorov-Smirnov test"),
        Option("-jarz", "bool", False, "Jarzynski estimation"),
        Option("-plot", "bool", True, "Plot work histograms"),
        Option("-nruns", "int", 100,
               "number of runs for bootstrapped BAR error"),
        ]

    files = [
        FileOption("-pa", "r/m", ["xvg"], "dgdl.xvg", "paths to 0->1 runs"),
        FileOption("-pb", "r/m", ["xvg"], "dgdl.xvg", "paths to 1->0 runs"),
        FileOption("-o", "w", ["dat"], "results.dat", "results"),
        FileOption("-cgi_plot", "w", ["png", "eps", "svg", "pdf"],
                   "cgi.png", "plot work histograms "),
        FileOption("-W_over_t", "w", ["png", "eps", "svg", "pdf"],
                   "W_over_t.png", "plot work over time "),
        FileOption("-i0", "r/m/o", ["dat"], "integ0.dat",
                   "read integrated W (0->1)"),
        FileOption("-i1", "r/m/o", ["dat"], "integ1.dat",
                   "read integrated W (1->0)"),
        FileOption("-o0", "w", ["dat"], "integ0.dat",
                   "write integrated W (0->1)"),
        FileOption("-o1", "w", ["dat"], "integ1.dat",
                   "write integrated W (1->0)"),
        ]

    help_text = ('Calculates free energies from fast growth  ',
                 'thermodynamic integration runs.',
                 'First method: Crooks-Gaussian Intersection (CGI)',
                 'Second method: Benett Acceptance Ratio (BAR)'
                 )

    cmdl = Commandline(argv, options=options, fileoptions=files,
                       program_desc=help_text, check_for_existing_files=False,
                       version=version)
    cmdl.argv = argv
    return cmdl


def main(cmdl):

    out = open(cmdl['-o'], 'w')
    print >>out, "# analyze_crooks.py, version = %s" % cmdl.version
    print >>out, "# pwd = %s" % os.getcwd()
    print >>out, "# %s (%s)" % (time.asctime(), os.environ.get('USER'))
    print >>out, "# command = %s" % ' '.join(cmdl.argv)
    print >>out, "#------------------------------------------------"

    if cmdl['-reverseB']:
        reverseB = True
    else:
        reverseB = False

    if not cmdl.opt['-i0'].is_set:
        run_ab = cmdl['-pa']
        run_ba = cmdl['-pb']
        run_ab = sort_file_list(run_ab)
        run_ba = sort_file_list(run_ba)
        res_ab, ab_data = work_from_crooks(run_ab, lambda0=0, reverse=False)
        res_ba, ba_data = work_from_crooks(run_ba, lambda0=1, reverse=reverseB)
        dump_integ_file(cmdl['-o0'], ab_data)
        dump_integ_file(cmdl['-o1'], ba_data)
    else:
        res_ab = []
        res_ba = []
        for fn in cmdl['-i0']:
            print '\t\tReading integrated values (0->1) from', fn
            res_ab.extend(data_from_file(fn))
        for fn in cmdl['-i1']:
            print '\t\tReading integrated values (1->0) from', fn
            res_ba.extend(data_from_file(fn))

    if cmdl['-integ_only']:
        print '\n    Integration done. Skipping analysis.'
        print '\n    ......done........\n'
        sys.exit(0)

    firstA = 0
    lastA = len(res_ab)
    firstB = 0
    lastB = len(res_ba)
    if cmdl.opt['-firstA'].is_set:
        firstA = cmdl['-firstA']
        tee(out, '   first trajectory to read from A: %d' % firstA)
    if cmdl.opt['-lastA'].is_set:
        lastA = cmdl['-lastA']
        tee(out, '   last trajectory to read from A : %d' % lastA)
    if cmdl.opt['-firstB'].is_set:
        firstB = cmdl['-firstB']
        tee(out, '   first trajectory to read from B: %d' % firstB)
    if cmdl.opt['-lastB'].is_set:
        lastB = cmdl['-lastB']
        tee(out, '   last trajectory to read from B : %d' % lastB)

    res_ab = res_ab[firstA:lastA]
    res_ba = res_ba[firstB:lastB]

    if cmdl.opt['-rand'].is_set:
        ntraj = cmdl['-rand']
        tee(out, ' select random subset of trajectories: %d' % ntraj)
        res_ab = select_random_subset(res_ab, ntraj)
        res_ba = select_random_subset(res_ba, ntraj)

    mf, devf, Af = data_to_gauss(res_ab)
    mb, devb, Ab = data_to_gauss(res_ba)
    tee(out, ' --------------------------------------------------------')
    tee(out, '             ANALYSIS: NUMBER OF TRAJECTORIES:')
    tee(out, '               0->1 : %d' % len(res_ab))
    tee(out, '               1->0 : %d' % len(res_ba))

    tee(out, ' --------------------------------------------------------')
    tee(out, '             ANALYSIS: Crooks-Gaussian Intersection     ')
    tee(out, ' --------------------------------------------------------')
    tee(out, '  Forward  : mean = %8.3f  std = %8.3f' % (mf, devf))
    tee(out, '  Backward : mean = %8.3f  std = %8.3f' % (mb, devb))

    if cmdl['-KS']:
        tee(out, '  Running KS-test ....')
        q0, lam00, check0, bOk0 = ks(res_ab)
        q1, lam01, check1, bOk1 = ks(res_ba)

        tee(out, '  Forward  : gaussian quality = %3.2f' % q0)
        if bOk0:
            tee(out, '             ---> KS-Test Ok')
        else:
            tee(out, '             ---> KS-Test Failed. sqrt(N)*Dmax = %4.2f, lambda0 = %4.2f' % (q0, check0))
        tee(out, '  Backward : gaussian quality = %3.2f' % q1)
        if bOk1:
            tee(out, '             ---> KS-Test Ok')
        else:
            tee(out, '             ---> KS-Test Failed. sqrt(N)*Dmax = %4.2f, lambda0 = %4.2f' % (q1, check1))

    tee(out, '  Calculating Intersection...')
    cgi_result = gauss_intersection([Af, mf, devf], [Ab, mb, devb])
    intersection = True
    if not cgi_result:
        tee(out, '\n  Gaussians too close for intersection calculation')
        tee(out, '   --> Taking difference of mean values')
        cgi_result = (mf+mb)*.5
        intersection = False
    tee(out, '  RESULT: dG ( CGI )  = %8.4f kJ/mol' % cgi_result)
    if intersection:
        cgi_err = cgi_error(1000, mf, devf, len(res_ab), mb, devb, len(res_ba))
    else:
        cgi_err = cgi_error_from_mean(1000, mf, devf, len(res_ab),
                                      mb, devb, len(res_ba))
    tee(out, '  RESULT: error_dG ( CGI ) = %8.4f kJ/mol' % cgi_err)
    tee(out, ' --------------------------------------------------------')
    tee(out, '             ANALYSIS: Bennett Acceptance Ratio     ')
    tee(out, ' --------------------------------------------------------')
    T = cmdl['-T']
    tee(out, '  Solving numerical equation with Nelder-Mead Simplex algorithm.. ')
    tee(out, '  Temperature used: %8.2f K' % T)

    bar = BAR(res_ab, res_ba, T=T, nboots=cmdl['-nruns'])
    # TODO: allow turning the bootstrapping off (e.g. nruns = 0 ):
    # it slows things down a lot but might not always be necessary
    tee(out, '  RESULT: dG (BAR ) = %8.4f kJ/mol' % bar.dg)
    tee(out, '  RESULT: error_dG_analyt (BAR ) = %8.4f kJ/mol' % bar.err)
    tee(out, '  RESULT: error_dG_bootstrap (BAR ) = %8.4f kJ/mol' % bar.err_boot)
    tee(out, '  RESULT: convergence (BAR ) = %8.4f' % bar.conv)
    tee(out, '  RESULT: convergence_error_bootstrap (BAR ) = %8.4f' % bar.conv_err_boot)
    tee(out, ' ------------------------------------------------------')
    diff = abs(cgi_result - bar.dg)
    mean = (cgi_result+bar.dg)*.5
    tee(out, '  Difference between BAR and CGI = %8.5f kJ/mol' % diff)
    tee(out, '  Mean of  BAR and CGI           = %8.5f kJ/mol' % mean)
    tee(out, ' ------------------------------------------------------')

    if cmdl['-jarz']:
        jarz = Jarz(wf=res_ab, wr=res_ba, T=T, nboots=cmdl['-nruns'])
        tee(out, ' --------------------------------------------------------')
        tee(out, '             ANALYSIS: Jarzynski estimator     ')
        tee(out, ' --------------------------------------------------------')
        tee(out, '  RESULT: dG_forward (Jarzynski) = %8.4f kJ/mol' % jarz.dg_for)
        tee(out, '  RESULT: dG_backward (Jarzynski) = %8.4f kJ/mol' % jarz.dg_rev)
        tee(out, '  RESULT: error_dG_bootstrap_forward (Jarzynski) = %8.4f kJ/mol' % jarz.err_boot_for)
        tee(out, '  RESULT: error_dG_bootstrap_backward (Jarzynski) = %8.4f kJ/mol' % jarz.err_boot_rev)
        tee(out, ' ------------------------------------------------------')
        tee(out, '  Mean of Jarzynski foward and backward = %8.5f kJ/mol' % jarz.dg_mean)
        tee(out, ' ------------------------------------------------------')

    if cmdl['-plot']:
        print '\n   Plotting histograms......'
        make_plot(cmdl['-cgi_plot'], res_ab, res_ba, cgi_result, cgi_err,
                  cmdl['-nbins'], cmdl['-dpi'])
        make_W_over_time_plot(cmdl['-W_over_t'], res_ab, res_ba, cgi_result,
                              cgi_err, cmdl['-nbins'], cmdl['-dpi'])

    tee(out, '\n   ......done...........\n')


if __name__ == '__main__':
    cmdl = parse_options(sys.argv)
    main(cmdl)
