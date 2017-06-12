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
from pmx.parser import *
from matplotlib import pyplot as plt
import numpy as np
from scipy.integrate import simps
from scipy.optimize import fmin
from scipy.special import erf
from random import gauss, randint, choice

# Constants
kb = 0.00831447215


# ==============================================================================
#                             Estimator Classes
# ==============================================================================
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
            mexp = mexp + np.exp(-beta*c*i)
            m = m + c*i
            m2 = m2 + i*i
        mexp = mexp/n
        m = m/n
        m2 = m2/n
        var = (m2-m*m)*(n/(n-1))
        # Jarzynski estimator
        dg = -kb*T*np.log(mexp)
        # Fluctuation-Dissipation estimator
        # FIXME: unused atm, remove or return?
        dg2 = m - beta*var/2.0

        return dg

    @staticmethod
    def calc_err_boot(w, T, c, nboots):
        out = []
        n = int(len(w))
        for k in range(nboots):
            sys.stdout.write('\r  Jarzynski error bootstrap: iteration %s/%s'
                             % (k+1, nboots))
            sys.stdout.flush()
            for i in range(n):
                val = [choice(w) for _ in xrange(n)]
            foo = -1.0 * Jarz.calc_dg(val, T, c)
            out.append(foo)
        sys.stdout.write('\n')
        err = np.std(out)

        return err


class Crooks(object):
    '''Crooks Gaussian Intersection (CGI) estimator. The forward and reverse work
    values are fitted to Gaussian functions and their intersection is taken
    as the free energy estimate. In some cases, when the two Gaussians are very
    close to each other, the intersection cannot be taken and the average of
    the two Gaussian means is taken as the free energy estimate insted. The
    standard error is by default calculated via bootstrap using 1000 samples.

    Parameters
    ----------
    wf : array_like
        array of forward work values.
    wr : array_like
        array of reverse work values.

    Examples
    --------
    >>> Crooks(wf, wr)
    >>> cgi_dg = Crooks.dg
    >>> cgi_err = Crooks.err_boot

    Attributes
    ----------
    dg : float
        the free energy estimate.
    err_boot : float
        standard error of the free energy estimate calculated via bootstrap.
    inters_bool : bool
        whether the interection could be taken. If False, the free energy
        estimate is the average of the two Gaussian means.
    Af : float
        height of the forward Gaussian.
    mf : float
        mean of the forward Gaussian.
    devf : float
        standard deviation of the forward Gaussian.
    Ar : float
        height of the reverse Gaussian.
    mr : float
        mean of the reverse Gaussian.
    devr : float
        standard deviation of the reverse Gaussian.
    '''

    def __init__(self, wf, wr):

        self.mf, self.devf, self.Af = data_to_gauss(wf)
        self.mr, self.devr, self.Ar = data_to_gauss(wr)

        # Calculate Crooks properties
        self.dg, self.inters_bool = self.calc_dg(A1=self.Af,
                                                 m1=self.mf,
                                                 s1=self.devf,
                                                 A2=self.Ar,
                                                 m2=self.mr,
                                                 s2=self.devr)
        self.err_boot = self.calc_err_boot(self.mf, self.devf, len(wf),
                                           self.mr, self.devr, len(wr),
                                           nboots=1000)

    @staticmethod
    def calc_dg(A1, m1, s1, A2, m2, s2):
        '''Calculates the free energy difference using the Crooks Gaussian
        Intersection method. It finds the intersection of two Gaussian
        functions. If the intersection cannot be computed, the average of
        the two Gaussian locations is returned.

        Parameters
        ----------
        A1 : float
            height of the forward Gaussian.
        m1 : float
            mean of the forward Gaussian.
        s1 : float
            standard deviation of the forward Gaussian.
        A2 : float
            height of the reverse Gaussian.
        m2 : float
            mean of the reverse Gaussian.
        s2 : float
            standard deviation of the reverse Gaussian.

        Returns
        -------
        float
            location of the intersection.
        bool
            whether the intersection could be calculated. If the intersection
            was calculated as expected a True value is returned.
            If the Gaussians are too close to each other, the intersection
            cannot be calculated and a False value is returned; in this case,
            the first float value retured is the average of the Gaussian means.
        '''

        p1 = m1/s1**2-m2/s2**2
        p2 = np.sqrt(1/(s1**2*s2**2)*(m1-m2)**2+2*(1/s1**2-1/s2**2)*np.log(s2/s1))
        p3 = 1/s1**2-1/s2**2
        x1 = (p1+p2)/p3
        x2 = (p1-p2)/p3

        # determine which solution to take
        if x1 > m1 and x1 < m2 or x1 > m2 and x1 < m1:
            dg = x1
            return dg, True
        elif x2 > m1 and x2 < m2 or x2 > m2 and x2 < m1:
            dg = x2
            return dg, True
        else:
            # we do not take the intersection but the average of the means
            dg = (m1 + m2) * 0.5
            return dg, False

    # Possible change of behaviour compared to the original script:
    # here it is not determined in advanced whether to take the intersection
    # or the mean, but for each bootstrap sample if the intersecion cannot
    # be taken, then the mean is used automatically.
    @staticmethod
    def calc_err_boot(m1, s1, n1, m2, s2, n2, nboots=1000):
        '''Calculates the standard error of the Crooks Gaussian Intersection
        via parametric bootstrap.

        Parameters
        ----------
        m1 : float
            mean of the forward Gaussian.
        s1 : float
            standard deviation of the forward Gaussian.
        n1 : int
            number of work values to which the first Gaussian was fit.
        m2 : float
            mean of the reverse Gaussian.
        s2 : float
            standard deviation of the reverse Gaussian.
        n2 : int
            number of work values to which the second Gaussian was fit.
        nboots: int
            number of bootstrap samples to use for the error estimate.
            Parametric bootstrap is used where work values are resampled from
            two Gaussians.
        intersection : bool
            True if the CGI free energy is to be taken from the Gaussians
            intersection, False if it is to be taken from the average of the
            means.

        Returns
        -------
        float
            standard error of the mean.
        '''

        iseq = []
        for k in range(nboots):
            g1 = []
            g2 = []
            for i in range(n1):
                g1.append(gauss(m1, s1))
            for i in range(n2):
                g2.append(gauss(m2, s2))

            m1, dev1, A1 = data_to_gauss(g1)
            m2, dev2, A2 = data_to_gauss(g2)
            iq, _ = Crooks.calc_dg(A1, m1, s1, A2, m2, s2)
            iseq.append(iq)
        err = np.std(iseq)
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
        self.wf = np.array(wf)
        self.wr = np.array(wr)
        self.T = float(T)
        self.nboots = nboots

        self.nf = len(wf)
        self.nr = len(wr)
        self.beta = 1./(kb*self.T)
        self.M = kb * self.T * np.log(float(self.nf) / float(self.nr))

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
        M = kb * T * np.log(nf/nr)

        def func(x, wf, wr):
            sf = 0
            for v in wf:
                sf += 1./(1+np.exp(beta*(M+v - x)))
                sr = 0
                for v in wr:
                    sr += 1./(1+np.exp(-beta*(M+v - x)))
            r = sf-sr
            return r**2

        avA = np.average(wf)
        avB = np.average(wr)
        x0 = (avA+avB)/2.
        dg = fmin(func, x0=x0, args=(wf, wr), disp=0)

        return float(dg)

    @staticmethod
    def calc_err(dg, wf, wr, T):
        '''Calculates the analytical error estimate.'''

        nf = float(len(wf))
        nr = float(len(wr))
        beta = 1./(kb*T)
        M = kb * T * np.log(nf/nr)

        err = 0
        for v in wf:
            err += 1./(2+2*np.cosh(beta * (M+v-dg)))
        for v in wr:
            err += 1./(2+2*np.cosh(beta * (M+v-dg)))
        N = nf + nr
        err /= float(N)
        tot = 1/(beta**2*N)*(1./err-(N/nf + N/nr))

        err = float(np.sqrt(tot))
        return err

    @staticmethod
    def calc_err_boot(wf, wr, nboots, T):
        '''Calculates the error by bootstrapping.'''

        nf = len(wf)
        nr = len(wr)
        res = []
        for k in range(nboots):
            sys.stdout.write('\r  BAR error bootstrap: iteration %s/%s'
                             % (k+1, nboots))
            sys.stdout.flush()
            for i in range(nf):
                valA = [choice(wf) for _ in xrange(nf)]
            for i in range(nr):
                valB = [choice(wr) for _ in xrange(nr)]
            foo = BAR.calc_dg(valA, valB, T)
            res.append(foo)
        sys.stdout.write('\n')
        err_boot = np.std(res)

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
    def calc_conv_err_boot(dg, wf, wr, nboots, T):
        nf = len(wf)
        nr = len(wr)
        res = []
        for k in range(nboots):
            sys.stdout.write('\r  Convergence error bootstrap: '
                             'iteration %s/%s' % (k+1, nboots))
            sys.stdout.flush()
            for i in range(nf):
                valA = [choice(wf) for _ in xrange(nf)]
            for i in range(nr):
                valB = [choice(wr) for _ in xrange(nr)]
            foo = BAR.calc_conv(dg, valA, valB, T)
            res.append(foo)
        sys.stdout.write('\n')
        err = np.std(res)
        return(err)


# ==============================================================================
#                               Functions
# ==============================================================================
def tee(fp, s):
    print >>fp, s
    print s


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
    '''Takes a one dimensional array and fits a Gaussian.

    Returns
    -------
    float
        mean of the distribution.
    float
        standard deviation of the distribution.
    float
        height of the curve's peak.
    '''
    m = np.average(data)
    dev = np.std(data)
    A = 1./(dev*np.sqrt(2*np.pi))
    return m, dev, A


def ks_norm_test(data, alpha=0.05, refks=None):
    '''Performs a Kolmogorov-Smirnov test of normality.

    Parameters
    ----------
    data : array_like
        a one-dimensional array of values. This is the distribution tested
        for normality.
    alpha : float
        significance level of the statistics. Default if 0.05.
    refks : ???
        ???

    Returns
    -------
    '''

    def ksref():
        f = 1
        potent = 10000
        lamb = np.arange(0.25, 2.5, 0.001)
        q = np.zeros(len(lamb), float)
        res = []
        for k in range(-potent, potent):
            q = q + f*np.exp(-2.0*(k**2)*(lamb**2))
            f = -f
        for i in range(len(lamb)):
            res.append((lamb[i], q[i]))
        return res

    def ksfunc(lamb):
        f = 1
        potent = 10000
        q = 0
        for k in range(-potent, potent):
            q = q + f*np.exp(-2.0*(k**2)*(lamb**2))
            f *= -1
        return q

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
        ndata = np.array(ndata)
        edf_ = np.array(edf_)
        return ndata, edf_

    def cdf(dg_data):
        data = deepcopy(dg_data)
        data.sort()
        mean = np.average(data)
        sig = np.std(data)
        cdf = 0.5*(1+erf((data-mean)/float(sig*np.sqrt(2))))
        return cdf

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
    check = np.sqrt(N)*dmax
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


def data_from_file(fn):
    data = read_and_format(fn, 'sf')
    return map(lambda a: a[1], data)


def dump_integ_file(fn, data):
    fp = open(fn, 'w')
    for fn, w in data:
        print >>fp, fn, w
    fp.close()


def gauss_func(A, mean, dev, x):
    x = np.array(x)
    y = A*np.exp(-(((x-mean)**2)/(2.0*(dev**2))))
    return y


# Is this needed? Seems redundandt since it's included in make_W_over_time_plot
def make_plot(fname, data1, data2, result, err, nbins, dpi=300):

    plt.figure(figsize=(8, 6))
    mf, devf, Af = data_to_gauss(data1)
    mb, devb, Ab = data_to_gauss(data2)

    maxi = max(data1+data2)
    mini = min(data1+data2)
    n1, bins1, patches1 = plt.hist(data1, range=(mini, maxi), bins=nbins,
                                   facecolor='blue', alpha=0.75, normed=True,
                                   label='0->1')
    n2, bins2, patches2 = plt.hist(data2, range=(mini, maxi), bins=nbins,
                                   facecolor='red', alpha=0.75, normed=True,
                                   label='1->0')
    plt.xlabel('W [kJ/mol]', fontsize=20)
    plt.ylabel('Probability', fontsize=20)
    plt.title(r'Work Distribution $\lambda$ 0->1 (blue) $\lambda$ 1->0 (red)')
    plt.grid(lw=2)
    loc, lab = plt.yticks()
    ll = []
    for i in range(len(lab)):
        ll.append("")
    plt.yticks(loc, ll)
    x = np.arange(mini, maxi, .5)
    y1 = gauss_func(Af, mf, devf, x)
    y2 = gauss_func(Ab, mb, devb, x)

    plt.plot(x, y1, 'b--', linewidth=2)
    plt.plot(x, y2, 'r--', linewidth=2)

    size = max([max(y1), max(y2)])
    res_x = [result, result]
    res_y = [0, size*1.2]
    plt.plot(res_x, res_y, 'k--', linewidth=2,
             label=r'$\Delta$G = %.2f $\pm$ %.2f kJ/mol' % (result, err))
    plt.legend(shadow=True, fancybox=True, prop={'size': 12})
    plt.ylim(0, size*1.2)
    xl = plt.gca()
    for val in xl.spines.values():
        val.set_lw(2)
    plt.savefig(fname, dpi=dpi)


def make_W_over_time_plot(fname, data1, data2, result, err, nbins, dpi=300):
    '''Plots work distributions and results for Crooks Gaussian Intersection'''

    def smooth(x, window_len=11, window='hanning'):

        if x.ndim != 1:
            raise ValueError("smooth only accepts 1 dimension arrays.")
        if x.size < window_len:
            raise ValueError("Input vector needs to be bigger than "
                             "window size.")
        if window_len < 3:
            return x
        if window not in ['flat', 'hanning', 'hamming',
                          'bartlett', 'blackman']:
            raise ValueError("Window is on of 'flat', 'hanning', 'hamming', "
                             "'bartlett', 'blackman'")
        s = np.r_[2*x[0]-x[window_len:1:-1], x, 2*x[-1]-x[-1:-window_len:-1]]
        # moving average
        if window == 'flat':
            w = ones(window_len, 'd')
        else:
            w = eval('np.' + window + '(window_len)')
        y = np.convolve(w/w.sum(), s, mode='same')
        return y[window_len-1:-window_len+1]

    plt.figure(figsize=(8, 6))
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

    sm1 = smooth(np.array(data1))
    sm2 = smooth(np.array(data2))
    plt.subplot(1, 2, 1)
    plt.plot(x1, data1, 'g-', linewidth=2, label="Forward (0->1)", alpha=.3)
    plt.plot(x1, sm1, 'g-', linewidth=3)
    plt.plot(x2, data2, 'b-', linewidth=2, label="Backward (1->0)", alpha=.3)
    plt.plot(x2, sm2, 'b-', linewidth=3)
    plt.legend(shadow=True, fancybox=True, loc='upper center',
               prop={'size': 12})
    plt.ylabel(r'W [kJ/mol]', fontsize=20)
    plt.xlabel(r'# Snapshot', fontsize=20)
    plt.grid(lw=2)
    plt.xlim(0, x[-1]+1)
    xl = plt.gca()
    for val in xl.spines.values():
        val.set_lw(2)
    plt.subplot(1, 2, 2)
    plt.hist(data1, bins=nbins, orientation='horizontal', facecolor='green',
             alpha=.75, normed=True)
    plt.hist(data2, bins=nbins, orientation='horizontal', facecolor='blue',
             alpha=.75, normed=True)

    x = np.arange(mini, maxi, .5)

    y1 = gauss_func(Af, mf, devf, x)
    y2 = gauss_func(Ab, mb, devb, x)

    plt.plot(y1, x, 'g--', linewidth=2)
    plt.plot(y2, x, 'b--', linewidth=2)
    size = max([max(y1), max(y2)])
    res_x = [result, result]
    res_y = [0, size*1.2]
    plt.plot(res_y, res_x, 'k--', linewidth=2,
             label=r'$\Delta$G = %.2f $\pm$ %.2f kJ/mol' % (result, err))
    plt.legend(shadow=True, fancybox=True, loc='upper center',
               prop={'size':12})
    plt.xticks([])
    plt.yticks([])
    xl = plt.gca()
    for val in xl.spines.values():
        val.set_lw(2)
    plt.subplots_adjust(wspace=0.0, hspace=0.1)
    plt.savefig(fname, dpi=dpi)


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

    from pmx import Option
    from pmx import FileOption
    from pmx import Commandline
    from pmx import __version__

    # TODO: option to choose units for output
    # TODO: choose which estimator to use as list. e.g. -est BAR CGI
    # TODO: pickle data and results?
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
        Option("-plot", "bool", False, "Plot work histograms"),
        Option("-nruns", "int", 0,
               "number of runs for bootstrapped BAR and Jarz errors. Default "
               "is 0, i.e. do not use bootstrap. For CGI, 1000 bootstrap "
               "samples are always used by default."),
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
                       version=__version__)
    cmdl.argv = argv
    return cmdl


def main(cmdl):

    out = open(cmdl['-o'], 'w')
    print >>out, "# analyze_crooks.py, pmx version = %s" % cmdl.version
    print >>out, "# pwd = %s" % os.getcwd()
    print >>out, "# %s (%s)" % (time.asctime(), os.environ.get('USER'))
    print >>out, "# command = %s" % ' '.join(cmdl.argv)
    print >>out, "\n\n"

    T = cmdl['-T']

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

    print('\n\n')
    tee(out, ' --------------------------------------------------------')
    tee(out, '                       ANALYSIS')
    tee(out, ' --------------------------------------------------------')
    tee(out, ' Number of forward (0->1) trajectories: %d' % len(res_ab))
    tee(out, ' Number of reverse (1->0) trajectories: %d' % len(res_ba))
    tee(out, ' Temperature : %.2f K' % T)

    # ============================
    # Crooks Gaussian Intersection
    # ============================
    tee(out, '\n --------------------------------------------------------')
    tee(out, '             Crooks Gaussian Intersection     ')
    tee(out, ' --------------------------------------------------------')

    print('  Calculating Intersection...')
    cgi = Crooks(wf=res_ab, wr=res_ba)

    tee(out, '  Forward  : mean = %8.3f  std = %8.3f' % (cgi.mf, cgi.devf))
    tee(out, '  Backward : mean = %8.3f  std = %8.3f' % (cgi.mr, cgi.devr))
    if cgi.inters_bool is False:
        tee(out, '\n  Gaussians too close for intersection calculation')
        tee(out, '   --> Taking difference of mean values')
    tee(out, '  RESULT: dG (CGI)  = %8.4f kJ/mol' % cgi.dg)
    tee(out, '  RESULT: error_dG (CGI) = %8.4f kJ/mol' % cgi.err_boot)

    # --------------
    # Normality test
    # --------------
    if cmdl['-KS']:
        print('\n  Running KS-test ....')
        q0, lam00, check0, bOk0 = ks_norm_test(res_ab)
        q1, lam01, check1, bOk1 = ks_norm_test(res_ba)

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

    # ========================
    # Bennett Acceptance Ratio
    # ========================
    tee(out, '\n --------------------------------------------------------')
    tee(out, '             Bennett Acceptance Ratio     ')
    tee(out, ' --------------------------------------------------------')

    bar = BAR(res_ab, res_ba, T=T, nboots=cmdl['-nruns'])

    print('  Solving numerical equation with Nelder-Mead Simplex algorithm... ')

    # TODO: allow turning the bootstrapping off (e.g. nruns = 0 ):
    # it slows things down a lot but might not always be necessary
    tee(out, '  RESULT: dG (BAR) = %8.4f kJ/mol' % bar.dg)
    tee(out, '  RESULT: error_dG_analyt (BAR) = %8.4f kJ/mol' % bar.err)
    if cmdl['-nruns'] > 0:
        tee(out, '  RESULT: error_dG_bootstrap (BAR) = %8.4f kJ/mol' % bar.err_boot)
    tee(out, '  RESULT: convergence (BAR) = %8.4f' % bar.conv)
    if cmdl['-nruns'] > 0:
        tee(out, '  RESULT: convergence_error_bootstrap (BAR) = %8.4f' % bar.conv_err_boot)
    tee(out, ' ------------------------------------------------------')
    # Should we make this optional?
    diff = abs(cgi.dg - bar.dg)
    mean = (cgi.dg + bar.dg) * 0.5
    tee(out, '  Difference between BAR and CGI = %8.5f kJ/mol' % diff)
    tee(out, '  Mean of  BAR and CGI           = %8.5f kJ/mol' % mean)
    tee(out, ' ------------------------------------------------------')

    # =========
    # Jarzynski
    # =========
    if cmdl['-jarz']:
        tee(out, '\n --------------------------------------------------------')
        tee(out, '             Jarzynski estimator     ')
        tee(out, ' --------------------------------------------------------')

        jarz = Jarz(wf=res_ab, wr=res_ba, T=T, nboots=cmdl['-nruns'])

        tee(out, '  RESULT: dG_forward (Jarzynski) = %8.4f kJ/mol' % jarz.dg_for)
        tee(out, '  RESULT: dG_backward (Jarzynski) = %8.4f kJ/mol' % jarz.dg_rev)
        if cmdl['-nruns'] > 0:
            tee(out, '  RESULT: error_dG_bootstrap_forward (Jarzynski) = %8.4f kJ/mol' % jarz.err_boot_for)
        if cmdl['-nruns'] > 0:
            tee(out, '  RESULT: error_dG_bootstrap_backward (Jarzynski) = %8.4f kJ/mol' % jarz.err_boot_rev)
        tee(out, ' ------------------------------------------------------')
        tee(out, '  Mean of Jarzynski foward and backward = %8.5f kJ/mol' % jarz.dg_mean)
        tee(out, ' ------------------------------------------------------')

    if cmdl['-plot']:
        print '\n   Plotting histograms......'
        make_plot(cmdl['-cgi_plot'], res_ab, res_ba, cgi.dg, cgi.err_boot,
                  cmdl['-nbins'], cmdl['-dpi'])
        make_W_over_time_plot(cmdl['-W_over_t'], res_ab, res_ba, cgi.dg,
                              cgi.err_boot, cmdl['-nbins'], cmdl['-dpi'])

    print('\n   ......done...........\n')


if __name__ == '__main__':
    cmdl = parse_options(sys.argv)
    main(cmdl)
