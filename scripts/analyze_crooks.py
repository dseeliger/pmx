#!/usr/bin/env python
# pmx  Copyright Notice
# ============================
#
# The pmx source code is copyrighted, but you can freely use and
# copy it as long as you don't change or remove any of the copyright
# notices.
#
# ----------------------------------------------------------------------
# pmx is Copyright (C) 2006-2017 by Daniel Seeliger
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

from __future__ import division
import sys
import os
import time
import re
from copy import deepcopy
from pmx.parser import read_and_format
from matplotlib import pyplot as plt
import numpy as np
from scipy.integrate import simps
from scipy.optimize import fmin
from scipy.special import erf
import scipy.stats
import pickle
import argparse

# Constants
kb = 0.00831447215   # kJ/(K*mol)


# ==============================================================================
#                             Estimator Classes
# ==============================================================================
# TODO: make estimators a separete module and do e.g.:
# from pmx.estimators import Jarz, BAR, Crooks
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

    def __init__(self, wf, wr, T, nboots=0, nblocks=1):
        self.wf = np.array(wf)
        self.wr = np.array(wr)
        self.T = float(T)
        self.nboots = nboots
        self.nblocks = nblocks

        # Calculate all Jarz properties available
        self.dg_for = self.calc_dg(w=self.wf, c=1.0, T=self.T)
        self.dg_rev = -1.0 * self.calc_dg(w=self.wr, c=-1.0, T=self.T)
        self.dg_mean = (self.dg_for + self.dg_rev) * 0.5

        if nboots > 0:
            self.err_boot_for = self.calc_err_boot(w=self.wf, T=self.T,
                                                   c=1.0, nboots=nboots)
            self.err_boot_rev = self.calc_err_boot(w=self.wr, T=self.T,
                                                   c=-1.0, nboots=nboots)

        if nblocks > 1:
            self.err_blocks_for = self.calc_err_blocks(w=self.wf, c=1.0,
                                                       T=self.T,
                                                       nblocks=nblocks)
            self.err_blocks_rev = self.calc_err_blocks(w=self.wr, c=-1.0,
                                                       T=self.T,
                                                       nblocks=nblocks)

    @staticmethod
    def calc_dg(w, T, c):
        '''to be filled
        '''
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
        '''Calculates the standard error via bootstrap. The work values are
        resampled randomly with replacement multiple (nboots) times,
        and the Jarzinski free energy recalculated for each bootstrap samples.
        The standard error of the estimate is returned as the standard d
        eviation of the bootstrapped free energies.

        Parameters
        ----------
        w : array_like
            work values.
        T : float
            temperature.
        c : [0,1]
            ???
        nboots: int
            number of bootstrap samples to use for the error estimate.

        Returns
        -------
        err : float
            standard error of the mean.
        '''
        dg_boots = []
        n = len(w)
        for k in range(nboots):
            sys.stdout.write('\r  Bootstrap (Std Err): iteration %s/%s'
                             % (k+1, nboots))
            sys.stdout.flush()

            boot = np.random.choice(w, size=n, replace=True)
            dg_boot = -1.0 * Jarz.calc_dg(boot, T, c)
            dg_boots.append(dg_boot)
        sys.stdout.write('\n')
        err = np.std(dg_boots)
        return err

    @staticmethod
    def calc_err_blocks(w, T, c, nblocks):
        '''Calculates the standard error based on a number of blocks the
        work values are divided into. It is useful when you run independent
        equilibrium simulations, so that you can then use their respective
        work values to compute the standard error based on the repeats.

        Parameters
        ----------
        w : array_like
            array of work values.
        T : float
            temperature.
        c : [0,1]
            ???
        nblocks: int
            number of blocks to divide the data into. This can be for
            instance the number of independent equilibrium simulations
            you ran.
        '''

        # get number of work values
        len_w = len(w)
        # here we assume/expect same lenght for wf and wr
        assert len_w % float(nblocks) == 0.0
        dg_blocks = []

        nf_block = int(len_w/nblocks)
        first = 0
        last = nf_block
        # calculate all dg
        for i in range(nblocks):
            w_block = np.array(w[first:last])
            dg_block = -1.0 * Jarz.calc_dg(w_block, T, c)
            dg_blocks.append(dg_block)
            first += nf_block
            last += nf_block

        # get std err
        err_blocks = scipy.stats.sem(dg_blocks, ddof=1)

        return err_blocks


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

    def __init__(self, wf, wr, nboots=0, nblocks=1):

        # inputs
        self.wf = np.array(wf)
        self.wr = np.array(wr)
        self.nboots = nboots
        self.nblocks = nblocks
        # params of the gaussians
        self.mf, self.devf, self.Af = data_to_gauss(wf)
        self.mr, self.devr, self.Ar = data_to_gauss(wr)

        # Calculate Crooks properties
        self.dg, self.inters_bool = self.calc_dg(wf=self.wf, wr=self.wr)

        self.err_boot1 = self.calc_err_boot1(m1=self.mf, s1=self.devf,
                                             n1=len(wf), m2=self.mr,
                                             s2=self.devr, n2=len(wr),
                                             nboots=1000)
        if nboots > 0:
            self.err_boot2 = self.calc_err_boot2(wf=self.wf, wr=self.wr,
                                                 nboots=nboots)

        if nblocks > 1:
            self.err_blocks = self.calc_err_blocks(self.wf, self.wr, nblocks)

    @staticmethod
    def calc_dg(wf, wr):
        '''Calculates the free energy difference using the Crooks Gaussian
        Intersection method. It finds the intersection of two Gaussian
        functions. If the intersection cannot be computed, the average of
        the two Gaussian locations is returned.

        Parameters
        ----------
        wf : array_like
            array of forward work values.
        wr : array_like
            array of reverse work values.

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

        m1, s1, A1 = data_to_gauss(wf)
        m2, s2, A2 = data_to_gauss(wr)

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
    def calc_err_boot1(m1, s1, n1, m2, s2, n2, nboots=1000):
        '''Calculates the standard error of the Crooks Gaussian Intersection
        via parametric bootstrap. Given the parameters of the forward and
        reverse Gaussian distributions, multiple (nboots) bootstrap samples
        are built by random sampling from these two Gaussian distributions.
        The CGI free energy is then calculated for each bootstrap sample
        (forward and reverse Gaussians). The standard error of the estimate
        is returned as the standard deviation of the bootstrapped free
        energies.

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

        Returns
        -------
        float
            standard error of the mean.
        '''

        dg_boots = []
        for k in range(nboots):
            bootA = np.random.normal(loc=m1, scale=s1, size=n1)
            bootB = np.random.normal(loc=m2, scale=s2, size=n2)

            dg_boot, _ = Crooks.calc_dg(bootA, bootB)
            dg_boots.append(dg_boot)
        err = np.std(dg_boots)
        return err

    @staticmethod
    def calc_err_boot2(wf, wr, nboots):
        '''Calculates the standard error of the Crooks Gaussian Intersection
        via non-parametric bootstrap. The work values are resampled randomly
        with replacement multiple (nboots) times, and the CGI free energy
        recalculated for each bootstrap samples. The standard error of
        the estimate is returned as the standard deviation of the bootstrapped
        free energies.

        Parameters
        ----------
        wf : array_like
            array of forward work values.
        wr : array_like
            array of reverse work values.
        nboots: int
            number of bootstrap samples to use for the error estimate.

        Returns
        -------
        err : float
            standard error of the mean.
        '''

        nf = len(wf)
        nr = len(wr)

        dg_boots = []
        for k in range(nboots):
            sys.stdout.write('\r  Bootstrap (Std Err): iteration %s/%s'
                             % (k+1, nboots))
            sys.stdout.flush()

            bootA = np.random.choice(wf, size=nf, replace=True)
            bootB = np.random.choice(wr, size=nr, replace=True)

            dg_boot, _ = Crooks.calc_dg(bootA, bootB)
            dg_boots.append(dg_boot)

        sys.stdout.write('\n')
        err = np.std(dg_boots)
        return err

    @staticmethod
    def calc_err_blocks(wf, wr, nblocks):
        '''Calculates the standard error based on a number of blocks the
        work values are divided into. It is useful when you run independent
        equilibrium simulations, so that you can then use their respective
        work values to compute the standard error based on the repeats.

        Parameters
        ----------
        wf : array_like
            array of forward work values.
        wr : array_like
            array of reverse work values.
        nblocks: int
            number of blocks to divide the data into. This can be for
            instance the number of independent equilibrium simulations
            you ran.
        '''

        # get number of work values
        len_wf = len(wf)
        len_wr = len(wr)
        # here we assume/expect same lenght for wf and wr
        assert len_wf == len_wr
        assert len_wf % float(nblocks) == 0.0
        dg_blocks = []

        nf_block = int(len_wf/nblocks)
        first = 0
        last = nf_block
        # calculate all dg
        for i in range(nblocks):
            wf_block = np.array(wf[first:last])
            wr_block = np.array(wr[first:last])
            dg_block, _ = Crooks.calc_dg(wf_block, wr_block)
            dg_blocks.append(dg_block)
            first += nf_block
            last += nf_block

        # get std err
        err_blocks = scipy.stats.sem(dg_blocks, ddof=1)

        return err_blocks


class BAR(object):
    '''Bennett acceptance ratio (BAR).

    Description...

    Parameters
    ----------

    Examples
    --------
    '''

    def __init__(self, wf, wr, T, nboots=0, nblocks=1):
        self.wf = np.array(wf)
        self.wr = np.array(wr)
        self.T = float(T)
        self.nboots = nboots
        self.nblocks = nblocks

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
        if nblocks > 1:
            self.err_blocks = self.calc_err_blocks(self.wf, self.wr, nblocks,
                                                   self.T)

    @staticmethod
    def calc_dg(wf, wr, T):
        '''Estimates and returns the free energy difference.

        Parameters
        ----------
        wf : array_like
            array of forward work values.
        wr : array_like
            array of reverse work values.
        T : float
            temperature

        Returns
        ----------
        dg : float
            the BAR free energy estimate.
        '''

        nf = float(len(wf))
        nr = float(len(wr))
        beta = 1./(kb*T)
        M = kb * T * np.log(nf/nr)

        def func(x, wf, wr):
            sf = 0
            for v in wf:
                sf += 1./(1+np.exp(beta*(M+v-x)))

            sr = 0
            for v in wr:
                sr += 1./(1+np.exp(-beta*(M+v-x)))

            r = sf-sr
            return r**2

        avA = np.average(wf)
        avB = np.average(wr)
        x0 = (avA+avB)/2.
        dg = fmin(func, x0=x0, args=(wf, wr), disp=0)

        return float(dg)

    @staticmethod
    def calc_err(dg, wf, wr, T):
        '''Calculates the analytical error estimate.

        Parameters
        ----------
        dg : float
            the BAR free energy estimate
        wf : array_like
            array of forward work values.
        wr : array_like
            array of reverse work values.
        T : float
            temperature
        '''

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
        '''Calculates the error by bootstrapping.

        Parameters
        ----------
        wf : array_like
            array of forward work values.
        wr : array_like
            array of reverse work values.
        T : float
            temperature
        nboots: int
            number of bootstrap samples.

        '''

        nf = len(wf)
        nr = len(wr)
        dg_boots = []
        for k in range(nboots):
            sys.stdout.write('\r  Bootstrap (Std Err): iteration %s/%s'
                             % (k+1, nboots))
            sys.stdout.flush()

            bootA = np.random.choice(wf, size=nf, replace=True)
            bootB = np.random.choice(wr, size=nr, replace=True)
            dg_boot = BAR.calc_dg(bootA, bootB, T)
            dg_boots.append(dg_boot)

        sys.stdout.write('\n')
        err_boot = np.std(dg_boots)

        return err_boot

    @staticmethod
    def calc_err_blocks(wf, wr, nblocks, T):
        '''Calculates the standard error based on a number of blocks the
        work values are divided into. It is useful when you run independent
        equilibrium simulations, so that you can then use their respective
        work values to compute the standard error based on the repeats.

        Parameters
        ----------
        wf : array_like
            array of forward work values.
        wr : array_like
            array of reverse work values.
        T : float
            temperature
        nblocks: int
            number of blocks to divide the data into. This can be for
            instance the number of independent equilibrium simulations
            you ran.
        '''

        # get number of work values
        len_wf = len(wf)
        len_wr = len(wr)
        # here we assume/expect same lenght for wf and wr
        assert len_wf == len_wr
        assert len_wf % float(nblocks) == 0.0
        dg_blocks = []

        nf_block = int(len_wf/nblocks)
        first = 0
        last = nf_block
        # calculate all dg
        for i in range(nblocks):
            wf_block = np.array(wf[first:last])
            wr_block = np.array(wr[first:last])
            dg_block = BAR.calc_dg(wf_block, wr_block, T)
            dg_blocks.append(dg_block)
            first += nf_block
            last += nf_block

        # get std err
        err_blocks = scipy.stats.sem(dg_blocks, ddof=1)

        return err_blocks

    @staticmethod
    def calc_conv(dg, wf, wr, T):
        '''Evaluates BAR convergence as described in Hahn & Then, Phys Rev E
        (2010), 81, 041117. Returns a value between -1 and 1: the closer this
        value to zero the better the BAR convergence.

        Parameters
        ----------
        dg : float
            the BAR free energy estimate
        wf : array_like
            array of forward work values.
        wr : array_like
            array of reverse work values.
        T : float
            temperature

        '''

        wf = np.array(wf)
        wr = np.array(wr)

        beta = 1./(kb*T)
        nf = len(wf)
        nr = len(wr)
        N = float(nf + nr)

        ratio_alpha = float(nf)/N
        ratio_beta = float(nr)/N
        bf = 1.0/(ratio_beta + ratio_alpha * np.exp(beta*(wf-dg)))
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
        conv_boots = []
        for k in range(nboots):
            sys.stdout.write('\r  Bootstrap (Conv): '
                             'iteration %s/%s' % (k+1, nboots))
            sys.stdout.flush()

            bootA = np.random.choice(wf, size=nf, replace=True)
            bootB = np.random.choice(wr, size=nr, replace=True)
            conv_boot = BAR.calc_conv(dg, bootA, bootB, T)
            conv_boots.append(conv_boot)

        sys.stdout.write('\n')
        err = np.std(conv_boots)
        return err


# ==============================================================================
#                               FUNCTIONS
# ==============================================================================
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
    Q : float
    lam0 : float
    check : float
    bOk : bool
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


def gauss_func(A, mean, dev, x):
    '''Given the parameters of a Gaussian and a range of the x-values, returns
    the y-values of the Gaussian function'''
    x = np.array(x)
    y = A*np.exp(-(((x-mean)**2.)/(2.0*(dev**2.))))
    return y


# -------------
# Files Parsing
# -------------
def parse_dgdl_files(lst, lambda0=0, invert_values=False):
    '''Takes a list of dgdl.xvg files and returns the integrated work values

    Parameters
    ----------
    lst : list
        list containing the paths to the dgdl.xvg files.
    lambda0 : [0,1]
        whether the simulations started from lambda 0 or 1. Default is 0.
    invert_values : bool
        whether to invert the sign of the returned work value.

    Returns
    -------
    w : array
        array of work values.
    lst : list
        sorted list of input dgdl.avg files corresponding to the work values
        in w.
    '''
    print '\nProcessing simulation data......'

    # check lambda0 is either 0 or 1
    assert lambda0 in [0, 1]

    _check_dgdl(lst[0], lambda0)
    first_w, ndata = integrate_dgdl(lst[0], lambda0=lambda0,
                                    invert_values=invert_values)
    w_list = [first_w]
    for idx, f in enumerate(lst[1:]):
        sys.stdout.write('\r------>  %s' % f)
        sys.stdout.flush()

        w, _ = integrate_dgdl(f, ndata=ndata, lambda0=lambda0,
                              invert_values=invert_values)
        if w is not None:
            w_list.append(w)

    return w_list


def integrate_dgdl(fn, ndata=-1, lambda0=0, invert_values=False):
    '''Integrates the data in a dgdl.xvg file.

    Parameters
    ----------
    fn : str
        the inpur dgdl.xvg file from Gromacs.
    ndata : int, optional
        number of datapoints in file. If -1, then ??? default is -1.
    lambda0 : [0,1]
        whether the simulations started from lambda 0 or 1. Default is 0.
    invert_values : bool
        whether to invert the sign of the returned work value.

    Returns
    -------
    integr : float
        result of the integration performed using Simpson's rule.
    ndata : int
        number of data points in the input file.
    '''

    # check lambda0 is either 0 or 1
    assert lambda0 in [0, 1]

    lines = open(fn).readlines()
    if not lines:
        return None, None

    # extract dgdl datapoints into r
    # TODO: we removed the check for file integrity. We could have an
    # optional files integrity check before calling this integration func

    lines = [l for l in lines if l[0] not in '#@&']
    r = map(lambda x: float(x.split()[1]), lines)

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

    # arrays for the integration
    # --------------------------
    # array of lambda values
    x = [lambda0+i*dlambda for i, dgdl in enumerate(r)]
    # array of dgdl
    y = r

    if lambda0 == 1:
        x.reverse()
        y.reverse()

    if invert_values is True:
        integr = simps(y, x) * (-1)
        return integr, ndata
    else:
        integr = simps(y, x)
        return integr, ndata


def _check_dgdl(fn, lambda0):
    '''Prints some info about a dgdl.xvg file.'''
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


def _dump_integ_file(outfn, f_lst, w_lst):
    with open(outfn, 'w') as f:
        for fn, w in zip(f_lst, w_lst):
            f.write('{dhdl} {work}\n'.format(dhdl=fn, work=w))


def _data_from_file(fn):
    data = read_and_format(fn, 'sf')
    return map(lambda a: a[1], data)


# ------------------
# Plotting functions
# ------------------
def make_cgi_plot(fname, data1, data2, result, err, nbins, dpi=300):
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
            w = np.ones(window_len, 'd')
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
               prop={'size': 12})
    plt.xticks([])
    plt.yticks([])
    xl = plt.gca()
    for val in xl.spines.values():
        val.set_lw(2)
    plt.subplots_adjust(wspace=0.0, hspace=0.1)
    plt.savefig(fname, dpi=dpi)


# ---------------------
# Some helper functions
# ---------------------
def _tee(fp, s):
    print >>fp, s
    print s


def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


def time_stats(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return h, m, s


# ==============================================================================
#                      COMMAND LINE OPTIONS AND MAIN
# ==============================================================================
def parse_options():

    parser = argparse.ArgumentParser(description='Calculates free energies '
            'from fast growth thermodynamic integration simulations. '
            'Available methods for free energy estimation: '
            'Crooks Gaussian Intersection (CGI); '
            'Benett Acceptance Ratio (BAR); '
            'Jarzinski equality (JARZ).')

    exclus = parser.add_mutually_exclusive_group()

    parser.add_argument('-fA',
                        metavar='dgdl',
                        dest='filesAB',
                        type=str,
                        help='dgdl.xvg files for the A->B simulations. Use '
                        'wildcard to select multiple xvg files: e.g. "-fa '
                        './forward_results/dgdl*.xvg"',
                        required=True,
                        nargs='+')
    parser.add_argument('-fB',
                        metavar='dgdl',
                        dest='filesBA',
                        type=str,
                        help='dgdl.xvg files for the B->A simulations Use '
                        'wildcard to select multiple xvg files: e.g. "-fb '
                        './backward_results/dgdl*.xvg"',
                        required=True,
                        nargs='+')
    parser.add_argument('-m',
                        metavar='method',
                        type=str.lower,
                        dest='methods',
                        help='Choose one or more estimators to use from the '
                        'available ones: CGI, BAR, JARZ. Default is all.',
                        default=['cgi', 'bar', 'jarz'],
                        nargs='+')
    parser.add_argument('-t',
                        metavar='temperature',
                        dest='temperature',
                        type=float,
                        help='Temperature in Kelvin. Default is 298.15.',
                        default=298.15)
    parser.add_argument('-o',
                        metavar='result file',
                        dest='outfn',
                        type=str,
                        help='Filename of output result file. Default is '
                        '"results.txt."',
                        default='results.txt')
    parser.add_argument('-b',
                        metavar='nboots',
                        dest='nboots',
                        type=int,
                        help='Number of bootstrap samples to use for the '
                        'bootstrap estimate of the standard errors. Default '
                        'is 0 (no bootstrap).',
                        default=0)
    parser.add_argument('-n',
                        metavar='nblocks',
                        dest='nblocks',
                        type=int,
                        help='Number of blocks to divide the data into for '
                        'an estimate of the standard error. You can use this '
                        'when multiple independent equilibrium simulations'
                        'have been run so to estimate the error from the '
                        'repeats. Default is 1 (i.e. no repeats). It assumes '
                        'the dgdl files for each repeat are read in order and '
                        'are contiguous, e.g. dgdl_0 to dgdl_9 is the first '
                        'repeat, dgdl_10 to dgdl_19 is the second one, etc.',
                        default=1)
    parser.add_argument('--integ_only',
                        dest='integ_only',
                        help='Whether to do integration only; the integrated '
                        'values are computed and saved, and the program '
                        'terminated. Default is False.',
                        default=False,
                        action='store_true')
    parser.add_argument('-iA',
                        metavar='work input',
                        dest='iA',
                        type=str,
                        help='Two-column dat file containing the list of input'
                        ' files and their respective integrated work values '
                        'for the forward (A->B) tranformation.')
    parser.add_argument('-iB',
                        metavar='work input',
                        dest='iB',
                        type=str,
                        help='Two-column dat file containing the list of input'
                        ' files and their respective integrated work values '
                        'for the reverse (B->A) tranformation.')
    parser.add_argument('-oA',
                        metavar='work output',
                        dest='oA',
                        type=str,
                        help='File where to save the list of input dgdl'
                        ' files and their respective integrated work values '
                        'for the forward (A->B) tranformation. Default is '
                        '"integA.dat"',
                        default='integA.dat')
    parser.add_argument('-oB',
                        metavar='work output',
                        dest='oB',
                        type=str,
                        help='File where to save the list of input dgdl'
                        ' files and their respective integrated work values '
                        'for the reverse (B->A) tranformation. Default is '
                        '"integB.dat"',
                        default='integB.dat')
    parser.add_argument('--reverseB',
                        dest='reverseB',
                        help='Whether to reverse the work values for the '
                        'backward (B->A) transformation. This is useful '
                        'when in Gromacs both forward and reverse simulations '
                        'were run from lambda zero to one.'
                        'Default is False.',
                        default=False,
                        action='store_true')
    # The following are mutually exclusive options
    exclus.add_argument('--skip',
                        metavar='',
                        dest='skip',
                        type=int,
                        help='Skip files, i.e. pick every nth work value. '
                        'Default is 1 (all); with 2, every other work value '
                        'is discarded, etc.',
                        default=1)
    exclus.add_argument('--slice',
                        metavar='',
                        dest='slice',
                        type=int,
                        help='Subset of trajectories to analyze.'
                        'Provide list slice, e.g. "10 50" will'
                        ' result in selecting dgdl_files[10:50].'
                        ' Default is all.',
                        default=None,
                        nargs=2)
    exclus.add_argument('--rand',
                        metavar='',
                        dest='rand',
                        type=int,
                        help='Take a random subset of trajectories. '
                        'Default is None (do not take random subset)',
                        default=None)
    exclus.add_argument('--index',
                        metavar='',
                        dest='index',
                        type=int,
                        help='Zero-based index of files to analyze (e.g.'
                        ' 0 10 20 50 60). It keeps '
                        'the dgdl.xvg files according to their position in the'
                        ' list, sorted according to the filenames. Default '
                        'is None (i.e. all dgdl are used).',
                        default=None,
                        nargs='+')
    parser.add_argument('--prec',
                        metavar='',
                        dest='precision',
                        type=int,
                        help='The decimal precision of the screen/file output.'
                        ' Default is 2.',
                        default=2)
    parser.add_argument('--units',
                        metavar='',
                        dest='units',
                        type=str.lower,
                        help='The units of the output. Choose from "kJ", '
                        '"kcal", "kT". Default is "kJ."',
                        default='kJ',
                        choices=['kj', 'kcal', 'kt'])
    parser.add_argument('--pickle',
                        dest='pickle',
                        help='Whether to save the free energy results from '
                        'the estimators in pickled files. Default is False.',
                        default=False,
                        action='store_true')
    parser.add_argument('--no_ks',
                        dest='do_ks_test',
                        help='Whether to do a Kolmogorov-Smirnov test '
                        'to check whether the Gaussian assumption for CGI '
                        'holds. Default is True; this flag turns it to False.',
                        default=True,
                        action='store_false')
    parser.add_argument('--cgi_plot',
                        metavar='',
                        dest='cgi_plot',
                        type=str,
                        help='Whether to plot the work histograms along with '
                        'the CGI results. If the flag is used, you also need'
                        'to specify a filename.',
                        default=None)
    parser.add_argument('--nbins',
                        metavar='',
                        dest='nbins',
                        type=int,
                        help='Number of histograms bins for the plot. '
                        'Default is 10.',
                        default=10)
    parser.add_argument('--dpi',
                        metavar='',
                        dest='dpi',
                        type=int,
                        help='Resolution of the plot. Default is 300.',
                        default=300)

    args = parser.parse_args()

    from pmx import __version__
    args.pmx_version = __version__

    return args


def main(args):
    """Run the main script.

    Parameters
    ----------
    args : argparse.Namespace
        The command line arguments
    """

    # start timing
    stime = time.time()

    # input arguments
    out = open(args.outfn, 'w')
    filesAB = natural_sort(args.filesAB)
    filesBA = natural_sort(args.filesBA)
    T = args.temperature
    skip = args.skip
    prec = args.precision
    methods = args.methods
    reverseB = args.reverseB
    integ_only = args.integ_only
    nboots = args.nboots
    nblocks = args.nblocks
    do_ks_test = args.do_ks_test

    # -------------------
    # Select output units
    # -------------------
    units = args.units
    if units.lower() == 'kj':
        # kJ is the input from GMX
        unit_fact = 1.
        units = 'kJ/mol'
    elif units == 'kcal':
        unit_fact = 1./4.184
        units = 'kcal/mol'
    elif units.lower() == 'kt':
        unit_fact = 1./(kb*T)
        units = 'kT'
    else:
        exit('No unit type \'%s\' available' % units)

    print >>out, "# analyze_crooks.py, pmx version = %s" % args.pmx_version
    print >>out, "# pwd = %s" % os.getcwd()
    print >>out, "# %s (%s)" % (time.asctime(), os.environ.get('USER'))
    print >>out, "# command = %s" % ' '.join(sys.argv)
    _tee(out, "\n\n")

    # ==========
    # Parse Data
    # ==========

    # If list of dgdl.xvg files are provided, parse dgdl
    if args.iA is None and args.iB is None:
        # If random selection is chosen, do this before reading files and
        # calculating the work values.
        if args.rand is not None:
            filesAB = np.random.choice(filesAB, size=args.rand, replace=False)
            filesBA = np.random.choice(filesBA, size=args.rand, replace=False)
            _tee(out, 'Selected random subset of %d trajectories.' % args.rand)

        # If slice values provided, select the files needed. Again before
        # reading files so speed up the process
        if args.slice is not None:
            first = args.slice[0]
            last = args.slice[1]
            _tee(out, 'First trajectories read: %s and %s'
                 % (filesAB[first], filesBA[first]))
            _tee(out, 'Last trajectories  read: %s and %s'
                 % (filesAB[last-1], filesBA[last-1]))
            filesAB = filesAB[first:last]
            filesBA = filesBA[first:last]

        # If index values provided, select the files needed
        if args.index is not None:
            filesAB = [filesAB[i] for i in args.index]
            filesBA = [filesBA[i] for i in args.index]

        # when skipping start count from end: in this way the last frame is
        # always included, and what can change is the first one
        filesAB = list(reversed(filesAB[::-skip]))
        filesBA = list(reversed(filesBA[::-skip]))

        # Now read in the data
        res_ab = parse_dgdl_files(filesAB, lambda0=0,
                                  invert_values=False)
        res_ba = parse_dgdl_files(filesBA, lambda0=1,
                                  invert_values=reverseB)

        _dump_integ_file(args.oA, filesAB, res_ab)
        _dump_integ_file(args.oB, filesBA, res_ba)

    # If work values are given as input instead, read those
    elif args.iA is not None and args.iB is not None:
        res_ab = []
        res_ba = []
        for fn in args.iA:
            print '\t\tReading integrated values (A->B) from', fn
            res_ab.extend(_data_from_file(fn))
        for fn in args.iB:
            print '\t\tReading integrated values (B->A) from', fn
            res_ba.extend(_data_from_file(fn))
    else:
        exit('\nERROR: you need to provide either none of both sets of '
             'integrated work values.')

    # If asked to only do the integration of dhdl.xvg, exit
    if integ_only:
        print('\n    Integration done. Skipping analysis.')
        print('\n    ......done........\n')
        sys.exit(0)

    # ==============
    # Begin Analysis
    # ==============
    print('\n\n')
    _tee(out, ' --------------------------------------------------------')
    _tee(out, '                       ANALYSIS')
    _tee(out, ' --------------------------------------------------------')
    _tee(out, '  Number of forward (0->1) trajectories: %d' % len(res_ab))
    _tee(out, '  Number of reverse (1->0) trajectories: %d' % len(res_ba))
    _tee(out, '  Temperature : %.2f K' % T)

    # ============================
    # Crooks Gaussian Intersection
    # ============================
    if 'cgi' in methods:
        _tee(out, '\n --------------------------------------------------------')
        _tee(out, '             Crooks Gaussian Intersection     ')
        _tee(out, ' --------------------------------------------------------')

        print('  Calculating Intersection...')
        cgi = Crooks(wf=res_ab, wr=res_ba, nboots=nboots, nblocks=nblocks)
        if args.pickle is True:
            pickle.dump(cgi, open("cgi_results.pkl", "wb"))

        _tee(out, '  CGI: Forward Gauss mean = {m:8.{p}f} {u} '
                  'std = {s:8.{p}f} {u}'.format(m=cgi.mf*unit_fact,
                                                s=cgi.devf*unit_fact,
                                                p=prec, u=units))
        _tee(out, '  CGI: Reverse Gauss mean = {m:8.{p}f} {u} '
                  'std = {s:8.{p}f} {u}'.format(m=cgi.mr*unit_fact,
                                                s=cgi.devr*unit_fact,
                                                p=prec, u=units))

        if cgi.inters_bool is False:
            _tee(out, '\n  Gaussians too close for intersection calculation')
            _tee(out, '   --> Taking difference of mean values')

        _tee(out, '  CGI: dG = {dg:8.{p}f} {u}'.format(dg=cgi.dg*unit_fact,
                                                       p=prec, u=units))
        _tee(out, '  CGI: Std Err (bootstrap:parametric) = {e:8.{p}f} {u}'.format(e=cgi.err_boot1*unit_fact,
                                                                                  p=prec, u=units))

        if nboots > 0:
            _tee(out, '  CGI: Std Err (bootstrap) = {e:8.{p}f} {u}'.format(e=cgi.err_boot2*unit_fact,
                                                                           p=prec, u=units))

        if nblocks > 1:
            _tee(out, '  CGI: Std Err (blocks) = {e:8.{p}f} {u}'.format(e=cgi.err_blocks*unit_fact,
                                                                        p=prec, u=units))

    # --------------
    # Normality test
    # --------------
    if do_ks_test:
        print('\n  Running KS-test...')
        q0, lam00, check0, bOk0 = ks_norm_test(res_ab)
        q1, lam01, check1, bOk1 = ks_norm_test(res_ba)

        _tee(out, '    Forward: gaussian quality = %3.2f' % q0)
        if bOk0:
            _tee(out, '             ---> KS-Test Ok')
        else:
            _tee(out, '             ---> KS-Test Failed. sqrt(N)*Dmax = %4.2f,'
                      ' lambda0 = %4.2f' % (q0, check0))
        _tee(out, '    Reverse: gaussian quality = %3.2f' % q1)
        if bOk1:
            _tee(out, '             ---> KS-Test Ok')
        else:
            _tee(out, '             ---> KS-Test Failed. sqrt(N)*Dmax = %4.2f,'
                      ' lambda0 = %4.2f' % (q1, check1))

    # ========================
    # Bennett Acceptance Ratio
    # ========================
    if 'bar' in methods:
        _tee(out, '\n --------------------------------------------------------')
        _tee(out, '             Bennett Acceptance Ratio     ')
        _tee(out, ' --------------------------------------------------------')

        print('  Running Nelder-Mead Simplex algorithm... ')

        bar = BAR(res_ab, res_ba, T=T, nboots=nboots, nblocks=nblocks)
        if args.pickle:
            pickle.dump(bar, open("bar_results.pkl", "wb"))

        _tee(out, '  BAR: dG = {dg:8.{p}f} {u}'.format(dg=bar.dg*unit_fact, p=prec, u=units))
        _tee(out, '  BAR: Std Err (analytical) = {e:8.{p}f} {u}'.format(e=bar.err*unit_fact, p=prec, u=units))

        if nboots > 0:
            _tee(out, '  BAR: Std Err (bootstrap)  = {e:8.{p}f} {u}'.format(e=bar.err_boot*unit_fact, p=prec, u=units))
        if nblocks > 1:
            _tee(out, '  BAR: Std Err (blocks)  = {e:8.{p}f} {u}'.format(e=bar.err_blocks*unit_fact, p=prec, u=units))

        _tee(out, '  BAR: Conv = %8.2f' % bar.conv)

        if nboots > 0:
            _tee(out, '  BAR: Conv Std Err (bootstrap) = %8.2f' % bar.conv_err_boot)
        _tee(out, ' ------------------------------------------------------')

    # =========
    # Jarzynski
    # =========
    if 'jarz' in methods:
        _tee(out, '\n --------------------------------------------------------')
        _tee(out, '             Jarzynski estimator     ')
        _tee(out, ' --------------------------------------------------------')

        jarz = Jarz(wf=res_ab, wr=res_ba, T=T, nboots=nboots, nblocks=nblocks)
        if args.pickle:
            pickle.dump(jarz, open("jarz_results.pkl", "wb"))

        _tee(out, '  JARZ: dG Forward = {dg:8.{p}f} {u}'.format(dg=jarz.dg_for*unit_fact,
                                                                p=prec, u=units))
        _tee(out, '  JARZ: dG Reverse = {dg:8.{p}f} {u}'.format(dg=jarz.dg_rev*unit_fact,
                                                                p=prec, u=units))
        _tee(out, '  JARZ: dG Mean    = {dg:8.{p}f} {u}'.format(dg=jarz.dg_mean*unit_fact,
                                                                p=prec, u=units))
        if nboots > 0:
            _tee(out, '  JARZ: Std Err Forward (bootstrap) = {e:8.{p}f} {u}'.format(e=jarz.err_boot_for*unit_fact,
                                                                                    p=prec, u=units))
            _tee(out, '  JARZ: Std Err Reverse (bootstrap) = {e:8.{p}f} {u}'.format(e=jarz.err_boot_rev*unit_fact,
                                                                                    p=prec, u=units))

        if nblocks > 1:
            _tee(out, '  JARZ: Std Err Forward (blocks) = {e:8.{p}f} {u}'.format(e=jarz.err_blocks_for*unit_fact,
                                                                                    p=prec, u=units))
            _tee(out, '  JARZ: Std Err Reverse (blocks) = {e:8.{p}f} {u}'.format(e=jarz.err_blocks_rev*unit_fact,
                                                                                    p=prec, u=units))
        _tee(out, ' ------------------------------------------------------')

    print '\n   Plotting histograms......'
    if 'cgi' in methods and args.cgi_plot is not None:
        make_cgi_plot(args.cgi_plot, res_ab, res_ba, cgi.dg, cgi.err_boot1,
                      args.nbins, args.dpi)

    print('\n   ......done...........\n')

    if args.pickle:
        print('   NOTE: units of results in pickled files are as in the\n'
              '   provided dgdl.xvg or integ.dat files. These are typically\n'
              '   in kJ/mol when using dgdl.xvg files from Gromacs.\n')
    # execution time
    etime = time.time()
    h, m, s = time_stats(etime-stime)
    print("   Execution time = %02d:%02d:%02d\n" % (h, m, s))


if __name__ == '__main__':
    args = parse_options()
    main(args)
