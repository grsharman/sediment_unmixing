# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 16:40:20 2017

With DZ data, people have come up with lots of comparison metrics to highlight differences in complex empirical distributions.
One challenge is that not all these lend themselves nicely to formal derivations of significance measures.
Fear not, numerical methods to the rescue! The permutation approach allows us to come up with a significance
measure for any choice of comparison function. If you are a person who likes 'Similarity' tests, ala Gehrels, 2000,
use that. If you're more into cross-correlation's of PDPs, ala Saylor and Sundell, 2016, that works too.

Conceptually the test works like this;

We have two groups of data, our null-hypothesis is that those data were drawn from the same distribution. If our
null-hypothesis were true, than the labelling we assigned to put them into those two groups was arbitrary. In other words,
the 'real' data labelling should result in no more of an extreme value of our comparison function (e.g. the D-max from
the KS-test) than any other arbitrary labelling. Slick.

This means that we can look at the distribution of our comparison metric with a bunch (e.g. 10,000+) randomly shuffled labels.
If the null hypothesis is true, our originial labelling won't be particularly different than most of these. If the null
hypothesis is false, our original labelling will be quite different than the randomly shuffled labels (because those
labels did define meaningful groups).

I will leave what determines 'particularly' different up to you, but for the classic p < 0.05 we would want our
comparison metric to fall outside the 5th or 95th percentile of shuffled label comparison metrics.
I say outside the 5th OR 95th, because some comparison metrics yield large values when two DZ datasets are similar
(e.g. cross correlation of PDFs) while others yield small values (e.g. Dmax).

As a warning, this implementation is particularly slow because for each permutation there is a considerable overhead. I
regenerate two detritalPopulation instances each iteration, and recalculate PDFs and CDFs for them. Not all comparison
functions require this, but because some do it was simpler to implement it this way. Folks who want to do a lot of these
tests for a specific metric may want to develop a less general approach that doesn't make all these (often unnecessary)
calculations.

This implements a class to compute and store two sample comparisons using a permutation approach. I came across this idea
on a post of Larry Wasserman's:

https://normaldeviate.wordpress.com/2012/07/14/modern-two-sample-tests/

"""

__author__ = "Sam Johnstone"
__copyright__ = "2016"
__credits__ = ["Sam Johnstone", "Glenn Sharman", "Lauren Shumaker"]
__license__ = "MIT"
__maintainer__ = "Sam Johnstone"
__email__ = "sjohnstone@usgs.gov"
__status__ = "Production"

# ==============================================================================
# Import helpful shtuff
# ==============================================================================

import numpy as np
import csv  # For reading files
from matplotlib import pyplot as plt  # plotting functions
from matplotlib import cm  # Colormaps
from matplotlib import patches
from detritalPopulation import population
from scipy.stats import kstwobign


# ==============================================================================
# define a class for computing probabilites
# ==============================================================================

class twoSampleProbability:
    ''' Class for assessing the probability that two populations where drawn from
    the same distribution using permutations
    '''

    def __init__(self, pop1, pop2, objFunc, areLargeValuesMoreDifferent, minAge=0.0, maxAge=4500.0, dt=1.0, nIters=5000):
        '''

        :param pop1: an instance of detrital population
        :param pop2: a different instance of detrital population to compare to pop1
        :param objFunc: one of the comparison function defined in populationMetrics.py (or a user defined function).
        objFunc needs to be a function of two variables (e.g. objFunc(pop1,pop2) that returns a single value comparing the
        two.
        :param areLargeValuesMoreDifferent: boolean, are large or small values 'good matches' between datasets?
        :param minAge: the minimum age axis to recompute pdfs/cdfs on for objFuncs that require those values
        :param maxAge: the maximum age axis to recompute pdfs/cdfs on for objFuncs that require those values
        :param dt: the age axis spacing to recompute pdfs/cdfs on for objFuncs that require those values
        :param nIters: the number of permutations to perform
        '''

        # Store the nunmber of iterations
        self.nIters = nIters

        # Store the directionalite of the objective function
        self.areLargeValuesMoreDifferent = areLargeValuesMoreDifferent

        # Concatenate ages and errors
        allAges = np.hstack((pop1.ages, pop2.ages))
        allErrors = np.hstack((pop1.errors, pop2.errors))

        # Create indexing variable for ages and errors
        indices = np.hstack((np.zeros_like(pop1.ages), np.ones_like(pop2.ages)))

        # Preallocate space for permutation results
        self._permVals = np.zeros(nIters)

        #What is the value of the true labelling?
        self._thisVal = objFunc(pop1, pop2)

        # Iterate through the requested number of permutations
        for i in range(nIters):
            # Shuffle the indices
            permIndices = np.random.permutation(len(indices))
            shuffledIndices = indices[permIndices]

            # Create new populations
            pop1_i = population(ages=allAges[shuffledIndices == 0], errors=allErrors[shuffledIndices == 0],tmin=minAge, tmax=maxAge, delt=dt)
            pop2_i = population(ages=allAges[shuffledIndices == 1], errors=allErrors[shuffledIndices == 1],tmin=minAge, tmax=maxAge, delt=dt)

            # # Make sure the axes match
            # pop1_i.calcDF(forceCalc=True, tmin=minAge, tmax=maxAge, delt=dt)
            # pop1_i.calcCDF(forceCalc=True, tmin=minAge, tmax=maxAge, delt=dt)
            #
            # pop2_i.calcDF(forceCalc=True, tmin=minAge, tmax=maxAge, delt=dt)
            # pop2_i.calcCDF(forceCalc=True, tmin=minAge, tmax=maxAge, delt=dt)

            # Calculate this objectiveFunctionValue
            self._permVals[i] = objFunc(pop1_i, pop2_i)

        # Calculate the probability of this value
        self.pVal = self.getPforVal(self._thisVal)

    def getPforVal(self, Val):
        ''' Lookup the probability associated with a specific value
        '''

        if self.areLargeValuesMoreDifferent:
            nMoreDiff = np.sum(self._permVals > Val)
        else:
            nMoreDiff = np.sum(self._permVals < Val)

        return float(nMoreDiff) / float(self.nIters)

    def getValForP(self, p):
        ''' What value is associated with some critical probability
        :param p: the confidence value of interese (e.g. 0.05 for 5% confidence)
        '''

        if self.areLargeValuesMoreDifferent:
            return np.percentile(self._permVals, (1.0-p) * 100.0)
        else:
            return np.percentile(self._permVals, p*100.0)

    def plot_frequencyValue(self, **kwargs):
        ''' Create a plot of the p value as a function of different objective function values
        '''

        psToPlot = np.logspace(-3, 0, 100)

        vals = np.array([self.getValForP(p) for p in psToPlot])

        plt.plot(vals, psToPlot, '-o', **kwargs)

        plt.xlabel('Function value', fontsize=14)
        plt.ylabel('$p$', fontsize=14)

    def plotCDF(self,**kwargs):
        '''
        Create a plot of the distribution of the permuted comparison values
        :param kwargs:
        :return:
        '''

        #Determine the x-axis
        minVal = np.min(self._permVals)
        maxVal = np.max(self._permVals)
        valsToPlot = np.linspace(minVal,maxVal,200)

        #Preallocate
        CDF = np.zeros_like(valsToPlot)

        #Loop through the values to calculate CDFs at
        for i,val in enumerate(valsToPlot):
            CDF[i] = np.float(np.sum(self._permVals >=val))/np.float(self.nIters)

        #plot the data
        plt.plot(valsToPlot,CDF,'-',**kwargs)
        plt.plot([self._thisVal, self._thisVal], [0, 1], '--k', label = r'$p =$ %.1e'%self.pVal)
        plt.ylabel('Fraction of larger values',fontsize = 13)
        plt.legend(loc = 'best')