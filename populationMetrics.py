#!/usr/bin/env python

'''
There are two classes of functions here. The mixture model class requires a function
that takes as input a daughter population, a list of parent populations, and a
list of mixture coefficients. This function returns some value that corresponds to
the quality of the fit (e.g. based on some statistical test or comparison function).
This ensures that all these tests are able to get all the inputs they need through the population class.

In some cases these functions will have optional additional arguments. In this case
these arguments can be fixed using pythons 'lambda' function.

Alternatively, I've also included some stripped down functions too that
just compare two sets of populations.

Dependencies:

 numpy

 scipy

 detritalPopulation


'''

__author__ = "Sam Johnstone"
__copyright__ = "2016"
__credits__ = ["Sam Johnstone", "Glenn Sharman", "Lauren Shumaker"]
__license__ = "MIT"
__maintainer__ = "Sam Johnstone"
__email__ = "sjohnstone@usgs.gov"
__status__ = "Production"

#==============================================================================
# Import helpful shtuff
#==============================================================================

import numpy as np
from detritalPopulation import population
from scipy.stats import ks_2samp


#==============================================================================
# Definition of utility functions for comparing two distributions
#==============================================================================

#==============================================================================
# Dmax of CDFs
#==============================================================================

def mixture_Dmax(daughterPop,parentPops,mixingCoefficients):
    '''Determine the Dmax value, the maximum difference between CDFs

    :param daughterPop: an instance of the population class describing a daughter population

    :param parentPops: a list of instances of the population class describing the parents

    :param mixingCoefficients: a list, or array, of fractional mixing coefficients that must sum to one
    
    ''' 
    mixedPop = population(parentPopulations = parentPops,mixingCoefficients = mixingCoefficients)
    
    return Dmax(daughterPop,mixedPop)
    
#==============================================================================
# Cross correlation coefficient (e.g. r^2) of pdfs
#==============================================================================
 
def mixture_correlationCoeff(daughterPop,parentPops,mixingCoefficients):
    
    '''Calculate the correlation coefficient (e.g. r^2) b/w the two pdfs.
    
    This method was proposed by Saylor and Sundell, 2016, Geosphere

    :param daughterPop: an instance of the population class describing a daughter population

    :param parentPops: a list of instances of the population class describing the parents

    :param mixingCoefficients: a list, or array, of fractional mixing coefficients that must sum to one
    '''
    
    mixedPop = population(parentPopulations = parentPops,mixingCoefficients = mixingCoefficients)

    return correlationCoeff(daughterPop,mixedPop)
#==============================================================================
# Similarity
#==============================================================================
def mixture_Similarity(daughterPop,parentPops,mixingCoefficients):
    '''Determine the \'similarity\' of the the mixed and daughtter DZ populations using the Gehrels, 2000 method
    
    This function is a wrapper to call calcSimilarity on the PDFs defined by 
    two population objects


    :param daughterPop: an instance of the population class describing a daughter population

    :param parentPops: a list of instances of the population class describing the parents

    :param mixingCoefficients: a list, or array, of fractional mixing coefficients that must sum to one
    '''
    
    mixedPop = population(parentPopulations = parentPops,mixingCoefficients = mixingCoefficients)
    
    return similarity(daughterPop,mixedPop)

#==============================================================================
# Likeness
#==============================================================================
def mixture_Likeness(daughterPop,parentPops,mixingCoefficients):
    '''Determine the \'likeness\' of the daughter and mixed parent populations.
    Here likeness is defined following Satkoski et al., 2013 as summarized in 
    Saylor and Sundell, 2016


    :param daughterPop: an instance of the population class describing a daughter population

    :param parentPops: a list of instances of the population class describing the parents

    :param mixingCoefficients: a list, or array, of fractional mixing coefficients that must sum to one
    
    ''' 
    mixedPop = population(parentPopulations = parentPops,mixingCoefficients = mixingCoefficients)
    
    return likeness(daughterPop,mixedPop)

#==============================================================================
# Vmax of Kuiper statistive
#==============================================================================

def mixture_Vmax(daughterPop,parentPops,mixingCoefficients):
    '''Calculate the V max value of the kuiper statistic as summarized.
    by Saylor and Sundell, 2016, Geosphere


    :param daughterPop: an instance of the population class describing a daughter population

    :param parentPops: a list of instances of the population class describing the parents

    :param mixingCoefficients: a list, or array, of fractional mixing coefficients that must sum to one
    '''
    mixedPop = population(parentPopulations = parentPops,mixingCoefficients = mixingCoefficients)
    return Vmax(daughterPop,mixedPop)
#==============================================================================
#
# Functions for direct comparisons of populations or distributions
#
#==============================================================================


def Dmax(pop1,pop2):
    '''return the maximum distance between the two CDFS contained within the
    populations

    :param pop1: an instance of the population class

    :param pop2: an instance of the population class

    '''
    if not isinstance(pop1,population):
        print('Error, both inputs must be population objects, first is not')
        pass
    
    if not isinstance(pop2,population):
        print('Error, both inputs must be population objects, second is not')
        pass

    if np.any(np.logical_not(pop1._tAxisCDF_ == pop2._tAxisCDF_)):
        print('Error, axes of CDFs do not match')
        pass
    
    return np.max(np.abs(pop1.CDF - pop2.CDF))


def likeness(pop1,pop2):
    
    '''Determine the \'likeness\' of the two populations.
    Here likeness is defined following Satkoski et al., 2013 as summarized in 
    Saylor and Sundell, 2016

    :param pop1: an instance of the population class

    :param pop2: an instance of the population class
    
    '''
    if not isinstance(pop1,population):
        print('Error, both inputs must be population objects, first is not')
        pass
    
    if not isinstance(pop2,population):
        print('Error, both inputs must be population objects, second is not')
        pass
    
    if np.any(np.logical_not(pop1._tAxisPDF_ == pop2._tAxisPDF_)):
        print('Error, age axes of PDFs do not match')
        pass
    
    return 1.0 - np.sum(np.abs(pop1.densFunc - pop2.densFunc))/2.0
    
def correlationCoeff(pop1,pop2):
    
    ''' Determine r^2 of the two PDFS. As proposed by Saylor and Sundell, 2016

    :param pop1: an instance of the population class

    :param pop2: an instance of the population class
    '''
        
    mean1 = np.mean(pop1.densFunc)
    mean2 = np.mean(pop2.densFunc)
    
    r_numerator = np.sum((pop1.densFunc - mean1)*(pop2.densFunc - mean2))
    r_denom = np.sqrt(np.sum((pop1.densFunc - mean1)**2))*np.sqrt(np.sum((pop2.densFunc - mean2)**2))
    
    return (r_numerator/r_denom)**2

def KSTest(pop1,pop2):
    '''Wrapper function to calculate the KS statistic from the CDFs defined by two populations.
    Following definition as provided in Saylor and Sundell 2016, Geosphere, Eqn. 2

    :param pop1: an instance of the population class

    :param pop2: an instance of the population class
    
   '''
    
    if not isinstance(pop1,population):
        print('Error, both inputs must be population objects, first is not')
        pass
    
    if not isinstance(pop2,population):
        print('Error, both inputs must be population objects, second is not')
        pass

    return ks_2samp(pop1.ages,pop2.ages).pvalue

def Vmax(pop1,pop2):
    '''Wrapper for calcVmax

    :param pop1: an instance of the population class

    :param pop2: an instance of the population class
    '''
    return calcVmax(pop1.CDF,pop2.CDF)

def calcVmax(CDF1,CDF2):
    '''Return the sum of the maximum distances (most positive & most negative)
    b/w the two CDFS, the V value of Saylor and Sundell 2016, Geosphere; Eqn. 6

    :param CDF1: a CDF

    :param CDF2: a different CDF... or is it the same?!
    '''
    return np.max(CDF1-CDF2) + np.max(CDF2-CDF1)


def similarity(pop1,pop2):
    '''Determine the \'similarity\' of the two DZ populations using the Gehrels, 2000 method
    
    This function is a wrapper to call calcSimilarity on the PDFs defined by 
    two population objects

    :param pop1: an instance of the population class

    :param pop2: an instance of the population class

    '''
    if not isinstance(pop1,population):
        print('Error, both inputs must be population objects, first is not')
        pass
    
    if not isinstance(pop2,population):
        print('Error, both inputs must be population objects, second is not')
        pass
    
    return calcSimilarity(pop1.densFunc,pop2.densFunc)
    
def calcSimilarity(PDF1,PDF2):
    '''Determine the \'similarity\' of the two DZ populations using the Gehrels, 2000 method
    
    This is the sum of the product to the two probability density functions

    :param PDF1: a PDF

    :param PDF2: a different PDF... or is it the same?!
    '''

    return np.sum(np.sqrt(PDF1*PDF2))