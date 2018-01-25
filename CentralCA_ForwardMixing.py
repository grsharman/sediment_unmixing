# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 09:31:32 2017

@author: sjohnstone, sharmang

Dependencies:

 numpy
 matplotlib
 os
 ternary : https://github.com/marcharper/python-ternary (For ternary plots)

 MixtureModel
 populationMetrics
 detritalPopulation

"""

import MixtureModel as mm #the class that performs mixture modelling
import populationMetrics as popMetrics #Functions for evaluating mixtures
import detritalPopulation as dp #contains 'population' the class describing detrital populations
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
import os

#==============================================================================
# Define the path to the data and the parameters for the mixture modeling
#==============================================================================

#Define the path to the data
filePath = os.path.join('TestData','Sickmann_etal_2016_Data.xlsx')

#Define the path to output the summary data
outfile = os.path.join('TestData','CA_mixingResults.txt')

#Name of the sheet that contains the data
excelSheet = 'U_Pb_Data'

#Name of the rows in the exceel sheet that contain information
ageHeader = 'BestAge' #For ages
errorHeader = 'BestAge_err_1s_Ma' #For 1-sigma errors
sampleID = 'Sample_ID' #For sample names

# Define parameters for the axis of CDFs and PDFs - NOTE THAT RESULTS ARE NOT INDEPENDENT OF THIS
minAge = 0.0
maxAge = 4501.0
tSpacing = 1.0

#Names of parents and daughters
ParentList =[['ANB','SNR','PAR'],'SAL', 'CAR']
ParentNames = ['AND,SNR,PAR','SAL','CAR']
DaughterList =['MAN','MLB','DR585','DR591','DR584','T755','MDP','T1135','GSB','T693','T978','T958'] 

#What function should be used to compare the data?
#objFunc = popMetrics.mixture_Likeness
#objFunc = popMetrics.mixture_Similarity
#objFunc = popMetrics.mixture_correlationCoeff #For the cross correlation of PDFs (e.g. Saylor and Sundell,2016)
objFunc = popMetrics.mixture_Dmax #For Dmax
#objFunc = popMetrics.mixture_Vmax #For Vmax

#Do 'good' mixtures get high values or low values with this objective function
mixOrder = 'normal' #Low values correspond to good mixtures (e.g. Dmax)
#mixOrder = 'reverse' #high values correspond to good mixtures (e.g. cross correlation)

#At what resolution to we want to create mixtures
mixSpacing = 0.01 #Keep coarse for testing (e.g. 0.1 (10%) ) to improve efficiency

# %%
#==============================================================================
# Load in the parents and daughters
#==============================================================================

parentSet = []
for parent in ParentList:
    parentSet.append(dp.population(excelFileName = filePath,excelSheetName = excelSheet,ageHeader = ageHeader,errorHeader = errorHeader,sampleIDfield=sampleID,sampleID = parent))
    
daughterSet = []
for daughter in DaughterList:
    daughterSet.append(dp.population(excelFileName = filePath,excelSheetName = excelSheet,ageHeader = ageHeader,errorHeader = errorHeader,sampleIDfield=sampleID,sampleID = daughter))

# %%
#==============================================================================
# Calculate the CDFs and PDFs for the parents and daughters
#==============================================================================
''' NOTE! when a population instance is created, the PDFs and CDFs are created by default
to highlight the data present in that population. This means that different populations
may have distributions calculated on different axes. Most comparisons of distributions
require the distributions to be calculated on matching axis - so we need to re-calculate distributions.

It is also important to recognize that the values of many comparison functions are sensitive to the axis they are calculated on.
That is, the same two sets of ages and errors could generate different values of a comparison/ objective function when calculated
on different time axes (e.g. differences in spacing or extent)
'''
#Plot the data
plt.figure()
plt.subplot(1,2,2)

#For different colormap options see the matplotlib colormap reference:
#https://matplotlib.org/examples/color/colormaps_reference.html
daughterCmap = 'Pastel1'
parentCmap = 'Pastel2'

cfun = cm.get_cmap(daughterCmap)
daughterColors = [cfun(i) for i in np.linspace(0,1,len(daughterSet))]

cfun = cm.get_cmap(parentCmap)
parentColors = [cfun(i) for i in np.linspace(0,1,len(daughterSet))]

for i,daughter in enumerate(daughterSet):
    daughter.calcDF(forceCalc = True, tmin = minAge, tmax = maxAge, delt = tSpacing)
    
    #The CDF also needs a method to know how to calculate that, either 'discrete' or 'integrated pdf'
    daughter.calcCDF(forceCalc = True, method = 'discrete', tmin = minAge, tmax = maxAge, delt = tSpacing)
    
    #plot this
    daughter.plotCDF(color = daughterColors[i],linewidth = 2, label = DaughterList[i])
    
plt.legend(loc = 'best')

plt.subplot(1,2,1)
for i,parent in enumerate(parentSet):
    parent.calcDF(forceCalc = True, tmin = minAge, tmax = maxAge, delt = tSpacing)
    
    #The CDF also needs a method to know how to calculate that, either 'discrete' or 'integrated pdf'
    parent.calcCDF(forceCalc = True, method = 'discrete', tmin = minAge, tmax = maxAge, delt = tSpacing)
    
    #plot this
    parent.plotCDF(color = parentColors[i],linewidth = 2, label = ParentNames[i])
    
plt.legend(loc = 'best')

# %%
#==============================================================================
# Create mixture models for all the daughters
#==============================================================================

mixSet = mm.mixtureModelSet(daughterSet,parentSet,objFunc,dFrac = mixSpacing,parentNames = ParentNames,daughterNames = DaughterList)
mixSet.sortMixtures(mixtureOrder = mixOrder)

#Write the results to a text files (will not overwrite existing files)
mixSet.outputResultsSummary(outfile)
# %%

############# PLOTTING THE RESULTS ##############

#==============================================================================
# Plot the results of the mixture models
#==============================================================================

mixSet.plotMixtureResults(plottingStyle='mixture value',cmap = 'RdBu')

#%%
#==============================================================================
# Plot the CDFs of the best mixtures
#==============================================================================

mixSet.plotMixtureResults(plottingStyle='best mixture')

#%%
#==============================================================================
# Plot the trends in the daughter populations, these will be drawn in the order they
# were loaded in above, so to highlight stratigraphic changes load in data in stratigraphic
# order
#==============================================================================

plt.figure()
mixSet.plotDaughterTrends()