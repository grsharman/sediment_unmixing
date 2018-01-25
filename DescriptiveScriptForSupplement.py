# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 09:58:56 2017
@author: sjohnstone
What follows is a script to walk through some different ways to manipulate detrital zircon data using the provided code.
Dependencies:
 numpy
 scipy
 matplotlib
 os
 populationMetrics
 detritalPopulation
 While not it this script, the mixture models will also require:
 ternary : https://github.com/marcharper/python-ternary (For ternary plots)
"""

#==============================================================================
# First, import the class to store detrital populations
#==============================================================================

from detritalPopulation import population #the class describing detrital populations
import numpy as np #For general data handling
import populationMetrics as popMetrics #For calculating some metrics comparing DZ data
from matplotlib import pyplot as plt #For plotting
import os
#%%
#==============================================================================
# Next we will explore various ways to initialize an instance of this class
# with data 
#==============================================================================

#First we will generate a population from an array of ages and errors,
#For this example we will create a bi-modal distribution of ages
#Create a function to generate ages at random 
def genRandomAges(n,mean,spread):
    ''' A function to return a selection of n random ages drawn from the population
    described by mean and spread
    '''
    return np.random.randn(int(n))*spread + mean

#Create a random population of ages with two components
#by concatenating two random arrays
n = 100 #number of ages
randAges = np.hstack((genRandomAges(n,30.0,5.0),genRandomAges(n,200.0,3.0)))
randErrors = randAges*0.10 #assign a constant relative error

#Create a population based on these ages and errors
popArray = population(ages = randAges, errors = randErrors)

#Next we will load in a population from a file
#This file must be a delimeted text file, must have headers titled 'age' and 'error'
# or header must be specified.
ageHdrName = 'ages'
errorHdrName = 'errors'
pathToFile = 'testData.csv'

#We will just create a text file to load in, or feel free to provide your own by commenting out the next few lines, and
#updating the above path to something of your own.

##### Starting file creation - comment out if providing your own.

#Again, we will just create some random distributions
n = 100 #number of ages
randAges = np.hstack((genRandomAges(n/2,75.0,5.0),genRandomAges(n/4,200.0,15.0),genRandomAges(n/4,250.0,15.0)))
randErrors = randAges*0.10 #assign a constant relative error

ageErrorArray = np.array([[age,error] for age,error in zip(randAges,randErrors)]) #Build the ages and errors into an n x 2 array

#Save this array as a delimited text file with a header, the comments = '' argument insures that nothing is appended before the header, default is a '#'
np.savetxt(pathToFile,ageErrorArray,fmt = '%.1f',delimiter=',',header= ageHdrName + ',' + errorHdrName,comments='')

##### Finished file creation - comment out if providing your own.


popLoad = population(filename = pathToFile, delimeter = ',', ageHeader = ageHdrName,
                     errorHeader = errorHdrName )

#Create a population by drawing grains from a PDF
#For illustration purposes, we will get the age axis and the
#actual density function from an existing population
popPDF = population(ageAxis = popLoad._tAxisDF_,PDF = popLoad.densFunc)
#Note that in a population the description of relative densities/probabilities
#is just referred to as a 'DF' or 'densFunc' - this can be created to be either a PDP or KDE

#Create a population by copying it from another population
popCopy = population(populationToCopy = popLoad)


#As a final example, we will create a population from a mixture of our previously
#created populations, before we do this it is important to first align the axis
#of the PDFs and CDFs
#Span of the age axis for density functinos
minAge = 0.0
maxAge = 2000.0
tSpacing = 1.0

#Force recalculation of density functions with same axis if making a mixture
popArray.calcDF(forceCalc = True, tmin = minAge, tmax = maxAge, delt = tSpacing)

#The `discrete' method specified a cdf calculated as a cumulative histogram, rather than
#an  `integrated pdf'
popArray.calcCDF(forceCalc = True, method = 'discrete',tmin = minAge,tmax = maxAge,delt = tSpacing)

popLoad.calcDF(forceCalc = True, tmin = minAge, tmax = maxAge, delt = tSpacing)
popLoad.calcCDF(forceCalc = True, method = 'discrete',tmin = minAge,tmax = maxAge,delt = tSpacing)

parentPops = (popArray,popLoad) #Create a tuple (or list) of the populations to be mixed
mixingCoeffs = (0.4,0.6) #create a tuple (or list) of the mixing coefficients
popMixture = population(parentPopulations = parentPops,mixingCoefficients = mixingCoeffs)

#%%
#==============================================================================
# Create some common plots
#==============================================================================

#plot the individual populations and mixed PDF
popArray.plotDF(color ='r', label ='Population 1')
popLoad.plotDF(color ='b', label ='Population 2')
popMixture.plotDF(color ='k', label ='Mixed population')
plt.xlim(0,500) #Adjust the x axis

plt.legend(loc = 'upper right')

#%%
#==============================================================================
#Alternatively, create KDE plots
#==============================================================================

plt.figure()
#In creating KDEs its important to specify an appropriate bandwidth.

'''
Optional you can also suppy a custom kernel as Kernel = fun(x), where fun(x) is
is the appropriate function of your choosing. Here x is expected to be:
 x = (ageAxis - observedAge)/bandwidth
Such that the KDE is then calculated as
    KDE = np.zeros(len(timeAxis))
    for age in self.ages:
       KDE += fun((timeAxis - age) / bandwidth)
'''
#Force the recalculation of density functions with a gaussian kernel
popLoad.calcDF(forceCalc = True, method = 'kde',bandwidth = 5.0,tmin = minAge,tmax = maxAge,delt = tSpacing)
popArray.calcDF(forceCalc= True, method = 'kde',bandwidth = 5.0,tmin = minAge,tmax = maxAge,delt = tSpacing)

#Create the mixture
parentPops = (popArray,popLoad) #Create a tuple (or list) of the populations to be mixed
mixingCoeffs = (0.25,0.75) #create a tuple (of list) of the mixing coefficients
popMixture = None #Clear the old mixture
popMixture = population(parentPopulations = parentPops,mixingCoefficients = mixingCoeffs)

#plot the individual populations and mixed PDF
popArray.plotDF(color ='r', label ='Population 1')
popLoad.plotDF(color ='b', label ='Population 2')
popMixture.plotDF(color ='k', label ='Mixed population')
plt.xlim(0,500) #Adjust the x axis

plt.legend(loc = 'upper right')

#%%
#==============================================================================
# A similar approach can be taken with excel files
#==============================================================================

#The path to an excel file
#If specifying a relative path, be sure that you are running the code in the
#directory required for that relative path
pathToExcelFile = os.path.join('TestData','Sickmann_etal_2016_Data.xlsx')

#What sheet stores the data in the spreadsheet?
excelSheet = 'U_Pb_Data'

#The name of the sheets where data is stored in that file
#Name of the rows in the exceel sheet that contain information
ageHeader = 'BestAge' #For ages
errorHeader = 'BestAge_err_1s_Ma' #For 1-sigma errors
sampleID = 'Sample_ID' #For sample names

#Names of parents and daughters
SampleNames =['ANB','SNR', 'PAR']

#Load in each of the files
allPops = [] #Create an empty list to store the populations in
for name in SampleNames:
 #Read in a set of data from the excel spreadsheet and add it to this list
 allPops.append(population(excelFileName = pathToExcelFile,excelSheetName = excelSheet,ageHeader = ageHeader,errorHeader = errorHeader,sampleIDfield=sampleID,sampleID = name))


# We can then iterate through this list and update properties of the samples
for i in range(len(allPops)):
    allPops[i].calcDF(forceCalc=True,method = 'kde',bandwidth = 10.0, tmin=minAge, tmax=maxAge, delt=tSpacing)

    # The CDF also needs a method to know how to calculate that, either 'discrete' or 'integrated pdf'
    # It is also necessary to calculate the CDF to perform some comparison tests
    allPops[i].calcCDF(forceCalc=True, method='discrete', tmin=minAge, tmax=maxAge, delt=tSpacing)

#Can also iterate through the data and plot them
for i in range(len(allPops)):
    allPops[i].plotDF(linewidth = 2, label = SampleNames[i])

plt.legend(loc = 'upper right')

#%%
#==============================================================================
# We can also utilize the functions within popmetrics to compare distributions
#==============================================================================

#Compare two samples with the KS test
pKS = popMetrics.KSTest(popArray,popLoad)
print('Acoording to the KS-test, the probability of the null hypothesis - that these populations where drawn from the same'
      ' distibution - is %.2e'%pKS)

#Compare two samples with r2 pdp
#With this metric its important to select the age axis appropriately and select the appropriate density function
popArray.calcDF(forceCalc = True,method = 'pdp', tmin = 0, tmax = 300, delt = 2.0)
popLoad.calcDF(forceCalc = True, method = 'pdp',tmin = 0, tmax = 300, delt = 2.0)

r2pdp = popMetrics.correlationCoeff(popArray,popLoad)
print('The correlation coefficient between the pdps is %.4f' %r2pdp)

#Or how about the 'similarity'?
sim = popMetrics.similarity(popArray,popLoad)
print('The \'similarity\' between the pdps is %.4f' %sim)

#%%
#==============================================================================
# Finally, we can also export these distributions
#==============================================================================

popArray.exportDensityFunctions('test_df.txt',type = 'df',delimiter=',')
popArray.exportDensityFunctions('test_cdf.txt',type = 'cdf',delimiter=',')
