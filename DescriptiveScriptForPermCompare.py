'''
This is a script to highlight the use of PermutationCompare to compare two distributions, uses data and code presented
in Sharman and Johnstone, 2017 EPSL and Sickmann et al., 2016.
'''

#Import necessary libraries
import detritalPopulation as detPop
import PermutationCompare as perm
import populationMetrics as popMetrics
from matplotlib import pyplot as plt
import os

#Path to data
pathToExcelFile = os.path.join('TestData','Sickmann_etal_2016_Data.xlsx')

#What sheet stores the data in the spreadsheet?
excelSheet = 'U_Pb_Data'

#The name of the sheets where data is stored in that file
#Name of the rows in the exceel sheet that contain information
ageHeader = 'BestAge' #For ages
errorHeader = 'BestAge_err_1s_Ma' #For 1-sigma errors
sampleID = 'Sample_ID' #For sample names

#How do we want to scale the PDFs?
minAge = 0.0
maxAge = 300.0
dt = 1.0

#Read in a few different samples, 'DR584'
DR584 = detPop.population(excelFileName = pathToExcelFile,excelSheetName = excelSheet,ageHeader = ageHeader,
                        errorHeader = errorHeader,sampleIDfield=sampleID,sampleID = 'DR584')

#Sample T693, which mixture modelling of Sharman and Johnstone, 2017 suggests is pretty similar to DR584
T693 = detPop.population(excelFileName = pathToExcelFile,excelSheetName = excelSheet,ageHeader = ageHeader,
                        errorHeader = errorHeader,sampleIDfield=sampleID,sampleID = 'T693')

#Sample T978, which mixture modelling of Sharman and Johnstone, 2017 suggests is pretty different to DR584
T978 = detPop.population(excelFileName = pathToExcelFile,excelSheetName = excelSheet,ageHeader = ageHeader,
                        errorHeader = errorHeader,sampleIDfield=sampleID,sampleID = 'T978')

#Force the calculation of the PDFs with the specified axis - important that we also pass these parameters to the
#permutation comparison function (so that the comparison function is calculated on equivalent data)
DR584.calcDF(forceCalc = True, tmin = minAge, tmax = maxAge, delt = dt)
T693.calcDF(forceCalc = True, tmin = minAge, tmax = maxAge, delt = dt)
T978.calcDF(forceCalc = True, tmin = minAge, tmax = maxAge, delt = dt)
DR584.calcCDF(forceCalc = True, tmin = minAge, tmax = maxAge, delt = dt)
T693.calcCDF(forceCalc = True, tmin = minAge, tmax = maxAge, delt = dt)
T978.calcCDF(forceCalc = True, tmin = minAge, tmax = maxAge, delt = dt)

#Plot these distributions
plt.subplot(2,1,1)
DR584.plotDF(color = 'dodgerblue',linewidth = 1.25,label = 'DR584')
T693.plotDF(color = 'gold',linewidth = 1.25, label = 'T693')
T978.plotDF(color = 'lightcoral',linewidth = 1.25, label = 'T978')

plt.legend(loc = 'best')

####### Do some comparisons

#First compare the ones we expect to be similar
comparisonfunction = popMetrics.Dmax #Here we will pass the actual function (so no ending parenthesis)
areLargeValuesMoreDifferent = True #For Dmax, large values are more different
comp1 = perm.twoSampleProbability(DR584,T693,comparisonfunction,areLargeValuesMoreDifferent,
                                  minAge=minAge, maxAge=maxAge,dt=dt,nIters = 10000)

plt.subplot(2,2,3)
comp1.plotCDF(color = 'mediumseagreen',linewidth = 2,label = 'DR584 - T693')
plt.xlabel('D-max',fontsize = 14)

comparisonfunction = popMetrics.correlationCoeff #Here we will pass the actual function (so no ending parenthesis)
areLargeValuesMoreDifferent = False #For correlation coefficient, large values are more similar
comp2 = perm.twoSampleProbability(DR584,T978,comparisonfunction,areLargeValuesMoreDifferent,
                                  minAge=minAge, maxAge=maxAge,dt=dt,nIters = 10000)

plt.subplot(2,2,4)
comp2.plotCDF(color = 'darkorchid',linewidth = 2, label = 'DR584 - T978')
plt.xlabel('Cross Correlation of PDPs',fontsize = 14)
plt.ylabel('')


