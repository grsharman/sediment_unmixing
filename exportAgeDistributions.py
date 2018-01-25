'''
Created on Mon Jun 26 09:32:14 2017

@author: sharmang, sjohnstone

Script for exporting age distributions to the format expected by Paterson and Heslop (2015)'s 'AnalySize'
'''

import detritalPopulation as dp #contains 'population' the class describing detrital populations
import os

#This example relies on an excel spreadsheet holding all the data, for alternatives for loading data
#consult the 'descriptiveScriptForSupplement.py'

#Define the path to the data
filePath = os.path.join('TestData','Sickmann_etal_2016_Data.xlsx')

densityFunctionMethod = 'pdp' #Options: 'kde' or 'pdp'
kdeBandwidth = 5.0 #Required if densityFunctionMethod = 'kde'

#Define the path to output the summary data
outfile = os.path.join('TestData','AnalysizeDistributions.txt')
outDelimiter = ','

agePrecision = 0 #Number of digits after decimal to report for ages
relProbPrecision = 6 #Number of digits afterdecimal to report for relative probabilities

#Create formating strings for numbers
ageFmt = '%.'+str(agePrecision)+'f'
pFmt = '%.'+str(relProbPrecision)+'f'

#Name of the sheet that contains the data
excelSheet = 'U_Pb_Data'

#Name of the rows in the Excel sheet that contain information
ageHeader = 'BestAge' #For ages
errorHeader = 'BestAge_err_1s_Ma' #For 1-sigma errors
sampleID = 'Sample_ID' #For sample names

# Define parameters for the axis of distributions - NOTE THAT RESULTS ARE NOT INDEPENDENT OF THIS
minAge = 0.0
maxAge = 3501.0
tSpacing = 2.0

#Samples to export
sampleList =['MAN','MLB','DR585','DR591','DR584','T755','MDP','T1135','GSB','T693','T978','T958']
with open(outfile,'w') as f:
    #Iterate through all the samples
    for i,sample in enumerate(sampleList):
        print(sample)
        pop = dp.population(excelFileName = filePath,excelSheetName = excelSheet,ageHeader = ageHeader,errorHeader = errorHeader,sampleIDfield=sampleID,sampleID = sample)
        pop.calcDF(forceCalc = True, tmin = minAge, tmax = maxAge, delt = tSpacing,method = densityFunctionMethod,bandwidth = kdeBandwidth)

        #Before first sample, write the header
        if i == 0:
            ts = pop._tAxisDF_
            tStrList = [' '+ageFmt%t for t in ts]
            header = 'Age' + outDelimiter + outDelimiter.join(tStrList) + '\n'
            f.write(header)

        #Convert density functions to strings and write them
        DF = pop.densFunc
        dfStrList = [' '+pFmt%p for p in DF]
        line = sample + outDelimiter + outDelimiter.join(dfStrList) + '\n'
        f.write(line)