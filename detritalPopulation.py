#!/usr/bin/env python

''' A class to store populations of detrital zircon U-Pb ages, and some helpful
functions for plotting them.


Dependencies:

 numpy
 csv
 pandas

'''

__author__ = "Sam Johnstone"
__copyright__ = "2016"
__credits__ = ["Sam Johnstone", "Glenn Sharman", "Lauren Shumaker"]
__license__ = "MIT"
__maintainer__ = "Sam Johnstone"
__email__ = "sjohnstone@usgs.gov"
__status__ = "Testing"

#==============================================================================
# Import helpful libraries
#==============================================================================

import numpy as np
import csv #For reading files
from matplotlib import pyplot as plt #plotting functions

import pandas as pd #Used for importing excel files


#==============================================================================
# Definitition of a class, population, that contains some useful methods
#==============================================================================

class population:
    '''Population is a class representing a suit of dz analysis that forms some
    cohesize unit. Really just a glorified wrapper for  numpy arrays of age and error, with some
    built in functions for plotting.

    A population may be initialized in a variety of ways:

    population(ages = X,errors = Y) , where ages and errors are numpy arrays

    population(filename = X), where filename specifies the path to a file of ages and errors.  This has a set of optional
    arguments and defaults: delimiter = ',',ageHeader = 'age',errorHeader = 'error'

    population(parentPopulations, mixingCoefficients), mixes a list of instances of parent populations based on the
    list of mixing coefficients.

    population(ageAxis = X,PDF = Y), where ageAxis and PDF are numpy arrays that describe the probabilities along the age axis

    population(populationToCopy = X), where populationToCopy is another population

    population(excelFileName = X, excelSheetName = Y, ageHeader = Z, errorHeader = A), where the arguments describe the
    path to an excel file and where the data of interest is stored within that file.  There are also two optional arguments
    sampleIDfield=None,sampleID = None . If all the age, error data is stored in a single sheet, an additional column can be
    used to specify the ID or the IDs of the samples of interest.


    '''
    
    required_inputs_and_actions = ((('ages','errors',) ,'_createFromArrays_'),
                                   (('filename',),'_loadFromFile_'),
                                (('ageAxis','PDF',), '_createFromPDF_'),
                                 (('parentPopulations','mixingCoefficients',), '_createFromParentMixture_'),
                                (('populationToCopy',), '_copyFrom_'),
                                (('excelFileName','excelSheetName','ageHeader','errorHeader',), '_loadFromExcelFile_'),)  
    
    def __init__(self,*args,**kwargs):
        '''

        Initialize a population.

        A population may be initialized in a variety of ways:

        population(ages = X,errors = Y) , where ages and errors are numpy arrays

        population(filename = X), where filename specifies the path to a file of ages and errors.  This has a set of optional
        arguments and defaults: delimiter = ',',ageHeader = 'age',errorHeader = 'error'

        population(ageAxis = X,PDF = Y), where ageAxis and PDF are numpy arrays that describe the probabilities along the age axis

        population(populationToCopy = X), where populationToCopy is another population

        population(excelFileName = X, excelSheetName = Y, ageHeader = Z, errorHeader = A), where the arguments describe the
        path to an excel file and where the data of interest is stored within that file.  There are also two optional arguments
        sampleIDfield=None,sampleID = None . If all the age, error data is stored in a single sheet, an additional column can be
        used to specify the ID or the IDs of the samples of interest.

        :param args:
        :param kwargs:
        '''
        
        self.ages = None
        self.errors = None
        self.n = None
        
        self._isDFCalcd_ = False
        self._isCDFCalcd_ = False
        self._CDFMethod_ = None
    
        if len(kwargs.keys()) == 0:
            return
        
        evaluative_action = self.__get_evaluative_action(*args, **kwargs)
        
        if evaluative_action == None:
            print('Input Error, Inputs not satisfied')
        
        eval('self.' + evaluative_action + '(*args, **kwargs)')
            
    def __get_evaluative_action(self, *args, **kwargs):
                
        for required_input_set, evaluative_action in self.required_inputs_and_actions:
            these_kw = set(kwargs.keys())
            required_kw = set(required_input_set)
            if required_kw.issubset(these_kw):
                return evaluative_action 
        
        return None
    
    def _createFromArrays_(self,ages,errors,**kwargs):
        '''Create the population by specifying ages and errors
        '''
        self.ages = ages
        self.errors = errors
        self.n = len(self.ages)  
        
        self.calcCDF(**kwargs)
        self.calcDF(**kwargs)
        self._isCDFCalcd_ = True
        self._isDFCalcd_ = True
        
    def _loadFromFile_(self,filename,delimiter = ',',ageHeader = 'age',errorHeader = 'error',**kwargs):
        '''Load the ages and errors from a text file
        '''
        ages = []
        errors = []
        
        #Open the file
        with open(filename,'r') as f:
            headers = f.readline().rstrip().split(delimiter)
            ageHeaderMatch = [hdr.lower()==ageHeader.lower() for hdr in headers]
            errorHeaderMatch = [hdr.lower()==errorHeader.lower() for hdr in headers]

            if (not any(ageHeaderMatch)):
                print('Whoopsy, couldn\'t find the appropriate age header')
                
            if (not any(errorHeaderMatch)):
                print('Whoopsy, couldn\'t find the appropriate error header')

            ageIdx = np.where(ageHeaderMatch)[0][0]
            errorIdx = np.where(errorHeaderMatch)[0][0]
            
            reader=csv.reader(f,delimiter=delimiter)
            for line in reader:
                errors.append(float(line[errorIdx]))
                ages.append(float(line[ageIdx]))
        
        self.ages = ages
        self.errors = errors
        self.n = len(ages)
        
        self.calcCDF(**kwargs)
        self.calcDF(**kwargs)
        self._isCDFCalcd_ = True
        self._isDFCalcd_ = True
        
    def _loadFromExcelFile_(self,excelFileName,excelSheetName = None,ageHeader = 'age',errorHeader = 'error',sampleIDfield=None,sampleID = None,**kwargs):
        ''' Load the ages and errors from an excel file
        '''
        
        #Load in the excel file using pandas
        xlFile = pd.read_excel(open(excelFileName,'rb'),sheetname = excelSheetName)
        
        #If requested, take a subset of the samples
        useSampleSubset = False
        if not(sampleIDfield is None) and not(sampleID is None):
            useSampleSubset = True
    
        #Index out just the requested data if necessary
        if useSampleSubset:
            ages = xlFile[ageHeader]
            errors = xlFile[errorHeader]            
            
            #If the user just specified a single string or value, convert it to a list
            if (isinstance(sampleID,str)) or (isinstance(sampleID,float)) or (isinstance(sampleID,int)):
                sampleID = [sampleID]
            
            #Get the field of the sample ids
            sampleNames = xlFile[sampleIDfield]            
            
            #Initialize a boolean index
            goodData = np.zeros_like(ages) == 1            
            for name in sampleID:
                goodData = goodData | (sampleNames == name)
            
            self.ages = np.array(ages[goodData])
            self.errors = np.array(errors[goodData])
        else:
            #Grab the specific data from the excel file
            self.ages = np.array(xlFile[ageHeader])
            self.errors = np.array(xlFile[errorHeader])

        # Clear the excel file from memory
        xlFile = None

        self.n = len(self.ages)
            
        #Calculate derivatives of the data
        self.calcCDF(**kwargs)
        self.calcDF(**kwargs)


        
        
    def _createFromPDF_(self,ageAxis,PDF,nGrains=100,percError = 0.10,**kwargs): 
        '''Sample ages from a PDF, assign them a constant percent Error
        
        Could improve this by introducing some improved functionality to assign errors'''
        
        PDF = PDF/np.sum(PDF) #Ensure PDF sums to 1    

        self.ages = np.random.choice(ageAxis,nGrains,p=PDF)
        self.errors = self.ages*percError
        self.n = nGrains
        
        self.calcCDF(**kwargs)
        self.calcDF(**kwargs)
        self._isCDFCalcd_ = True
        self._isDFCalcd_ = True
        
    
    def _createFromParentMixture_(self,parentPopulations,mixingCoefficients):
        '''Create a population from a mixture of parent populations
        '''

        #Copy details from other parents
        self._copyFrom_(parentPopulations[0])
        
        #This population shouldn't have any 
        self.n = None
        self.ages = None
        self.errors = None
        
        mixedCDF = np.zeros_like(parentPopulations[0].CDF)
        mixedDF = np.zeros_like(parentPopulations[0].densFunc)
    
        #Mixing coefficients may not sum perfectly to 1, must correct
        #for this in PDFs
        mixCoeffSum = 0 
        for i,m in enumerate(mixingCoefficients):
            mixCoeffSum+=m
            mixedCDF += m*parentPopulations[i].CDF
            mixedDF += m*parentPopulations[i].densFunc
         
        #Assign these to the daughter
        self.CDF = mixedCDF/mixCoeffSum
        self.densFunc = mixedDF / mixCoeffSum
        self._isDFCalcd_ = True
        self._isCDFCalcd = True

    def _copyFrom_(self,populationToCopy):
        '''Copies data from an existing population
        '''
       
        self.__dict__.update(populationToCopy.__dict__)
    
    def plot(self,plotStyle = 'cdf',**kwargs):
        if plotStyle == 'cdf':
            self.plotCDF(**kwargs)
        elif plotStyle == 'pdf' or plotStyle == 'pdp' or plotStyle == 'kde':
            self.plotDF(**kwargs)
        else:
            print('Error: Unknown plotting style, options are: pdf,pdp,kde cdf')
    

    def plotCDF(self,**kwargs):
        '''Plot a cumulative probability distribution of the ages and errors associated with this population'''
        if not self._isCDFCalcd_:
            self.calcCDF(**kwargs)        
            
        plt.plot(self._tAxisCDF_,self.CDF,**kwargs)
        plt.xlabel('Age [Ma]', fontsize=14)

    def plotDF(self, **kwargs):
        '''Plot a probability distribution of the ages and errors associated with this population'''
        if not self._isDFCalcd_:
            self.calcDF(**kwargs)
        
        plt.plot(self._tAxisDF_, self.densFunc, **kwargs)
        plt.xlabel('Age [Ma]',fontsize = 14)
        plt.ylabel(r'$P$',fontsize = 14)
    
    def calcDF(self, forceCalc = False, method = 'pdp', **kwargs):
        '''
        Estimate the relative distribution of ages in the sample. Either through summing the error distributions of
        individual grains or constructing a kernel density estimate
        
        
        :param forceCalc: True or False, default is False. Forces the re-calculation of the density function
        :param method: either 'pdp' or 'kde'

        Optional arguments depend on methods and whether time axis of Density function is to be recalculated.

        Time axis can be manipulated with the following optional arguments:

        tmin= None, tmax = None, delt = None

        if method == \'kde\', the following arguments can be used to manipulate the kernel density estimate:

        :param bandwidth: The scaling term for the kernel in Myrs, default = 10. Larger numbers generate smoother estimates
        :param Kernel: optional, a function that operates on an array and returns an array of the same length


        '''       
        if self._isDFCalcd_ and (not forceCalc):
            pass
        else:    
            
            self._isDFCalcd_ = True
            
            self._getTimeAxisOfDF_(**kwargs)
            
            if method.lower() == 'pdp':
                self.densFunc = self._getPDP_()
            elif method.lower() == 'kde':
                self.densFunc = self._getKDE_(**kwargs)
            else:
                print('Whoopsy, valid method not detected. Please set method to either \'kde\' or \'pdp\'')

    def _getPDP_(self):
        '''
        
        Sum the individual grain probabilities and normalize
        
        :return: The Probability density plot as defined by the sum of the individual gaussian's associacted with each grain
        '''

        pdp = np.zeros_like(self._tAxisDF_)

        for i in range(self.n):  # iterate through the age, error data
            pdp += self.normDistribution(i)

        pdp = pdp / (np.sum(pdp * self._deltDF_))  # Ensure PDF integrates to 1
        
        return pdp

    def _getKDE_(self,**kwargs):
        '''

        Calculate the Kernel density estimate for a fixed bandwidth. Can optionally specify a kernel, but otherwise
        a gaussian kernel will be used.

        Optional parameters:
        :param bandwidth: The scaling term for the kernel in Myrs. Larger numbers generate smoother estimates
        :param Kernel: optional, a function that operates on an array and returns an array of the same length

        :return: The kernel density estimate
        '''

        if 'bandwidth' in kwargs:
            bandwidth = kwargs.get('bandwidth')
        else:
            bandwidth = 10.0
            print('No bandwidth selected, defaulting to %.0f'%bandwidth)

        if 'Kernel' in kwargs:
            Kernel = kwargs.get('Kernel')
        else:
            Kernel = lambda t: (1/np.sqrt(2.0*np.pi))*np.exp((-t**2)/2)

        KDE = np.zeros_like(self._tAxisDF_)

        #Iterate through the ages
        for age in self.ages:
            KDE+= Kernel((self._tAxisDF_ - age) / bandwidth)

        return (1/(self.n*bandwidth))*KDE

    def _getTimeAxisOfDF_(self, **kwargs):
        '''Get the time axis of the PDF'''

        if 'tmin' in kwargs:
            tmin = kwargs.get('tmin')
        else: #if minimum axis is not specified, get it
            minIdx = np.argmin(self.ages)
            tmin = self.ages[minIdx]-3*self.errors[minIdx]
        self._tminPDF_ = tmin
        
        if 'tmax' in kwargs:
            tmax = kwargs.get('tmax')

        else: #if maximum of axis is not specified, get it
            maxIdx = np.argmax(self.ages)
            tmax = self.ages[maxIdx]+3*self.errors[maxIdx]
        self._tmaxPDF_ = tmax
        
        if 'delt' in kwargs:
            delt = kwargs.get('delt')
        else:
            delt = (tmax-tmin)/1000.0
            
        self._deltDF_ = delt
        
        self._tAxisDF_ = np.arange(tmin, tmax, delt)
        
    
    def calcCDF(self,forceCalc = False, method = 'discrete',**kwargs):
        ''' calculate the cumulative density function of the distribution.
        
        Allows for one of two methods. Either method = 'discrete' for essentially
        a cumulative sum of a histogram, or method = 'integrated pdf' for a cumulative
        sum of the PDF
        Time axis can be manipulated with the following optional arguments:

        tmin= None, tmax = None, delt = None

        '''        
        method = method.lower()
        
        if not (method == 'discrete' or method == 'integrated pdf'):
            print('Error, method must be either \'discrete\' or \'integrated pdf\'')
            pass
        if self._isCDFCalcd_ and (not forceCalc):
            #Have we already caculated the CDF with this method?
            if self._CDFMethod_ == method:
                pass
        else:
            self._isCDFCalcd_ = True
            self._CDFMethod_ = method
            self._getTimeAxisOfCDF_(**kwargs)
            
            
            self.CDF = np.zeros_like(self._tAxisCDF_)
        
            if method == 'discrete':
                #CDF based on discrete ages
                for i in range(len(self._tAxisCDF_)):
                    self.CDF[i] = np.sum(self.ages <= self._tAxisCDF_[i])
                
                #Normalized CDF to sum to 1
                self.CDF = self.CDF*1.0/len(self.ages)
            elif method == 'integrated pdf':
                
                #CDF based on probabilities of  ages
                self.calcDF()
                self.CDF = np.cumsum(self.densFunc * self._deltDF_)
    
    def _getTimeAxisOfCDF_(self,**kwargs):
        ''' Get the time axis of the CDF, currently dividing this by using regular steps...'''

        if 'tmin' in kwargs:
            tmin = kwargs.get('tmin')
        else:  # if minimum axis is not specified, get it
            minIdx = np.argmin(self.ages)
            tmin = self.ages[minIdx] - 3 * self.errors[minIdx]
        self._tminPDF_ = tmin

        if 'tmax' in kwargs:
            tmax = kwargs.get('tmax')

        else:  # if maximum of axis is not specified, get it
            maxIdx = np.argmax(self.ages)
            tmax = self.ages[maxIdx] + 3 * self.errors[maxIdx]
        self._tmaxPDF_ = tmax

        if 'delt' in kwargs:
            delt = kwargs.get('delt')
        else:
            delt = (tmax - tmin) / 1000.0
        self._deltCDF_ = delt
        
        self._tAxisCDF_ = np.arange(tmin,tmax,delt)
    
    def normDistribution(self,i):
        '''Calculate the probabilities along the axis t, for the mean
        and error of the ith sample'''
        mean = self.ages[i]
        error = self.errors[i]
        
        return (1.0/np.sqrt(2.0*error**2*np.pi))*np.exp(-(self._tAxisDF_ - mean) ** 2 / (2.0 * error ** 2))

    def exportDensityFunctions(self,outputFileName,type = 'df',delimiter = ','):
        '''

        Write the density function and cumulative density function as a text file

        :param type: either 'df' for density function or 'cdf' for cumulative
        :param outputFileName: the name with which to save the file
        :param delimiter: default comma, the seperator for the age axis and probability rows
        '''


        if type.lower() == 'df':
            header = delimiter.join(('Age','Prob'))

            with open(outputFileName,'w') as f:
                f.write(header+'\n')

                for i in range(len(self.densFunc)):
                    line = delimiter.join((str(self._tAxisDF_[i]),str(self.densFunc[i])))
                    f.write(line+'\n')

        elif type.lower() == 'cdf':
            header = delimiter.join(('Age', 'Cum Prob'))

            with open(outputFileName, 'w') as f:
                f.write(header + '\n')

                for i in range(len(self.CDF)):
                    line = delimiter.join((str(self._tAxisCDF_[i]), str(self.CDF[i])))
                    f.write(line + '\n')

        else:
            print('Whoopsy, unrecognized type')