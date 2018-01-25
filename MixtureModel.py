#!/usr/bin/env python

''' A class to perform mixture modelling of detrital zircon data. Constructs forward
mixtures based on defined mixing coefficients. Uses special 'mixing functions' to 
evaluate the goodness of fit of these mixtures. Examples of these 'mixing functions'
 are provided in 'populationMetrics'.
 Dependencies:
 numpy
 matplotlib
 os
 ternary : https://github.com/marcharper/python-ternary (For ternary plots)
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
import ternary
from matplotlib import pyplot as plt #plotting functions
from matplotlib import cm #Colormaps
from matplotlib import patches
from detritalPopulation import population
import os

#==============================================================================
# Create helper functinos
#==============================================================================
def getMixes(nParents, nSteps):
    """
    Computes mixture coefficients at regular intervals.
    
    Parameters
    ----------
    nParents : scalar, the number of parents     
    nSteps : Step interval over which to create the mixture
        (e.g., 1001 steps = 0, 0.001, 0.002, . . . 1.000)
        
    Returns
    -------
    mixes : An array of the mixture coefficients
    
    Notes
    -----
    Created by Jonathan Sharman
    """
    mixes = getMixesHelper(nParents, nSteps, nSteps)
    return np.asarray(mixes)

def getMixesHelper(nParents, nSteps, maxSteps):
    """
    Helper function for getMixes().
    
    Parameters
    ----------
    nParents : scalar, the number of parents        
    nSteps : Step interval over which to create the mixture
        (e.g., 1001 steps = 0, 0.001, 0.002, . . . 1.000)
    maxSteps : Internal parameter used in the recursive function
    
    Returns
    -------
    mixes : An array of the mixture coefficients
    
    Notes
    -----
    Created by Jonathan Sharman
    """
    if (nParents == 1):
        return [[(maxSteps - 1.) / (nSteps -1.)]]
    else:
        mixes = []
        for i in range(maxSteps):
            submixes = getMixesHelper(nParents - 1, nSteps, maxSteps - i)
            for j in range(len(submixes)):
                mixes.append([i / (nSteps - 1.)] + submixes[j])
    return mixes  

#==============================================================================
# Functions for mixing models
#==============================================================================

class MixtureModel:
    def __init__(self,daughterPop,parentPops,objFunc,dFrac = 0.05,parentNames = None,daughterName = 'Observed daughter',mixtureOrder = 'normal'):
        '''
        A class that contains the data and performs the operations necessary to conduct forward, or 'top-down' mixtures
        of empirically defined parents, comparing these to an empirically defined daughter population.
        :param daughterPop: an instance of the population class describing a daughter population
        :param parentPops: a list of instances of the population class describing the parents
        :param objFunc: A function of three variables. The first two are an instance of the population class that describes
        the daughter and a list of instances of population classes describing the daughters. The second is a list or array of
        floats describing the mixing coefficients that sum to 1.
        :param dFrac: A float describing the precision/ spacing with which to construct mixtures. E.g. 0.05 would be looking
        at all mixtures at 5% intervals. Reduction in this term will exponentially increase the computation time.
        :param parentNames: Optionally, a list of strings of the parent names. Used in plotting results
        :param daughterName: Optionally, the name of the daughter population. Used in plotting results.
        :param mixtureOrder: Determines whether the objective function defines good mixtures as high values or low values.
        If mixtureOrder = 'normal', small values of the objective function define good mixtures (e.g. \'Dmax\').
        If mixtureOrder = 'reverse', high values are good mixtures (e.g. \'similarity\' and \'r2 pdp\')
        '''
        #First things first - check that distribution axes match
        unmatchedAxis = False        
        for parent in parentPops:
            unmatchedAxis = unmatchedAxis or (np.any(daughterPop._tAxisCDF_
                != parent._tAxisCDF_)) or (np.any(daughterPop._tAxisDF_
                != parent._tAxisDF_))
                
        if unmatchedAxis:
            print('Whoopsy, the distributions need to be calculated with matching axes.')
            return None
        else:
            self.dFrac = dFrac
            self.objFunc = objFunc
            self.areMixingCoeffsSorted = False
            self.nParents = len(parentPops)
            if parentNames is None:
                self.parentNames = ['P'+str(x) for x in range(self.nParents)]
            elif len(parentNames) < self.nParents:
                print('Whoops, fewer parents names than parents - proceeding with integer names \n')
                self.parentNames = ['P'+str(x) for x in range(self.nParents)]
            else:
                self.parentNames = parentNames
                
            self.daughterName = daughterName
            
            self.simulateMixtures(daughterPop,parentPops,objFunc,dFrac,mixtureOrder)
            
            self.parentPops = parentPops
            self.daughterPop = daughterPop

    def simulateMixtures(self,daughterPop,parentPops,objFunc, dFrac = 0.05,mixtureOrder = 'normal'):
        '''Compares the daughter population to combinations of the parent populations
        present in the tuple of populations 'parentPops' using objFunc.
        
        objFunc is a function that must take as input the daughter population,
        parent populations, a mixing coefficients for the parents. It must return
        some value representing the goodness of fit
        
        Mixture coefficients will be spaced dFrac apart.
        '''
        
        self.getMixingCoefficients(dFrac) #Construct mixture coefficients
        self.mixtureFunctionValue = np.zeros(self.nTotalMixes) #Preallocate space for mixture vales
        
        for i in range(self.nTotalMixes): #Iterate through each mixture
            self.mixtureFunctionValue[i] = objFunc(daughterPop,parentPops,self.mixingCoeffs[i])
            
        self.sortMixtures(mixtureOrder)
        
    def exportResults(self,filename):
        ''' Exports the values of the objective function for each of the mixture values tested as a comma-
        delimited file
        '''
        if os.path.isfile(filename):
            print('Whoopsy, that filename already exists, please choose another')
        else:
            #Open the file for writing
            with open(filename,'w') as f:
                #instead of mixtureFunctionValue could use self.objFunc.func_name
                #Create a header line for the file and write that line
                header = ','.join(self.parentNames)+', mixtureFunctionValue \n'
                f.write(header)
                #Iterate through the values tested
                for i,coeffs in enumerate(self.mixingCoeffs):
                    #Write values to a string, and those to the file
                    strVals = ['%.3f'%self.mixingCoeffs[i][j] for j in self.nParents]
                    line =  ','.join(strVals) + ',' + self.mixtureFunctionValue[i]
                    f.write(line)
                    
    def sortMixtures(self,mixtureOrder):
        ''' Sort the mixture function values and mixture coefficients
        mixtureOder may be either \'normal\' or \'reverse\' if good matches
        have high values or low values        
        '''
        sortIdcs = np.argsort(self.mixtureFunctionValue)
        
        if mixtureOrder.lower() == 'reverse':
            sortIdcs = sortIdcs[::-1]
        
        self.mixtureFunctionValue = self.mixtureFunctionValue[sortIdcs]
        self.mixingCoeffs = self.mixingCoeffs[sortIdcs]
            
        
    def getMixedPopulation(self,mixtureNumber,parentPops,**kwargs):
        '''Return a population produced from the mixture definined by the ith mixture coefficients
        '''
        return population(parentPopulations = parentPops,mixingCoefficients = self.mixingCoeffs[mixtureNumber],**kwargs)        
        
    def getMixingCoefficients(self,dFrac):
        ''' Constructs a m x n array of all possible mixing coefficients, where n
        is the number of populations being mixed, and m is the number of allowable
        mixtures.
        '''
        self.dFrac = dFrac
        self.mixingCoeffs = getMixes(self.nParents, np.int(1.0/self.dFrac)+1)
        self.nTotalMixes = len(self.mixingCoeffs[:,1])
    
    def plot(self,plottingStyle = 'mixture value',**kwargs):
        '''Plot the result of the mixture model. Value of objective function for individual
        mixtures can be displayed (\'mixture value\') or the best mixture can be plotted (\'best mixture\'),
        or everything can be collapsed to a single axis (\'single axis\')
        Different plotting styles have different optional arguments.
        for plottingStyle = 'mixture value':
             if there are 2 parents, optional arguments are:
                **kwargs
                Where **kwargs modifiy the color, line arguments, etc of a line plot.
                pyplot.plot(percentOfParent1,ObjectiveFunctionValues,'-',**kwargs)
             if there are 3 parents, optional arguments are:
                cmap = 'viridis' , specifies the colormap that will be used for ternary plots
             if there are more than this, optional arguments are:
                numMixturesToPlot = 10, colormap = 'Accent'
                Where colormap specifies the matplotlib colormap used to color different parents. NumMixturesToPlot specifies
                the number of mixtures to plot, starting from the best mixture and moving towards poorer matches.
        for plottingStyle = 'best mixture':
             optional arguments are:
                colormap = 'Accent',numMixturesToPlot = 1
                These describe the matplotlib colormap to use to draw line colors, and the number of mixtures to plot.
                If numMixturesToPlot = 10, this will plot the 10 best fitting mixtures.
        for plottingStyle = 'single axis':
            optional arguments are:
                colormap = 'Accent', specifies the colors that are used to plot the different colors of parents
        '''
        
        if plottingStyle.lower() == 'mixture value':
            if self.nParents == 2:
                self._plotBinaryMixture_(**kwargs)
            elif self.nParents == 3:
                self._plotTernaryMixture_(**kwargs)
            else:
                self._plotMixturesYAxis_(**kwargs)
        elif plottingStyle.lower() == 'best mixture':
            self._plotBestMixtureCDF_(**kwargs)
        elif plottingStyle.lower() == 'single axis':
            self._plotCollapsedMixturesSingleAxis(**kwargs)
        else:
            print('Whoopsy, an unrecognized plottingStyle was given. Please use either \'mixture value\', \'best mixture\', or \'single axis\'')
    
    def _plotBestMixtureCDF_(self,colormap = 'Accent',numMixturesToPlot = 1, **kwargs):
        '''Plot the best mixture
        '''
        plt.figure()
        ##Plot the parents
        plt.subplot(1,2,1)
        
        cmap = cm.get_cmap(colormap)  
        colors = [cmap(x) for x in np.linspace(0,1,self.nParents)]
        
        for i in range(self.nParents):
            plt.plot(self.parentPops[i]._tAxisCDF_,self.parentPops[i].CDF,linewidth = 1.5, color = colors[i],label = self.parentNames[i]+', %.1f%%'%(self.mixingCoeffs[0,i]*100))
    
        plt.legend(loc = 'lower right')
        plt.xlabel('Age (Ma)',fontsize = 12)
        plt.ylabel(r'$P$',fontsize = 12)
        
        plt.subplot(1,2,2)
        plt.plot(self.daughterPop._tAxisCDF_,self.daughterPop.CDF,'-k',linewidth = 1.5,label = self.daughterName)
        
        if numMixturesToPlot > 1:
            #Use colormap to denote mixture quality? Is this necessary since you are lookng at the PDFS themselves...
        
            for i in range(1,numMixturesToPlot-1):
                mixedDaughter = self.getMixedPopulation(i,self.parentPops)
                plt.plot(mixedDaughter._tAxisCDF_,mixedDaughter.CDF,'--k',linewidth = 1, alpha = 0.5)
        
            #Just label on of these plots
            mixedDaughter = self.getMixedPopulation(i,self.parentPops)
            #Get the string describing the range of mixing coefs
            label = ''
            for i in range(self.nParents):
                label+= self.parentNames[i]+', %.1f %% - %.1f%%,'%(100.0*np.min(self.mixingCoeffs[i,:numMixturesToPlot]),100.0*np.max(self.mixingCoeffs[i,:numMixturesToPlot]))
                
            plt.plot(mixedDaughter._tAxisCDF_,mixedDaughter.CDF,'--k',linewidth = 1, alpha = 0.5,label = label)
            
        mixedDaughter = self.getMixedPopulation(0,self.parentPops) #Get the best mixture
        plt.plot(mixedDaughter._tAxisCDF_,mixedDaughter.CDF,'--r',linewidth = 2,label = 'Best mixture')                       
        
        plt.xlabel('Age (Ma)',fontsize = 12)
        plt.legend(loc = 'lower right')
    
    def _plotTernaryMixture_(self,cmap = 'viridis'):
        ''' When there are three parents, plot the value of the comparison function
        in a ternary plot '''
        
        #Create a ternary figure
        figure, tax = ternary.figure(scale = 100.0)
        
        ## Boundary and Gridlines
    
        # Draw Boundary and Gridlines
        tax.boundary(linewidth=2.0)
        
        # Set Axis labels and Title
        tax.left_axis_label('% '+self.parentNames[2], fontsize=13)
        tax.right_axis_label('% '+self.parentNames[1], fontsize=13)
        tax.bottom_axis_label('% '+self.parentNames[0], fontsize=13)
        
        # Set ticks
        tax.ticks(axis='lbr', linewidth=1, multiple=10)
        
        # Remove default Matplotlib Axes
        tax.clear_matplotlib_ticks()
        
        d = dict()
        for i,mixCoeff in enumerate(self.mixingCoeffs):
            d[(round(mixCoeff[0],3)*100.0,round(mixCoeff[1],3)*100.0,round(mixCoeff[2],3)*100.0)] = self.mixtureFunctionValue[i]
        
        tax.heatmap(d, style='hexagonal',cmap = cmap)
        bestFit = self.mixingCoeffs[0]
        tax.scatter([(round(bestFit[0],3)*100.0,round(bestFit[1],3)*100.0,round(bestFit[2],3)*100.0)], marker='o', lw=3, color='black')
    
        return tax

    def _plotBinaryMixture_(self,**kwargs):
        ''' When there are only two parents, plot the value of the comparison function on the y axis
        and the proportion of those daughters on the x axis
        '''
        
        #Get the proportions as a vector, sort that vector so that data can
        #be plotted as a line
        parent1Coeffs = self.mixingCoeffs[:,0]
        vals = self.mixtureFunctionValue
        sortIdcs = np.argsort(parent1Coeffs)
        parent1Coeffs = parent1Coeffs[sortIdcs]
        vals = vals[sortIdcs]
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        ax.plot(100.0*parent1Coeffs,vals,'-',**kwargs)
        plt.xlabel('Proportion of '+self.parentNames[0]+' (%)',fontsize = 14)
        plt.ylabel('Func. value',fontsize = 12)
    
    def _plotMixturesYAxis_(self,numMixturesToPlot = 10, colormap = 'Accent'):
        '''
        Plot the results of the mixture model, here each parents contribution is represented as a vertically oriented
        polygon - whose x-axis extent describes the proportion of that parent, and the y-axis describes the value of the objective
        function at that mixing value.
        '''
        fig = plt.figure()
        ax = fig.add_subplot(111)        
        
        xVals = self._getMixtureCoeffsXAxis_()[:numMixturesToPlot]
        yVals = self.mixtureFunctionValue[:numMixturesToPlot]
        
        
        cmap = cm.get_cmap(colormap)  
        colors = [cmap(x) for x in np.linspace(0,1,self.nParents)]
        
        for i in range(self.nParents):
            if i==0:
                polygonX = np.hstack((np.zeros_like(yVals), xVals[:,i][::-1], 0))
            
            else:
                polygonX = np.hstack((xVals[:,i-1], xVals[:,i][::-1], xVals[0,i-1]))
            
            polygonY = np.hstack((yVals, np.flipud(yVals), [yVals[0]]))

            xys = np.vstack((polygonX,polygonY)).T
            
            patch = patches.Polygon(xys,facecolor = colors[i],edgecolor = None,label = self.parentNames[i])
            
            ax.add_patch(patch)
        
        plt.ylim(self.mixtureFunctionValue[0],self.mixtureFunctionValue[numMixturesToPlot])
        plt.xlabel('Parent proportion',fontsize = 12)
        plt.ylabel(self.objFunc.func_name,fontsize = 12)
        plt.legend(loc='upper center')
        
    def _plotCollapsedMixturesSingleAxis(self,colormap = 'Accent'):
        
        cmap = cm.get_cmap(name = colormap)
        colors = [cmap(i) for i in np.linspace(0,1.0,self.nParents)]
        
        for i in range(self.nParents):
            theseCoeffs = np.array([phis[i] for phis in self.mixingCoeffs])
            coeffAxis = np.sort(np.unique(theseCoeffs))
            
            y_mean = np.zeros_like(coeffAxis)
            y_max = np.zeros_like(coeffAxis)
            y_min = np.zeros_like(coeffAxis)
            for j,phi in enumerate(coeffAxis):
                idcs = theseCoeffs == phi
                y_mean[j] = np.sum(self.mixtureFunctionValue[idcs])/np.sum(idcs)
                y_min[j] = np.min(self.mixtureFunctionValue[idcs])
                y_max[j] = np.max(self.mixtureFunctionValue[idcs])

            plt.plot(coeffAxis,y_mean,'-',linewidth = 2,color = colors[i],label = self.parentNames[i])
            plt.plot(coeffAxis,y_min,'--',color = colors[i],label = self.parentNames[i])
            plt.plot(coeffAxis,y_max,'--',color = colors[i],label = self.parentNames[i])
        
        plt.xlabel('Mixing coefficient')
        plt.legend(loc = 'best')
        
        #Loop seperately post legend addition to not add this to legend
        bestCoeffs = self.mixingCoeffs[0]
        ax = plt.gca()
        Ylim = [y for y in ax.get_ylim()]
        for i in range(self.nParents):
            plt.plot([bestCoeffs[i], bestCoeffs[i]],Ylim,'--',color = colors[i])
    
    def _mixFunValueToAlpha_(self, numMixturesToPlot,i):
        '''Convert the ith value of the mixture function to a transparency value
        '''
        
        minAlpha = 0.2
        maxAlpha = 1.0

        minVal = np.min(self.mixtureFunctionValue[:numMixturesToPlot])
        maxVal = np.max(self.mixtureFunctionValue[:numMixturesToPlot])
        
        thisVal = self.mixtureFunctionValue[i]
        
        return minAlpha + ((thisVal - minVal)/(maxVal - minVal))*(maxAlpha - minAlpha)
    
    def _getMixtureCoeffsXAxis_(self):
        ''' Get the x-axis coordinates of the boundaries of each mixture
        '''
        
        return np.cumsum(self.mixingCoeffs,axis=1)
        
#==============================================================================
# Class for multiple mixture models (e.g. when dealing with a set of daughters)
#==============================================================================

class mixtureModelSet:

    ''' A class for handling sets of mixtures models, that is when the user wants to perform mixture modelling on more than a single
    daughter.
        This class is largely just a list of individual mixtures models, but also provides utilities for calling the plotting functions
        on the results.
        Most notably, this adds the ability to plot trends in the best fitting mixture coefficients (.plotDaughterTrends)
        and save summaries (.outputResultsSummary)
    '''

    def __init__(self,daughterSet,parentPops,objFunc,dFrac = 0.05,parentNames = None,daughterNames = None,mixtureOrder = 'normal'):
        '''
        Initialize an instance of the mixtureModelSet class.
        :param daughterSet: a list of instances of the population class describing a daughter population
        :param parentPops: a list of instances of the population class describing the parents
        :param objFunc: A function of three variables. The first two are an instance of the population class that describes
        the daughter and a list of instances of population classes describing the daughters. The second is a list or array of
        floats describing the mixing coefficients that sum to 1.
        :param dFrac: A float describing the precision/ spacing with which to construct mixtures. E.g. 0.05 would be looking
        at all mixtures at 5% intervals. Reduction in this term will exponentially increase the computation time.
        :param parentNames: Optionally, a list of strings of the parent names. Used in plotting results
        :param daughterName: Optionally, the name of the daughter population. Used in plotting results.
        :param mixtureOrder: Determines whether the objective function defines good mixtures as high values or low values.
        If mixtureOrder = 'normal', small values of the objective function define good mixtures (e.g. \'Dmax\').
        If mixtureOrder = 'reverse', high values are good mixtures (e.g. \'similarity\' and \'r2 pdp\')
        '''
        #First things first - check that distribution axes match
        unmatchedAxis = False
        
        for daughterPop in daughterSet:           
            for parent in parentPops:
                unmatchedAxis = unmatchedAxis or (np.any(daughterPop._tAxisCDF_
                    != parent._tAxisCDF_)) or (np.any(daughterPop._tAxisDF_
                    != parent._tAxisDF_))
                
        if unmatchedAxis:
            print('Whoopsy, the distributions need to be calculated with matching axes.')
            return None
        else:
            
            self.dFrac = dFrac
            self.parentNames = parentNames
            self.objFuncName = objFunc.__name__
            self.nDaughters = len(daughterSet) #Number of daughters
            self.nParents = len(parentPops)
            self.nParentGrains = [parent.n for parent in parentPops]
            self.nDaughterGrains = [daughter.n for daughter in daughterSet]
            
            #If the daughters are unnamed give them some sort of name
            if daughterNames is None:
                daughterNames = ['D%.0f'%n for n in range(self.nDaughters)] 
            elif len(daughterNames) < self.nDaughters:
                print('Whoopsy, fewer daughter names than daughters, processding with integer names')
                daughterNames = ['D%.0f'%n for n in range(self.nDaughters)]
                
            self.daughterNames = daughterNames
            
            #Do the same things with the parent names            
            if parentNames is None:
                self.parentNames = ['P'+str(x) for x in range(self.nParents)]
            elif len(parentNames) < self.nParents:
                print('Whoops, fewer parents names than parents - proceeding with integer names \n')
                self.parentNames = ['P'+str(x) for x in range(self.nParents)]
            else:
                self.parentNames = parentNames
            
            #Create a list of mixture models for each of the daughters, an improved version of this would
            #not store all the repeat information this creates. In particular there
            #is a significant amount of data duplicated with the mixture coeffecients
            self.mixModels = []
            
            for i,daughter in enumerate(daughterSet):
                self.mixModels.append(MixtureModel(daughter,parentPops,objFunc,dFrac, parentNames,daughterNames[i],mixtureOrder))

            self.areMixturesSorted = False

    def sortMixtures(self,mixtureOrder):
        ''' Sort the mixture functino values and mixture coefficients
        mixtureOder may be either \'normal\' or \'reverse\' if good matches
        have high values or low values        
        '''
        self.areMixturesSorted = True
        for mixModel in self.mixModels:
            mixModel.sortMixtures(mixtureOrder = mixtureOrder)
    
    def plotMixtureResults(self,plottingStyle = 'mixture value',**kwargs):
        '''Plot the individual results of the mixture models. To adjust the type of plot
        specifying an optional argument to \'plottingStyle\'. The value of objective function for individual
        mixtures can be displayed (\'mixture value\') or the best mixture can be plotted (\'best mixture\')
        '''
        
        for i in range(self.nDaughters):
            self.mixModels[i].plot(plottingStyle,**kwargs)
            plt.title(self.daughterNames[i])

    def plotSingleMixture(self,daughterName,plottingStyle = 'mixture value',**kwargs):
        '''
        Plot one of the mixture model results
        :param daughterName:
        :return:
        '''

        #Did we find the requested daughter name yet?
        foundDaughter = False

        #Iterate through each of the daughter names
        for i in range(self.nDaughters):

            #If this is the requested daughter name
            if self.daughterNames[i] == daughterName:

                foundDaughter = True #Hooray! A name was correctly specified
                self.mixModels[i].plot(plottingStyle, **kwargs)
                plt.title(self.daughterNames[i])

                #No need to keep looking
                break

        #No daughter of the specified name found
        if not foundDaughter:
            print('Whoopsy, could not find a daughter with that name, the available names are:\n')
            print(self.daughterNames)



    def plotDaughterTrends(self,cmap = 'Accent'):
        ''' Plot trends in the mixture results
        '''
        
        #Only make this plot if the mixtures have already been sorted
        if self.areMixturesSorted:
            allBestFits = np.zeros((self.nDaughters,self.nParents))
            allFitQualities = np.zeros(self.nDaughters)
            #Iterate through the mixtures
            for i in range(self.nDaughters):
                allBestFits[i,:] = self.mixModels[i].mixingCoeffs[0][:]
                allFitQualities[i] = self.mixModels[i].mixtureFunctionValue[0]
    
            colorFun = cm.get_cmap(cmap)
            colors = [colorFun(i) for i in np.linspace(0,1,self.nParents)]
            #Iterate through the parents and plot
            xAxis = range(self.nDaughters)
            for i in range(self.nParents):
                plt.plot(xAxis,allBestFits[:,i],'--ok',color = colors[i],label = self.parentNames[i])
            
            plt.xticks(xAxis, self.daughterNames,rotation = 45) #Replace the x axis labels with the daughter names
            plt.xlabel('Daughters',fontsize = 13)
            plt.ylabel('Parent proportions',fontsize = 13)
            plt.ylim(0,1)
            plt.legend(loc = 'best')
            
        else:
            print('Whoopsy, must sort the mixtures first')
            
    def outputResultsSummary(self,filename,delimiter='\t'):
        ''' Write the results to an ouput file
        '''
        
        if os.path.isfile(filename):
            print('Whoopsy, that filename already exists, please choose another')
        else:
            #Open the file for writing
            with open(filename,'w') as f:
                
                #Write out some information about the parameters
                f.write('Mixture coefficient spacing:'+delimiter+'%.3f \n'%self.dFrac)
                f.write('Objective function name:'+delimiter+self.objFuncName+'\n')

                #Write out information about the parents  
                f.write('Parent name' + delimiter + 'N\n')
                for i in range(self.nParents):
                    line = self.parentNames[i] + delimiter + '%.0f\n'%self.nParentGrains[i]
                    f.write(line)
                
                #Write out information about the daughters
                mixHeader = 'Daughter name' + delimiter + 'N' + delimiter
                mixHeader+= delimiter.join(self.parentNames)
                mixHeader+= delimiter + 'Mix. fun value \n'
                f.write(mixHeader)
                for i in range(self.nDaughters):
                    line = self.daughterNames[i] + delimiter + '%.0f'%self.nDaughterGrains[i] + delimiter
                    strVals = ['%.3f'%self.mixModels[i].mixingCoeffs[0][j] for j in range(self.nParents)] #Get the values of the best mixing coeffs as strings
                    line+= delimiter.join(strVals)
                    line+= delimiter + '%.3e\n'%self.mixModels[i].mixtureFunctionValue[0]
                    f.write(line)
