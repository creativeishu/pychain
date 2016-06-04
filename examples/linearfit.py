"""
Example demonstrating the usage for the linear-fit scenario.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'../src')
from mcmc import MCMC

#----------------------------------------------------------

class LinearFit(MCMC):
	"""
	Extends the original MCMC class to sample the parameters of a linear model.

	Parameters
	-----------
	MCMC (Class): Parent MCMC class.
	m (Float): Feducial value of the slope of the linear data.
	c (Float): Feducial value of the intercept of the linear data.
	RedStd (Float): Feducial value of the standard deviation of the linear data.
	"""
	def __init__(self, m=5.0, c=25.0, RedStd=15.0):
		"""
		Instantiates the class by synthetically generating data.
		"""
		MCMC.__init__(self, NumberOfSteps=100000, \
				NumberOfParams=2, Mins=[0.0,20.0], Maxs=[10.0,30.0], SDs=[1.0,2.0], alpha=0.2,\
				write2file=True, outputfilename='chain.mcmc', randomseed=250192)		

		self.X=np.linspace(-10, 10, 25)
		self.delta = np.random.uniform(low=-1*RedStd, high=RedStd, size=len(self.X))
		self.Y = (m*self.X + c) + self.delta

#----------------------------------------------------------

	def FittingFunction(self, Params):
		"""
		Parametric form of the model.

		Parameters
		----------
		Params (1d array): Numpy array containing values of the parameters. 

		Returns
		-------
		model values (y = mx + c)
		"""
		return Params[0]*self.X + Params[1]

#----------------------------------------------------------

	def chisquare(self, Params):
		"""
		Computes Chi-square.

		Parameters
		----------
		Params (1d array): Numpy array containing values of the parameters. 

		Returns
		-------
		chi square.
		"""
		kisquare = ((self.Y-self.FittingFunction(Params))/self.delta)**2
		return np.sum(kisquare)

#==============================================================================

if __name__=="__main__":
	co = LinearFit()
	print co.MainChain()