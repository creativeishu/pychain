"""
Example demonstrating the usage for the quadratic-fit scenario.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'../src')
from mcmc import MCMC

#----------------------------------------------------------

class QuadraticFit(MCMC):
	"""
	Extends the original MCMC class to sample the parameters of a quadratic model (y = ax^2 + bx + c)

	Parameters
	-----------
	MCMC (Class): Parent MCMC class.
	a (Float): Feducial value of the a parameter.
	b (Float): Feducial value of the b parameter.
	c (Float): Feducial value of the c parameter.
	RedStd (Float): Feducial value of the standard deviation of the linear data.
	"""
	def __init__(self, a=1.0, b=10.0, c=25.0, RedStd=15.0):

		MCMC.__init__(self, TargetAcceptedPoints=10000, \
				NumberOfParams=3, Mins=[-5.0, 5.0, 20], Maxs=[5.0, 15.0, 30.0], SDs=[0.7,4.0,7.0], alpha=0.01,\
				write2file=True, outputfilename='chain.mcmc', randomseed=250192)		

		self.X=np.linspace(0, 10, 25)
		self.delta = np.random.uniform(low=-1*RedStd, high=RedStd, size=len(self.X))
		self.Y = (a*self.X**2 + b*self.X + c) + self.delta

#----------------------------------------------------------

	def FittingFunction(self, Params):
		"""
		Parametric form of the model.

		Parameters
		----------
		Params (1d array): Numpy array containing values of the parameters. 

		Returns
		-------
		model values (y = ax^2 + bx + c)
		"""
		return Params[0]*self.X**2 + Params[1]*self.X + Params[2]

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
	co = QuadraticFit()
	print co.MainChain()