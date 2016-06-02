import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'../src')
from mcmc import MCMC

class LinearFit(MCMC):
	'''

	'''
	def __init__(self, m=5.0, c=25.0, RedStd=15.0):

		MCMC.__init__(self, NumberOfSteps=100000, \
				NumberOfParams=2, Mins=[0.0,20.0], Maxs=[10.0,30.0], SDs=[1.0,2.0], alpha=0.2,\
				write2file=True, outputfilename='chain.mcmc', randomseed=250192)		

		self.X=np.linspace(-10, 10, 25)
		self.delta = np.random.uniform(low=-1*RedStd, high=RedStd, size=len(self.X))
		self.Y = (m*self.X + c) + self.delta
		# plt.errorbar(self.X, self.Y, self.delta)
		# plt.show()

	def FittingFunction(self, Params):
		return Params[0]*self.X + Params[1]

	def chisquare(self, Params):
		kisquare = ((self.Y-self.FittingFunction(Params))/self.delta)**2
		# print np.sum(kisquare)
		return np.sum(kisquare)

if __name__=="__main__":
	co = LinearFit()
	print co.MainChain()