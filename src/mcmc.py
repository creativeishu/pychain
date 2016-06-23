"""
Markov Chain Monte Carlo (MCMC) method using Metropolis-Hastings algorithm.

Copyright: MarkovianLabs 
"""

import numpy as np 
import matplotlib.pyplot as plt
from sys import exit

#----------------------------------------------------------
__author__ = ("Irshad Mohammed <creativeishu@gmail.com>")
#----------------------------------------------------------

class MCMC(object):
	"""
	Implements single MCMC using MH algorithm.

	Parameters
	----------
	NumberOfSteps (Integer): Maximum number of steps.
	NumberOfParams (Integer): Total number of model parameters. 
	Mins (1d array): Array containing the minimum values of each model parameter.
	Maxs (1d array): Array containing the maximum values of each model parameter.
	SDs (1d array): Array containing the expected standard deviation values of each model parameter.
	alpha (Float): Scaling for chain steps.
	write2file (Boolean): Flag indication whether to write the output to a file.
							Default is False.
	outputfilename (String): Name of the output file.
	randomseed (Integer): Random Seed.
	"""
	def __init__(self, NumberOfSteps=10000, \
				NumberOfParams=2, Mins=[0.0,-1.0], Maxs=[2.0,1.0], SDs=[1.0,1.0], alpha=1.0,\
				write2file=False, outputfilename='chain.mcmc', randomseed=250192, debug=False):
		"""
		Instantiates the class.
		"""
		np.random.seed(randomseed)

		if not (NumberOfParams == len(Mins) and \
			NumberOfParams==len(Maxs) and NumberOfParams==len(SDs)):
			print "Length of Mins, Maxs and SDs should be same as NumberOfParams"
			exit()
		
		self.write2file=write2file
		self.outputfilename=outputfilename

		self.NumberOfSteps = NumberOfSteps
		self.NumberOfParams = NumberOfParams
		self.mins = np.array(Mins)
		self.maxs = np.array(Maxs)
		self.SD = np.array(SDs)

		self.alpha = alpha
		self.CovMat = 100.0*self.alpha*np.diag(self.SD**2)

		self.debug = debug

#----------------------------------------------------------
		
	def FirstStep(self):
		"""
		Initiates the chain.

		Returns
		-------
		A numpy array containing random initial value for each parameter. 
		"""
		return self.mins + \
				np.random.uniform(size=self.NumberOfParams)*\
				(self.maxs - self.mins)

#----------------------------------------------------------

	def NextStep(self,Oldstep):
		"""
		Generates the next step for the chain.

		Parameters
		----------
		Oldstep (1d array): A numpy array containing the current values of the parameters.

		Returns
		-------
		A numpy array containing the next values of the parameters.
		"""
		NS = np.random.multivariate_normal(Oldstep,self.CovMat)
		while np.any(NS<self.mins) or np.any(NS>self.maxs):
			NS = np.random.multivariate_normal(Oldstep,self.CovMat)
		return NS

#----------------------------------------------------------

	def MetropolisHastings(self,Oldchi2,Newchi2):
		"""
		Determines the acceptance of new step based upon current step.

		Parameters
		----------
		Oldchi2 (Float): Chi-square of the current step.
		Newchi2 (Float): Chi-square of the new step.

		Returns
		-------
		True if the new step is accepted, False otherwise.
		"""
		likelihoodratio = np.exp(-(Newchi2-Oldchi2)/2)
		if likelihoodratio < np.random.uniform():
			return False
		else:
			return True


#----------------------------------------------------------

	def chisquare(self, Params):
		"""
		Computes Chi-square of the parameters.

		Parameters
		----------
		Params (1d array): Numpy array containing values of the parameters.

		Returns
		-------
		Value of the Chi-square.
		"""
		return np.random.chisquare(df=len(Params))

#----------------------------------------------------------

	def MainChain(self):
		"""
		Runs the chain.

		Returns
		-------
		Acceptance rate.
		"""
		# Initialising the chain
		OldStep = self.FirstStep()
		Oldchi2 = self.chisquare(OldStep)
		Bestchi2 = Oldchi2

		# Preparing output file
		if self.write2file:
			outfile = open(self.outputfilename,'w')
			writestring = '%1.6f \t'*self.NumberOfParams

		# Initialising multiplicity and accepted number of points.
		multiplicity = 0
		acceptedpoints = 0

		# Chain starts here...
		for i in range(self.NumberOfSteps):
			multiplicity += 1

			# Generating next step and its chi-square
			NewStep = self.NextStep(OldStep)
			Newchi2 = self.chisquare(NewStep)

			if self.debug:
				strFormat = 7 * '{:10f} '
				print "Step Number: %i \t Accepted Points: %i"%(i, acceptedpoints)
				print 'Old: ', Oldchi2, strFormat.format(*OldStep)
				print 'New: ', Newchi2, strFormat.format(*NewStep)
				print

			# Checking if it is to be accepted.
			GoodPoint = self.MetropolisHastings(Oldchi2,Newchi2)

			# Updating step scale using a threshold chi-square.
			if Newchi2<350:
						self.CovMat = self.alpha*np.diag(self.SD**2)

			if GoodPoint:
				# Updating best chi-square so far in the chain.
				if Newchi2<Bestchi2:
					strFormat = 7 * '{:10f} '
					Bestchi2=Newchi2
					print i, acceptedpoints, Bestchi2, strFormat.format(*NewStep)

				# Writing accepted steps into the output file
				if self.write2file:
					print >>outfile, '%i \t'%i, '%1.6f \t'%Newchi2,'%i \t'%multiplicity,\
								writestring%tuple(NewStep)

				# Updating number of accepted points.
				acceptedpoints += 1
				multiplicity = 0

				# Updating the old step. 
				OldStep = NewStep
				Oldchi2 = Newchi2
			else:
				continue
		# Writing Best chi-square of the full chain and the acceptance ratio.
		print "Best chi square: %1.5f"%Bestchi2
		print "Acceptance Ratio: %1.5f"%(float(acceptedpoints)/i) 

		return float(acceptedpoints)/i

#==============================================================================

if __name__=="__main__":
	print "Class: Markov Chain Monte Carlo (MCMC) method using Metropolis-Hastings algorithm."

	

