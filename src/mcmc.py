import numpy as np 
import matplotlib.pyplot as plt
from sys import exit

#----------------------------------------------------------
__author__ = ("Irshad Mohammed <creativeishu@gmail.com>")
#----------------------------------------------------------

class MCMC(object):
	'''

	'''
	def __init__(self, NumberOfSteps=10000, \
				NumberOfParams=2, Mins=[0.0,-1.0], Maxs=[2.0,1.0], SDs=[1.0,1.0], alpha=1.0,\
				write2file=False, outputfilename='chain.mcmc', randomseed=250192):
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
		self.CovMat = 10.0*self.alpha*np.diag(self.SD**2)

#----------------------------------------------------------
		
	def FirstStep(self):
		return self.mins + \
				np.random.uniform(size=self.NumberOfParams)*\
				(self.maxs - self.mins)

#----------------------------------------------------------

	def NextStep(self,Oldstep):
		NS = np.random.multivariate_normal(Oldstep,self.CovMat)
		while np.any(NS<self.mins) or np.any(NS>self.maxs):
			NS = np.random.multivariate_normal(Oldstep,self.CovMat)
		return NS

#----------------------------------------------------------

	def MetropolisHastings(self,Oldchi2,Newchi2):
		likelihoodratio = np.exp(-(Newchi2-Oldchi2)/2)
		if likelihoodratio < np.random.uniform():
			return False
		else:
			return True

#----------------------------------------------------------

	def MainChain(self):
		OldStep = self.FirstStep()
		Oldchi2 = self.chisquare(OldStep)
		Bestchi2 = Oldchi2

		if self.write2file:
			outfile = open(self.outputfilename,'w')
			writestring = '%1.6f \t'*self.NumberOfParams

		multiplicity = 0
		acceptedpoints = 0

		for i in range(self.NumberOfSteps):
			multiplicity += 1
			NewStep = self.NextStep(OldStep)
			Newchi2 = self.chisquare(NewStep)
			GoodPoint = self.MetropolisHastings(Oldchi2,Newchi2)
			if Newchi2<2*len(self.Y):
						self.CovMat = self.alpha*np.diag(self.SD**2)

			# print Oldchi2, Newchi2, GoodPoint, OldStep, NewStep
			if GoodPoint:
				if Newchi2<Bestchi2:
					Bestchi2=Newchi2
					print Bestchi2, NewStep

				if self.write2file:
					print >>outfile, '%1.6f \t'%Newchi2,'%i \t'%multiplicity,\
								writestring%tuple(NewStep)
				acceptedpoints += 1
				multiplicity = 0
				Oldstep = NewStep
				Oldchi2 = Newchi2
				plt.plot(NewStep[0], NewStep[1],'.k')
			else:
				continue
		print "Best chi square: %1.5f"%Bestchi2
		print "Acceptance Ratio: %1.5f"%(float(acceptedpoints)/i) 
		# plt.xlim(0,10)
		# plt.ylim(20,30)
		plt.show()
		return float(acceptedpoints)/i

#----------------------------------------------------------

	def chisquare(self, Params):
		return np.random.chisquare(df=len(Params))

#==============================================================================

if __name__=="__main__":
	obj = MCMC()
	print obj.MainChain()

	

