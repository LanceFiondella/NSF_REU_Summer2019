import sys, numpy as np
from math import *
import scipy, time, warnings
import scipy.optimize

model = sys.argv[1]

warnings.filterwarnings("ignore")



#Weibull model
data = [3,33,146,227,342,351,353,444,556,571,709,759,836,860,968,1056,1726,1846,1872,1986,2311,2366,2608,2676,3098,3278,3288,4434,5034,5049,5085,5089,5089,5097,5324,5389,5565,5623,6080,6380,6477,6740,7192,7447,7644,7837,7843,7922,8738,10089,10237,10258,10491,10625,10982,11175,11411,11442,11811,12559,12559,12791,13121,13486,14708,15251,15261,15277,15806,16185,16229,16358,17168,17458,17758,18287,18568,18728,19556,20567,21012,21308,23063,24127,25910,26770,27753,28460,28493,29361,30085,32408,35338,36799,37642,37654,37915,39715,40580,42015,42045,42188,42296,42296,45406,46653,47596,48296,49171,49416,50145,52042,52489,52875,53321,53443,54433,55381,56463,56485,56560,57042,62551,62651,62661,63732,64103,64893,71043,74364,75409,76057,81542,82702,84566,88682]
data = np.array([float(j) for j in data])
#Best solution  - 966.0803348790324
def RLLWei(x): 
	n = len(data)
	tn = data[n-1] 
	b,c=x
	aHat =  n / (1-np.exp(-b*(tn**c)))
	r = -(-n +sum(np.log(aHat*b*c*np.power(data,(c-1))*np.exp(-b*np.power(data,c)))))
	if np.isnan(r):
		return float('inf')
	return r

def RLLWeiRand(pop_size):
	out_pop = []
	for i in range(pop_size):
		b = np.random.uniform(0, 0.1)
		c = np.random.uniform(0, 5)
		out_pop.append([b, c])
		
	return out_pop

def RLLWeiEst(pop_size):
	out_pop = []
	alpha = 2

	b0 = float(len(data)) / sum(data)
	c0 = 1

	for i in range(pop_size):
		b = np.random.uniform(b0/alpha, b0*alpha)
		c = np.random.uniform(c0/alpha, c0*alpha)
		out_pop.append([b, c])
		
	return out_pop

#--- Covariate model ------------------------------------------------------------+
########### Best solution  - -23.0067  ###########################################

kVec = np.array([2, 11, 2, 4, 3, 1, 1, 2, 4, 0, 4, 1, 3, 0])
EVec = np.array([0.05, 1, 0.19, 0.41, 0.32, 0.61, 0.32, 1.83, 3.01, 1.79, 3.17,3.4, 4.2, 1.2])
FVec = np.array([1.3, 17.8, 5.0, 1.5, 1.5, 3.0, 3.0, 8, 30, 9, 25, 15, 15, 2])
CVec = np.array([0.5, 2.8, 1, 0.5, 0.5, 1, 0.5, 2.5, 3.0, 3.0, 6, 4, 4, 1])

#--- COST/OBJECTIVE FUNCTION ------------------------------------------------------------+
def RLLCV(x):
	n = len(kVec)
	b, b1, b2, b3 = x
	second = []
	prodlist = []
	for i in range(n):
		sum1=1
		sum2=1
		try:
			sum1=1-((1-b)**(exp(EVec[i]*b1)*exp(FVec[i]*b2)*exp(CVec[i]*b3)))
		except OverflowError:
			return float('inf')
		for k in range(i):
			sum2 = sum2*((1-b)**(exp(EVec[k]*b1)*exp(FVec[k]*b2)*exp(CVec[k]*b3)))
		second.append(sum2)
		prodlist.append(sum1*sum2)
	
	firstTerm = -sum(kVec) #Verified

	secondTerm = sum(kVec)*np.log(sum(kVec)/sum(prodlist))
	
	logTerm = [] #Verified
	for i in range(n):
		logTerm.append(kVec[i]*np.log(prodlist[i]))
	thirdTerm = sum(logTerm)
	
	factTerm = [] #Verified
	for i in range(n):
		factTerm.append(np.log(factorial(kVec[i])))
	fourthTerm = sum(factTerm)
	cv = -(firstTerm + secondTerm + thirdTerm - fourthTerm)
	if np.isnan(cv):
		return float('inf')
	return cv

#--- multivariable optimization

def Sphere(vector):
	return np.sum(np.square(vector)) + 1

def Rastrigin(vec):
	A = 10
	n = len(vec)
	return A*n + np.sum([x**2 - A*np.cos(2*np.pi*x) for x in vec]) + 1

def Rastrigin8Rand(pop_size):
	out_pop = []
	for i in range(pop_size):
		out_pop.append([np.random.uniform(-5.12, 5.12) for i in range(8)])
		
	return out_pop


def Ackley(vec):
	rsq = np.sum(np.square(vector))
	cos_sum = np.sum([2*np.pi*x for x in vec])
	return -20 * np.exp(-0.2 * np.sqrt(0.5 * rsq) ) - np.exp(0.5 * cos_sum) + np.exp(1) + 20 + 1

def Beale(vec):
	x,y = vec
	return (1.5 - x + (x*y))**2 + (2.25 - x + (x*y**2))**2 + (2.625 - x + (x*y**3))**2 + 1

def GoldsteinPrice(vec):
	x, y = vec
	return (1 + (x + y + 1)**2 * (19 - 14*x + 3*x**2 - 14*y + 6*x*y + 3*y**2) ) + \
			(30 + (2*x - 3*y)**2 * (18 - 32*x + 12*x**2 + 48*y - 36*x*y + 27*y**2) )

def Booth(vec):
	x, y = vec
	return (x + 2*y - 7)**2 + (2*x + y - 5)**2 + 1


#--- EM calculation -----------------------------------------------------------------

def calcMLEs(x):
	n = len(data)
	tn = data[n-1] 
	b, c = x 
	aMLE = n / (1 - np.exp(-b * np.power(tn,c)))
	bHat = -aMLE*np.power(tn,c)*np.exp(-b *np.power(tn,c))+ sum((1-b*np.power(data,c))/(b))
	cHat = -aMLE*b*np.power(tn,c)*np.log(tn)*np.exp(-b*np.power(tn,c)) + (n/c) + sum(np.log(data)-b*np.log(data)*np.power(data,c))
	if np.isnan(bHat) or np.isnan(cHat):
		return [b, c]
	return [bHat, cHat]

def calcMLEsSecondorder(x): #Verified so be sure to use this code throughout
	n = len(data)
	tn = data[n-1] 
	b, c = x 
	aMLE = n / (1 - np.exp(-b * np.power(tn,c)))
	bHat2 = -(n/np.power(b,2)) + ((np.exp(b*np.power(tn,c))*n*np.power(tn,2*c))/(np.power((1-np.exp(b*np.power(tn,c))),2)))
	cFirstTerm = aMLE*b*np.power(tn,c)*np.exp(-b*np.power(tn,c))*(-1+b*np.power(tn,c))*np.power(np.log(tn),2)
	cSecondTerm = -n/np.power(c,2)
	cThirdTerm = np.sum(-b*np.power(np.log(data),2)*np.power(data,c))
	cHat2 = cFirstTerm + cSecondTerm +cThirdTerm
	return [bHat2, cHat2]

def NM(estimates):
	'''
	Newton's method
	'''
	bMLEinit, cMLEinit = estimates
	try:
		result = scipy.optimize.newton(calcMLEs,x0=(bMLEinit,cMLEinit), fprime=calcMLEsSecondorder, tol=1e-10, maxiter=10000)
		#print(time.time() - ts)
		return result.root, all(result.converged)
	except RuntimeError:	#scipy throws an error if it fails to converge - catch and shortcircuit
		#print(time.time() - ts)
		return estimates, False

# ------------ SCIPY METHODS -----------------------------------------------------

def nelder_mead(x):
	r = scipy.optimize.minimize(models[model]['objective'], x, method='Nelder-Mead', tol=1e-10, options={'maxiter':100})
	return r.x, r.success

def powell(x):
	r = scipy.optimize.minimize(models[model]['objective'], x, method='Powell', tol=1e-10, options={'maxiter':100})
	return r.x, r.success

def cg(x):
	r = scipy.optimize.minimize(models[model]['objective'], x, method='CG', tol=1e-10, options={'maxiter':100})
	return r.x, r.success

def bfgs(x):
	r = scipy.optimize.minimize(models[model]['objective'], x, method='BFGS', tol=1e-10, options={'maxiter':100})
	return r.x, r.success

def lbfgsb(x):
	r = scipy.optimize.minimize(models[model]['objective'], x, method='L-BFGS-B', tol=1e-10, options={'maxiter':100})
	return r.x, r.success

def tnc(x):
	r = scipy.optimize.minimize(models[model]['objective'], x, method='TNC', tol=1e-10, options={'maxiter':100})
	return r.x, r.success

def cobyla(x):
	r = scipy.optimize.minimize(models[model]['objective'], x, method='COBYLA', tol=1e-10, options={'maxiter':100})
	return r.x, r.success

def slsqp(x):
	r = scipy.optimize.minimize(models[model]['objective'], x, method='SLSQP', tol=1e-10, options={'maxiter':100})
	return r.x, r.success

def dogleg(x):
	r = scipy.optimize.minimize(models[model]['objective'], x, method='dogleg', tol=1e-10, options={'maxiter':100})
	return r.x, r.success

def trustncg(x):
	r = scipy.optimize.minimize(models[model]['objective'], x, method='trust-ncg', tol=1e-10, options={'maxiter':100})
	return r.x, r.success






# ------------ FORMATTED / ESTIMATES -----------------------------------------------------

models = {
	"Weibull":{
		"objective":	RLLWei,
		"dimensions":	2,
		"search_space":	[0, 1],
		"rand":			RLLWeiRand,
		"estimates":	RLLWeiEst,#[ [0.1, 0.1],[0.9, 0.1] ],
		"result":		966.0803348790324
	},
	"Covariate":{
		"objective":	RLLCV,
		"dimensions":	4,
		"search_space":	[0.00000001,1-0.00000001],
		"estimates":	[ [0.2, 0.2], [0.2, 0.2], [0.2, 0.2], [0.2, 0.2] ],
		"result":		23.0067
	},
	"Sphere":{
		"objective":	Sphere,
		"dimensions":	2,
		"search_space":	[-1,1],
		"estimates":	[[ 0, 1] for i in range(2)],
		"result":		1
	},

	"Rastrigin8":{
			"objective":    Rastrigin,
			"dimensions":   8,
			"search_space": [-5.12,5.12],
			"rand":			Rastrigin8Rand,
			"estimates":    Rastrigin8Rand,
			"result":               1# added 1 to eval to calculate error
	}

}
