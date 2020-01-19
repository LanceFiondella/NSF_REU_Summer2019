import sys, numpy as np
from math import *
import scipy, time, warnings
import scipy.optimize

model = sys.argv[1] if len(sys.argv) > 1 else 'Covariate'

warnings.filterwarnings("ignore")

#Weibull model
data = [3,33,146,227,342,351,353,444,556,571,709,759,836,860,968,1056,1726,1846,1872,1986,2311,2366,2608,2676,3098,3278,3288,4434,5034,5049,5085,5089,5089,5097,5324,5389,5565,5623,6080,6380,6477,6740,7192,7447,7644,7837,7843,7922,8738,10089,10237,10258,10491,10625,10982,11175,11411,11442,11811,12559,12559,12791,13121,13486,14708,15251,15261,15277,15806,16185,16229,16358,17168,17458,17758,18287,18568,18728,19556,20567,21012,21308,23063,24127,25910,26770,27753,28460,28493,29361,30085,32408,35338,36799,37642,37654,37915,39715,40580,42015,42045,42188,42296,42296,45406,46653,47596,48296,49171,49416,50145,52042,52489,52875,53321,53443,54433,55381,56463,56485,56560,57042,62551,62651,62661,63732,64103,64893,71043,74364,75409,76057,81542,82702,84566,88682]
data = np.array([float(j) for j in data])
#--- COST/OBJECTIVE FUNCTION ------------------------------------------------------------+
########### Best solution  - -966.0803348790324  ###########################################

# function we are attempting to optimize (minimize), 2-variable
def RLLWei(x): 
	n = len(data)
	tn = data[n-1] 
	b,c=x
	aHat =  n / (1-np.exp(-b*(tn**c)))
	r = -(-n +sum(np.log(aHat*b*c*np.power(data,(c-1))*np.exp(-b*np.power(data,c)))))
	if np.isnan(r):
		return float('inf')
	return r

# ZTP model ----------------------------------------------------------------------

tVec = np.array([3,33,146,227,342,351,353,444,556,571,709,759,836,860,968,1056,1726,1846,1872,1986,2311,2366,2608,2676,3098,3278,3288,4434,5034,5049,5085,5089,5089,5097,5324,5389,5565,5623,6080,6380,6477,6740,7192,7447,7644,7837,7843,7922,8738,10089,10237,10258,10491,10625,10982,11175,11411,11442,11811,12559,12559,12791,13121,13486,14708,15251,15261,15277,15806,16185,16229,16358,17168,17458,17758,18287,18568,18728,19556,20567,21012,21308,23063,24127,25910,26770,27753,28460,28493,29361,30085,32408,35338,36799,37642,37654,37915,39715,40580,42015,42045,42188,42296,42296,45406,46653,47596,48296,49171,49416,50145,52042,52489,52875,53321,53443,54433,55381,56463,56485,56560,57042,62551,62651,62661,63732,64103,64893,71043,74364,75409,76057,81542,82702,84566,88682])

def LLInitiaPDFZTP(x): 

	n = len(tVec)
	tn = tVec[-1]

	b, a1, b1, c,p = x
	Second = []
	for i in range(n):
		FIFirst = c*exp(b*tVec[i])*(1+a1)
		FISecond = (-1+exp(b*tVec[i]))/(a1+exp(b*tVec[i]))
		FIPower = -1+((c*(p-b1))/b)
		FIDenom = (a1+exp(b*tVec[i]))**2
		Second.append(log((FIFirst*(FISecond**FIPower))/FIDenom))
	SecondTerm = sum(Second)

def RLLZTP(x):

	n = len(tVec)
	tn = tVec[-1]

	b, a1, b1, c, p = x
	print(b, a1, b1, c, p)
	#aMLE
	firstT = (exp(-b*tn)*(1+a1))/(1+exp(-b*tn)*a1)
	SecondTPower =  (-c*(p-b1))/b
	ThirdT = (p-b1)
	aMLE = n *((1-firstT)**SecondTPower)*ThirdT

	#MVF
	mtFirst = aMLE/(p-b1)
	mtSecond = ((1+a1)*exp(-b*tn))/(1+a1*exp(-b*tn))
	mtThree = (c*(p-b1))/b

	FirstTerm = mtFirst*((1-mtSecond)**mtThree)
	
	#FI
	Second = []
	for i in range(n):
		FIFirst = aMLE*c*exp(b*tVec[i])*(1+a1)
		FISecond = (-1+exp(b*tVec[i]))/(a1+exp(b*tVec[i]))
		FIPower = -1+((c*(p-b1))/b)
		FIDenom = (a1+exp(b*tVec[i]))**2
		sum1 = (FIFirst*(FISecond**FIPower))/(FIDenom)
		Second.append(log(sum1))
	return -(-FirstTerm+sum(Second))

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
		result = scipy.optimize.newton(calcMLEs,x0=(bMLEinit,cMLEinit), fprime=calcMLEsSecondorder, tol=1e-10, maxiter=10000, full_output=True)
		#print(time.time() - ts)
		return result.root, all(result.converged)
	except RuntimeError:	#scipy throws an error if it fails to converge - catch and shortcircuit
		#print(time.time() - ts)
		return estimates, False


#--- ECM calculation ----------------------------------------------------------------

def logL(b,c):
	n = len(data)
	tn = data[n-1]
	aHat = n / (1-np.exp(-b*(tn**c)))
	return (-n +sum(np.log(aHat*b*c*np.power(data,(c-1))*np.exp(-b*np.power(data,c)))))

def bMLE(b, c):
	n = len(data)
	tn = data[n-1]
	aMLE = n/ (1 - np.exp(-b * np.power(tn,c)))
	return -aMLE*np.power(tn,c)*np.exp(-b *np.power(tn,c))+ sum((1-b*np.power(data,c))/(b))

def cMLE(c, b):
	n = len(data)
	tn = data[n-1]
	aMLE = n/ (1 - np.exp(-b * np.power(tn,c)))
	return -aMLE*b*np.power(tn,c)*np.log(tn)*np.exp(-b*np.power(tn,c)) + (n/c) + sum(np.log(data)-b*np.log(data)*np.power(data,c))

def ECM(estimates):
	ts = time.time()
	b_est, c_est = estimates                            # initial passed estimates for b and c
	ll_list = [logL(b_est, c_est), logL(b_est, c_est)]  # used to compare current error to previous error in loop
	ll_error = 1                                        # initial error
	tolerance = 1e-10
	max_iterations = 1000

	iterations = 0
	while (ll_error > tolerance) and (iterations < max_iterations):
		b_est = scipy.optimize.fsolve(bMLE, x0 = b_est, args=(c_est))   # solve for b using c estimate
		c_est = scipy.optimize.fsolve(cMLE, x0 = c_est, args=(b_est))   # solve for c using b estimate
		ll_list[1] = (logL(b_est, c_est))   # log likelihood of new estimates
		ll_error = ll_list[1] - ll_list[0]  # calculate ll error, new ll - old ll
		ll_list[0] = ll_list[1]             # calculated ll becomes old ll for next iteration
		iterations += 1
	roots = np.array([b_est[0], c_est[0]])
	if (RLLWei(roots) / models["Weibull"]["result"] <= tolerance + 1) and (iterations < max_iterations):
		converged = True
	else:
		#if iterations == 100:
			#print('eeeeeee')
		converged = False
	#print(time.time() - ts, "ECM")
	return roots, converged


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
		"search_space":	[0.000001,1-0.000001],
		"estimates":	[ [0.1, 0.1],[0.9, 0.1] ],
		"result":		966.0803348790324
	},
	"Covariate":{
		"objective":	RLLCV,
		"dimensions":	4,
		"search_space":	[0.00000001,1-0.00000001],
		"estimates":	[ [0.2, 0.2], [0.2, 0.2], [0.2, 0.2], [0.2, 0.2] ],
		"result":		23.0067
	},
	"ZTP":{
		"objective":	RLLZTP,
		"dimensions":	5,
		"search_space":	[-1, 1],
		"estimates":	[ [0.05, 0.2], [0.05, 0.2], [0.05, 0.2], [0.05, 0.2], [0.05, 0.2] ],
		"result":		966.102
	},
	"Sphere":{
		"objective":	Sphere,
		"dimensions":	2,
		"search_space":	[-1,1],
		"estimates":	[[ 0, 1] for i in range(2)],
		"result":		1
	},
	"Rastrigin2":{
			"objective":    Rastrigin,
			"dimensions":   2,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(2)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin3":{
			"objective":    Rastrigin,
			"dimensions":   3,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(3)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin4":{
			"objective":    Rastrigin,
			"dimensions":   4,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(4)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin5":{
			"objective":    Rastrigin,
			"dimensions":   5,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(5)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin6":{
			"objective":    Rastrigin,
			"dimensions":   6,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(6)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin7":{
			"objective":    Rastrigin,
			"dimensions":   7,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(7)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin8":{
			"objective":    Rastrigin,
			"dimensions":   8,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(8)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin9":{
			"objective":    Rastrigin,
			"dimensions":   9,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(9)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin10":{
			"objective":    Rastrigin,
			"dimensions":   10,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(10)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin11":{
			"objective":    Rastrigin,
			"dimensions":   11,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(11)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin12":{
			"objective":    Rastrigin,
			"dimensions":   12,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(12)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin13":{
			"objective":    Rastrigin,
			"dimensions":   13,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(13)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin14":{
			"objective":    Rastrigin,
			"dimensions":   14,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(14)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin15":{
			"objective":    Rastrigin,
			"dimensions":   15,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(15)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin16":{
			"objective":    Rastrigin,
			"dimensions":   16,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(16)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin17":{
			"objective":    Rastrigin,
			"dimensions":   17,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(17)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin18":{
			"objective":    Rastrigin,
			"dimensions":   18,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(18)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin19":{
			"objective":    Rastrigin,
			"dimensions":   19,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(19)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin20":{
			"objective":    Rastrigin,
			"dimensions":   20,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(20)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin21":{
			"objective":    Rastrigin,
			"dimensions":   21,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(21)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin22":{
			"objective":    Rastrigin,
			"dimensions":   22,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(22)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin23":{
			"objective":    Rastrigin,
			"dimensions":   23,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(23)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin24":{
			"objective":    Rastrigin,
			"dimensions":   24,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(24)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin25":{
			"objective":    Rastrigin,
			"dimensions":   25,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(25)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin26":{
			"objective":    Rastrigin,
			"dimensions":   26,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(26)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin27":{
			"objective":    Rastrigin,
			"dimensions":   27,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(27)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin28":{
			"objective":    Rastrigin,
			"dimensions":   28,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(28)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin29":{
			"objective":    Rastrigin,
			"dimensions":   29,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(29)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin30":{
			"objective":    Rastrigin,
			"dimensions":   30,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(30)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin31":{
			"objective":    Rastrigin,
			"dimensions":   31,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(31)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin32":{
			"objective":    Rastrigin,
			"dimensions":   32,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(32)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin33":{
			"objective":    Rastrigin,
			"dimensions":   33,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(33)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin34":{
			"objective":    Rastrigin,
			"dimensions":   34,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(34)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin35":{
			"objective":    Rastrigin,
			"dimensions":   35,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(35)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin36":{
			"objective":    Rastrigin,
			"dimensions":   36,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(36)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin37":{
			"objective":    Rastrigin,
			"dimensions":   37,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(37)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin38":{
			"objective":    Rastrigin,
			"dimensions":   38,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(38)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin39":{
			"objective":    Rastrigin,
			"dimensions":   39,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(39)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin40":{
			"objective":    Rastrigin,
			"dimensions":   40,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(40)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin41":{
			"objective":    Rastrigin,
			"dimensions":   41,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(41)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin42":{
			"objective":    Rastrigin,
			"dimensions":   42,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(42)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin43":{
			"objective":    Rastrigin,
			"dimensions":   43,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(43)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin44":{
			"objective":    Rastrigin,
			"dimensions":   44,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(44)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin45":{
			"objective":    Rastrigin,
			"dimensions":   45,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(45)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin46":{
			"objective":    Rastrigin,
			"dimensions":   46,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(46)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin47":{
			"objective":    Rastrigin,
			"dimensions":   47,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(47)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin48":{
			"objective":    Rastrigin,
			"dimensions":   48,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(48)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin49":{
			"objective":    Rastrigin,
			"dimensions":   49,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(49)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin50":{
			"objective":    Rastrigin,
			"dimensions":   50,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(50)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin51":{
			"objective":    Rastrigin,
			"dimensions":   51,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(51)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin52":{
			"objective":    Rastrigin,
			"dimensions":   52,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(52)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin53":{
			"objective":    Rastrigin,
			"dimensions":   53,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(53)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin54":{
			"objective":    Rastrigin,
			"dimensions":   54,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(54)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin55":{
			"objective":    Rastrigin,
			"dimensions":   55,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(55)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin56":{
			"objective":    Rastrigin,
			"dimensions":   56,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(56)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin57":{
			"objective":    Rastrigin,
			"dimensions":   57,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(57)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin58":{
			"objective":    Rastrigin,
			"dimensions":   58,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(58)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin59":{
			"objective":    Rastrigin,
			"dimensions":   59,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(59)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin60":{
			"objective":    Rastrigin,
			"dimensions":   60,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(60)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin61":{
			"objective":    Rastrigin,
			"dimensions":   61,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(61)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin62":{
			"objective":    Rastrigin,
			"dimensions":   62,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(62)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin63":{
			"objective":    Rastrigin,
			"dimensions":   63,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(63)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin64":{
			"objective":    Rastrigin,
			"dimensions":   64,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(64)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin65":{
			"objective":    Rastrigin,
			"dimensions":   65,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(65)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin66":{
			"objective":    Rastrigin,
			"dimensions":   66,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(66)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin67":{
			"objective":    Rastrigin,
			"dimensions":   67,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(67)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin68":{
			"objective":    Rastrigin,
			"dimensions":   68,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(68)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin69":{
			"objective":    Rastrigin,
			"dimensions":   69,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(69)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin70":{
			"objective":    Rastrigin,
			"dimensions":   70,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(70)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin71":{
			"objective":    Rastrigin,
			"dimensions":   71,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(71)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin72":{
			"objective":    Rastrigin,
			"dimensions":   72,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(72)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin73":{
			"objective":    Rastrigin,
			"dimensions":   73,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(73)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin74":{
			"objective":    Rastrigin,
			"dimensions":   74,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(74)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin75":{
			"objective":    Rastrigin,
			"dimensions":   75,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(75)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin76":{
			"objective":    Rastrigin,
			"dimensions":   76,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(76)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin77":{
			"objective":    Rastrigin,
			"dimensions":   77,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(77)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin78":{
			"objective":    Rastrigin,
			"dimensions":   78,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(78)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin79":{
			"objective":    Rastrigin,
			"dimensions":   79,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(79)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin80":{
			"objective":    Rastrigin,
			"dimensions":   80,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(80)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin81":{
			"objective":    Rastrigin,
			"dimensions":   81,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(81)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin82":{
			"objective":    Rastrigin,
			"dimensions":   82,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(82)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin83":{
			"objective":    Rastrigin,
			"dimensions":   83,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(83)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin84":{
			"objective":    Rastrigin,
			"dimensions":   84,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(84)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin85":{
			"objective":    Rastrigin,
			"dimensions":   85,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(85)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin86":{
			"objective":    Rastrigin,
			"dimensions":   86,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(86)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin87":{
			"objective":    Rastrigin,
			"dimensions":   87,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(87)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin88":{
			"objective":    Rastrigin,
			"dimensions":   88,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(88)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin89":{
			"objective":    Rastrigin,
			"dimensions":   89,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(89)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin90":{
			"objective":    Rastrigin,
			"dimensions":   90,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(90)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin91":{
			"objective":    Rastrigin,
			"dimensions":   91,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(91)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin92":{
			"objective":    Rastrigin,
			"dimensions":   92,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(92)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin93":{
			"objective":    Rastrigin,
			"dimensions":   93,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(93)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin94":{
			"objective":    Rastrigin,
			"dimensions":   94,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(94)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin95":{
			"objective":    Rastrigin,
			"dimensions":   95,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(95)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin96":{
			"objective":    Rastrigin,
			"dimensions":   96,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(96)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin97":{
			"objective":    Rastrigin,
			"dimensions":   97,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(97)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin98":{
			"objective":    Rastrigin,
			"dimensions":   98,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(98)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin99":{
			"objective":    Rastrigin,
			"dimensions":   99,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(99)],
			"result":               1# added 1 to eval to calculate error
	},
	"Rastrigin100":{
			"objective":    Rastrigin,
			"dimensions":   100,
			"search_space": [-5,5],
			"estimates":    [[ 0, 5] for i in range(100)],
			"result":               1# added 1 to eval to calculate error
	},

	"Ackley":{
		"objective":	Ackley,
		"dimensions":	2,
		"search_space":	[-5,5],
		"estimates":	[[ 0, 5] for i in range(2)],
		"result":		1	# added 1 to eval to calculate error
	},
	"Beale":{
		"objective":	Beale,
		"dimensions":	2,
		"search_space":	[-5,5],
		"estimates":	[[ 0, 5] for i in range(2)],
		"result":		1	# added 1 to eval to calculate error
	},
	"GoldsteinPrice":{
		"objective":	GoldsteinPrice,
		"dimensions":	2,
		"search_space":	[-5,5],
		"estimates":	[[ 0, 5] for i in range(2)],
		"result":		3	# added 1 to eval to calculate error
	},
	"Booth":{
		"objective":	Booth,
		"dimensions":	2,
		"search_space":	[-5,5],
		"estimates":	[[ 0, 5] for i in range(2)],
		"result":		1	# added 1 to eval to calculate error
	}
	# TODO - MATYAS FNS DOWN

}
