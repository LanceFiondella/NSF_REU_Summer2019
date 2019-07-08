import numpy as np
from math import *
import scipy, time, warnings

warnings.filterwarnings("ignore")

#Weibull model
data = [479,745,1022,1576,2610,2859,3552,4149,4266,4436,4553,5827,6296,7470,8163,10071,10206,10483,11079,11836,12273,14503,14940,15280,15685,16220,16497,16860,17382,17995,18272,19572,20393,20606,22226,23827,24125,24999,25617,28257,28262,28411,29445,31886,32346,32911,34030,34467,35394,39856,40570,40751,42236,42993,46147,48262,49146,51183,52664,53223,53713,54306,56075,56160,58996,59209,61075,61565,63052,67374,68792,69815,75305,76825,80106,82822,84997,88502,89227,91190,95169,96259,96504,97698,98692,102594]
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
	if (r < 0):
		print(b,c,r)
	return r

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
		sum1=1-((1-b)**(exp(EVec[i]*b1)*exp(FVec[i]*b2)*exp(CVec[i]*b3)))
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
	
	return -(firstTerm + secondTerm + thirdTerm - fourthTerm)

#--- ECM calculation ----------------------------------------------------------------
def logL(b,c):
	aHat = n / (1-np.exp(-b*(tn**c)))
	return (-n +sum(np.log(aHat*b*c*np.power(data.FT,(c-1))*np.exp(-b*np.power(data.FT,c)))))

def bMLE(b):
	aMLE = n/ (1 - np.exp(-b * np.power(tn,c)))
	return -aMLE*np.power(tn,c)*np.exp(-b *np.power(tn,c))+ sum((1-b*np.power(data.FT,c))/(b))

def cMLE(c):
	aMLE = n/ (1 - np.exp(-b * np.power(tn,c)))
	return -aMLE*b*np.power(tn,c)*np.log(tn)*np.exp(-b*np.power(tn,c)) + (n/c) + sum(np.log(data.FT)-b*np.log(data.FT)*np.power(data.FT,c))

def ECM():

	TimeList = []
	LLList = []
	bList = []
	cList = []
	jList = []

	for i in range(100):
		stime = time.perf_counter()
		brule = [n/np.sum(data)]
		crule = [1.0]

		ll_list = [logL(brule[0], crule[0])]
		ll_error_list = []
		ll_error = 1
		j = 0

		while ll_error > 1e-10:
			c = crule[j]
			b_est = scipy.optimize.fsolve(bMLE, x0 = brule[j])
			brule.append(b_est)
			del c

			b = brule[j+1]
			c_est = scipy.optimize.fsolve(cMLE, x0 = crule[j])
			crule.append(c_est)
			del b

			ll_list.append(logL(b_est, c_est))
			j += 1
			ll_error = ll_list[j] - ll_list[j-1]
			ll_error_list.append(ll_error)

		TimeList.append(time.perf_counter() - stime)
		LLList.append(ll_list[-1])
		bList.append(brule[-1])
		cList.append(crule[-1])
		jList.append(j)

		return 0	#??????

def Sphere(vector):
	return np.sum(np.square(vector)) + 1

models = {
	"Weibull":{
		"objective":	RLLWei,
		"dimensions":	2,
		"search_space":	[-1,1],
		"estimates":	[
							[0.1, 0.1],
							[0.9, 0.1]#[404/1000000, 1],
						],
		"result":		686.34571
	},
	"Covariate":{
		"objective":	RLLCV,
		"dimensions":	4,
		"search_space":	[0,1],
		"estimates":	[
							[0.2, 0.2],
							[0.2, 0.2],
							[0.2, 0.2],
							[0.2, 0.2]
						],#[0.1,0,0,0],
		"result":		23.0067
	},
	"Sphere10":{
		"objective":	Sphere,
		"dimensions":	2,
		"search_space":	[-1,1],
		"estimates":	[[ 0, 1] for i in range(2)],
		"result":		1
	}
}


if __name__ == "__main__":

	print(RLLCV([0.05775846, 0.14350439,0.30174806, 0.32566799]))