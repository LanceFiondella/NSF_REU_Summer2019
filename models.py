import sys, numpy as np
from math import *
import scipy, time, warnings
import scipy.optimize

model = sys.argv[1]

warnings.filterwarnings("ignore")



#Weibull model
sets = {
	"sys1": [3,33,146,227,342,351,353,444,556,571,709,759,836,860,968,1056,1726,1846,1872,1986,2311,2366,2608,2676,3098,3278,3288,4434,5034,5049,5085,5089,5089,5097,5324,5389,5565,5623,6080,6380,6477,6740,7192,7447,7644,7837,7843,7922,8738,10089,10237,10258,10491,10625,10982,11175,11411,11442,11811,12559,12559,12791,13121,13486,14708,15251,15261,15277,15806,16185,16229,16358,17168,17458,17758,18287,18568,18728,19556,20567,21012,21308,23063,24127,25910,26770,27753,28460,28493,29361,30085,32408,35338,36799,37642,37654,37915,39715,40580,42015,42045,42188,42296,42296,45406,46653,47596,48296,49171,49416,50145,52042,52489,52875,53321,53443,54433,55381,56463,56485,56560,57042,62551,62651,62661,63732,64103,64893,71043,74364,75409,76057,81542,82702,84566,88682],
	"sys2": [479,745,1022,1576,2610,2859,3552,4149,4266,4436,4553,5827,6296,7470,8163,10071,10206,10483,11079,11836,12273,14503,14940,15280,15685,16220,16497,16860,17382,17995,18272,19572,20393,20606,22226,23827,24125,24999,25617,28257,28262,28411,29445,31886,32346,32911,34030,34467,35394,39856,40570,40751,42236,42993,46147,48262,49146,51183,52664,53223,53713,54306,56075,56160,58996,59209,61075,61565,63052,67374,68792,69815,75305,76825,80106,82822,84997,88502,89227,91190,95169,96259,96504,97698,98692,102594],
	"sys3": [39,49,53,89,93,98,102,193,242,243,268,269,273,303,345,354,403,447,479,482,560,561,591,796,801,930,1033,1257,1443,1496,1510,1519,1521,1531,1532,1566,1736,1865,1869,1873,1908,1913,1918,1940,1976,2011,2132,2155,2188,2236,2268,2289,2293,2316,2325,2338,2503,2517,2539,2580,2592,2730,2825,2874,2936,2938,2973,3062,3152,3221,3243,3258,3277,3319,3333,3344,3385,3595,3611,3641,3678,3744,3753,3769,3783,3807,3819,3978,4067,4185,4214,4235,4253,4255,4369,4406,4452,4469,4470,4620,5002,5162,5228,5434,5443,5469,5531,5770,5783,5787,5872,5957,6197,6375,6409,6511,6520,6666,6725,6773,6798,6823,6934,6939,6970,7021,7027,7220,7247,7272,7368,7394,7424,7454,7471,7791,7869,7908,7921,7934,7953,8081,8115,8199,8239,8416,8765,9039,9121,9179,9210,9324,9363,9451,9535,9767,9875,9913,9999,10006,10028,10108,10347,10350,10389,10452,10604,10667,10747,10992,11188,11234,11386,11488,11497,11725,11945,12153,12231,12234,12317,12323,12535,12626,12629,12639,12811,12832,13005,13376,13416,13464,13590,13680,13829,13859,14176,14676,15349,15781,15847,16015,16081,16147,16275,16324,16656],
	"s2":   [191,413,693,983,1273,1658,2228,2838,3203,3593,3868,4228,5028,6238,6645,6695,7355,8862,9487,10399,11037,11330,12542,13154,13829,15044,17759,21310,22110,26020,32920,36220,37730,37925,39881,40016,40677,40727,41456,42356,42536,46761,62361,62361,62361,62661,71682,74201,81091,84439,87189,93864,100809,108708],
	"s27":  [20336,32112,73045,107839,124975,273421,281416,283052,298882,320814,323299,334299,337179,398361,403161,441166,457366,463366,464366,474366,474586,510166,591166,1234261,1282118,1436918,1607378,1715918,1789718,1791578,2128178,2396318,2471198,2757398,2782718,2789798,2849618,2937518,3013718,3102998,4312598],
	"ss3":  [107.4,124.62,124.8,157.68,158.64,184.74,228.9,562.62,580.44,621.3,640.08,641.04,642,721.86,722.1,722.22,724.02,724.5,725.28,762.54,764.64,836.7,1095.404,1095.884,1117.784,1596.404,1677.164,1678.364,1759.064,2447.924,2450.144,3209.024,3375.644,3383.924,4335.278,4336.598,4351.298,4354.718,4357.238,4519.718,5040.038,5136.758,5554.958,5989.718,6533.498,6542.318,7030.598,7031.078,7031.618,7033.838,7034.918,7172.258,7264.118,7286.918,7309.838,7783.178,8138.079,8507.559,8887.779,9736.419,9736.539,9739.955,9814.115,10076.615,10955.915,10956.275,10964.435,10964.615,11202.535,11202.655,11273.455,11286.415,11286.715,11286.835,11845.375,12033.415,12089.695,12090.115,12504.579,12745.359,12951.999,12956.739,12966.879,12967.179,12971.319,13443.399,13443.699,13531.299,13579.539,13621.479,14198.091,14269.911,14353.011,14353.911,14594.211,14667.951,14837.751,14837.752,15140.032,15143.392,15145.732,15227.992,15787.912,15788.692,15799.432,15799.612,16230.472,16397.212,16397.812,16773.952,16779.052,17328.592,17329.132,17330.032,17851.284,17851.704,18370.344,18371.364,18375.504,18375.984,18376.164,18376.764,18430.524,18512.964,18513.144,18786.144,18846.024,18846.864,18854.004,18930.324,19079.004,19316.844,19321.404,19323.324,19340.184,19417.224,19491.984,20230.164,20377.164,20453.844,20524.644,20590.824,20618.364,20673.384,20673.504,20970.3,21060.48,21785.04,21952.14,22058.34,22058.82,22176.18,22182.66,22182.72,22280.58,22679.16,23070.54,23070.72,23070.9,23071.14,23071.68,23408.58,23673.06,24520.14,24546.6,24895.92,24900,24964.68,24965.52,24966.06,25556.04,25888.32,25982.46,26222.52,26225.22,26226.12,26227.2,26238.78,26240.94,26433.66,26521.5,26605.86,26983.98,27042.48,27126.36,27285,27285.66,27288.84,27290.4,27293.58,27299.28,27525.84,27535.68,27604.74,27673.62,27739.08,28141.98,28217.46,28597.68,29302.648,29808.328,29862.748,30181.768,30276.988,30282.088,30288.328,30337.768,30338.188,31005.508,31005.628,31012.828,31081.768,31108.588,31557.208,31896.628,31897.108,32939.788,33719.368,33727.408,34885.648,35792.788,35851.288,36235.228,38274.688,38796.928,38862.928,38906.428,38908.468,38909.068,39135.388,39462.988,39664.288,39891.268,40444.708,40445.728,40446.688,40959.448,41778.688,42580.348,42740.728,42812.368,43176.358,43185.448,43413.418,43430.608,44028.508,44717.908,44729.428,44753.278,44829.148,44952.178,44978.188,45053.428,45121.558,45932.608,46430.968,47054.248,47057.578,47064.868,47112.028,48440.428,48550.228,48894.118,50509.978,50524.918,51205.678,51231.898,51608.008,51789.898,51854.218,52322.398,53890.978,54224.698,54224.878,54225.688,54547.798,54569.758,54933.358],
	"csr1": [60,90,630,697,737,760,765,818,822,838,932,947,952,1042,1119,1187,1202,1339,1362,1363,1467,1483,1492,1502,1514,1518,1528,1610,1616,1670,1695,1738,1750,1798,1821,1827,1837,1849,1863,1896,1905,1909,1975,1975.5,1993.5,2008.5,2083.5,2113.5,2229.5,2243.5,2258.5,2299.5,2300.5,2399.5,2408.5,2453.5,2521.5,2557.5,2607.5,2688.5,2777.5,2862.5,2916.5,2919.5,2934.5,2940.5,2948.5,2984.5,3082.5,3114.5,3150.5,3155.5,3164.5,3224.5,3258.5,3274.5,3438.5,3561.5,3580.5,3599.5,3725.5,3761.5,3815.5,3830.5,3846.5,4000.5,4084.5,4176.5,4423.5,4667.5,4727.5,4732.5,4734.5,4739.5,4869.5,4870,5103,5153,5207,5259,5316,5317,5319,5324,5325,5342,5356,5358,5445,5464,5493,5493.5,5554.5,5672.5,5692.5,5695.5,5706.5,5793.5,5798.5,6047.5,6075.5,6119.5,6150.5,6153.5,6163.5,6166.5,6174.5,6191.5,6194.5,6249.5,6256.5,6268.5,6274.5,6278.5,6447.5,6477.5,6481.5,6519.5,6526.5,6530.5,6534.5,6547.5,6557.5,6597.5,6654.5,6676.5,6713.5,6840.5,6852.5,6860.5,6861.5,6882.5,6986.5,6994.5,7017.5,7018.5,7024.5,7165.5,7168.5,7189.5,7192.5,7210.5,7211,7286,7378,7417,7446,7449,7451,7609,7639,7663,7668,7760,7767,7800,7820,7836,8128,8131,8140,8152,8170,8179,8254,8269,8315,8324,8418,8443,8618,8623,8635,8653,8723,8726,8736,8850,9063,9275,9279,9284,9390,9654,9923,10199,10200,10403,10520,10521,10566,10571,10681,10699,10709,10888,10954,10955,11061,11063,11082,11199,11229,11359,11390,11418,11422,11724,12086,12091,12154,12196,12282,12540,12834,13090,13208,13221,13268,13360,13703,13831,14223,14313,14429,14464,14635,14774,14884,14982,15042,15132,15214,15219,15249,15284,15515,15577,15735,17357,17710,17743,17813,17848,17964,18773,20483,21228,21578,22048,22170,22414,24798,25047,25654,25737,25739,25765,26351,26703,27376,27706,28355,28478,30267,31555,31666,31741,31815,32148,32435,32436,33317,33330,34644,35116,35479,35485,35489,35544,35953,35989,36004,37408,37425,37496,37530,37964,38024,38043,38063,38142,38166,38706,39746,39784,39862,39903,41660,41865,43960,44748,44749,47417,47887,47897,47917,48255,48477,48505,48561,49122,49187,49287,50187,50399,50686,50739,50742,55715,55912,57086,57869,59215,59274,59372,60966,60991,61089,61811,62039,62117,62150,62603,63623,67950,68875,69177,69826,69869,70054,70211,70241,72012,73100,73656,73711,78603,78684,78745,79221,79284,79287,79290,79352,79353,79397,80633,82039,82148,83619,85416,87165,89261,89337,91504,93563,95740,97633,97831,101157,104257,104843,107529,107653,107882,108890],
	"csr2": [760,1518,1821,1827,1849,1863,1905,1909,1993,2008,2229,2243,2258,2299,2300,2453,2862,2916,2940,2984,3164,3561,3580,3725,3761,3815,5152,5315,5323,5324,5341,5357,5444,5463,5492,5492,5497,5857,5867,5878,5978,6230,6690,6869,6872,6896,7149,7312,7366,7503,7831,7834,7843,7855,7873,7882,7957,7972,8338,8766,8978,9093,9357,9626,9902,9903,10902,10932,11427,11899,12243,12793,12924,12971,13063,13926,14917,14952,24501,24750,25357,25440,26054,26406,27079,31258,31369,31444,31851,32139,33033,34347,35192,35247,35656,35692,35707,37667,37727,37746,37766,37845,37869,39606,47590,47600,47620,47958,48208,49890,50102,50389,50445,55418,58918,58977,59075,61514,63326,69529,69914,73414,78306,78993,79055,81851,85119,88964,89040],
	"csr3": [33,42,46,112,112.5,130.5,279.5,293.5,308.5,358.5,439.5,473.5,558.5,612.5,615.5,630.5,636.5,644.5,774.5,793.5,812.5,924.5,939.5,955.5,1109.5,1159.5,1169.5,1171.5,1193.5,1246.5,1265.5,1323.5,1343.5,1346.5,1438.5,1443.5,1509.5,1798.5,1801.5,1810.5,1822.5,1840.5,1849.5,1924.5,1939.5,2230.5,2442.5,2446.5,2451.5,2759.5,3028.5,3304.5,3305.5,3705.5,3999.5,4226.5,4344.5,4357.5,4404.5,4493.5,4735.5,4834.5,5441.5,5524.5,5526.5,5552.5,6138.5,6846.5,6852.5,6856.5,6911.5,7320.5,7356.5,7371.5,7944.5,8527.5,8587.5,8606.5,8626.5,8705.5,8729.5,9269.5,9321.5,10917.5,11231.5,11232.5,11995.5,12005.5,12025.5,12169.5,12197.5,12253.5,12729.5,12794.5,12892.5,13776.5,13988.5,14275.5,14328.5,14331.5,15162.5,15205.5,15260.5,15369.5],
	#"ss4":  [242460,248760,502440,526560,564180,1492560,2026740,2218440,2230200,2389440,2408580,3004560,3480780,3539760,3545580,3646200,3652500,3980640,3985200,3991980,4061400,4640520,5217600,5849520,5849640,5997840,6925980,6969840,6973440,6973740,7634367,7802247,7989087,8268131,8688251,8756051,9079511,9191531,9196811,9831072,9832992,10735332,10737432,10741632,10742112,10919172,10961532,10998552,11053092,11388612,11388912,11540952,11652312,12144072,12224412,12542472,12545592,12550992,12735732,12751152,12852012,13087148,13088768,13710428,14270656,14775076,14781856,15027676,15509964,15939624,16935160,17010280,17470300,17502460,17704960,18022780,18288040,18324220,18498400,18518140,18526540,18566800,18569740,18578440,18585520,18649540,18649900,18734740,18734920,19643260,19650940,19709800,19732360,19740640,20309736,20379396,20673276,20680176,20939436,20943216,21254992,21366028,21779968,22051588,22116148,22415548,22448428,23406388,24404168,24404228,24517568,24984548,25495748,25771748,25956908,26280848,26599268,27377228,27381428,27954568,27956788,29194588,29417548,29525848,29584828,29604928,29687248,30049408,30050008,30059008,31304668,31321468,32322508,32322656,32322772,32323012,32323492,32324332,32743672,32743792,32803852,32890732,32896912,32897092,32995192,32998492,33340732,34452052,34460452,34552132,35173912,36049612,36684232,36684652,37064032,37082032,37144132,37152952,37656532,37932352,38280532,38285632,38763532,39399952,39975292,40197652,40284232,40286992,40289812,40291372,40291792,40292992,40294432,40401472,40432912,40737172,40796272,41064232,41085952,41154712,41417092,41421652,41752672,43474972,43522492,43823662,43824682,43992682,44105992,45719542,45811732,46258072,46323322,48223942,48644002,48889282,50236822]
}	# ss4 omitted: does not converge

setsres = {
	"sys1": 966.0803348790863,
	"sys2": 686.2916747961536,
	"sys3": 1092.9987779219996,
	"s2":	446.22654018382457,
	"s27":	1302.971448570725, 
	"ss3":	1731.3035150212263, 
	"csr1":	2379.348044059513, 
	"csr2":	920.2931069655324, 
	"csr3":	600.9963885102122
}

if "-s" in sys.argv:
	setchoice = sys.argv[sys.argv.index("-s")+1]
	print(setchoice)
else:
	setchoice = "sys1"

data = np.array([float(j) for j in sets[setchoice]])
data_res = setsres[setchoice]

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
		"result":		data_res
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
			"result":               1# added 1 to eval; does not affect actual run, just error calculation
	}

}
