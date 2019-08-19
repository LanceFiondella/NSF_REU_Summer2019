# TODO - explore variants of best ones (bee, bat, pso, etc)
# TODO - separate each model in file into model folder

# TODO - use both CV datasets
# TODO - add estimate fn for wei (a = n, b = n/sum(bi), c = 1)
# TODO - before illustrations, talk about methods then covariates
# TODO - table of NSGA variables


# TODO - population density plots but for rastr with d = 10, 20, 30


import sys, os, models, csv, multiprocessing as mp, matplotlib.pyplot as plt

sys.path.insert(0, f"algorithms")
import bat, bee, cuckoo, firefly, fish, pollination, pso, wolf

#---- NSGA SETTINGS ------------------------------

nsga_max_gens = 64				# number of generations in NSGA-II
nsga_pop_size = 128 				# how many combinations there are, must be even
nsga_p_cross = 0.98				# mutation crossover probability
nsga_fn_evals = 16				# how many evaluations to average
nsga_bpp = 32 					# bits / precision to use for each parameter
nsga_cross_breed = False		# allow algorithms to change during the process
nsga_free_range = False			# sets all parameter ranges to 0.5 +- 0.5 instead of predefined (needs more gens)
nsga_cross_head_tail = True		# uses head/tail crossover instead of random bits

model = models.models[sys.argv[1] if len(sys.argv) > 1 else 'Covariate']
model_pop_count = [22, 16]		# pop size for method
model_generations = [22, 16]	# generations used in method


stage1 = [lambda x: x]*2
stage2 = [	{	
				"algo":	firefly.search,
				"params":[
					[0.95, 0.05],
					[0.96, 0.04]
				]
			},
			{						# algo points to search function in swarm file
				"algo":	pso.search,
				"params":[			# contains all constant parameters in terms of 
					[0.5, 0.15],	# [A, B] where the passed value is A + B*f(bits)
					[0.1, 0.05],	
					[0.1, 0.05]
				]
			},
			{
				"algo": cuckoo.search,
				"params":[
					[0.97, 0.03],
					[0.97, 0.03],
					[0.25, 0.1]
				]
			},
			{
				"algo": bat.search,
				"params":[
					[0.1, 0.1],
					[0.9, 0.1],
					[0.9, 0.08],
					[0.8, 0.15]

				]
			},
			{
				"algo": bee.search,
				"params":[
					[0.33, 0.10],
					[0.33, 0.10],
					[0.50, 0.25],
					[0.50, 0.25],
					[0.20, 0.10]
				]
			},
			{
				"algo":	pollination.search,
				"params":[
					[0.7, 0.2]
				]
			},
			{
				"algo": fish.search,
				"params":[
					[0.125, 0.075],
					[4, 2],
					[0.125, 0.075],
					[2, 1],
				]
			},
			{
				"algo": wolf.search,
				"params":[
					[0.5, 0.1],
					[0.85, 0.15],
					[0.5, 0.1],
					[0.9, 0.075]
				]
			},
			{
				"algo": None,
				"params":[]
			}
			]
stage3 = [
			models.nelder_mead, 
			models.powell, 
			models.cg,
			models.bfgs, 
			models.lbfgsb, 
			models.tnc, 
			#models.cobyla, 
			models.slsqp 
			#models.dogleg 
			#models.trustncg
		]

#---------------begin NSGA-II---------------------


import numpy as np
import time


def search(max_gens, pop_size, p_cross, result_block = None):
	'''
	Main routine - creates population of algorithms,
	runs initial sorting and new population, iteratively
	runs re-evaluation until convergence
	'''
	param_count = len(max(stage2, key = lambda x: len(x["params"]))["params"])

	total_bit_count = nsga_bpp * (param_count + 5)
											# get total number of bits for bitstring

	pop = [	{"bitstring":random_bitstring(total_bit_count)} 
			for i in range(pop_size)]		# bitstring consists of algo bits, parameter bits, pop size bits, and gen count bits
										
	print(" > begin initial objectives & sort")
	calculate_measures(pop)					# prelim. results pt 1: get initial pop scoring
	fast_nondominated_sort(pop)				# sort initial pop by domination

	selected = [better(pop[np.random.randint(pop_size)], pop[np.random.randint(pop_size)]) for i in range(pop_size)]
	children = reproduce(selected, pop_size, p_cross)
											# select a random well-fit parent pop and generate a new child pop with it

	print(" > objectives of first child pop")
	calculate_measures(children)
	print(" > starting generations")		# prelim. results pt 2: score child pop

	for gen in range(max_gens):				# begin generational aspect

		
		t = {'NONE':0}								# code to visualize population percentages of each algorithm
		for i in [ x[:-3] for x in os.listdir('algorithms') if x[-3:] == '.py']:
			t[i] = 0
		for p in pop:
			_,b,_,_ = decode(p["bitstring"])
			m = b.__module__ if b != None else 'NONE'
			t[m] += 1
		algs.append(t)
		

		union = pop + children				# sort union of parent/child pops, create new pop based on least dominated
		fronts = fast_nondominated_sort(union)
		parents = select_parents(fronts, pop_size)
										
		selected = [better(parents[np.random.randint(pop_size)], parents[np.random.randint(pop_size)])
					for i in range(pop_size)]
											# select well-fit parents at random
		
		pop = children.copy()
		children = reproduce(selected, pop_size, p_cross)

		calculate_measures(children, gen, max_gens)
		snapshots.append([(t['bitstring'], t['objectives']) for t in children])
		res.append([t['objectives'] for t in children])
		pops.append([(len(decode(t['bitstring'])[3][3]),np.prod(t['objectives'])) for t in children])
		gens.append([(decode(t['bitstring'])[3][2],np.prod(t['objectives'])) for t in children])


	union = pop + children
	fronts = fast_nondominated_sort(union)
	parents = select_parents(fronts, pop_size)
	if result_block == None:
		return parents
	else:
		result_block.put(parents)
		return


def random_bitstring(num_bits):
	'''
	Generates a binary string of some length,
		no limit on size 
	'''
	return "".join([str(np.random.randint(2)) for i in range(num_bits)]).zfill(num_bits)


def calculate_measures(pop, gen=None, maxgen=None):
	'''
	Given a population member, decodes and evaluates the
	time and accuracy of its algorithms 
	'''
	for index, p in enumerate(pop):			# find fitness of each member of a population in order to find pareto-fitness
		
		alg_1, alg_2, alg_3, params = decode(p["bitstring"])

		pop_size = len(params[3]) 			# get size of population to make new one
		expanded = list(params) 			# convert it to a list to replace it

		runtime, error = 0, 0


		for i in range(nsga_fn_evals):		# get average of N runs


			expanded[3] = [						# generate population, pseudorandom
				[t[0] + t[1]*np.random.uniform(-1, 1) for t in model["estimates"]]
				for x in range(pop_size)]

			
			params = tuple(expanded) 		# convert back to tuple

			col = int(os.popen('stty size', 'r').read().split()[1])-5
			print(' ' * col, end='\r')
			st = f" {str(gen+1).zfill(len(str(maxgen)))}/{maxgen}:" if gen != None else ""
			st = f"{st}\t > {str(index+1).zfill(len(str(len(pop))))}/{len(pop)} ({str(i+1).zfill(len(str(nsga_fn_evals)))}/{nsga_fn_evals})\t{alg_2.__module__ if (alg_2 != None) else 'NONE'}: {alg_3.__name__} {[round(x,3) for x in params[4:]]}"
			print(f" {st}", end='\r')

			stime = time.time()
			lst =  alg_2(*params) if alg_2 != None else params[3]
											# run swarm estimator if enabled, otherwise just take random values

			candidate = min(lst, key = model['objective'])
											# get most-fit particle, check if it converges

			result, conv = alg_3(candidate)	# do newton's method or ECM on result for runtime
											# don't care about the actual accuracy of that convergence
			runtime += time.time() - stime

			error += (model['objective'](result) / model['result'])
			
		p["objectives"] = [runtime/nsga_fn_evals, error/nsga_fn_evals]
											# metrics: total runtime, and accuracy of swarm before convergence
	col = int(os.popen('stty size', 'r').read().split()[1])
	print(' '*col, end='\r')


def decode(bitstring):
	'''
	Takes some bit-pattern (some member of the population),
	gives back the corresponding algorithms (phase 1, 2, 3)
	as well as the parameters passed to ph1
	'''
	groups = [bitstring[x*nsga_bpp:][:nsga_bpp] for x in range( int( len(bitstring) / nsga_bpp ) )]

	index_1 = round(bin_signum(groups[0], (len(stage1)-1)/2, (len(stage1)-1)/2))
	index_2 = round(bin_signum(groups[1], (len(stage2)-1)/2, (len(stage2)-1)/2))
	index_3 = round(bin_signum(groups[2], (len(stage3)-1)/2, (len(stage3)-1)/2))
											# get indexes of each stage, by rounding
											# stage_num / 2 +- stage_num / 2

	pop_size = round(bin_signum(groups[3], model_pop_count[0], model_pop_count[1]))
	gen_count = round(bin_signum(groups[4], model_generations[0], model_generations[1]))
											# get bits corresponding to model gens and size
											# and generate pop size / generations based on values

	model_search_space = [model["search_space"] for i in range(model["dimensions"])]

	params = (model["objective"], model_search_space, gen_count, [0]*pop_size)
											# get parameters for swarm algorithms
	param_bits = groups[5:]
	for idx, prm in enumerate(stage2[index_2]["params"]):
		this_bits = param_bits[idx]
											# get part of bitstring corresponding to parameter
		if nsga_free_range:
			params += (bin_signum(this_bits, 0.5, 0.5),)
		else:
			params += (bin_signum(this_bits, prm[0], prm[1]),)

	return stage1[index_1], stage2[index_2]["algo"], stage3[index_3], params	# send back algorithms & parameters for timing


def bin_signum(bitstring, lmbda, mu, bpp = nsga_bpp):
	'''
	Use a bitstring to get a number within lambda +- mu
	'''
	uniform = int(bitstring, 2) / (2 ** bpp - 1)
	uniform = 2 * uniform - 1
	return lmbda + uniform * mu


def fast_nondominated_sort(pop):
	'''
	Takes a population and sorts it into fronts, based on
	which candidates dominate other candidates. 
	'''
	fronts = [[]]                       # generate a list of fronts for the next population
	for p1 in pop:                      # lower "rank" value indicates dominated by less and therefore better
		p1["dom_count"], p1["dom_set"] = 0, []
		for p2 in pop:
			if dominates(p1, p2):		# get list of dominated pop members & count for each
				p1["dom_set"].append(p2)
			elif dominates(p2, p1):
				p1["dom_count"] += 1
		if p1["dom_count"] == 0:
			p1["rank"] = 0
			fronts[0].append(p1)
	curr = 0
	while True:
		next_front = []
		for p1 in fronts[curr]:
			for p2 in p1["dom_set"]:	# create fronts based on rank
				p2["dom_count"] -= 1
				if p2["dom_count"] == 0:
					p2["rank"] = (curr + 1)
					next_front.append(p2)
		curr += 1
		if len(next_front) > 0:
			fronts.append(next_front)
		if not (curr < len(fronts)):
			break
	return fronts


def dominates(p1, p2):
	'''
	Decides if a candidate p1 is dominant, i.e.
	all objectives are no worse than the other p2
	'''
	for i in range(len(p1["objectives"])): # all traits are no worse than other
		if p1["objectives"][i] > p2["objectives"][i]:	
			return False
	return True


def point_mutation(bitstring, rate=None):
	'''
	Mutates a bitstring according to some rate,
	(flips bits at random) - will not change the
	algorithm-choice bits if not enabled
	'''
	if rate == None:					# flip bits at random according to rate
		rate = 1/len(bitstring)
	child = ""
	for i in range(len(bitstring)):
		bit = bitstring[i]				# only change fn bits if crossbreed is enable
		if nsga_cross_breed or not (nsga_bpp < i < 2*nsga_bpp):
			child += str(1-int(bit)) if (np.random.random_sample()<rate) else bit
		else:
			child += bit
	return child


def reproduce(selected, pop_size, p_cross):
	'''
	Generates new sub-population based on a selection
	of an initial population with mutations and crossover
	of random members - if cross-breed is disabled, 
	parents will only have same algorithms
	'''
	children = []                       # generate new children population based off of parent population

	for i, p1 in enumerate(selected):   # with crossovers and mutation

		if(not nsga_cross_breed):
			count = 1
			while True:					# if no cross breed, get two parents with same algorithms
				p2 = selected[(i+count)%len(selected)]
				if p2["bitstring"][nsga_bpp:][:nsga_bpp] == p1["bitstring"][nsga_bpp:][:nsga_bpp]:
					break				# compare similar middle bitstrings
				count += 1
		else:
			p2 = selected[0] if (i == len(selected) - 1) else \
					selected[i+1] if i%2 == 0 else selected[i-1]

		child = {}						# do crossover and mutations to new child
		child["bitstring"] = crossover(p1["bitstring"], p2["bitstring"], p_cross)
		child["bitstring"] = point_mutation(child["bitstring"])
		children.append(child)
		if len(children) >= pop_size:
			break
	return children


def crossover(parent1, parent2, rate):
	'''
	Creates hybrid algorithm based on two parents,
	picks random bits from either, may rarely just
	use one parent entirely
	'''
	if np.random.random_sample() >= rate:
		return ""+parent1

	if nsga_cross_head_tail:
		centerpt = np.random.randint(len(parent1))
		return parent1[:centerpt] + parent2[centerpt:]

	child = ""	# otherwise use random bits method
	for i in range(len(parent1)):		# get random bits from either parent
		child += parent1[i] if np.random.random_sample() < 0.5 else parent2[i]
	return child


def calculate_crowding_distance(pop):
	'''
	Sets crowding distance of each population member,
	that is, the distance of its neighbors' objectives
	in order to preserve diverse options within the population.
	'''
	for p in pop:
		p["dist"] = 0
	num_obs = len(pop[0]["objectives"])
	for i in range(num_obs):
		mn = min(pop, key = lambda x: x["objectives"][i])
		mx = max(pop, key = lambda x: x["objectives"][i])
		rge = mx["objectives"][i] - mn["objectives"][i]
		pop[0]["dist"], pop[-1]["dist"] = float('inf'), float('inf')
		if rge == 0:
			continue
		for j in range(1, len(pop)-1):
			pop[j]["dist"] += (pop[j+1]["objectives"][i] - pop[j-1]["objectives"][i]) / rge


def better(x, y):
	'''
	Given two candidates, gives the better of the two
	sorted by rank and then crowding distance
	'''
	if ("dist" in x.keys()) and (x["rank"] == y["rank"]):
		return x if (x["dist"] > y["dist"]) else y
	return x if (x["rank"] < y["rank"]) else y


def select_parents(fronts, pop_size):
	'''
	Generates new parent population based off of fronts
	from sorting, in order to fill the rest of the pop
	with better candidates
	'''
	for f in fronts:
		calculate_crowding_distance(f)	# generate new parent pop and sub-pop for new pop
	offspring, last_front = [], 0
	for front in fronts:
		if (len(offspring) + len(front)) > pop_size:
			break
		for p in front:
			offspring.append(p)
		last_front += 1
	remaining = pop_size - len(offspring)
	if remaining > 0:
		fronts[last_front].sort(key = lambda x: (x["rank"], x["dist"]))
		offspring += fronts[last_front][0:remaining]
	return offspring

# ---------- DONE, RUN SEARCH ALGORITHM -------------------------------------------------------------

''' MULTITHREAD SEARCH TO FIND FRONTRUNNERS
qu = mp.Queue(100)
threads = []
t = time.time()

for i in range(100):
	thread = mp.Process(target = search, args = (nsga_max_gens, nsga_pop_size, nsga_p_cross, qu))
	threads.append(thread)
	thread.start()

populations = [qu.get() for i in range(100)]

for thread in threads:
	thread.join()

print(len(populations))
print('done')

frontrunners = []
for p in populations:
	frontrunners.append(min(p, key = lambda x: (x['objectives'][1], x['objectives'][0]) ).copy())
pop = populations[0]
#finals = []
#res = []
#algs = []
#pops = []
#gens = []
#pop = search(nsga_max_gens, nsga_pop_size, nsga_p_cross)
#	frontrunners.append(min(pop, key = lambda x: (x['objectives'][1], x['objectives'][0]) ).copy())
#	print(f"done iteration {i+1}")
#finals.append(algs[-1])
#algs = finals
print("done!\n")
'''

''' 
times = 0
for i in range(100):
	pop = search(nsga_max_gens, nsga_pop_size, nsga_p_cross)
	best = min(pop, key = lambda x: (x["objectives"][1], x["objectives"][0]))
	species = [decode(x['bitstring'])[1] for x in pop]
	most_frq = max(set(species), key = species.count)
	_, spc, _, _ = decode(best['bitstring'])
	if spc == most_frq:
		times += 1
print(times/100)
'''
res = []
pops = []
gens = []
algs = []
snapshots = []

t = time.time()
pop = search(nsga_max_gens, nsga_pop_size, nsga_p_cross)
t = time.time() - t
print(f'done! {str(int(np.floor(t/3600))).zfill(2)}:{str(int(np.floor(t/60) % 60)).zfill(2)}:{str(int(np.floor(t) % 60)).zfill(2)}')
# ----------- VISUALIZATION -------------------------------------------------------------------------

pop.sort(key = lambda x: (x["objectives"][1], x["objectives"][0]))	# sort by conv then time

colors = {
	'bat':		'#020202',
	'bee':		'#bbaf22',
	'cuckoo':	'#ec994a',
	'firefly': 	'#ff0000',
	'fish':		'#3d7bff',
	'pollination':'#ff3dff',
	'pso':		'#3dd6ff',
	'wolf':		'#e03e3e',
	'NONE':		'#00aa00'
}

sep = "\t"

print('\033[4m' + sep.join(['runtime', 'avg conv', 'algo', 'conv method', 'gens', 'pop', 'algo params', 'init estimate score']) + "\\\\" + '\033[0m')
for p in pop:
	r, r2, r3, params = decode(p["bitstring"])
	n = r2.__module__ if r2 != None else "NONE"
	c = colors[n]

	print(sep.join([str(round(x,12)).ljust(12, ' ') for x in p["objectives"]]), end=sep)
	print(f"{n[:6] if r2 != None else 'NONE'}{sep}{r3.__name__[:6] }{sep}{params[2] if r2 != None else 'n/a'}{sep}{len(params[3])if r2 != None else 'n/a'}", end=sep)

	if r2 != None:
		for i in params[5:]:
			print(round(i,3), end=sep)
	print("\\\\")

	#plt.title("1+ep vs time")
	#plt.plot(p["objectives"][0],p["objectives"][1],c+"o")


''' NONE TYPES GRAPH
types = {}
for p in pop:
	ptype, conv = decode(p['bitstring'])[1:][:2]
	if ptype == None:
		conv = conv.__name__
		if conv in types:
			types[conv] += 1
		else:
			types[conv] = 1
plt.bar(list(range(len(types))),
		[types[x] for x in types],
		tick_label = [str(x) for x in types.keys()])	
plt.show()
'''


''' PLOT RUNTIME CHANGES OVER DIFFERENT RUNS
for idx, pop in enumerate(pops):
	for p in pop:
		plt.plot(idx+2, p['objectives'][0], 'r*')
	plt.plot(idx+2, sum([x['objectives'][0] for x in pop])/nsga_pop_size, 'b')
	print(len(pop))

plt.show()
'''


''' PLOT FRONTRUNNERS OF 100 RUNS
names = [x[:-3] for x in os.listdir('algorithms') if x[-3:] == '.py']
names.append('NONE')
prog = {}
for n in names:
	prog[n] = 0
for p in frontrunners:
	_, s2, _, _ = decode(p["bitstring"])
	s2 = s2.__module__ if s2 != None else "NONE"
	if s2 != None:
		prog[s2] += 1
	else:
		prog['NONE'] += 1
klist = list(prog.keys()).copy()
for t in klist:
	if prog[t] == 0:
		del prog[t]
keys = prog.keys()
values = [prog[x] for x in keys]
plt.bar(list(range(len(keys))), values, align='center', alpha = 0.5)
plt.xticks(list(range(len(keys))), keys)
plt.ylabel('Frontrunner algorithm type (100 gens)')
plt.show()
'''


''' 			PLOT SIZE AND GENERATIONS VS ITERATIONS
plt.plot([sum([x[0] for x in l])/len(l) for l in pops] , 'b', label='pop size')
plt.plot([sum([x[0] for x in l])/len(l) for l in gens], 'r', label='gen count')
plt.title('Population Size & Swarm Generations vs Iterations (Average)')
plt.legend(loc='best')
plt.show()
'''


'''				PLOT TIME/EPSILON VS ITERATIONS
ers = []
tss = []
for iteration, pt in enumerate(res):
	er = 0
	ts = 0
	c = 1
	for memb in pt:
		plt.plot(iteration, memb[0],'r*')
		ts += memb[0]
		if (not np.isnan(memb[1])) and (0 < memb[1] < 3):
			plt.plot(iteration, memb[1]-1,'b*')
			er += memb[1]-1
			c += 1
	er /= c
	ts /= len(pt)
	ers.append(er)
	tss.append(ts)
plt.plot(tss, 'r',label='time')
plt.plot(ers, 'b', label='epsilon')
plt.title('Time & Error vs Iterations')
plt.legend(loc='best')
plt.show()
'''


names = [x[:-3] for x in os.listdir('algorithms') if x[-3:] == '.py']
names.append('NONE')
prog = {}
for n in names:
	prog[n] = [x[n] for x in algs]
for n in names:
	plt.plot([l/nsga_pop_size for l in prog[n]],colors[n],label=n)
plt.title('Algorithm percentage of population vs Iterations')
plt.legend(loc='best')
plt.show()


for p in pop:
	plt.plot(p['objectives'][0], p['objectives'][1], 'b*')
plt.title('Error vs Runtime Pareto Front')
plt.xlabel('Runtime (s)')
plt.ylabel('1+epsilon error')
plt.show()

with open("bitstrings.txt","w") as f:
	f.write(str(snapshots))

				# WRITE OUTPUTS TO CSV
with open('output_populations.csv','w') as csvfile:
	writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
	writer.writerow([f"NSGA GENS: {nsga_max_gens}",f"NSGA POP: {nsga_pop_size}",f"NSGA CROSSOVER: {nsga_p_cross}",f"NSGA AVG EVALS: {nsga_fn_evals}",f"NSGA CROSSBREED: {nsga_cross_breed}"])
	writer.writerow([f"MODEL: {[x for x in models.models if models.models[x]==model][0]}", f"MODEL POP: {model_pop_count}", f"MODEL GENS: {model_generations}"])
	writer.writerow([])
	writer.writerow(["algorithm", "convergence mtd", "runtime", "AVG error", "BEST error", "best candidate's model parameters", "score of best"])
	for p in pop:
		r, r2, r3, params = decode(p["bitstring"])
		pop_size = len(params[3])
		rep = list(params)
		rep[3] = [						# generate population, pseudorandom
				[t[0] + t[1]*np.random.uniform(-1, 1) for t in model["estimates"]]
				for x in range(pop_size)]

		params = tuple(rep)

		newl = r2(*params) if r2 != None else rep[3]
		bst = min(newl, key = model['objective'])
		rs, con = r3(bst)
		sc = model['objective'](rs) / model['result']
		writer.writerow([r2.__module__ if r2 != None else 'NONE', r3.__name__, p["objectives"][0],p["objectives"][1], sc, rs, model["objective"](rs)])