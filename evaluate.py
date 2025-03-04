
import sys, os, models

sys.path.insert(0, "algorithms")
import bat, bee, cuckoo, firefly, fish, pollination, pso, wolf

predef = [	# 3 pre-defined candidates for weibull
	"10010111000101111000101010110101101001001010010110100101100001111011101010011000101011100011111000000000011000011010000011001000000000001100001000011001111101100010010000010100001011110001111111010011011010100000110100011101101101101100100111110101101011011111100111110010100010010110100001110100110101111100001001010000",
	"1000011111110011011000100000111100000110110110001011000110001110011001000101101111010010011010100000001000111000010110110010111010000110101001100010101010111010",
	"1100011011011011011110100000111100000110110110001011000110001110011001000101100111110010011110100000011000111011110010110011110110000100101001100100101000111010"
]

predef = [	# 2 (really 3) pre-defined covariate candidates, comment out section if using weibull

	#"111111100101000101110110001110101000010001100010110010001001010010111010111011100010101000110101000000110101100100011000000110010000001001010110011111111001110100011011110111001111110001110011100101100000100000001011111110000000100000111110111110110111010110010100001001111110111000011011",
	#"111111100101000101110110001110101000010001100010110010001001010010111010001011011010100000010111000000010010110000000000100110110100001001010010111111111001010000011011110111001110110011010011100111000100110000011010111001000000100000011111010111110110010110010100001001111110111000010001",
	"111111100101000101111110001101001000010001100010110010001001010011111010111011100010101000110101000000110101110100011000000110010000001001010010011111111001110100111011110101001110110011110011100111010010101110011001111110000000100000111110001111111110010001011010010011111111111000010001",	#candidate 4 (bat conv)
	#"1100101000010001001100000011001100001111100010111000100010000100101111100011010101001000110100100000000011110100001111001001100111111111011100101001101000011110",
	#"1100000001011001001101000001100100001111100010111000100010000100001011110110110111011110000100100000000011110100101011101000110111011111100000011001101100001110",
	"1000100000010001001100000001100100001111100010111000100010000100111011000011010111011110000100100000000011111100101001001010101011001111001100101001101001000001"	# candidate 5 (numeric best)
]


#---- NSGA SETTINGS ------------------------------

model = models.models[sys.argv[1]]	# pick model from model file using given name

nsga_max_gens = 128				# number of generations in NSGA-II
nsga_pop_size = 128 			# how many combinations there are, must be even
nsga_p_cross = 0.98				# mutation crossover probability
nsga_fn_evals = 63				# how many evaluations to take median of (take n/2 + 1 entry)
nsga_bpp = 32 					# bits / precision to use for each parameter
nsga_cross_breed = False		# allow algorithms to change during the process
nsga_free_range = False			# sets all parameter ranges to 0.5 +- 0.5 instead of predefined (needs more gens)
nsga_cross_head_tail = True		# uses head/tail crossover instead of random bits

model_pop_count = [22, 16]		# pop size for method
model_generations = [22, 16]	# generations used in method

stage1 = [
			model['rand'],		# two initial values, being either random in search space
			model['estimates']	# or EM-based initial estimates
		 ]
stage2 = [	
			{
				"algo": None,	# numeric only
				"params":[]
			},
			{	
				"algo":	firefly.search,	# list swarm types as well as ranges for input params
				"params":[ [0.95, 0.05], [0.96, 0.04] ]
			},
			{						# contains all constant parameters in terms of 
				"algo":	pso.search, # [A, B] where the passed value is A + B*f(bits)
				"params":[ [0.5, 0.15], [0.1, 0.05], [0.1, 0.05] ]
			},
			{
				"algo": cuckoo.search,
				"params":[ [0.97, 0.03], [0.97, 0.03], [0.25, 0.1] ]
			},
			{
				"algo": bat.search,
				"params":[ [0.1, 0.1], [0.9, 0.1], [0.9, 0.08], [0.8, 0.15] ]
			},
			{
				"algo": bee.search,
				"params":[ [0.33, 0.10], [0.33, 0.10], [0.50, 0.25], [0.50, 0.25], [0.20, 0.10] ]
			},
			{
				"algo":	pollination.search,
				"params":[ [0.7, 0.2] ]
			},
			{
				"algo": fish.search,
				"params":[ [0.125, 0.075], [4, 2], [0.125, 0.075], [2, 1]  ]
			},
			{
				"algo": wolf.search,
				"params":[ [0.5, 0.1], [0.85, 0.15], [0.5, 0.1], [0.9, 0.075] ]
			}
			]
stage3 = [
			models.nelder_mead, # numeric methods to finalize optimizations
			models.powell, 
			models.cg,
			models.bfgs, 
			models.lbfgsb, 
			models.tnc, 
			models.slsqp 
			#models.cobyla, 	# derivative-based methods
			#models.dogleg, 
			#models.trustncg
		]

if sys.argv[1] == "Weibull":	# newtons method not applicable to covariate
	stage3.append(models.NM)

#---------------begin NSGA-II---------------------


import numpy as np
import time


def search(max_gens, pop_size, p_cross):
	'''
	Main routine - creates population of algorithms,
	runs initial sorting and new population, iteratively
	runs re-evaluation until convergence
	'''
	param_count = len(max(stage2, key = lambda x: len(x["params"]))["params"])

	total_bit_count = nsga_bpp * (param_count + 5)
											# get total number of bits for individual chromosome

	pop = [	{"bitstring":random_bitstring(total_bit_count)} 	# generate initial random population
			for i in range(pop_size)]		# bitstring consists of algo bits, parameter bits, pop size bits, and gen count bits
	
	for i, p in enumerate(pop):							# truncate extra bits in bitstring to ensure proper crossover
		a1, a2, a3, prms = decode(p['bitstring'])
		numparams = len(prms) - 4
		needed_bits = (5 + numparams) * nsga_bpp
		p['bitstring'] = p['bitstring'][:needed_bits]

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

		if not verbose:
			print("\t{}/{}     ".format(gen+1, max_gens),end="\r")

		union = pop + children				# sort union of parent/child pops, create new pop based on least dominated
		fronts = fast_nondominated_sort(union)
		parents = select_parents(fronts, pop_size)
										
		selected = [better(parents[np.random.randint(pop_size)], parents[np.random.randint(pop_size)])
					for i in range(pop_size)]
											# select well-fit parents at random
		
		pop = children.copy()
		children = reproduce(selected, pop_size, p_cross)

		calculate_measures(children)

		for candidate in children:	# candidate is tuple of bitstring, array of results
			outfile.write( " ".join( [candidate['bitstring'], str(candidate['objectives'][0]), str(candidate['objectives'][1])] ) )
			outfile.write(" ")
		outfile.write("\n")
											# take a snapshot of bitstrings and objective values for plotting

	union = pop + children
	fronts = fast_nondominated_sort(union)
	parents = select_parents(fronts, pop_size)	# finalize population


def evaluate_predef():
	pop = [{"bitstring":x} for x in predef]
	calculate_measures(pop)		# evaluate chromosome rather than letting nsga tune it

	for candidate in pop:	# candidate is tuple of bitstring, array of results
		outfile.write( " ".join( [candidate['bitstring'], str(candidate['objectives'][0]), str(candidate['objectives'][1])] ) )
		outfile.write(" ")
	outfile.write("\n")



def random_bitstring(num_bits):
	'''
	Generates a binary string of some length,
		no limit on size 
	'''
	return "".join([str(np.random.randint(2)) for i in range(num_bits)]).zfill(num_bits)


def calculate_measures(pop):
	'''
	Given a population member, decodes and evaluates the
	time and accuracy of its algorithms 
	'''
	for index, p in enumerate(pop):			# find fitness of each member of a population in order to find pareto-fitness
		
		alg_1, alg_2, alg_3, params = decode(p["bitstring"])

		pop_size = len(params[3]) 			# get size of population to make new one
		expanded = list(params) 			# convert it to a list to replace it

		runtime, error = 0, 0
		runtimes = []
		errors = []

		#print(p['bitstring'])

		for i in range(nsga_fn_evals):		# get average of N runs

			expanded[3] = alg_1(pop_size)# generate a population with size pop_size

			
			params = tuple(expanded) 		# convert back to tuple

			stime = time.time()				# stopwatch swarm alg
			lst =  alg_2(*params) if alg_2 != None else params[3]
											# run swarm estimator if enabled, otherwise just take pop as is

			candidate = min(lst, key = model['objective'])
											# get most-fit particle, check if it converges

			result, conv = alg_3(candidate)	# do newton's method or ECM, etc on result for runtime
											# don't care about the actual accuracy of that convergence

			runtimes.append(time.time() - stime)
			errors.append((model['objective'](result) / model['result']))	# store runtimes and errors in order to take median
		

		runtimes.sort()
		errors.sort()
		if "-h" in sys.argv:
			print('r-----------------------')
			print(runtimes)
			print('e-----------------------')
			print(errors)
			print()
		p["objectives"] = [runtimes[int(nsga_fn_evals/2 + 1)], errors[int(nsga_fn_evals/2 + 1)]]
											# take median runtime and errors
	if verbose:
		col = int(os.popen('stty size', 'r').read().split()[1])
		print(' '*col, end='\r')


def decode(bitstring, model = None):
	'''
	Takes some bit-pattern (some member of the population),
	gives back the corresponding algorithms (phase 1, 2, 3)
	as well as the parameters passed to ph1
	'''

	if model == None:
		model = models.models[sys.argv[1]]

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
	return fronts	# send back groups of algorithms of different tiers


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
	bitstring_sw = bitstring[nsga_bpp:2*nsga_bpp]
	if rate == None:					# flip bits at random according to rate
		rate = 1/len(bitstring)
	child = ""
	for i in range(len(bitstring)):
		bit = bitstring[i]				# only change fn bits if crossbreed is enable
		child += str(1-int(bit)) if (np.random.random_sample()<rate) else bit

	if (not nsga_cross_breed):	# if swarm doesnt mutate
		child = child[:nsga_bpp] + bitstring_sw + child[2*nsga_bpp:]

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
if __name__ == "__main__":

	verbose = "-v" in sys.argv
	outfile = open(sys.argv[2] ,"a+")

	if "-p" in sys.argv:
		evaluate_predef()
	else:
		t = time.time()
		pop = search(nsga_max_gens, nsga_pop_size, nsga_p_cross)
		t = time.time() - t

	outfile.close()

	print("done!")
	