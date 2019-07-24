import numpy as np
import time, models, sys

model = models.models[sys.argv[1]]


def set_global(g_vars):
	for key in g_vars.keys():
		if not (key[:2] == "__" and key[-2:] == "__"):
			globals()[key] = g_vars[key]


def search(fn_evals, max_gens, pop_size, p_cross):
	'''
	Main routine - creates population of algorithms,
	runs initial sorting and new population, iteratively
	runs re-evaluation until convergence
	'''
	param_count = len(max(stage2, key = lambda x: len(x["params"]))["params"])

	total_stage_bits = stage_bits(0)			# get total number of bits for bitstring

	total_bit_count = total_stage_bits + (2 + param_count) * nsga_bpp
												# bitstring consists of algo bits, parameter bits
												# pop size bits, and gen count bits
	pop = [	{"bitstring":random_bitstring(total_bit_count)} 
			for i in range(pop_size)]
										
	print(" > begin initial objectives & sort")
	calculate_measures(pop)					# prelim. results pt 1: get initial pop scoring
	fast_nondominated_sort(pop)				# sort initial pop by domination

	parent_pop = [better(pop[np.random.randint(pop_size)], pop[np.random.randint(pop_size)]) for i in range(pop_size)]
	children = reproduce(parent_pop, pop_size, p_cross)
											# select a random well-fit parent pop and generate a new child pop with it

	print(" > objectives of first child pop")
	calculate_measures(children)
	print(" > starting generations")		# prelim. results pt 2: score child pop

	for gen in range(max_gens):				# begin generational aspect

		starting_time = time.time()

		'''
		t = {}								# code to visualize population percentages of each algorithm
		for i in [ x[:-3] for x in os.listdir('algorithms') if x[-3:] == '.py']:
			t[i] = 0
		for p in pop:
			m = decode(p["bitstring"])[1].__module__
			t[m] += 1
		algs.append(t)
		'''

		union = pop + children				# sort union of parent/child pops, create new pop based on least dominated
		fronts = fast_nondominated_sort(union)
		parents = select_parents(fronts, pop_size)
										
		parent_pop = [better(parents[np.random.randint(pop_size)], parents[np.random.randint(pop_size)])
					for i in range(pop_size)]
											# select well-fit parents at random
		
		pop = children.copy()
		children = reproduce(parent_pop, pop_size, p_cross)

		calculate_measures(children)
		print(f"    > gen = {gen+1} / {max_gens}\tfronts = {len(fronts)}\tdelta={round(time.time()-starting_time, 1)}")


	union = pop + children
	fronts = fast_nondominated_sort(union)
	parents = select_parents(fronts, pop_size)
	return parents


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

		results = []						

		for i in range(nsga_fn_evals):		# get average of N runs


			expanded[3] = [						# generate population, pseudorandom
				[t[0] + t[1]*np.random.uniform(-1, 1) for t in model["estimates"]]
				for x in range(pop_size)]

			
			params = tuple(expanded) 		# convert back to tuple

			print(' '*55, end='\r')
			print(f"       > {index+1}/{len(pop)} ({i+1}/{nsga_fn_evals})\t{alg_2.__module__}: {alg_3.__name__}", end='\r')


			stime = time.time()
			lst =  alg_2(*params)			# get time of swarm algorithm

			candidate = min(lst, key = model['objective'])
											# get most-fit particle, check if it converges

			result, conv = alg_3(candidate)	# do newton's method or ECM on result for runtime
											# don't care about the actual accuracy of that convergence
			runtime += time.time() - stime

			error += (model['objective'](result) / model['result'])
			
		p["objectives"] = [runtime/nsga_fn_evals, error/nsga_fn_evals]
											# metrics: total runtime, and accuracy of swarm before convergence

	print(' '*75, end= '\r')


def decode(bitstring):
	'''
	Takes some bit-pattern (some member of the population),
	gives back the corresponding algorithms (phase 1, 2, 3)
	as well as the parameters passed to ph1
	'''
	s1, s2, s3 = stage_bits(-1)				# get bits per stage for indexing
	s1_bits, s2_bits, s3_bits = bitstring[ : s1], bitstring[s1 : s1+s2], bitstring[ s1+s2 : s1+s2+s3 ]
											# get corresponding bits for each stage
	index_1 = round(bin_signum(s1_bits, (len(stage1)-1)/2, (len(stage1)-1)/2, s1))
	index_2 = round(bin_signum(s2_bits, (len(stage2)-1)/2, (len(stage2)-1)/2, s2))
	index_3 = round(bin_signum(s3_bits, (len(stage3)-1)/2, (len(stage3)-1)/2, s3))
											# get indexes of each stage, by rounding
											# stage_num / 2 +- stage_num / 2

	pop_gen_bits = bitstring[s1+s2+s3:][:nsga_bpp*2]
	pop_size = round(bin_signum(pop_gen_bits[:nsga_bpp], model_pop_count[0], model_pop_count[1]))
	gen_count = round(bin_signum(pop_gen_bits[nsga_bpp:], model_generations[0], model_generations[1]))
											# get bits corresponding to model gens and size
											# and generate pop size / generations based on values

	model_objective = model["objective"]	# get specific parameters for model
	model_search_space = [model["search_space"] for i in range(model["dimensions"])]


	params = (model_objective, model_search_space, gen_count, [0]*pop_size)
											# get parameters for swarm algorithms
	param_bits = bitstring[s1+s2+s3+len(pop_gen_bits):]
	for idx, prm in enumerate(stage2[index_2]["params"]):
		this_bits = param_bits[idx * nsga_bpp :][:nsga_bpp]
											# get part of bitstring corresponding to parameter
		if nsga_free_range:
			params += (bin_signum(this_bits, 0.5, 0.5),)
		else:
			params += (bin_signum(this_bits, prm[0], prm[1]),)

	return stage1[index_1], stage2[index_2]["algo"], stage3[index_3], params 	# send back algorithms & parameters for timing

def stage_bits(status):
	'''
	Get the bits per stage, or the total;
	function has no special purpose other
	than to save repeated usage of logs
	'''
	if status == -1:						# return individual bit counts for all stages
		return stage_bits(1), stage_bits(2), stage_bits(3)
	elif status == 0:						# get total bit count
		return sum(stage_bits(-1))
	elif status == 1:						# get bits for stage 1
		return int(np.ceil(np.log2(len(stage1))))
	elif status == 2:						# stage 2 bits
		return int(np.ceil(np.log2(len(stage2))))
	elif status == 3:						# stage 3 bits
		return int(np.ceil(np.log2(len(stage3))))

def bin_signum(bitstring, lmbda, mu, bpp = None):
	'''
	Use a bitstring to get a number within lambda +- mu
	'''
	if bpp == None:
		bpp = nsga_bpp
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
	s1, s2, s3 = stage_bits(-1)

	if rate == None:					# flip bits at random according to rate
		rate = 1/len(bitstring)
	child = ""
	for i in range(len(bitstring)):
		bit = bitstring[i]				# only change fn bits if crossbreed is enable
		if nsga_cross_breed or not (s1 < i < (s1 + s2)):
			child += str(1-int(bit)) if (np.random.random_sample()<rate) else bit
		else:
			child += bit
	return child


def reproduce(parent_pop, pop_size, p_cross):
	'''
	Generates new sub-population based on a selection
	of an initial population with mutations and crossover
	of random members - if cross-breed is disabled, 
	parents will only have same algorithms
	'''
	children = []                       # generate new children population based off of parent population
	s1, s2, s3 = stage_bits(-1)

	for i, p1 in enumerate(parent_pop):   # with crossovers and mutation

		if(not nsga_cross_breed):
			count = 1
			while True:					# if no cross breed, get two parents with same algorithms
				p2 = parent_pop[(i+count)%len(parent_pop)]
				if p2["bitstring"][s1:s1+s2] == p1["bitstring"][s1:s1+s2]:
					break				# compare similar middle bitstrings
				count += 1
		else:
			p2 = parent_pop[0] if (i == len(parent_pop) - 1) else \
					parent_pop[i+1] if i%2 == 0 else parent_pop[i-1]

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