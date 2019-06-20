
import sys, os 				# import all algorithms in sub-directory
for directory in os.listdir("algorithms"):
	sys.path.insert(0, f"algorithms/{directory}")
import firefly, pso, pollination
import models

#---- NSGA SETTINGS ------------------------------

nsga_max_gens = 50			# number of generations in NSGA-II
nsga_pop_size = 100 		# how many combinations there are
nsga_p_cross = 0.98			# mutation crossover probability
nsga_fn_evals = 25 			# how many evaluations to average
nsga_bits_per_param = 8 	# bits / precision to use for each parameter

model_objective = models.RLLWei
model_dimensions = 2
model_search_space = [[0, 1] for x in range(model_dimensions)]	# 3 dimensional co-variate search space
model_pop_count = 25											# pop size for method
model_generations = 25											# generations used in method
model_init_pop = [[404/1000000, 1] for x in range(model_pop_count)]	# None 		# set to none if randomized

stage1 = [	{	"algo":	firefly.search,	#formatted by algorithm then variable parameters base +- diff
				"params":[
					[0.95, 0.05],
					[0.96, 0.04]
				]
			},
			{
				"algo":	pso.search,
				"params":[
					[0.5, 0.15],
					[0.1, 0.05],
					[0.1, 0.05]
				]
			},
			{
				"algo":	pollination.search,
				"params":[
					[0.7, 0.2]
				]
			}]
stage2 = [lambda x: x]*2				# s2 algos; return a list of vectors
stage3 = [lambda x: min(model_objective(y) for y in x)]*2		# s3 algos; return a margin of error


#---------------begin NSGA-II---------------------

import numpy as np
import time

def decode(bitstring, bpp):
	params = (model_objective, model_search_space, model_generations, model_init_pop, model_pop_count)

	s1, s2, s3 = int(np.log2(len(stage1))), int(np.log2(len(stage2))), int(np.log2(len(stage3)))
	i1 = int(bitstring[0:s1], 2)
	i2 = int(bitstring[s1:s1+s2], 2)
	i3 = int(bitstring[s1+s2:s1+s2+s3], 2)

	param_bits = bitstring[s1+s2+s3:]
	for idx, prm in enumerate(stage1[i1]["params"]):
		this_bits = param_bits[idx * bpp :][:8]	# get part of bitstring corresponding to parameter

		rng = int(this_bits, 2) / (2 ** bpp - 1)		# get range of 0-1 that 000... -> 111... is
		diff = 2 * rng - 1								# transform it from -1 <-> 1
		params += (prm[0] + prm[1]*diff,)				# append param +- diff to list of parameters

			# this returns a function that evals s1, hands it to s2, then s3, which returns margin of error
	return (lambda x: stage3[i3]( stage2[i2]( stage1[i1]["algo"]( *x ) ) )), params

def name(bitstring, bpp):

	params = ()

	s1, s2, s3 = int(np.log2(len(stage1))), int(np.log2(len(stage2))), int(np.log2(len(stage3)))
	i1 = int(bitstring[0:s1], 2)
	i2 = int(bitstring[s1:s1+s2], 2)
	i3 = int(bitstring[s1+s2:s1+s2+s3], 2)

	param_bits = bitstring[s1+s2+s3:]
	for idx, prm in enumerate(stage1[i1]["params"]):
		this_bits = param_bits[idx * bpp :][:8]	# get part of bitstring corresponding to parameter

		rng = int(this_bits, 2) / (2 ** bpp - 1)		# get range of 0-1 that 000... -> 111... is
		diff = 2 * rng - 1								# transform it from -1 <-> 1
		params += (prm[0] + prm[1]*diff,)				# append param +- diff to list of parameters

	return [stage1[i1]["algo"].__module__, params]

def calculate_objectives(pop, fn_count, bpp):
	for p in pop:                       # find fitness of each member of a population in order to find pareto-fitness
		fn, params = decode(p["bitstring"], bpp)
		runtime = 0
		errorsum = 0
		for i in range(fn_count):
			stime = time.time()
			errorsum += fn(params)
			runtime += time.time() - stime
											# objectives are average function runtime and error
											# future: consider making error a list and checking std dev
			p["objectives"] = [runtime/fn_count, errorsum/fn_count]


def random_bitstring(num_bits):       	# generate some n-length string of random bits
	return str(bin(np.random.randint(2**num_bits)))[2:].zfill(num_bits)


def point_mutation(bitstring, rate=None):
	if rate == None:                    # basic genetic mutation, common to all gen methods
		rate = 1/len(bitstring)
	child = ""
	for i in range(len(bitstring)):
		bit = bitstring[i]
		child += str(1-int(bit)) if (np.random.random_sample()<rate) else bit
	return child


def crossover(parent1, parent2, rate):
	if np.random.random_sample() >= rate:                # basic crossover common to all gen methods
		return ""+parent1
	child = ""
	for i in range(len(parent1)):
		child += parent1[i] if np.random.random_sample() < 0.5 else parent2[i]
	return child


def reproduce(selected, pop_size, p_cross):
	children = []                       # generate new children population based off of parent population
	for i, p1 in enumerate(selected):   # with crossovers and mutation
		p2 = selected[i+1] if i%2 == 0 else selected[i-1]
		if i == (len(selected) - 1):
			p2 = selected[0]
		child = {}
		child["bitstring"] = crossover(p1["bitstring"], p2["bitstring"], p_cross)
		child["bitstring"] = point_mutation(child["bitstring"])
		children.append(child)
		if len(children) >= pop_size:
			break
	return children


def dominates(p1, p2):                # used to find whether one population is more dominant
	for i in range(len(p1["objectives"])):
		if p1["objectives"][i] > p2["objectives"][i]:
			return False
	return True


def fast_nondominated_sort(pop):
	fronts = [[]]                       # generate a list of fronts for the next population
	for p1 in pop:                      # lower "rank" value indicates dominated by less and therefore better
		p1["dom_count"], p1["dom_set"] = 0, []
		for p2 in pop:
			if dominates(p1, p2):
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
			for p2 in p1["dom_set"]:
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


def calculate_crowding_distance(pop):
	for p in pop:                       # use crowding distance to preserve diverse options
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


def better(x, y):                     # pick by rank then by crowding distance
	if ("dist" in x.keys()) and (x["rank"] == y["rank"]):
		return x if (x["dist"] > y["dist"]) else y
	return x if (x["rank"] < y["rank"]) else y


def select_parents(fronts, pop_size):
	for f in fronts:
		calculate_crowding_distance(f)
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


def search(fn_evals, max_gens, pop_size, p_cross, bpp):

	paramcount = len(max(stage1, key = lambda x: len(x["params"]))["params"])

	fnbits = int(np.log2(len(stage1)) + np.log2(len(stage2)) + np.log2(len(stage3)))	# calculate total bits used for picking fn
	pop = [{"bitstring":random_bitstring(fnbits + bpp * paramcount)} for i in range(pop_size)]

	calculate_objectives(pop, fn_evals, bpp)
	fast_nondominated_sort(pop)

	selected = [better(pop[np.random.randint(pop_size)], pop[np.random.randint(pop_size)]) for i in range(pop_size)]
	children = reproduce(selected, pop_size, p_cross)

	calculate_objectives(children, fn_evals, bpp)
	print(" > starting generations")

	for gen in range(max_gens):

		union = pop + children
		fronts = fast_nondominated_sort(union)
		parents = select_parents(fronts, pop_size)

		selected = [better(parents[np.random.randint(pop_size)], parents[np.random.randint(pop_size)]) for i in range(pop_size)]
		pop = children
		children = reproduce(selected, pop_size, p_cross)

		calculate_objectives(children, fn_evals, bpp)

		print(" > gen = {}, fronts = {}".format(gen+1, len(fronts)))

	union = pop + children
	fronts = fast_nondominated_sort(union)
	parents = select_parents(fronts, pop_size)
	return parents

pop = search(nsga_fn_evals, nsga_max_gens, nsga_pop_size, nsga_p_cross, nsga_bits_per_param)
print("done!")

pop.sort(key = lambda x: 1500*x["objectives"][0] + x["objectives"][1])
for p in pop:
	print("{} {}".format('\t'.join([str(x) for x in p["objectives"]]), name(p["bitstring"], nsga_bits_per_param)))
