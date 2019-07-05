# TODO - run cross-validation on other datasets after finding final convergent model parameters

import sys, os, csv			# import all algorithms in sub-directory
import matplotlib.pyplot as plt
from models import models

for directory in os.listdir("algorithms"):
	sys.path.insert(0, f"algorithms/{directory}")
import firefly, pso, pollination, bat, cuckoo, bee, fish



#---- NSGA SETTINGS ------------------------------

nsga_max_gens = 30			# number of generations in NSGA-II
nsga_pop_size = 64 			# how many combinations there are, must be even
nsga_p_cross = 0.98			# mutation crossover probability
nsga_fn_evals = 16 			# how many evaluations to average
nsga_bits_per_param = 8 	# bits / precision to use for each parameter
nsga_cross_breed = False	# allow algorithms to change during the process

model = models[sys.argv[1]]
model_pop_count = 18		# pop size for method
model_generations = 20		# generations used in method


stage1 = [	#{	"algo":	firefly.search,	#formatted by algorithm then variable parameters base +- diff
			#	"params":[
			#		[0.95, 0.05],
			#		[0.96, 0.04]
			#	]
			#},
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
			}]
stage2 = [lambda x: x]*2	# todo: change phase 1 to phase 2
stage3 = [lambda lst:min([model["objective"](x)/model["result"] for x in lst])]*2
							# stage 3 calculates margin of error used as stat

#---------------begin NSGA-II---------------------

import numpy as np
import time

def decode(bitstring, bpp):

	model_objective = model["objective"]				# get specific parameters for model
	model_search_space = [model["search_space"] for i in range(model["dimensions"])]
	model_init_pop = [									# generate population, pseudorandom
		[t[0] + t[1]*np.random.uniform()*np.sign(np.random.uniform()-0.5) for t in model["estimates"]]
		for x in range(model_pop_count)]
														# turn params into tuple
	params = (model_objective, model_search_space, model_generations, model_init_pop, model_pop_count)

	s1, s2, s3 = int(np.ceil(np.log2(len(stage1)))), int(np.ceil(np.log2(len(stage2)))), int(np.ceil(np.log2(len(stage3))))
	i1 = int(bitstring[0:s1], 2) % len(stage1)			# get bits for each stage in bitstring
	i2 = int(bitstring[s1:s1+s2], 2) % len(stage2)
	i3 = int(bitstring[s1+s2:s1+s2+s3], 2) % len(stage3)

	param_bits = bitstring[s1+s2+s3:]					# separate bits for variables
	for idx, prm in enumerate(stage1[i1]["params"]):
		this_bits = param_bits[idx * bpp :][:nsga_bits_per_param]	# get part of bitstring corresponding to parameter

		rng = int(this_bits, 2) / (2 ** bpp - 1)		# get range of 0-1 that 000... -> 111... is
		diff = 2 * rng - 1								# transform it from -1 <-> 1
		params += (prm[0] + prm[1]*diff,)				# append param +- diff to list of parameters

	return stage1[i1]["algo"], stage2[i2], stage3[i3], params 	# send back algorithms & parameters for timing

def calculate_objectives(pop, fn_count, bpp):
	for index, p in enumerate(pop):                       # find fitness of each member of a population in order to find pareto-fitness
		f1, f2, f3, params = decode(p["bitstring"], bpp)

		runtime = 0
		errorsum = 0
		for i in range(fn_count):			# get average of N runs for timing
			stime = time.time()
			lst =  f1(*params)				# get time of swarm algorithm as well as values for error calc
			runtime += time.time() - stime

			newerror = f3(lst)
			if np.isnan(newerror):
				newerror = float('inf')		# possible hotfix for NaN results, please investigate and fix
			errorsum += newerror

			p["objectives"] = [runtime/fn_count, errorsum/fn_count]
											# average runtime and error becomes objectives
		print(f"   > {f1.__module__}: {index+1} / {len(pop)}                  " , end = '\r')
	print('                              ', end= '\r')

def random_bitstring(num_bits):       	# generate some n-length string of random bits
	return "".join([str(np.random.randint(2)) for i in range(num_bits)]).zfill(num_bits)
										# cannot use bin(random) as that has 64-bit max


def point_mutation(bitstring, rate=None):

	fnbits = int(np.ceil(np.log2(len(stage1))))

	if rate == None:					# flip bits at random according to rate
		rate = 1/len(bitstring)
	child = ""
	for i in range(len(bitstring)):
		bit = bitstring[i]				# only change fn bits if crossbreed is enable
		if nsga_cross_breed or (i > fnbits):
			child += str(1-int(bit)) if (np.random.random_sample()<rate) else bit
		else:
			child += bit
	return child


def crossover(parent1, parent2, rate):
	if np.random.random_sample() >= rate:
		return ""+parent1
	child = ""
	for i in range(len(parent1)):		# get random bits from either parent
		child += parent1[i] if np.random.random_sample() < 0.5 else parent2[i]
	return child


def reproduce(selected, pop_size, p_cross, bit_count):
	children = []                       # generate new children population based off of parent population
	for i, p1 in enumerate(selected):   # with crossovers and mutation

		if(not nsga_cross_breed):
			count = 1
			while True:					# if no cross breed, get two parents with same algorithms
				p2 = selected[(i+count)%len(selected)]
				if p2["bitstring"][:bit_count] == p1["bitstring"][:bit_count]:
					break
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


def dominates(p1, p2):                # used to find whether one population is more dominant
	for i in range(len(p1["objectives"])): # all traits are no worse than other
		if p1["objectives"][i] > p2["objectives"][i]:	
			return False
	return True


def fast_nondominated_sort(pop):
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


def better(x, y):                     # decide which of two is best, for parent selection
	if ("dist" in x.keys()) and (x["rank"] == y["rank"]):
		return x if (x["dist"] > y["dist"]) else y
	return x if (x["rank"] < y["rank"]) else y


def select_parents(fronts, pop_size):
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


def search(fn_evals, max_gens, pop_size, p_cross, bpp):
	starting_time = time.time()
	paramcount = len(max(stage1, key = lambda x: len(x["params"]))["params"])

	algbits = int(np.ceil(np.log2(len(stage1))))
	fnbits = int(np.ceil(np.log2(len(stage1))) + np.ceil(np.log2(len(stage2))) + np.ceil(np.log2(len(stage3))))
						# calculate total bits used for picking fn
	pop = [{"bitstring":random_bitstring(fnbits + bpp * paramcount)} for i in range(pop_size)]

	print(" > begin initial objectives & sort")
	calculate_objectives(pop, fn_evals, bpp)
	fast_nondominated_sort(pop)

	selected = [better(pop[np.random.randint(pop_size)], pop[np.random.randint(pop_size)]) for i in range(pop_size)]
	children = reproduce(selected, pop_size, p_cross, algbits)

	print(" > objectives of first child pop")
	calculate_objectives(children, fn_evals, bpp)
	print(" > starting generations")

	for gen in range(max_gens):

		union = pop + children
		fronts = fast_nondominated_sort(union)
		parents = select_parents(fronts, pop_size)

		selected = [better(parents[np.random.randint(pop_size)], parents[np.random.randint(pop_size)]) for i in range(pop_size)]
		pop = children
		children = reproduce(selected, pop_size, p_cross, algbits)

		calculate_objectives(children, fn_evals, bpp)

		print(f"    > gen = {gen+1}\tfronts = {len(fronts)}\telapsed={round(time.time()-starting_time, 1)}")

	union = pop + children
	fronts = fast_nondominated_sort(union)
	parents = select_parents(fronts, pop_size)
	return parents

pop = search(nsga_fn_evals, nsga_max_gens, nsga_pop_size, nsga_p_cross, nsga_bits_per_param)
print("done!\n")




# ----------- VISUALIZATION -------------------------------------------------------------------------


pop.sort(key = lambda x: 1+x["objectives"][0] + x["objectives"][1])
print('\033[4m' + "runtime         error margin    algo    params(algo-specific)" + '\033[0m')
for p in pop:
	r, r2, r3, params = decode(p["bitstring"], nsga_bits_per_param)
	n = r.__module__
	c = "b" if n == "bat" else "m" if n == "pso" else "r" if n == "firefly" else "k" if n == "cuckoo" else "y" if n == "bee" else ""

	print('\t'.join([str(round(x,8)) for x in p["objectives"]]).zfill(10), end="\t")
	print(n, end="\t")
	#outparams = r(*params)

	#print(outparams, end="\t")
	for i in params[5:]:
		print(round(i,3), end="\t")
	print()
	plt.plot(p["objectives"][0],p["objectives"][1],c+"o")

with open('output_populations.csv','w') as csvfile:
	writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
	writer.writerow([f"NSGA GENS: {nsga_max_gens}",f"NSGA POP: {nsga_pop_size}",f"NSGA CROSSOVER: {nsga_p_cross}",f"NSGA AVG EVALS: {nsga_fn_evals}",f"NSGA CROSSBREED: {nsga_cross_breed}"])
	writer.writerow([f"MODEL: {[x for x in models if models[x]==model][0]}", f"MODEL POP: {model_pop_count}", f"MODEL GENS: {model_generations}"])
	writer.writerow([])
	writer.writerow(["algorithm", "runtime", "1+epsilon (error)", "best candidate's model parameters"])
	for p in pop:
		r, _, _, params = decode(p["bitstring"], nsga_bits_per_param)
		rs = r(*params)
		objs = [model["objective"](x) for x in rs]
		best = objs.index(min(objs))
		writer.writerow([r.__module__, p["objectives"][0],p["objectives"][1], rs[best]])

plt.title("1+ep vs time")
plt.show()