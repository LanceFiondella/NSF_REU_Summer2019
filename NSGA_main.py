from random import randrange, random
from time import time
import numpy as np
# also import other algorithm files

stage1 = [] # list of references to stage 1 algorithms
stage2 = [] # stage 2 algorithms
stage3 = [] # stage 3 if more than 1

# stage 1 and 2 algorithms will need to return a reference to the points they operate on
# stage 3 algorithms must return some measure of error
# needs at least 1 for each

def weighted_sum(x):	# function that outright chooses a "best" candidate
						# by default minimizes sum of error & runtime, change to value one over other
	return sum(x["objectives"])

def calculate_objectives(pop, fn_count):
	for p in pop:                       # find fitness of each member of a population in order to find pareto-fitness
		fn, params = decode(p["bitstring"])
		runtime = 0
		errorsum = 0
		for i in range(fn_count):
			start = time()

			errorsum += fn(params)
			runtime += (time() - start)
										# objectives are average function runtime and error
										# future: consider making error a list and checking std dev
		p["objectives"] = [1,1]#[runtime/fn_count, errorsum/fn_count]

def decode(bitstring):

	s1, s2, s3 = int(np.log2(len(stage1))), int(np.log2(len(stage2))), int(np.log2(len(stage3)))
	i1 = 0 if s1 == 0 else int(bitstring[0:s1], 2)
	i2 = 0 if s1 == 0 else int(bitstring[s1:s1+s2], 2)
	i3 = 0 if s1 == 0 else int(bitstring[s1+s2:s1+s2+s3], 2)

			# future: change none to a decoded list of parameters based on bitstring
			# this returns a function that evals s3 on s2's results, s2's results are based on s1's
	return (lambda x: stage3[i3]( stage2[i2]( stage1[i1]( x ) ) )), None

def random_bitstring(num_bits):       	# generate some n-length string of random bits
	return str(bin(randrange(2**num_bits)))[2:].zfill(num_bits)

def point_mutation(bitstring, rate=None):
	if rate == None:                    # basic genetic mutation, common to all gen methods
		rate = 1/len(bitstring)
	child = ""
	for i in range(len(bitstring)):
		bit = bitstring[i]
		child += str(1-int(bit)) if (random()<rate) else bit
	return child

def crossover(parent1, parent2, rate):
	if random() >= rate:                # basic crossover common to all gen methods
		return ""+parent1
	child = ""
	for i in range(len(parent1)):
		child += parent1[i] if random() < 0.5 else parent2[i]
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

def search(fn_evals, max_gens, pop_size, p_cross):
	fnbits = int(np.log2(len(stage1)) + np.log2(len(stage2)) + np.log2(len(stage3)))	# calculate total bits used for picking fn
	pop = [{"bitstring":random_bitstring(fnbits)} for i in range(pop_size)]
	print(pop)
	calculate_objectives(pop, fn_evals)
	fast_nondominated_sort(pop)
	selected = [better(pop[randrange(pop_size)], pop[randrange(pop_size)]) for i in range(pop_size)]
	children = reproduce(selected, pop_size, p_cross)

	calculate_objectives(children, fn_evals)
	for gen in range(max_gens):
		union = pop + children
		fronts = fast_nondominated_sort(union)
		parents = select_parents(fronts, pop_size)
		selected = [better(parents[randrange(pop_size)], parents[randrange(pop_size)]) for i in range(pop_size)]
		pop = children
		children = reproduce(selected, pop_size, p_cross)
		calculate_objectives(children, fn_evals)
		best = min(parents, key = weighted_sum)

		print(" > gen = {}, fronts = {}".format(gen+1, len(fronts)))

	union = pop + children
	fronts = fast_nondominated_sort(union)
	parents = select_parents(fronts, pop_size)
	return parents


max_gens = 50
pop_size = 100
p_cross = 0.98
fn_evals = 100

pop = search(fn_evals, max_gens, pop_size, p_cross)
print("done!")
