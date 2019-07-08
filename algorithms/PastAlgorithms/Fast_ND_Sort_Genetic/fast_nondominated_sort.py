from random import randrange, random
import matplotlib.pyplot as plt
from math import exp, sqrt, cos, pi

def objective1(vector):               # one of two objective functions, maximizes at (0, 0, ..., 0)
	#return 0
	return sum([x**2 for x in vector])

def objective2(vector):               # other function, maximizes at (2, 2, ..., 2)
	#return sum([-20*exp(-.2 * sqrt(.5*x**2)) - exp(.5 * cos(2 * pi * x)) + 20 + exp(1) for x in vector])
	return sum([(x-2)**2 for x in vector])

def decode(bitstring, search_space, bits_per_param):
	vector = []													# convert some bitstring to a form that can be placed into objectives
	for i, bounds in enumerate(search_space):
		off, sm = i*bits_per_param, 0
		param = bitstring[off:(off+bits_per_param)][::-1]
		for j in range(len(param)):
			sm += ((1 if param[j] == '1' else 0) * (2 ** j))
		mn, mx = bounds
		vector.append( mn + ((mx-mn)/((2 ** bits_per_param) - 1)) * sm )
	return vector

def random_bitstring(num_bits):       # generate some n-length string of random bits
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

def calculate_objectives(pop, search_space, bits_per_param):
	for p in pop:                       # find fitness of each member of a population in order to find pareto-fitness
		p["vector"] = decode(p["bitstring"], search_space, bits_per_param)
		p["objectives"] = [objective1(p["vector"]), objective2(p["vector"])]

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

def weighted_sum(x):
	return sum(x["objectives"])

def search(search_space, max_gens, pop_size, p_cross, bits_per_param = 16):
	pop = [{"bitstring":random_bitstring(len(search_space)*bits_per_param)} for i in range(pop_size)]
	calculate_objectives(pop, search_space, bits_per_param)
	fast_nondominated_sort(pop)
	selected = [better(pop[randrange(pop_size)], pop[randrange(pop_size)]) for i in range(pop_size)]
	children = reproduce(selected, pop_size, p_cross)
	calculate_objectives(children, search_space, bits_per_param)
	for gen in range(max_gens):
		union = pop + children
		fronts = fast_nondominated_sort(union)
		parents = select_parents(fronts, pop_size)
		selected = [better(parents[randrange(pop_size)], parents[randrange(pop_size)]) for i in range(pop_size)]
		pop = children
		children = reproduce(selected, pop_size, p_cross)
		calculate_objectives(children, search_space, bits_per_param)
		parents.sort(key = weighted_sum)
		best = parents[0]
		best_s = "[x={}, objs={}]".format(best["vector"], best["objectives"])

		obj1s.append(best["objectives"][0])
		obj2s.append(best["objectives"][1])
		distances.append(best["dist"])

		print(" > gen = {}, fronts = {}, best = {}".format(gen+1, len(fronts), best_s))
	union = pop + children
	fronts = fast_nondominated_sort(union)
	parents = select_parents(fronts, pop_size)
	return parents


problem_size = 1
search_space = [[-10, 10] for i in range(problem_size)]

max_gens = 50
pop_size = 100
p_cross = 0.98

distances = []
obj1s = []
obj2s = []

pop = search(search_space, max_gens, pop_size, p_cross)
print("done!")







convfig = plt.figure()
plt.plot(distances)
plt.plot(obj1s)
plt.plot(obj2s)
plt.title("Objectives & Crowding Distance vs Iteration")
plt.xlabel('Iteration Count')
plt.legend(["Crowding Distance", "Obj 1 (Sphere 0,0)", "Obj 2 (Sphere 2, 2)"])

plt.show()
