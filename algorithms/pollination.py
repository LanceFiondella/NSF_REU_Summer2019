from math import sqrt, pi, exp, sin
from scipy.stats import gamma, levy
import numpy as np

def random_solution(search_space):
	return np.array([x[0] + random() * (x[1]-x[0]) for x in search_space])

def search(objective, search_space, max_gen, pop = None, switch_prob = 0.8):
	if pop == None:									# create/fit new population
		pop = [{'vector':random_solution(search_space)} for x in range(pop_count)]
	else:
		pop = [{'vector':np.array(x)} for x in pop]

	for p in pop:
		p['objective'] = objective(p['vector'])		# calculate initial objectives of population

	best = min(pop, key = lambda x: x['objective'])	# keep track of best for local search

	l = levy()										# levy generator for random walk

	for gen in range(max_gen):
		for f in pop:
			if np.random.uniform() < switch_prob:
				step = np.array([l.ppf(np.random.uniform()) for x in f["vector"]])	# create levy walk
				diff = best['vector'] - f['vector']
				f['vector'] = f['vector'] + step * diff			# do levy walk in direction of best
			else:
				ep = np.random.uniform()						# global - pick two points and move using their difference
				j, k = np.random.choice(pop), np.random.choice(pop)
				diff = k['vector'] - j['vector']
				f['vector'] = f['vector'] + ep * diff 
				
			for i in range(len(f["vector"])):		# keep it inside search space
				f["vector"][i] = min(f["vector"][i], search_space[i][1])
				f["vector"][i] = max(f["vector"][i], search_space[i][0])
						
			f['objective'] = objective(f['vector'])	
				
		best = min(pop, key = lambda x: x['objective'])
	return [x['vector'] for x in pop]

if __name__ == "__main__":
	obj = lambda x: x[0]**2 + x[1]**2 + x[2]**2
	ss = [[-3,3],[-3,3],[-3,3]]
	pop = [[np.random.uniform(),np.random.uniform(),np.random.uniform()] for i in range(25)]
	print(search(obj, ss, 30, pop))