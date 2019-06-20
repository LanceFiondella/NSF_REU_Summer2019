from random import random, randrange, choice
from math import sqrt, pi, exp, sin
from scipy.stats import gamma, levy
import numpy as np

def random_solution(search_space):
	return [x[0] + random() * (x[1]-x[0]) for x in search_space]

def search(objective, search_space, max_gen, pop = None, switch_prob = 0.8):
	if pop == None:
		pop = [random_solution(search_space) for x in range(pop_count)]
	best = min(pop, key = objective)
	l = levy()
	for gen in range(max_gen):
		for f in pop:
			if random() < switch_prob:
				step = [l.ppf(random()) for x in f]
				diff = np.subtract(best, f)
				term = np.multiply(step, diff)
				f = np.add(f, term)
			else:
				ep = random()
				j, k = choice(pop), choice(pop)
				diff = np.subtract(k, j)
				term = np.multiply(ep, diff)
				f = np.add(f, term) 
				
		best = min(pop, key = objective)
	print(best)
	return pop