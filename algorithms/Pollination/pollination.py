from random import random, randrange
from math import sqrt, pi, exp
import numpy as np
import scipy

def objective(vector):
	return sum([x**2 for x in vector])

def random_solution(search_space):
	return [x[0] + random() * (x[1]-x[0]) for x in search_space]

def levy(d):
	beta = 3/2
	sigma = (scipy.special.gamma(1 + beta) * np.sin(np.pi * beta / 2) / (scipy.special.gamma((1+beta)/2)*beta * 2**((beta-1)/2)))**(1/beta)
	u = (1 + random()*(d-1))*sigma
	v = 1 + random()*(d-1)
	return 0.01 * (u/np.abs(v)**(1/beta))

def search(search_space, pop_count, max_gen, switch_prob):
	pop = [random_solution(search_space) for x in range(pop_count)]
	best = min(pop, key = objective)
	for gen in range(max_gen):
		for f in pop:
			if random() < switch_prob:
				levy_vector = [levy(search_space[0][1]) for x in search_space]
				levy_vector = np.multiply(levy_vector, np.subtract(best, f))
				f = list(np.add(f, levy_vector) )
			else:
				sm = np.multiply(random(), np.subtract(xj, xk))
				f = list(np.add(f, sm) )
		best = min(pop, key = objective)
	print(best)
	return pop

search_space = [[-10, 10] for x in range(1)]
pcount = 50
max_gen = 30
switch_prob = 0.5

search(search_space, pcount, max_gen, switch_prob)