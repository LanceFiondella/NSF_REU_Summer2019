import math
from random import randint
import numpy as np
import random

'''
Objective Function: Caculates fitness value

Takes: 
	vector:        List of coordinates

Returns: 
	x**2 + y**2 : Float that is the fitness value
'''
def objective_function(vector):
	x = vector[0]
	y = vector[1]
	return x**2 + y**2          #Sphere Function

'''
Random Vector: Caculates random vector value between minmax

Takes: 
	minmax:        List of lists containing min and max values of each dimension
	
Returns: 
	rand_vec: List of coordinates for single point in space
'''
def random_vector(minmax):
	rand_vec = []
	for i in range(0, len(minmax)):
		rand_vec.append(minmax[i][0] + ((minmax[i][1] - minmax[i][0]) * random.random()))
	return rand_vec

'''
Iinitialize Population: Creates random population of size nests

Takes:
	minmax:         List of lists containing min and max values of each dimension
	nests:          Integer of size of population

Returns: 
	pop: List of dictionaries containing the randomized population
'''
def init_population(minmax, nests):
	pop = [dict() for x in range(0, nests)]
	for i in range(0, len(pop)):
		pop[i]['vector'] = random_vector(minmax)
	for i in pop:
		i['fitness']= objective_function(i['vector'])
	return pop

# This Levy Functions Made By
# https://github.com/YutaUme/CS/blob/master/individual.py
# Adapted
'''
Levy Flight: Creates new Levy Solution

Takes: 
	objective:          objective_function()
	search_space:       List of lists containing min and max values of each dimension
	rc:                 Single dictionary for the random candidate that the solution will be derived from
	ld:                 Integer lambda value
	ps:                 Integer value for how many dimensions
	alpha:              Float value for alpha multiplier to step size
	
Returns: 
	levy_solution: Single dictionary for the created Levy Solution
'''
def levy_flight(objective, search_space ,rc, ld,  ps, alpha):
	levy_solution = {'vector':[]}
	sigma1 = np.power((math.gamma(1 + ld) * np.sin((np.pi * ld) / 2)) \
					  / math.gamma((1 + ld) / 2) * np.power(2, (ld - 1) / 2), 1 / ld)
	sigma2 = 1
	u = np.random.normal(0, sigma1, size=ps)
	v = np.random.normal(0, sigma2, size=ps)
	step = u / np.power(np.fabs(v), 1 / ld)

	#Keeping the resultant vector contained
	for i in range(0, ps):
		x = step[i] + alpha *rc['vector'][i]
		if x > search_space[i][1]:
			levy_solution['vector'].append(search_space[i][1])
		elif x < search_space[i][0]:
			levy_solution['vector'].append(search_space[i][0])
		else:
			levy_solution['vector'].append(x)
	levy_solution['fitness'] = objective(levy_solution['vector'])

	return levy_solution


'''
Search: Main search routine for Cuckoo Search

Takes: 
	objective: objective_function()
	search_space:           List of lists containing min and max values of each dimension
	max_generations:        Integer that holds how many iterations 
	population:             List of dictionaries that contain solutions
	nests:                  Integer that holds the population size
	problem_size:           Integer that holds number of dimensional spaces
	ld:                     Float that holds the lambda value
	alpha:                  Float that holds the alpha value
	pa:                     Float that holds the percentage of population thrown out each iteration 

Returns: 
	final_positions: List of lists that have the final vectors for each in the population during the last iteration
'''
def search (objective, search_space, max_generations, population = None, nests = 100, ld = 1.0, alpha =1.0, pa = 0.25):

	if population == None:
		population = init_population(search_space, nests)
	else:
		population = [{'vector':x} for x in population]
		for j in population:
			j['fitness'] = objective_function(j['vector'])

	problem_size = len(search_space)
	throw = int(pa * len(population))
	keep =len(population) - throw
	for gen in range(0, max_generations):
		for test in range(0, nests):
			rand_candidate = population[randint(0,len(population)-1)]
			levy_solution = levy_flight(objective,search_space, rand_candidate,ld, problem_size,alpha)
			rand2 = population[randint(0, len(population) - 1)]
			if rand2['fitness'] > levy_solution['fitness']:
				population.remove(rand2)
				population.append(levy_solution)
		population = sorted(population, key=lambda i: i['fitness'])
		partial_nests = population[0:keep]
		new_nests = init_population(search_space, throw)
		population = partial_nests + new_nests

	final_positions = []
	for i in population:
		final_positions.append(i['vector'])
	return final_positions


def main():
	problem_size = 2                                                #Dimensional Space
	search_space = []
	for i in range(0, problem_size): search_space.append([0, 5])
	max_generations = 100
	nests = 100
	ld = 1
	alpha = 1
	pa = 0.25
	population = None
	best = search(objective_function, search_space, max_generations, population, nests, problem_size, ld, alpha, pa)
	for i in best:
		print(i)

if __name__== "__main__":
	main()
