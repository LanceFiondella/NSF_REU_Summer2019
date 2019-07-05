import math
import numpy as np
import random

'''
Euclidean Distance: Calculates distance between two points in space

Takes: 
	vector:        Two lists of coordinates

Returns: 
	x**2 + y**2 : Float that is the Euclidean Distance 
'''

def euc_2d(c1,c2):
	return round(math.sqrt((c1[0] - c2[0])**2.0 + (c1[1] - c2[1])**2.0))




'''
Objective Function: Calculates fitness value

Takes: 
	vector:        List of coordinates

Returns: 
	x**2 + y**2 : Float that is the fitness value
'''

def objective_function(vector):
	x = vector[0]
	y = vector[1]
	return x ** 2 + y ** 2  # Sphere Function




'''
Random Vector: Caculates random vector value between minmax

Takes: 
	minmax:        List of lists containing min and max values of each dimension

Returns: 
	rand_vec: List of coordinates for single point in space
'''

def random_vector(minmax):
	rand_vec = np.array([])
	for i in range(0, len(minmax)):
		rand_vec = np.append(rand_vec, minmax[i][0] + ((minmax[i][1] - minmax[i][0]) * random.random()))
	return rand_vec




'''
Initialize Population: Creates random population of size s

Takes:
	minmax:         List of lists containing min and max values of each dimension
	nests:          Integer of size of population

Returns: 
	pop: List of dictionaries containing the randomized population
'''

def init_population(p):
	pop = [{'vector': x} for x in p]
	for i in pop:
		i['fitness'] = objective_function(i['vector'])
	return pop


'''
Generate Visual Range: Creates a the bounds for visual range for a wolf solution

Takes:
	wolf:           Dictionary of main solution which is the center for the visual range
	r:              Float of radius of search area
	search_space:   List of Lists of bounds for main search area 

Returns: 
	visual_range = List of List of bounds for a solution's search area
'''

def gen_visual_range(wolf, r , search_space):
	visual_range = ([])
	for i in wolf['vector']:
		visual_range.append([i-r, i+r])
	for i in range(len(visual_range)):
		if visual_range[i][0] < search_space[i][0]: visual_range[i][0] = search_space[i][0]
		if visual_range[i][0] > search_space[i][1]: visual_range[i][0] = search_space[i][1]
		if visual_range[i][1] < search_space[i][0]: visual_range[i][1] = search_space[i][0]
		if visual_range[i][1] > search_space[i][1]: visual_range[i][1] = search_space[i][1]
	return visual_range



'''
Generate Within Range: Finds all solutions within the visual range of a wolf

Takes:
	current:        Dictionary solution of the searching wolf
	pop:            List of dictionaries of whole population
	visual_range:   List of lists of bounds for a solution's search area

Returns: 
	visual_range = List of dictionaries of all wolves in search space
'''
def gen_within_range(current, pop, visual_range):
	true = 0;
	cand_within_range=([])
	for i in range(len(pop)):
		for j in range(len(pop[i]['vector'])):
			if (pop[i]['vector'][j] > visual_range[j][1]) | (pop[i]['vector'][j] < visual_range[j][0]) | np.array_equal(current['vector'], pop[i]['vector']):
				true = 0
				break
			else:
				true = true + 1
		if true > 0:
			cand_within_range.append(pop[i])
		true = 0;
	return cand_within_range




'''
Brownian Motion: Generates random walk for a solution

Takes:
	obective:       Function of objective 
	wolf:           Solution that will preform the walk
	alpha:          Float that serves as a step multiplier 
	r:              Float that is the radius of the visual range
	search_space:   Lists of Lists that shows bounds of the problem in each dimension 

Returns: 
	new_pos = Dictionary of new solution after the walk
'''

def bm_motion(objective, wolf, alpha, r ,search_space):
	new_pos = {}
	vector = np.array([])

	for i in range(len(wolf['vector'])):
		v = wolf['vector'][i] + (alpha * r * np.random.normal())
		if v < search_space[i][0]: v = search_space[i][0]
		if v > search_space[i][1]: v = search_space[i][1]
		vector = np.append(vector, v)
	new_pos['vector'] = vector
	new_pos['fitness'] = objective(new_pos['vector'])
	return new_pos




'''
Escape Local: A step towards the best solution in the visual range

Takes:
	obective:       Function of objective 
	best:           Dictionary of best solution defined by the fitness in the visual range
	cand:           Dictionary of solution that is doing the moving
	visual_range:   List of Lists that shows bounds of a wolf's visual range
	search_space:   List of Lists that shows bounds of the problem in each dimension 

Returns: 
	new_pos = Dictionary of new solution after the walk
'''

def escape_local(objective, best, cand, visual_range, search_space):
	new_pos = {}
	s = []
	for i in visual_range:
		s.append([0,i[0]])
	term1 = np.multiply(best['vector'], np.exp(-euc_2d(best['vector'],cand['vector'])**2))
	term2 = np.subtract(best['vector'],cand['vector'])
	term3 = random_vector(s)
	new_pos['vector'] = cand['vector'] + np.multiply(term1,term2) + term3

	#Check Within search_space
	for i in range(len(new_pos['vector'])):
		if new_pos['vector'][i] < search_space[i][0]: new_pos['vector'][i] = search_space[i][0]
		if new_pos['vector'][i] > search_space[i][1]: new_pos['vector'][i] = search_space[i][1]
	new_pos['fitness'] = objective(new_pos['vector'])
	return new_pos




'''
Escape Global: A leap outside the visual range but within the search_space

Takes:
	obective:       Function of objective 
	cand:           Dictionary of solution that is doing the moving
	alpha:          Float of alpha value that is a constant multiplier of step size
	step:           Float of step value that is a constant multiplier of step size
	visual_range:   List of Lists that shows bounds of a wolf's visual range
	search_space:   List of Lists that shows bounds of the problem in each dimension 

Returns: 
	new_pos = Dictionary of new solution after the jump
'''


def escape_global(objective, cand, alpha, step, visual_range, search_space):
	new_pos = {}
	s = []
	for i in range(len(search_space)):
		mx = visual_range[i][0]
		mn =(search_space[i][1])/2
		if mn > mx:
			mx, mn = mn, mx
		s.append([mn, mx])
	steps = np.multiply(alpha * step, random_vector(s))
	new_pos['vector'] = steps + cand['vector']
	for i in range(len(new_pos['vector'])):
		if new_pos['vector'][i] < search_space[i][0]: new_pos['vector'][i] = search_space[i][0]
		if new_pos['vector'][i] > search_space[i][1]: new_pos['vector'][i] = search_space[i][1]
	new_pos['fitness'] = objective(new_pos['vector'])
	return new_pos





'''
Search: Main search routine for Wolf Search

Takes: 
	objective: objective_function()
	search_space:           List of lists containing min and max values of each dimension
	max_generations:        Integer that holds how many iterations 
	population:             List of dictionaries that contain solutions
	wolves:                 Integer that holds the population size
	r:                      Float that holds the visual range radius
	alpha:                  Float that holds the alpha step multiplier
	step:                   Float that holds the step multiplier used in the global escape
	pa:                     Float that holds the percentage of a global escape to occur for each wolf

Returns: 
	final_positions: List of lists that have the final vectors for each in the population during the last iteration
'''
def search(objective, search_space, max_generations, population=None, wolves=100, r=0.5, alpha=1.0, step = 0.5, pa=0.90):
	#Generates random pop if none is passed
	population = init_population(population)

	#Main loop that happens for max_generations times
	for gen in range(max_generations):

		#Iterate through each solution each iteration
		for wolf in range(len(population)):
			#First preform the BM to walk the current solution in a random direction. If better solution is found, replace it with the original
			cand = bm_motion(objective, population[wolf], alpha, r, search_space)
			if cand['fitness'] < population[wolf]['fitness']:
				population[wolf] = cand
			#Generate a visual range for the current candidate wolf and select all in the population that are in that visual range
			visual_range = gen_visual_range(population[wolf], r, search_space)
			cand_within_range = gen_within_range(population[wolf], population, visual_range)
			#If there are no other solution in range, do the BM again
			if not cand_within_range:
				cand = bm_motion(objective, population[wolf], alpha, r, search_space)
				if cand['fitness'] < population[wolf]['fitness']:
					population[wolf] = cand
			#Else move towards the best solution in the visual range
			else:
				best_in_range = min(cand_within_range, key = lambda x: x["fitness"])
				cand = escape_local(objective, best_in_range, population[wolf], visual_range, search_space)
				if cand['fitness'] < population[wolf]['fitness']:
					population[wolf] = cand
			#If the random occurence passes, jump out of the visual range to do a global search
			if random.random() > pa:
				population[wolf] = escape_global(objective, population[wolf], alpha, step, visual_range, search_space)

	final_positions = []
	for i in population: final_positions.append(i['vector'])
	return final_positions




if __name__ == "__main__":
	objective = lambda x: x[0]**2 + x[1]**2
	search_space = [[-5,5],[-5,5]]
	max_gen = 10
	pop = [[np.random.uniform(-5,5),np.random.uniform(-5,5)] for x in range(75)]
	npop = search(objective, search_space, max_gen, pop)
	print("done!")
	print(min(npop, key = lambda x: objective(x)))
