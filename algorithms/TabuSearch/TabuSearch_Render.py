#!/usr/bin/python

# http://www.cleveralgorithms.com/nature-inspired/stochastic/tabu_search.html
# Python port of Jason Brownlee's Tabu-Search salesman implementation
# Josh Steakelum

import sys
from datasets import sets
from random import randrange
from matplotlib import pyplot as plt
from matplotlib import animation
from time import time

def euc_2d(c1, c2):								# euclidean distance from two x,y pairs
	return round( ((c1[0] - c2[0])**2 + (c1[1] - c2[1])**2)**(1/2) )

def cost(permutation, cities):					# get total distance from each point to neighbor
	distance = 0								# does NOT use random permutations
	for index, value in enumerate(permutation):
		neighbor = permutation[0] if (index == len(permutation)-1) else permutation[index+1]
		distance += euc_2d(cities[value], cities[neighbor])
	return distance


def random_permutation(cities):					#generate list of random unique numbers 0<=n<len(cities)
	perm = list(range(len(cities)))
	for i in range(len(perm)):
		r = randrange(len(perm)-i) + i 			#rand number between i and size of list
		perm[r], perm[i] = perm[i], perm[r]		#swap list elements
	return perm


def stochastic_two(parent):
	perm = parent[:]
	permsize = len(perm)

	c1, c2 = randrange(permsize), randrange(permsize)
	exclude = [c1]								#exclude some point and its neighbors
	exclude.append(permsize-1 if c1==0 else c1-1)
	exclude.append(0 if c1==permsize-1 else c1+1)

	while c2 in exclude:						#ensure c2 isn't 
		c2 = randrange(permsize)

	if(c2 < c1):								# swap order if out of order
		c1, c2 = c2, c1

	perm[c1:c2] = reversed(perm[c1:c2])			# reverse that subsection of list 
	return perm, [[parent[c1-1], parent[c1]], [parent[c2-1], parent[c2]]]


def is_tabu(permutation, tabu_list):
	for index, value in enumerate(permutation):
	 	c2 = permutation[0] if index==len(permutation)-1 else permutation[index+1]
	 	for forbidden in tabu_list:
	 		if(forbidden == [value, c2]):
	 			return True
	return False


def generate_candidate(best, tabu_list, cities):
	permutation, edges = None, None
	lastp = []
	while True:									# find a new candidate as well as its reversed edges
		permutation, edges = stochastic_two(best["vector"])
		if(not is_tabu(permutation, tabu_list)):
			break								# (repeatedly) until a unique candidate is found
												# store its "tabu" edge by indicating which values begin and end the reversal
	candidate = {"vector": permutation}
	candidate["cost"] = cost(candidate["vector"], cities)
	return candidate, edges


def search(cities, tabu_list_size, max_candidates, max_iterations):
	stime = time()
	current = {"vector": random_permutation(cities)}
	current["cost"] = cost(current["vector"], cities)
	best = current
	tabu_list = []
	for iteration in range(max_iterations):
		candidates = [generate_candidate(current, tabu_list, cities) for x in range(max_candidates)]
		candidates.sort(key=(lambda x: x[0]["cost"]))
		best_candidate = candidates[0][0]
		best_candidate_edges = candidates[0][1]
		if best_candidate["cost"] < current["cost"]:
			current = best_candidate
			if best_candidate["cost"] < best["cost"]:
				best = best_candidate
			for edge in best_candidate_edges:
				tabu_list.append(edge)
			while len(tabu_list) > tabu_list_size:
				tabu_list.pop()

		time_between.append(time() - stime)
		iteration_set.append(best["vector"])
		costs.append(best["cost"])
		if verbose:
			print(" > iteration {}, best = {}".format(iteration+1, best["cost"]))
	return best

# configure algorithm
dataset = sets[sys.argv[1]]
max_iterations = int(sys.argv[2])
tabu_list_size = 15
max_candidates = 50

#set up plot
iteration_set = []
costs = []
time_between = []

# execute search

verbose = bool(int(sys.argv[3]))
start_time = time()
best = search(dataset, tabu_list_size, max_candidates, max_iterations)
timetaken = time() - start_time
print("Search completed in {} seconds - best solution: c = {}".format( timetaken, best["cost"]))



# draw plot animation
if bool(int(sys.argv[3])):
	def animate(i):
		if(i >= len(iteration_set)):		#add pause to end of animation
			return lines
		x, y = [dataset[x][0] for x in iteration_set[i]], [dataset[x][1] for x in iteration_set[i]]
		for line in lines:
			line.set_data(x,y)
		return lines

	def findaxes(dset):
		xs = [dset[x][0] for x in range(len(dset))]
		ys = [dset[x][1] for x in range(len(dset))]
		xs.sort()
		ys.sort()
		return [int(1.1*xs[0] - 0.1*xs[-1]), int(1.1*xs[-1] - 0.1*xs[0]),
				int(1.1*ys[0] - 0.1*ys[-1]), int(1.1*ys[-1] - 0.1*ys[0])]

	fig, ax = plt.subplots()
	ax = plt.axis(findaxes(dataset))

	x, y = [dataset[x][0] for x in iteration_set[0]], [dataset[x][1] for x in iteration_set[0]]
	dots = plt.plot(x,y,'o')
	lines = plt.plot(x,y)

	roundtime = round(timetaken * 10000) / 10000
	plt.title("TSP Path, took {} seconds".format(roundtime))
	plt.xlabel("X Coord (units)")
	plt.ylabel("Y Coord (units)")

	anim = animation.FuncAnimation(fig, animate, frames=range(int(max_iterations*1.5)), interval = 10, repeat=False)
	#anim.save('traveling_salesman.gif', fps=30, writer='imagemagick')

	# plot the changing cost of the TSP algorithm
	convfig = plt.figure(2)
	plt.plot(costs)
	plt.title("Cost (distance) vs iteration")
	plt.xlabel('Iteration Count')
	plt.ylabel('Cost (distance)')

	plt.show()
