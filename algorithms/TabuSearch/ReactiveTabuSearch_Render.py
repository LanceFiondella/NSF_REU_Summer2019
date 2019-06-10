#!/usr/bin/python

# http://www.cleveralgorithms.com/nature-inspired/stochastic/reactive_tabu_search.html
# Python port of Jason Brownlee's Reactive-Tabu-Search salesman implementation
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


def is_tabu(edge, tabu_list, iteration, prohib_period):
	for entry in tabu_list:
		if entry["edge"] == edge:
			return entry["iter"] >= iteration - prohib_period
	return False

def make_tabu(tabu_list, edge, iteration):
	for entry in tabu_list:
		if entry["edge"] == edge:
			entry["iter"] = iteration
			return entry
	entry = {"edge": edge, "iter": iteration}
	tabu_list.append(entry)
	return entry

def to_edge_list(permutation):
	lst = []
	for index, c1 in enumerate(permutation):
		c2 = permutation[0] if index == len(permutation)-1 else permutation[index+1]
		if c1 > c2:
			c1, c2 = c2, c1
		lst.append([c1, c2])
	return lst

def equivalent(el1, el2):
	for e in el1:
		if not e in el2:
			return False
	return True


def generate_candidate(best, cities):
	candidate = {}
	candidate["vector"], edges = stochastic_two(best["vector"])
	candidate["cost"] = cost(candidate["vector"], cities)
	return candidate, edges

def get_candidate_entry(visited_list, permutation):
	edges = to_edge_list(permutation)
	for entry in visited_list:
		if equivalent(edges, entry["edgelist"]):
			return entry
	return None

def store_permutation(visited_list, permutation, iteration):
	entry = {	"edgelist": to_edge_list(permutation),
				"iter": iteration,
				"visits": 1	}
	visited_list.append(entry)
	return entry

def sort_neighborhood(candidates, tabu_list, prohib_period, iteration):
	tabu, admissable = [], []
	for a in candidates:
		if (is_tabu(a[1][0], tabu_list, iteration, prohib_period) or is_tabu(a[1][1], tabu_list, iteration, prohib_period)):
			tabu.append(a)
		else:
			admissable.append(a)
	return tabu, admissable

def search(cities, max_cand, max_iter, increase, decrease):
	stime = time()
	current = {"vector": random_permutation(cities)}
	current["cost"] = cost(current["vector"], cities)
	best = current
	tabu_list, prohib_period = [], 1
	visited_list, avg_size, last_change = [], 1, 0
	for iteration in range(max_iter):
		candidate_entry = get_candidate_entry(visited_list, current["vector"])
		if candidate_entry != None:
			repetition_interval = iteration - candidate_entry["iter"]
			candidate_entry["iter"] = iteration
			candidate_entry["visits"] += 1
			if repetition_interval < 2*(len(cities) - 1):
				avg_size = 0.1 * (iteration - candidate_entry["iter"]) + 0.9 * avg_size
				prohib_period = float(prohib_period) * increase
				last_change = iteration
		else:
			store_permutation(visited_list, current["vector"], iteration)
		if (iteration - last_change) > avg_size:
			prohib_period = max(1, prohib_period*decrease)
			last_change = iteration
		candidates = [generate_candidate(current, cities) for x in range(max_cand)]
		candidates.sort(key = lambda x: x[0]["cost"])
		tabu, admis = sort_neighborhood(candidates, tabu_list, prohib_period, iteration)
		if len(admis) < 2:
			prohib_period = len(cities) - 2
			last_change = iteration
		current, best_move_edges = tabu[0] if len(admis) == 0 else admis[0]
		if len(tabu) > 0:
			tf = tabu[0][0]
			if (tf["cost"] < best["cost"]) and (tf["cost"] < current["cost"]):
				current, best_move_edges = tabu[0]
		for edge in best_move_edges:
			make_tabu(tabu_list, edge, iteration)
		if candidates[0][0]["cost"] < best["cost"]:
			best = candidates[0][0]

		time_between.append(time() - stime)
		iteration_set.append(best["vector"])
		costs.append(best["cost"])
		if verbose:
			print(" > iteration {}, tenure = {}, best = {}".format(iteration+1, round(prohib_period), best["cost"]))

	return best


# configure algorithm
dataset = sets[sys.argv[1]]
max_iterations = int(sys.argv[2])
max_candidates = 50
increase = 1.3
decrease = 0.9

#set up plot
iteration_set = []
costs = []
time_between = []

# execute search

verbose = bool(int(sys.argv[3]))
start_time = time()
best = search(dataset, max_candidates, max_iterations, increase, decrease)
timetaken = time() - start_time

print("Search completed in {} seconds - best solution: c = {}".format(timetaken, best["cost"]))



# draw plot animation
if verbose:
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
