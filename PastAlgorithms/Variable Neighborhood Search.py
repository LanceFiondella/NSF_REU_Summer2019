# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from mpl_toolkits import mplot3d
from matplotlib import animation, rc
from IPython.display import HTML, Image
from matplotlib.animation import FuncAnimation,FFMpegFileWriter
#! /usr/bin/env python3
from matplotlib import animation
from common import path_cost, random_permutation
import random
import matplotlib.pyplot as plt
"""
2.7
Variable Neighborhood Search involves iterative exploration of larger and larger
neighborhoods for a given local optima until an improvement is located after
which time the search across expanding neighborhoods is repeated. The strategy
is motivated by three principles: 1) a local minimum for one neighborhood
structure may not be a local minimum for another neighborhood structure, 2) a
global minimum is local minimum for all possible neighborhood structures, and
3) local minima are relatively close to global minima for many problem classes.
"""

def stochastic_two_opt(perm):
	randlimit = len(perm) - 1
	c1, c2 = random.randint(0, randlimit), random.randint(0, randlimit)
	exclude = [c1]
	exclude.append(randlimit if c1 == 0 else c1 -1)
	exclude.append(0 if c1 == randlimit else c1 + 1)

	while c2 in exclude:
		c2 = random.randint(0, randlimit)

	c1, c2 = c2, c1 if c2 < c1 else None
	perm[c1:c2] = perm[c1:c2][::-1]
	return perm

def local_search(best, cities, max_no_improv, neighborhood_size):
	count = 0

	while count < max_no_improv:
		candidate = {}
		candidate["vector"] = [v for v in best["vector"]]

		for _ in range(neighborhood_size):
			stochastic_two_opt(candidate["vector"])

		candidate["cost"] = path_cost(candidate["vector"], cities)

		if candidate["cost"] < best["cost"]:
			count, best = 0, candidate
		else:
			count += 1

	return best

def search(cities, neigborhoods, max_no_improv, max_no_improv_ls):
	best = {}
	best["vector"] = random_permutation(cities)
	best["cost"] = path_cost(best["vector"], cities)
	iter_, count = 0, 0

	while count < max_no_improv:
		for neigh in neighborhoods:
			candidate = {}
			candidate["vector"] = [v for v in best["vector"]]

			for _ in range(neigh):
				stochastic_two_opt(candidate["vector"])

			candidate["cost"] = path_cost(candidate["vector"], cities)
			candidate = local_search(candidate, cities, max_no_improv_ls, neigh)
			framevalues.append(candidate["vector"])
			print("> iteration #%s, neigh=%s, best=%s" % (iter_ + 1, neigh, best["cost"]))
			iter_ += 1

			if candidate["cost"] < best["cost"]:
				best, count = candidate, 0
				print("New best, restarting neighborhood search")
				break
			else:
				count += 1
			
	return best

if __name__ == "__main__":
	from common import berlin52
	max_no_improv = 50
	max_no_improv_ls = 70
	neighborhoods = list(range(20))
	framevalues = []
	best = search(berlin52, neighborhoods, max_no_improv, max_no_improv_ls)
	print("Done. Best Solution: c=%s, v=%s" % (best["cost"], best["vector"]))

#    plt.plot(berlin52)
#    plt.ylabel('some numbers')
#    plt.show()


def animate(i):							# function is ran every frame with the frame number as the argument
	if(i >= len(best)):		# do nothing if the number is over the number of frame data lists
		return lines
	x, y = [berlin52[x][0] for x in framevalues[i]], [berlin52[x][1] for x in framevalues[i]]
	for line in lines:					# look up the new line positioning in the saved frame data
		line.set_data(x,y) 				# assign the values to the lines in the graph
	return lines

def findaxes(dset): 					# supply this with the set of actual data values
	xs = [dset[x][0] for x in range(len(dset))] # all it does is calculate the min/max and return some fitting rectangle around it
	ys = [dset[x][1] for x in range(len(dset))]
	xs.sort()
	ys.sort()
	return [	int(1.1*xs[0] - 0.1*xs[-1]), int(1.1*xs[-1] - 0.1*xs[0]),
			int(1.1*ys[0] - 0.1*ys[-1]), int(1.1*ys[-1] - 0.1*ys[0])]


fig, ax = plt.subplots()			# create plot objects
ax = plt.axis(findaxes(berlin52)) 	# set the viewing area to include all data points

x, y = [berlin52[x][0] for x in framevalues[0]], [berlin52[x][1] for x in framevalues[0]] # find initial dot/line positions
dots = plt.plot(x,y,'o')			# plot them - dots don't move, so this one function call is sufficient
lines = plt.plot(x,y)

# create an animation with 1.5*the amount of frames to give a delay at the end, 10ms between frames (can be changed), repeating
anim = animation.FuncAnimation(fig, animate, frames=range(int(len(framevalues)*1.5)), interval = 10, repeat=True)
#anim.save('traveling_salesman.gif', fps=30, writer='imagemagick')	# save the animation as a gif, however this doesn't work on windows (need imagemagick)
mywriter = FFMpegFileWriter(fps=25,codec="libx264")
anim.save("test.mp4", writer=mywriter)




