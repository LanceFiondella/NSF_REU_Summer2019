#  The following code is a python adaptation of Simulated Annealing originally written in the Ruby language by Jason
#  Brownlee
#  Uses the Travelling Salesperson Problem
from math import *
from random import *
import matplotlib.pyplot as plt
from matplotlib import animation
iteration_set = []


def euc_2d(c1, c2):
    return round( ((c1[0] - c2[0])**2 + (c1[1] - c2[1])**2)**(1/2) )

def cost(permutation, cities):
    distance = 0
    for index, value in enumerate(permutation):
        neighbor = permutation[0] if (index == len(permutation)-1) else permutation[index+1]
        distance += euc_2d(cities[value], cities[neighbor])
    return distance


def random_permutation(cities):
    citiessize=len(cities)
    perm = [x for x in range(citiessize)]
    for i in range(citiessize):
        r = randrange(citiessize - i) + i
        perm[r], perm[i] = perm[i], perm[r]
    print(perm)
    return perm



def stochastic_two_opt(perm):
    permsize = len(perm)

    c1, c2 = randrange(permsize), randrange(permsize)
   #  print(c1, c2)
    exclude = [c1]  # exclude some point and its neighbors
    exclude.append(permsize - 1 if c1 == 0 else c1 - 1)
    exclude.append(0 if c1 == permsize - 1 else c1 + 1)

    while c2 in exclude:  # ensure c2 isn't
        c2 = randrange(permsize)

    if (c2 < c1):  # swap order if out of order
        c1, c2 = c2, c1

    perm[c1:c2] = perm[c1:c2][::-1]  # reverse that subsection of list
    return perm


def create_neighbor(current, cities):
    candidate = {}
    candidate["vector"] = current["vector"][:]
    candidate["vector"] = stochastic_two_opt(candidate["vector"])
    candidate["cost"]=cost(candidate["vector"],cities)
    return candidate


def should_accept(candidate, current, temp):
    if candidate["cost"] <= current["cost"]:
        return True
    return exp((current["cost"] - candidate["cost"])/temp) > random()


def search(cities, max_iter, max_temp, temp_change):
    current = {"vector": random_permutation(cities)}
    current["cost"] = cost(current["vector"], cities)
    temp, best = max_temp, current
    for i in range(max_iter):
        candidate = create_neighbor(current, cities)
        # print(candidate["vector"])
        temp = temp * temp_change
        if should_accept(candidate, current, temp):
            current = candidate
        if candidate["cost"] < best["cost"]:
            best = candidate
        if (i + 1) % 10 == 0:
            # print(f"{best['vector']}")
            print(f" > iteration {i+1}, temp={temp}, best={best['cost']}")
            #pass
        iteration_set.append(best["vector"].copy())
    return best


# main
berlin52=[[565,575],[25,185],[345,750],[945,685],[845,655],
   [880,660],[25,230],[525,1000],[580,1175],[650,1130],[1605,620],
   [1220,580],[1465,200],[1530,5],[845,680],[725,370],[145,665],
   [415,635],[510,875],[560,365],[300,465],[520,585],[480,415],
   [835,625],[975,580],[1215,245],[1320,315],[1250,400],[660,180],
   [410,250],[420,555],[575,665],[1150,1160],[700,580],[685,595],
   [685,610],[770,610],[795,645],[720,635],[760,650],[475,960],
   [95,260],[875,920],[700,500],[555,815],[830,485],[1170,65],
   [830,610],[605,625],[595,360],[1340,725],[1740,245]]
# algorithm configuration
max_iterations = 2000
max_temp = 100000.0
temp_change = .98
# execute the algorithm
best = search(berlin52, max_iterations, max_temp, temp_change)
print(f"Done, Best Solution: c={best['cost']}, v={best['vector']}")


#  Josh TSP Render

def animate(i):							# function is ran every frame with the frame number as the argument
	if(i >= len(iteration_set)):		# do nothing if the number is over the number of frame data lists
		return lines
	x, y = [berlin52[x][0] for x in iteration_set[i]], [berlin52[x][1] for x in iteration_set[i]]
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

x, y = [berlin52[x][0] for x in iteration_set[0]], [berlin52[x][1] for x in iteration_set[0]] # find initial dot/line positions
dots = plt.plot(x,y,'o')			# plot them - dots don't move, so this one function call is sufficient
lines = plt.plot(x,y)


#  create an animation with 1.5*the amount of frames to give a delay at the end, 10ms between frames (can be changed), repeating
anim = animation.FuncAnimation(fig, animate, frames=range(int(max_iterations*1.5)), interval = 10, repeat=True)
anim.save('traveling_salesman.gif', fps=30, writer='imagemagick')	# save the animation as a gif, however this doesn't work on windows (need imagemagick)
print("done saving")
plt.show()
